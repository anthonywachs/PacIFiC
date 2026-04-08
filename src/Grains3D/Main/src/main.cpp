

#include <mpi.h> 

#include "Grains.hh"
#include "GrainsBuilderFactory.hh"
#include "ReaderXML.hh"
// #include "Version.hh"

#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct MpiSession {
  MpiSession(int &argc, char **&argv) { MPI_Init(&argc, &argv); }
  ~MpiSession() { MPI_Finalize(); }
  MpiSession(const MpiSession &) = delete;
  MpiSession &operator=(const MpiSession &) = delete;
};

struct CliOptions {
  std::string input_xml;
  bool keep_temp = false;
  bool dry_run = false;
  bool help = false;
  bool version = false;
};

static void print_usage(std::ostream &os, const char *prog) {
  os << "Usage:\n"
     << "  " << prog << " --input <file.xml> [--keep-temp] [--dry-run]\n"
     << "  " << prog << " <file.xml>\n"
     << "  " << prog << " --version\n\n"
     << "Options:\n"
     << "  -i, --input <file.xml>   Input XML file\n"
     << "      --keep-temp          Do not delete the generated temporary XML\n"
     << "      --dry-run            Parse/build, print info, exit before "
        "Simulation\n"
     << "  -V, --version            Print version information and exit\n"
     << "  -h, --help               Show this help\n";
}

static bool ends_with(const std::string &s, const std::string &suffix) {
  return s.size() >= suffix.size() &&
         s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}
static std::optional<CliOptions> parse_args(int argc, char *argv[],
                                            std::string &err) {
  CliOptions opt;

  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];

    if (a == "-h" || a == "--help") {
      opt.help = true;
      return opt;
    } else if (a == "-V" || a == "--version") {
      opt.version = true;
      return opt;
    } else if (a == "-i" || a == "--input") {
      if (i + 1 >= argc) {
        err = "Missing value after " + a;
        return std::nullopt;
      }
      opt.input_xml = argv[++i];
    } else if (a == "--keep-temp") {
      opt.keep_temp = true;
    } else if (a == "--dry-run") {
      opt.dry_run = true;
    } else if (!a.empty() && a[0] == '-') {
      err = "Unknown option: " + a;
      return std::nullopt;
    } else {
      // positional
      if (!opt.input_xml.empty()) {
        err = "Multiple input files provided (use only one).";
        return std::nullopt;
      }
      opt.input_xml = a;
    }
  }

  // Only require input file if not help/version
  if (opt.input_xml.empty() && !opt.help && !opt.version) {
    err = "No input file provided.";
    return std::nullopt;
  }

  return opt;
}

int main(int argc, char *argv[]) {
  MpiSession mpi(argc, argv);

  int rank = 0, nprocs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  auto rank0 = (rank == 0);

  std::string parse_err;
  auto opt_maybe = parse_args(argc, argv, parse_err);

  if (!opt_maybe) {
    if (rank0) {
      std::cerr << "ERROR: " << parse_err << "\n\n";
      print_usage(std::cerr, argv[0]);
    }
    return 2;
  }

  CliOptions opt = *opt_maybe;

  if (opt.version) {
    if (rank0) {
      // std::cout << "PacIFiC " << Version::semver << " (" << Version::full
      //           << ")\n"
      //           << "Type:   " << Version::build_type << "\n";
    }
    return 0;
  }

  if (opt.help) {
    if (rank0)
      print_usage(std::cout, argv[0]);
    return 0;
  }

  // Validate input
  if (!ends_with(opt.input_xml, ".xml")) {
    if (rank0)
      std::cerr << "ERROR: input file must have .xml extension\n";
    return 2;
  }

  if (!fs::exists(opt.input_xml)) {
    if (rank0)
      std::cerr << "ERROR: input file does not exist: " << opt.input_xml
                << "\n";
    return 2;
  }

  std::string filename_exe;

  try {
    // Create a temporary input file with the proper XML header
    filename_exe = GrainsBuilderFactory::init(opt.input_xml, rank, nprocs);

    // Creates the Grains application
    ReaderXML::initialize();
    DOMElement *rootNode = ReaderXML::getRoot(filename_exe);
    string option = ReaderXML::getNodeAttr_String(rootNode, "Type");

    Grains *grains = NULL;
    grains = GrainsBuilderFactory::create(rootNode);

    if (!grains) {
      if (rank0)
        std::cerr << "ERROR: GrainsBuilderFactory::create returned null\n";
      // Cleanup temp below
      if (rank0 && !opt.keep_temp && !filename_exe.empty())
        fs::remove(filename_exe);
      return 3;
    }

    // Initial output message
    grains->initialOutputMessage();

    // Tasks to perform before time-stepping
    grains->do_before_time_stepping(rootNode);
    ReaderXML::terminate();

    if (rank0 && !opt.keep_temp && !filename_exe.empty()) {
      std::error_code ec;
      fs::remove(filename_exe, ec);
      if (ec) {
        std::cerr << "Warning: failed to remove temp file " << filename_exe
                  << ": " << ec.message() << "\n";
      }
    }

    if (opt.dry_run) {
      if (rank0)
        std::cout << "Dry run: skipping simulation\n";
      delete grains;
      return 0;
    }

    // Run the simulation
    grains->Simulation();

    // Tasks to perform after time-stepping
    grains->do_after_time_stepping();

    delete grains;

    if (Component::getNbCreatedComponents()) {
      std::cerr << "Warning: " << Component::getNbCreatedComponents()
                << " component(s) not properly destroyed on proc " << rank
                << "\n";
    }

  } catch (const std::exception &e) {
    if (rank0)
      std::cerr << "ERROR: exception: " << e.what() << "\n";
    // Best effort cleanup of temp file
    if (rank0 && !opt.keep_temp && !filename_exe.empty()) {
      std::error_code ec;
      fs::remove(filename_exe, ec);
    }
    return 3;
  } catch (...) {
    if (rank0)
      std::cerr << "ERROR: unknown exception\n";
    if (rank0 && !opt.keep_temp && !filename_exe.empty()) {
      std::error_code ec;
      fs::remove(filename_exe, ec);
    }
    return 3;
  }

  return 0;
}
