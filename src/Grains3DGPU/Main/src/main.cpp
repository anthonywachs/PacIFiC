#include "Grains.hh"
#include "GrainsFactory.hh"
#include "ReaderXML.hh"

#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>

namespace fs = std::filesystem;

struct CliOptions
{
    std::string input_xml;
    bool        keep_temp = false;
    bool        dry_run   = false;
    bool        help      = false;
};

static void print_usage(std::ostream& os, const char* prog)
{
    os << "Usage:\n"
       << "  " << prog << " --input <file.xml> [--keep-temp] [--dry-run]\n"
       << "  " << prog << " <file.xml>\n\n"
       << "Options:\n"
       << "  -i, --input <file.xml>   Input XML file\n"
       << "      --keep-temp          Keep temporary XML generated for execution\n"
       << "      --dry-run            Parse/build only, skip simulate()\n"
       << "  -h, --help               Show this help\n";
}

static bool ends_with(const std::string& s, const std::string& suffix)
{
    return s.size() >= suffix.size()
           && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::optional<CliOptions> parse_args(int argc, char* argv[], std::string& err)
{
    CliOptions opt;

    for(int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if(a == "-h" || a == "--help")
        {
            opt.help = true;
            return opt;
        }
        else if(a == "-i" || a == "--input")
        {
            if(i + 1 >= argc)
            {
                err = "Missing value after " + a;
                return std::nullopt;
            }
            opt.input_xml = argv[++i];
        }
        else if(a == "--keep-temp")
        {
            opt.keep_temp = true;
        }
        else if(a == "--dry-run")
        {
            opt.dry_run = true;
        }
        else if(!a.empty() && a[0] == '-')
        {
            err = "Unknown option: " + a;
            return std::nullopt;
        }
        else
        {
            if(!opt.input_xml.empty())
            {
                err = "Multiple input files provided (use only one).";
                return std::nullopt;
            }
            opt.input_xml = a;
        }
    }

    if(opt.input_xml.empty() && !opt.help)
    {
        err = "No input file provided.";
        return std::nullopt;
    }
    return opt;
}

template <typename T>
static int run_grains(DOMElement* rootNode, bool dry_run)
{
    std::unique_ptr<Grains<T>> grains(GrainsFactory<T>::create(rootNode));
    if(!grains)
    {
        std::cerr << "ERROR: GrainsFactory::create returned null\n";
        return 3;
    }

    grains->initialize(rootNode);
    if(dry_run)
    {
        std::cout << "Dry run: skipping simulation\n";
        return 0;
    }

    grains->simulate();
    grains->finalize();
    return 0;
}

int main(int argc, char* argv[])
{
    std::string parse_err;
    auto        opt_maybe = parse_args(argc, argv, parse_err);

    if(!opt_maybe)
    {
        std::cerr << "ERROR: " << parse_err << "\n\n";
        print_usage(std::cerr, argv[0]);
        return 2;
    }

    CliOptions opt = *opt_maybe;
    if(opt.help)
    {
        print_usage(std::cout, argv[0]);
        return 0;
    }

    if(!ends_with(opt.input_xml, ".xml"))
    {
        std::cerr << "ERROR: input file must have .xml extension\n";
        return 2;
    }
    if(!fs::exists(opt.input_xml))
    {
        std::cerr << "ERROR: input file does not exist: " << opt.input_xml << "\n";
        return 2;
    }

    std::string filename_exe;
    DOMElement*  rootNode         = nullptr;
    bool         reader_initialized = false;
    try
    {
        filename_exe = GrainsFactory<double>::init(opt.input_xml);
        ReaderXML::initialize();
        reader_initialized = true;
        rootNode           = ReaderXML::getRoot(filename_exe);
        std::string prc      = ReaderXML::getNodeAttr_String(rootNode, "Precision");

        int rc = 0;
        if(prc == "Single")
        {
            rc = run_grains<float>(rootNode, opt.dry_run);
        }
        else
        {
            if(prc != "Double")
            {
                std::cout << "Invalid precision; defaulting to Double precision.\n";
            }
            rc = run_grains<double>(rootNode, opt.dry_run);
        }

        rootNode->getOwnerDocument()->release();
        rootNode = nullptr;
        ReaderXML::terminate();
        reader_initialized = false;

        if(!opt.keep_temp && !filename_exe.empty())
        {
            std::error_code ec;
            fs::remove(filename_exe, ec);
            if(ec)
            {
                std::cerr << "Warning: failed to remove temp file " << filename_exe << ": "
                          << ec.message() << "\n";
            }
        }
        return rc;
    }
    catch(const std::exception& e)
    {
        std::cerr << "ERROR: exception: " << e.what() << "\n";
    }
    catch(...)
    {
        std::cerr << "ERROR: unknown exception\n";
    }

    if(rootNode)
    {
        rootNode->getOwnerDocument()->release();
    }
    if(reader_initialized)
    {
        ReaderXML::terminate();
    }

    if(!opt.keep_temp && !filename_exe.empty())
    {
        std::error_code ec;
        fs::remove(filename_exe, ec);
    }
    return 3;
}
