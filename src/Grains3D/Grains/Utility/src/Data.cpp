#include <dlfcn.h>
#include <filesystem>
#include <vector>

std::filesystem::path get_library_dir() {
  Dl_info info;
  dladdr((void *)&get_library_dir, &info); // find our own symbol
  std::filesystem::path lib_path =
      std::filesystem::path(info.dli_fname).parent_path();
  return lib_path;
}

std::filesystem::path get_dtd_dir() {
  std::filesystem::path libdir = get_library_dir();
  std::vector<std::filesystem::path> candidates = {
      libdir / "./Dtd",                  // in-tree build
      libdir / "../share/Grains/Dtd",    // install prefix
      libdir / "../../../share/Grains/Dtd", // mpi-install prefix
      libdir / "../../usr/share/Grains/Dtd" // debian multi-arch 
  };

  for (auto &p : candidates) {
    std::error_code ec;
    if (std::filesystem::exists(p, ec))
      return std::filesystem::canonical(p, ec);
  }

  throw std::runtime_error("Cannot find DTD directory near " + libdir.string());
}

std::filesystem::path get_tools_dir() {
  std::filesystem::path libdir = get_library_dir();
  std::vector<std::filesystem::path> candidates = {
      libdir / "./Tools",                  // in-tree build
      libdir / "../share/Grains/Tools",    // install prefix
      libdir / "../../../share/Grains/Tools", // mpi-install prefix (/usr/lib64/openmpi/lib/ + ../../../share/Grains/Tools)
      libdir / "../../usr/share/Grains/Tools" // debian multi-arch 
  };

  for (auto &p : candidates) {
    std::error_code ec;
    if (std::filesystem::exists(p, ec))
      return std::filesystem::canonical(p, ec);
  }

  throw std::runtime_error("Cannot find Tools directory near " + libdir.string());
}
