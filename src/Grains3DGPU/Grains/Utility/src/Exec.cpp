#include "Exec.hh"

#include <array>
#include <system_error>

namespace PacIFiC::Grains {

namespace fs = std::filesystem;

static bool matches_prefix_suffix(const std::string &name,
                                  const std::string &prefix,
                                  const std::string &suffix) {
  if (name.size() < prefix.size() + suffix.size())
    return false;
  if (name.rfind(prefix, 0) != 0)
    return false; // doesn't start with prefix
  if (!suffix.empty() &&
      name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
    return false;
  return true;
}

static bool remove_if_exists(const std::filesystem::path &p) {
  std::error_code ec;
  if (!std::filesystem::exists(p, ec) || ec)
    return false;
  if (!std::filesystem::is_regular_file(p, ec) || ec)
    return false;
  return std::filesystem::remove(p, ec) && !ec;
}

static std::size_t remove_matching_files_in_dir(const fs::path &dir,
                                                const std::string &prefix,
                                                const std::string &suffix,
                                                bool recursive) {
  std::size_t removed_count = 0;

  std::error_code ec;
  if (!fs::exists(dir, ec) || ec || !fs::is_directory(dir, ec) || ec) {
    return 0; // rm -f style: do nothing if dir missing/unreadable
  }

  auto try_remove = [&](const fs::path &p) {
    std::error_code st_ec;
    if (!fs::is_regular_file(p, st_ec) || st_ec)
      return;

    std::error_code rm_ec;
    if (fs::remove(p, rm_ec) && !rm_ec) {
      ++removed_count;
    }
  };

  const auto opts = fs::directory_options::skip_permission_denied;

  if (recursive) {
    for (fs::recursive_directory_iterator it(dir, opts, ec), end;
         it != end && !ec; it.increment(ec)) {
      const fs::path p = it->path();
      const std::string name = p.filename().string();
      if (matches_prefix_suffix(name, prefix, suffix)) {
        try_remove(p);
      }
    }
  } else {
    for (fs::directory_iterator it(dir, opts, ec), end; it != end && !ec;
         it.increment(ec)) {
      const fs::path p = it->path();
      const std::string name = p.filename().string();
      if (matches_prefix_suffix(name, prefix, suffix)) {
        try_remove(p);
      }
    }
  }

  return removed_count;
}

std::size_t cleanup_paraview_outputs(const fs::path &dir,
                                     const std::string &pv_prefix,
                                     bool recursive) {
  std::size_t removed = 0;

  // Paraview files: $1/$2*.{pvd,vtu,pvtu,vtp,pvtp}
  removed += remove_matching_files_in_dir(dir, pv_prefix, ".pvd", recursive);
  removed += remove_matching_files_in_dir(dir, pv_prefix, ".vtu", recursive);
  removed += remove_matching_files_in_dir(dir, pv_prefix, ".pvtu", recursive);
  removed += remove_matching_files_in_dir(dir, pv_prefix, ".vtp", recursive);
  removed += remove_matching_files_in_dir(dir, pv_prefix, ".pvtp", recursive);

  // Other diagnostics: $1/ErreurContact*, etc.
  removed += remove_matching_files_in_dir(dir, "ErreurContact", "", recursive);
  removed +=
      remove_matching_files_in_dir(dir, "ErreurImposedMotion", "", recursive);
  removed += remove_matching_files_in_dir(dir, "LinkedCell", "", recursive);

  return removed;
}

std::size_t cleanup_text_outputs(const fs::path &dir,
                                 const std::string &filerootname) {
  static constexpr std::array<const char *, 29> kTextFiles = {
      {"position_x.dat",
       "position_y.dat",
       "position_z.dat",
       "velocity_x.dat",
       "velocity_y.dat",
       "velocity_z.dat",
       "rotation_x.dat",
       "rotation_y.dat",
       "rotation_z.dat",
       "coordinationNumber.dat",
       "totalForce_x.dat",
       "totalForce_y.dat",
       "totalForce_z.dat",
       "cumulatedContactForce_x.dat",
       "cumulatedContactForce_y.dat",
       "cumulatedContactForce_z.dat",
       "instantaneousContactForce_x.dat",
       "instantaneousContactForce_y.dat",
       "instantaneousContactForce_z.dat",
       "hydroForce_x.dat",
       "hydroForce_y.dat",
       "hydroForce_z.dat",
       "slipVelocity_x.dat",
       "slipVelocity_y.dat",
       "slipVelocity_z.dat",
       "particleTemperature.dat",
       "cumulatedLubriForce_x.dat",
       "cumulatedLubriForce_y.dat",
       "cumulatedLubriForce_z.dat"}};

  std::size_t removed = 0;
  const fs::path base = dir.empty() ? fs::current_path() : dir;

  for (const char *name : kTextFiles) {
    fs::path p = base / (filerootname + "_" + name);
    if (remove_if_exists(p)) {
      ++removed;
    }
  }
  return removed;
}

bool cleanup_force_stats_outputs(const fs::path &dir) {
  const fs::path base = dir.empty() ? fs::current_path() : dir;
  fs::path p = base / "ForceStats.res";
  return remove_if_exists(p);
}

} // namespace PacIFiC::Grains
