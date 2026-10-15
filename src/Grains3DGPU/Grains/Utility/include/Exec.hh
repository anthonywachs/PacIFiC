 #pragma once

#include <cstddef>
#include <filesystem>
#include <string>

namespace PacIFiC::Grains {

/**
 * Remove Paraview and diagnostic output files in `dir`.
 *
 * Mirrors the bash script behavior:
 *   - Removes: <prefix>*.{pvd,vtu,pvtu,vtp,pvtp}
 *   - Removes: ErreurContact*, ErreurImposedMotion*, LinkedCell*
 *
 * @param dir        Directory containing output files.
 * @param pv_prefix  Prefix used for Paraview outputs (corresponds to "$2" in script).
 * @param recursive  If true, search recursively under `dir` (default: false).
 * @return Number of files removed.
 *
 * Notes:
 *   - Missing directories / missing files are not treated as errors (rm -f style).
 *   - Only regular files are removed.
 */
std::size_t cleanup_paraview_outputs(const std::filesystem::path& dir,
                            const std::string& pv_prefix,
                            bool recursive = false);


// Deletes <root>_<name>.dat files listed by the legacy Text_clear.exec script.
// If `dir` is empty, it uses current_path().
std::size_t cleanup_text_outputs(const std::filesystem::path& dir,
                                 const std::string& filerootname);

bool cleanup_force_stats_outputs(const std::filesystem::path& dir);

} // namespace pacific::grains
