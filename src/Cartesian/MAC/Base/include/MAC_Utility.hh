#pragma once

#include <cstddef>
#include <string>

// Delete all regular files directly contained in both directories.
// Mirrors the legacy FLUID Pacific_clear_exec script.
int pacific_clear_exec( std::string const& results_directory,
                        std::string const& savings_directory ) ;

// Delete current_*.tmp, *_time.xml, *.bin, *.pbin, *.m, and *.pm files
// directly contained in directory.
int FVMatlab_clear_exec( std::string const& directory ) ;

// Delete *.pvd, *.pvtr, and *.vtr files directly contained in directory.
int FVParaview_clear_exec( std::string const& directory ) ;
