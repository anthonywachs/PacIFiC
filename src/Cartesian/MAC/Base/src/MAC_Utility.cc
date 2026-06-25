#include <MAC_Utility.hh>

#include <filesystem>
#include <initializer_list>
#include <utility>

static bool name_matches(std::string const &name, std::string const &prefix,
                         std::string const &suffix) {
  return name.size() >= prefix.size() + suffix.size() &&
         name.compare(0, prefix.size(), prefix) == 0 &&
         name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static int remove_matching_files(
    std::string const &directory,
    std::initializer_list<std::pair<std::string, std::string>> patterns) {
  namespace fs = std::filesystem;

  try {
    fs::path const dir(directory);
    if (!fs::is_directory(dir)) {
      return 0;
    }

    for (fs::directory_entry const &entry : fs::directory_iterator(dir)) {
      if (!entry.is_regular_file()) {
        continue;
      }

      std::string const name = entry.path().filename().string();
      for (auto const &pattern : patterns) {
        if (name_matches(name, pattern.first, pattern.second)) {
          fs::remove(entry.path());
          break;
        }
      }
    }
  } catch (fs::filesystem_error const &) {
    return 1;
  }

  return 0;
}

//---------------------------------------------------------------------------
int pacific_clear_exec(std::string const &results_directory,
                       std::string const &savings_directory)
//---------------------------------------------------------------------------
{
  /*
  #!/bin/bash
  for file in $1/*
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  for file in $2/*
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done
  */

  int result = 0;
  result |= remove_matching_files(results_directory, {{"", ""}});
  result |= remove_matching_files(savings_directory, {{"", ""}});
  return result;
}

//---------------------------------------------------------------------------
int FVMatlab_clear_exec(std::string const &directory)
//---------------------------------------------------------------------------
{
  /*
  #!/bin/bash

  # This script deletes all current_*.tmp, *_time.xml, .bin, .pbin, .m
  # and .pm files in directory $1

  # current_*.tmp file
  for file in $1/current_*.tmp
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  # *_time.xml files
  for file in $1/*_time.xml
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  # BIN files
  for file in $1/*.bin
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  # PBIN files
  for file in $1/*.pbin
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  # M files
  for file in $1/*.m
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  # PM files
  for file in $1/*.pm
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done
  */
  
  return remove_matching_files(directory, {{"current_", ".tmp"},
                                           {"", "_time.xml"},
                                           {"", ".bin"},
                                           {"", ".pbin"},
                                           {"", ".m"},
                                           {"", ".pm"}});
}

//---------------------------------------------------------------------------
int FVParaview_clear_exec(std::string const &directory)
//---------------------------------------------------------------------------
{
  /*
  #!/bin/bash
  # This script deletes all .pvd, .pvtr and .vtr files in directory $1

  # PVD file
  for file in $1/*.pvd
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  # PVTR files
  for file in $1/*.pvtr
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done

  # VTR files
  for file in $1/*.vtr
    do
      if [ -f $file ]
        then
          rm -f $file
      fi
    done
  */
  return remove_matching_files(directory,
                               {{"", ".pvd"}, {"", ".pvtr"}, {"", ".vtr"}});
}
