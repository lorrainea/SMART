# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.4

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/lorraine/Documents/smart

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lorraine/Documents/smart

# Utility rule file for uninstall.

# Include the progress variables for this target.
include ext/libdivsufsort/CMakeFiles/uninstall.dir/progress.make

ext/libdivsufsort/CMakeFiles/uninstall:
	cd /home/lorraine/Documents/smart/ext/libdivsufsort && /usr/local/bin/cmake -P /home/lorraine/Documents/smart/ext/libdivsufsort/CMakeModules/cmake_uninstall.cmake

uninstall: ext/libdivsufsort/CMakeFiles/uninstall
uninstall: ext/libdivsufsort/CMakeFiles/uninstall.dir/build.make

.PHONY : uninstall

# Rule to build all files generated by this target.
ext/libdivsufsort/CMakeFiles/uninstall.dir/build: uninstall

.PHONY : ext/libdivsufsort/CMakeFiles/uninstall.dir/build

ext/libdivsufsort/CMakeFiles/uninstall.dir/clean:
	cd /home/lorraine/Documents/smart/ext/libdivsufsort && $(CMAKE_COMMAND) -P CMakeFiles/uninstall.dir/cmake_clean.cmake
.PHONY : ext/libdivsufsort/CMakeFiles/uninstall.dir/clean

ext/libdivsufsort/CMakeFiles/uninstall.dir/depend:
	cd /home/lorraine/Documents/smart && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lorraine/Documents/smart /home/lorraine/Documents/smart/ext/libdivsufsort /home/lorraine/Documents/smart /home/lorraine/Documents/smart/ext/libdivsufsort /home/lorraine/Documents/smart/ext/libdivsufsort/CMakeFiles/uninstall.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/libdivsufsort/CMakeFiles/uninstall.dir/depend
