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
CMAKE_SOURCE_DIR = /home/lorraine/Documents/SMART

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lorraine/Documents/SMART

# Include any dependencies generated for this target.
include ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/depend.make

# Include the progress variables for this target.
include ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/progress.make

# Include the compile flags for this target's objects.
include ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/flags.make

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/flags.make
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o: ext/libdivsufsort/lib/divsufsort.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lorraine/Documents/SMART/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/divsufsort.o   -c /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/divsufsort.c

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/divsufsort.i"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/divsufsort.c > CMakeFiles/divsufsort64.dir/divsufsort.i

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/divsufsort.s"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/divsufsort.c -o CMakeFiles/divsufsort64.dir/divsufsort.s

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires:

.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires
	$(MAKE) -f ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/build.make ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides.build
.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.provides.build: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o


ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/flags.make
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o: ext/libdivsufsort/lib/sssort.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lorraine/Documents/SMART/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/sssort.o   -c /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/sssort.c

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/sssort.i"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/sssort.c > CMakeFiles/divsufsort64.dir/sssort.i

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/sssort.s"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/sssort.c -o CMakeFiles/divsufsort64.dir/sssort.s

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.requires:

.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.requires

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.provides: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.requires
	$(MAKE) -f ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/build.make ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.provides.build
.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.provides

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.provides.build: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o


ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/flags.make
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o: ext/libdivsufsort/lib/trsort.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lorraine/Documents/SMART/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/trsort.o   -c /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/trsort.c

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/trsort.i"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/trsort.c > CMakeFiles/divsufsort64.dir/trsort.i

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/trsort.s"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/trsort.c -o CMakeFiles/divsufsort64.dir/trsort.s

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.requires:

.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.requires

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.provides: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.requires
	$(MAKE) -f ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/build.make ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.provides.build
.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.provides

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.provides.build: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o


ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/flags.make
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o: ext/libdivsufsort/lib/utils.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lorraine/Documents/SMART/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/divsufsort64.dir/utils.o   -c /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/utils.c

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/divsufsort64.dir/utils.i"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/utils.c > CMakeFiles/divsufsort64.dir/utils.i

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/divsufsort64.dir/utils.s"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && /usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/utils.c -o CMakeFiles/divsufsort64.dir/utils.s

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.requires:

.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.requires

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.provides: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.requires
	$(MAKE) -f ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/build.make ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.provides.build
.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.provides

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.provides.build: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o


# Object files for target divsufsort64
divsufsort64_OBJECTS = \
"CMakeFiles/divsufsort64.dir/divsufsort.o" \
"CMakeFiles/divsufsort64.dir/sssort.o" \
"CMakeFiles/divsufsort64.dir/trsort.o" \
"CMakeFiles/divsufsort64.dir/utils.o"

# External object files for target divsufsort64
divsufsort64_EXTERNAL_OBJECTS =

ext/libdivsufsort/lib/libdivsufsort64.a: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o
ext/libdivsufsort/lib/libdivsufsort64.a: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o
ext/libdivsufsort/lib/libdivsufsort64.a: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o
ext/libdivsufsort/lib/libdivsufsort64.a: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o
ext/libdivsufsort/lib/libdivsufsort64.a: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/build.make
ext/libdivsufsort/lib/libdivsufsort64.a: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lorraine/Documents/SMART/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C static library libdivsufsort64.a"
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && $(CMAKE_COMMAND) -P CMakeFiles/divsufsort64.dir/cmake_clean_target.cmake
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/divsufsort64.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/build: ext/libdivsufsort/lib/libdivsufsort64.a

.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/build

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/requires: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/divsufsort.o.requires
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/requires: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/sssort.o.requires
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/requires: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/trsort.o.requires
ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/requires: ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/utils.o.requires

.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/requires

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/clean:
	cd /home/lorraine/Documents/SMART/ext/libdivsufsort/lib && $(CMAKE_COMMAND) -P CMakeFiles/divsufsort64.dir/cmake_clean.cmake
.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/clean

ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/depend:
	cd /home/lorraine/Documents/SMART && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lorraine/Documents/SMART /home/lorraine/Documents/SMART/ext/libdivsufsort/lib /home/lorraine/Documents/SMART /home/lorraine/Documents/SMART/ext/libdivsufsort/lib /home/lorraine/Documents/SMART/ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/libdivsufsort/lib/CMakeFiles/divsufsort64.dir/depend

