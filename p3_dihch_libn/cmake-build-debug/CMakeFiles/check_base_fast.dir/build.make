# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /tmp/tmp.2zHWrVoTcO

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /tmp/tmp.2zHWrVoTcO/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/check_base_fast.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/check_base_fast.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/check_base_fast.dir/flags.make

CMakeFiles/check_base_fast.dir/test_accuracy.c.o: CMakeFiles/check_base_fast.dir/flags.make
CMakeFiles/check_base_fast.dir/test_accuracy.c.o: ../test_accuracy.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.2zHWrVoTcO/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/check_base_fast.dir/test_accuracy.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/check_base_fast.dir/test_accuracy.c.o   -c /tmp/tmp.2zHWrVoTcO/test_accuracy.c

CMakeFiles/check_base_fast.dir/test_accuracy.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/check_base_fast.dir/test_accuracy.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /tmp/tmp.2zHWrVoTcO/test_accuracy.c > CMakeFiles/check_base_fast.dir/test_accuracy.c.i

CMakeFiles/check_base_fast.dir/test_accuracy.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/check_base_fast.dir/test_accuracy.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /tmp/tmp.2zHWrVoTcO/test_accuracy.c -o CMakeFiles/check_base_fast.dir/test_accuracy.c.s

# Object files for target check_base_fast
check_base_fast_OBJECTS = \
"CMakeFiles/check_base_fast.dir/test_accuracy.c.o"

# External object files for target check_base_fast
check_base_fast_EXTERNAL_OBJECTS =

check_base_fast: CMakeFiles/check_base_fast.dir/test_accuracy.c.o
check_base_fast: CMakeFiles/check_base_fast.dir/build.make
check_base_fast: CMakeFiles/check_base_fast.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/tmp/tmp.2zHWrVoTcO/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable check_base_fast"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/check_base_fast.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/check_base_fast.dir/build: check_base_fast

.PHONY : CMakeFiles/check_base_fast.dir/build

CMakeFiles/check_base_fast.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/check_base_fast.dir/cmake_clean.cmake
.PHONY : CMakeFiles/check_base_fast.dir/clean

CMakeFiles/check_base_fast.dir/depend:
	cd /tmp/tmp.2zHWrVoTcO/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /tmp/tmp.2zHWrVoTcO /tmp/tmp.2zHWrVoTcO /tmp/tmp.2zHWrVoTcO/cmake-build-debug /tmp/tmp.2zHWrVoTcO/cmake-build-debug /tmp/tmp.2zHWrVoTcO/cmake-build-debug/CMakeFiles/check_base_fast.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/check_base_fast.dir/depend

