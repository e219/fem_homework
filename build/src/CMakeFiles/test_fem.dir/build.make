# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dpt/cpp-project/读书报告

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dpt/cpp-project/读书报告/build

# Include any dependencies generated for this target.
include src/CMakeFiles/test_fem.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/test_fem.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/test_fem.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/test_fem.dir/flags.make

src/CMakeFiles/test_fem.dir/MG.cpp.o: src/CMakeFiles/test_fem.dir/flags.make
src/CMakeFiles/test_fem.dir/MG.cpp.o: ../src/MG.cpp
src/CMakeFiles/test_fem.dir/MG.cpp.o: src/CMakeFiles/test_fem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dpt/cpp-project/读书报告/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/test_fem.dir/MG.cpp.o"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/test_fem.dir/MG.cpp.o -MF CMakeFiles/test_fem.dir/MG.cpp.o.d -o CMakeFiles/test_fem.dir/MG.cpp.o -c /home/dpt/cpp-project/读书报告/src/MG.cpp

src/CMakeFiles/test_fem.dir/MG.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_fem.dir/MG.cpp.i"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dpt/cpp-project/读书报告/src/MG.cpp > CMakeFiles/test_fem.dir/MG.cpp.i

src/CMakeFiles/test_fem.dir/MG.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_fem.dir/MG.cpp.s"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dpt/cpp-project/读书报告/src/MG.cpp -o CMakeFiles/test_fem.dir/MG.cpp.s

src/CMakeFiles/test_fem.dir/P1Fem.cpp.o: src/CMakeFiles/test_fem.dir/flags.make
src/CMakeFiles/test_fem.dir/P1Fem.cpp.o: ../src/P1Fem.cpp
src/CMakeFiles/test_fem.dir/P1Fem.cpp.o: src/CMakeFiles/test_fem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dpt/cpp-project/读书报告/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/test_fem.dir/P1Fem.cpp.o"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/test_fem.dir/P1Fem.cpp.o -MF CMakeFiles/test_fem.dir/P1Fem.cpp.o.d -o CMakeFiles/test_fem.dir/P1Fem.cpp.o -c /home/dpt/cpp-project/读书报告/src/P1Fem.cpp

src/CMakeFiles/test_fem.dir/P1Fem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_fem.dir/P1Fem.cpp.i"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dpt/cpp-project/读书报告/src/P1Fem.cpp > CMakeFiles/test_fem.dir/P1Fem.cpp.i

src/CMakeFiles/test_fem.dir/P1Fem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_fem.dir/P1Fem.cpp.s"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dpt/cpp-project/读书报告/src/P1Fem.cpp -o CMakeFiles/test_fem.dir/P1Fem.cpp.s

src/CMakeFiles/test_fem.dir/test_fem.cpp.o: src/CMakeFiles/test_fem.dir/flags.make
src/CMakeFiles/test_fem.dir/test_fem.cpp.o: ../src/test_fem.cpp
src/CMakeFiles/test_fem.dir/test_fem.cpp.o: src/CMakeFiles/test_fem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dpt/cpp-project/读书报告/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/test_fem.dir/test_fem.cpp.o"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/test_fem.dir/test_fem.cpp.o -MF CMakeFiles/test_fem.dir/test_fem.cpp.o.d -o CMakeFiles/test_fem.dir/test_fem.cpp.o -c /home/dpt/cpp-project/读书报告/src/test_fem.cpp

src/CMakeFiles/test_fem.dir/test_fem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_fem.dir/test_fem.cpp.i"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/dpt/cpp-project/读书报告/src/test_fem.cpp > CMakeFiles/test_fem.dir/test_fem.cpp.i

src/CMakeFiles/test_fem.dir/test_fem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_fem.dir/test_fem.cpp.s"
	cd /home/dpt/cpp-project/读书报告/build/src && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/dpt/cpp-project/读书报告/src/test_fem.cpp -o CMakeFiles/test_fem.dir/test_fem.cpp.s

# Object files for target test_fem
test_fem_OBJECTS = \
"CMakeFiles/test_fem.dir/MG.cpp.o" \
"CMakeFiles/test_fem.dir/P1Fem.cpp.o" \
"CMakeFiles/test_fem.dir/test_fem.cpp.o"

# External object files for target test_fem
test_fem_EXTERNAL_OBJECTS =

../bin/test_fem: src/CMakeFiles/test_fem.dir/MG.cpp.o
../bin/test_fem: src/CMakeFiles/test_fem.dir/P1Fem.cpp.o
../bin/test_fem: src/CMakeFiles/test_fem.dir/test_fem.cpp.o
../bin/test_fem: src/CMakeFiles/test_fem.dir/build.make
../bin/test_fem: src/CMakeFiles/test_fem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dpt/cpp-project/读书报告/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable ../../bin/test_fem"
	cd /home/dpt/cpp-project/读书报告/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_fem.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/test_fem.dir/build: ../bin/test_fem
.PHONY : src/CMakeFiles/test_fem.dir/build

src/CMakeFiles/test_fem.dir/clean:
	cd /home/dpt/cpp-project/读书报告/build/src && $(CMAKE_COMMAND) -P CMakeFiles/test_fem.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/test_fem.dir/clean

src/CMakeFiles/test_fem.dir/depend:
	cd /home/dpt/cpp-project/读书报告/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dpt/cpp-project/读书报告 /home/dpt/cpp-project/读书报告/src /home/dpt/cpp-project/读书报告/build /home/dpt/cpp-project/读书报告/build/src /home/dpt/cpp-project/读书报告/build/src/CMakeFiles/test_fem.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/test_fem.dir/depend

