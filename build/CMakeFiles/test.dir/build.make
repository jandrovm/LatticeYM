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
CMAKE_SOURCE_DIR = /home/ale/Documentos/lattice

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ale/Documentos/lattice/build

# Include any dependencies generated for this target.
include CMakeFiles/test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/test.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test.dir/flags.make

CMakeFiles/test.dir/test.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/test.cc.o: ../test.cc
CMakeFiles/test.dir/test.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test.dir/test.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/test.cc.o -MF CMakeFiles/test.dir/test.cc.o.d -o CMakeFiles/test.dir/test.cc.o -c /home/ale/Documentos/lattice/test.cc

CMakeFiles/test.dir/test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/test.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/test.cc > CMakeFiles/test.dir/test.cc.i

CMakeFiles/test.dir/test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/test.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/test.cc -o CMakeFiles/test.dir/test.cc.s

CMakeFiles/test.dir/src/Lattice.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/Lattice.cc.o: ../src/Lattice.cc
CMakeFiles/test.dir/src/Lattice.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/test.dir/src/Lattice.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/Lattice.cc.o -MF CMakeFiles/test.dir/src/Lattice.cc.o.d -o CMakeFiles/test.dir/src/Lattice.cc.o -c /home/ale/Documentos/lattice/src/Lattice.cc

CMakeFiles/test.dir/src/Lattice.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/src/Lattice.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/src/Lattice.cc > CMakeFiles/test.dir/src/Lattice.cc.i

CMakeFiles/test.dir/src/Lattice.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/Lattice.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/src/Lattice.cc -o CMakeFiles/test.dir/src/Lattice.cc.s

CMakeFiles/test.dir/src/Link.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/Link.cc.o: ../src/Link.cc
CMakeFiles/test.dir/src/Link.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/test.dir/src/Link.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/Link.cc.o -MF CMakeFiles/test.dir/src/Link.cc.o.d -o CMakeFiles/test.dir/src/Link.cc.o -c /home/ale/Documentos/lattice/src/Link.cc

CMakeFiles/test.dir/src/Link.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/src/Link.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/src/Link.cc > CMakeFiles/test.dir/src/Link.cc.i

CMakeFiles/test.dir/src/Link.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/Link.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/src/Link.cc -o CMakeFiles/test.dir/src/Link.cc.s

CMakeFiles/test.dir/src/SU3.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/SU3.cc.o: ../src/SU3.cc
CMakeFiles/test.dir/src/SU3.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/test.dir/src/SU3.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/SU3.cc.o -MF CMakeFiles/test.dir/src/SU3.cc.o.d -o CMakeFiles/test.dir/src/SU3.cc.o -c /home/ale/Documentos/lattice/src/SU3.cc

CMakeFiles/test.dir/src/SU3.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/src/SU3.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/src/SU3.cc > CMakeFiles/test.dir/src/SU3.cc.i

CMakeFiles/test.dir/src/SU3.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/SU3.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/src/SU3.cc -o CMakeFiles/test.dir/src/SU3.cc.s

CMakeFiles/test.dir/src/Simulation.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/Simulation.cc.o: ../src/Simulation.cc
CMakeFiles/test.dir/src/Simulation.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/test.dir/src/Simulation.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/Simulation.cc.o -MF CMakeFiles/test.dir/src/Simulation.cc.o.d -o CMakeFiles/test.dir/src/Simulation.cc.o -c /home/ale/Documentos/lattice/src/Simulation.cc

CMakeFiles/test.dir/src/Simulation.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/src/Simulation.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/src/Simulation.cc > CMakeFiles/test.dir/src/Simulation.cc.i

CMakeFiles/test.dir/src/Simulation.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/Simulation.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/src/Simulation.cc -o CMakeFiles/test.dir/src/Simulation.cc.s

CMakeFiles/test.dir/src/SimulationMT.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/SimulationMT.cc.o: ../src/SimulationMT.cc
CMakeFiles/test.dir/src/SimulationMT.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/test.dir/src/SimulationMT.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/SimulationMT.cc.o -MF CMakeFiles/test.dir/src/SimulationMT.cc.o.d -o CMakeFiles/test.dir/src/SimulationMT.cc.o -c /home/ale/Documentos/lattice/src/SimulationMT.cc

CMakeFiles/test.dir/src/SimulationMT.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/src/SimulationMT.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/src/SimulationMT.cc > CMakeFiles/test.dir/src/SimulationMT.cc.i

CMakeFiles/test.dir/src/SimulationMT.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/SimulationMT.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/src/SimulationMT.cc -o CMakeFiles/test.dir/src/SimulationMT.cc.s

CMakeFiles/test.dir/src/ThreadPool.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/ThreadPool.cc.o: ../src/ThreadPool.cc
CMakeFiles/test.dir/src/ThreadPool.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/test.dir/src/ThreadPool.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/ThreadPool.cc.o -MF CMakeFiles/test.dir/src/ThreadPool.cc.o.d -o CMakeFiles/test.dir/src/ThreadPool.cc.o -c /home/ale/Documentos/lattice/src/ThreadPool.cc

CMakeFiles/test.dir/src/ThreadPool.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/src/ThreadPool.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/src/ThreadPool.cc > CMakeFiles/test.dir/src/ThreadPool.cc.i

CMakeFiles/test.dir/src/ThreadPool.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/ThreadPool.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/src/ThreadPool.cc -o CMakeFiles/test.dir/src/ThreadPool.cc.s

CMakeFiles/test.dir/src/utils.cc.o: CMakeFiles/test.dir/flags.make
CMakeFiles/test.dir/src/utils.cc.o: ../src/utils.cc
CMakeFiles/test.dir/src/utils.cc.o: CMakeFiles/test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/test.dir/src/utils.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test.dir/src/utils.cc.o -MF CMakeFiles/test.dir/src/utils.cc.o.d -o CMakeFiles/test.dir/src/utils.cc.o -c /home/ale/Documentos/lattice/src/utils.cc

CMakeFiles/test.dir/src/utils.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test.dir/src/utils.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ale/Documentos/lattice/src/utils.cc > CMakeFiles/test.dir/src/utils.cc.i

CMakeFiles/test.dir/src/utils.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test.dir/src/utils.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ale/Documentos/lattice/src/utils.cc -o CMakeFiles/test.dir/src/utils.cc.s

# Object files for target test
test_OBJECTS = \
"CMakeFiles/test.dir/test.cc.o" \
"CMakeFiles/test.dir/src/Lattice.cc.o" \
"CMakeFiles/test.dir/src/Link.cc.o" \
"CMakeFiles/test.dir/src/SU3.cc.o" \
"CMakeFiles/test.dir/src/Simulation.cc.o" \
"CMakeFiles/test.dir/src/SimulationMT.cc.o" \
"CMakeFiles/test.dir/src/ThreadPool.cc.o" \
"CMakeFiles/test.dir/src/utils.cc.o"

# External object files for target test
test_EXTERNAL_OBJECTS =

test: CMakeFiles/test.dir/test.cc.o
test: CMakeFiles/test.dir/src/Lattice.cc.o
test: CMakeFiles/test.dir/src/Link.cc.o
test: CMakeFiles/test.dir/src/SU3.cc.o
test: CMakeFiles/test.dir/src/Simulation.cc.o
test: CMakeFiles/test.dir/src/SimulationMT.cc.o
test: CMakeFiles/test.dir/src/ThreadPool.cc.o
test: CMakeFiles/test.dir/src/utils.cc.o
test: CMakeFiles/test.dir/build.make
test: CMakeFiles/test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ale/Documentos/lattice/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test.dir/build: test
.PHONY : CMakeFiles/test.dir/build

CMakeFiles/test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test.dir/clean

CMakeFiles/test.dir/depend:
	cd /home/ale/Documentos/lattice/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ale/Documentos/lattice /home/ale/Documentos/lattice /home/ale/Documentos/lattice/build /home/ale/Documentos/lattice/build /home/ale/Documentos/lattice/build/CMakeFiles/test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test.dir/depend
