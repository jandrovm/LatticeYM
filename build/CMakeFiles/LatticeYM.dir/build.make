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

# Utility rule file for LatticeYM.

# Include any custom commands dependencies for this target.
include CMakeFiles/LatticeYM.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/LatticeYM.dir/progress.make

CMakeFiles/LatticeYM: latticeYM

LatticeYM: CMakeFiles/LatticeYM
LatticeYM: CMakeFiles/LatticeYM.dir/build.make
.PHONY : LatticeYM

# Rule to build all files generated by this target.
CMakeFiles/LatticeYM.dir/build: LatticeYM
.PHONY : CMakeFiles/LatticeYM.dir/build

CMakeFiles/LatticeYM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LatticeYM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LatticeYM.dir/clean

CMakeFiles/LatticeYM.dir/depend:
	cd /home/ale/Documentos/lattice/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ale/Documentos/lattice /home/ale/Documentos/lattice /home/ale/Documentos/lattice/build /home/ale/Documentos/lattice/build /home/ale/Documentos/lattice/build/CMakeFiles/LatticeYM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LatticeYM.dir/depend

