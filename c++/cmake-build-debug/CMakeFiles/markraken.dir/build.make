# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /tmp/tmp.dY6zw5puj1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /tmp/tmp.dY6zw5puj1/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/markraken.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/markraken.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/markraken.dir/flags.make

CMakeFiles/markraken.dir/main.cpp.o: CMakeFiles/markraken.dir/flags.make
CMakeFiles/markraken.dir/main.cpp.o: ../main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /tmp/tmp.dY6zw5puj1/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/markraken.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/markraken.dir/main.cpp.o -c /tmp/tmp.dY6zw5puj1/main.cpp

CMakeFiles/markraken.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/markraken.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /tmp/tmp.dY6zw5puj1/main.cpp > CMakeFiles/markraken.dir/main.cpp.i

CMakeFiles/markraken.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/markraken.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /tmp/tmp.dY6zw5puj1/main.cpp -o CMakeFiles/markraken.dir/main.cpp.s

CMakeFiles/markraken.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/markraken.dir/main.cpp.o.requires

CMakeFiles/markraken.dir/main.cpp.o.provides: CMakeFiles/markraken.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/markraken.dir/build.make CMakeFiles/markraken.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/markraken.dir/main.cpp.o.provides

CMakeFiles/markraken.dir/main.cpp.o.provides.build: CMakeFiles/markraken.dir/main.cpp.o

CMakeFiles/markraken.dir/HPC.cpp.o: CMakeFiles/markraken.dir/flags.make
CMakeFiles/markraken.dir/HPC.cpp.o: ../HPC.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /tmp/tmp.dY6zw5puj1/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/markraken.dir/HPC.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/markraken.dir/HPC.cpp.o -c /tmp/tmp.dY6zw5puj1/HPC.cpp

CMakeFiles/markraken.dir/HPC.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/markraken.dir/HPC.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /tmp/tmp.dY6zw5puj1/HPC.cpp > CMakeFiles/markraken.dir/HPC.cpp.i

CMakeFiles/markraken.dir/HPC.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/markraken.dir/HPC.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /tmp/tmp.dY6zw5puj1/HPC.cpp -o CMakeFiles/markraken.dir/HPC.cpp.s

CMakeFiles/markraken.dir/HPC.cpp.o.requires:
.PHONY : CMakeFiles/markraken.dir/HPC.cpp.o.requires

CMakeFiles/markraken.dir/HPC.cpp.o.provides: CMakeFiles/markraken.dir/HPC.cpp.o.requires
	$(MAKE) -f CMakeFiles/markraken.dir/build.make CMakeFiles/markraken.dir/HPC.cpp.o.provides.build
.PHONY : CMakeFiles/markraken.dir/HPC.cpp.o.provides

CMakeFiles/markraken.dir/HPC.cpp.o.provides.build: CMakeFiles/markraken.dir/HPC.cpp.o

CMakeFiles/markraken.dir/include/FastaTools.cpp.o: CMakeFiles/markraken.dir/flags.make
CMakeFiles/markraken.dir/include/FastaTools.cpp.o: ../include/FastaTools.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /tmp/tmp.dY6zw5puj1/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/markraken.dir/include/FastaTools.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/markraken.dir/include/FastaTools.cpp.o -c /tmp/tmp.dY6zw5puj1/include/FastaTools.cpp

CMakeFiles/markraken.dir/include/FastaTools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/markraken.dir/include/FastaTools.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /tmp/tmp.dY6zw5puj1/include/FastaTools.cpp > CMakeFiles/markraken.dir/include/FastaTools.cpp.i

CMakeFiles/markraken.dir/include/FastaTools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/markraken.dir/include/FastaTools.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /tmp/tmp.dY6zw5puj1/include/FastaTools.cpp -o CMakeFiles/markraken.dir/include/FastaTools.cpp.s

CMakeFiles/markraken.dir/include/FastaTools.cpp.o.requires:
.PHONY : CMakeFiles/markraken.dir/include/FastaTools.cpp.o.requires

CMakeFiles/markraken.dir/include/FastaTools.cpp.o.provides: CMakeFiles/markraken.dir/include/FastaTools.cpp.o.requires
	$(MAKE) -f CMakeFiles/markraken.dir/build.make CMakeFiles/markraken.dir/include/FastaTools.cpp.o.provides.build
.PHONY : CMakeFiles/markraken.dir/include/FastaTools.cpp.o.provides

CMakeFiles/markraken.dir/include/FastaTools.cpp.o.provides.build: CMakeFiles/markraken.dir/include/FastaTools.cpp.o

CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o: CMakeFiles/markraken.dir/flags.make
CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o: ../include/NcbiTaxonomy.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /tmp/tmp.dY6zw5puj1/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o -c /tmp/tmp.dY6zw5puj1/include/NcbiTaxonomy.cpp

CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /tmp/tmp.dY6zw5puj1/include/NcbiTaxonomy.cpp > CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.i

CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /tmp/tmp.dY6zw5puj1/include/NcbiTaxonomy.cpp -o CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.s

CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.requires:
.PHONY : CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.requires

CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.provides: CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.requires
	$(MAKE) -f CMakeFiles/markraken.dir/build.make CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.provides.build
.PHONY : CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.provides

CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.provides.build: CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o

CMakeFiles/markraken.dir/markerizer.cpp.o: CMakeFiles/markraken.dir/flags.make
CMakeFiles/markraken.dir/markerizer.cpp.o: ../markerizer.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /tmp/tmp.dY6zw5puj1/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/markraken.dir/markerizer.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/markraken.dir/markerizer.cpp.o -c /tmp/tmp.dY6zw5puj1/markerizer.cpp

CMakeFiles/markraken.dir/markerizer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/markraken.dir/markerizer.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /tmp/tmp.dY6zw5puj1/markerizer.cpp > CMakeFiles/markraken.dir/markerizer.cpp.i

CMakeFiles/markraken.dir/markerizer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/markraken.dir/markerizer.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /tmp/tmp.dY6zw5puj1/markerizer.cpp -o CMakeFiles/markraken.dir/markerizer.cpp.s

CMakeFiles/markraken.dir/markerizer.cpp.o.requires:
.PHONY : CMakeFiles/markraken.dir/markerizer.cpp.o.requires

CMakeFiles/markraken.dir/markerizer.cpp.o.provides: CMakeFiles/markraken.dir/markerizer.cpp.o.requires
	$(MAKE) -f CMakeFiles/markraken.dir/build.make CMakeFiles/markraken.dir/markerizer.cpp.o.provides.build
.PHONY : CMakeFiles/markraken.dir/markerizer.cpp.o.provides

CMakeFiles/markraken.dir/markerizer.cpp.o.provides.build: CMakeFiles/markraken.dir/markerizer.cpp.o

CMakeFiles/markraken.dir/hasher.cpp.o: CMakeFiles/markraken.dir/flags.make
CMakeFiles/markraken.dir/hasher.cpp.o: ../hasher.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /tmp/tmp.dY6zw5puj1/cmake-build-debug/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/markraken.dir/hasher.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/markraken.dir/hasher.cpp.o -c /tmp/tmp.dY6zw5puj1/hasher.cpp

CMakeFiles/markraken.dir/hasher.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/markraken.dir/hasher.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /tmp/tmp.dY6zw5puj1/hasher.cpp > CMakeFiles/markraken.dir/hasher.cpp.i

CMakeFiles/markraken.dir/hasher.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/markraken.dir/hasher.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /tmp/tmp.dY6zw5puj1/hasher.cpp -o CMakeFiles/markraken.dir/hasher.cpp.s

CMakeFiles/markraken.dir/hasher.cpp.o.requires:
.PHONY : CMakeFiles/markraken.dir/hasher.cpp.o.requires

CMakeFiles/markraken.dir/hasher.cpp.o.provides: CMakeFiles/markraken.dir/hasher.cpp.o.requires
	$(MAKE) -f CMakeFiles/markraken.dir/build.make CMakeFiles/markraken.dir/hasher.cpp.o.provides.build
.PHONY : CMakeFiles/markraken.dir/hasher.cpp.o.provides

CMakeFiles/markraken.dir/hasher.cpp.o.provides.build: CMakeFiles/markraken.dir/hasher.cpp.o

# Object files for target markraken
markraken_OBJECTS = \
"CMakeFiles/markraken.dir/main.cpp.o" \
"CMakeFiles/markraken.dir/HPC.cpp.o" \
"CMakeFiles/markraken.dir/include/FastaTools.cpp.o" \
"CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o" \
"CMakeFiles/markraken.dir/markerizer.cpp.o" \
"CMakeFiles/markraken.dir/hasher.cpp.o"

# External object files for target markraken
markraken_EXTERNAL_OBJECTS =

markraken: CMakeFiles/markraken.dir/main.cpp.o
markraken: CMakeFiles/markraken.dir/HPC.cpp.o
markraken: CMakeFiles/markraken.dir/include/FastaTools.cpp.o
markraken: CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o
markraken: CMakeFiles/markraken.dir/markerizer.cpp.o
markraken: CMakeFiles/markraken.dir/hasher.cpp.o
markraken: CMakeFiles/markraken.dir/build.make
markraken: CMakeFiles/markraken.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable markraken"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/markraken.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/markraken.dir/build: markraken
.PHONY : CMakeFiles/markraken.dir/build

CMakeFiles/markraken.dir/requires: CMakeFiles/markraken.dir/main.cpp.o.requires
CMakeFiles/markraken.dir/requires: CMakeFiles/markraken.dir/HPC.cpp.o.requires
CMakeFiles/markraken.dir/requires: CMakeFiles/markraken.dir/include/FastaTools.cpp.o.requires
CMakeFiles/markraken.dir/requires: CMakeFiles/markraken.dir/include/NcbiTaxonomy.cpp.o.requires
CMakeFiles/markraken.dir/requires: CMakeFiles/markraken.dir/markerizer.cpp.o.requires
CMakeFiles/markraken.dir/requires: CMakeFiles/markraken.dir/hasher.cpp.o.requires
.PHONY : CMakeFiles/markraken.dir/requires

CMakeFiles/markraken.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/markraken.dir/cmake_clean.cmake
.PHONY : CMakeFiles/markraken.dir/clean

CMakeFiles/markraken.dir/depend:
	cd /tmp/tmp.dY6zw5puj1/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /tmp/tmp.dY6zw5puj1 /tmp/tmp.dY6zw5puj1 /tmp/tmp.dY6zw5puj1/cmake-build-debug /tmp/tmp.dY6zw5puj1/cmake-build-debug /tmp/tmp.dY6zw5puj1/cmake-build-debug/CMakeFiles/markraken.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/markraken.dir/depend

