# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = "C:/Program Files/CMake/bin/cmake.exe"

# The command to remove a file.
RM = "C:/Program Files/CMake/bin/cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:/Users/remof/source/repos/IMT_FAST/IMT_FAST

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:/Users/remof/source/repos/IMT_FAST/IMT_FAST

# Include any dependencies generated for this target.
include CMakeFiles/MyExecutable.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MyExecutable.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MyExecutable.dir/flags.make

CMakeFiles/MyExecutable.dir/binned_conv.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/binned_conv.obj: binned_conv.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/MyExecutable.dir/binned_conv.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/binned_conv.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/binned_conv.c

CMakeFiles/MyExecutable.dir/binned_conv.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/binned_conv.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/binned_conv.c > CMakeFiles/MyExecutable.dir/binned_conv.i

CMakeFiles/MyExecutable.dir/binned_conv.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/binned_conv.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/binned_conv.c -o CMakeFiles/MyExecutable.dir/binned_conv.s

CMakeFiles/MyExecutable.dir/conv.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/conv.obj: conv.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/MyExecutable.dir/conv.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/conv.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/conv.c

CMakeFiles/MyExecutable.dir/conv.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/conv.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/conv.c > CMakeFiles/MyExecutable.dir/conv.i

CMakeFiles/MyExecutable.dir/conv.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/conv.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/conv.c -o CMakeFiles/MyExecutable.dir/conv.s

CMakeFiles/MyExecutable.dir/emgpdf.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/emgpdf.obj: emgpdf.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/MyExecutable.dir/emgpdf.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/emgpdf.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/emgpdf.c

CMakeFiles/MyExecutable.dir/emgpdf.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/emgpdf.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/emgpdf.c > CMakeFiles/MyExecutable.dir/emgpdf.i

CMakeFiles/MyExecutable.dir/emgpdf.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/emgpdf.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/emgpdf.c -o CMakeFiles/MyExecutable.dir/emgpdf.s

CMakeFiles/MyExecutable.dir/gp_max.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/gp_max.obj: gp_max.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/MyExecutable.dir/gp_max.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/gp_max.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/gp_max.c

CMakeFiles/MyExecutable.dir/gp_max.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/gp_max.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/gp_max.c > CMakeFiles/MyExecutable.dir/gp_max.i

CMakeFiles/MyExecutable.dir/gp_max.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/gp_max.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/gp_max.c -o CMakeFiles/MyExecutable.dir/gp_max.s

CMakeFiles/MyExecutable.dir/imt_analysis.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/imt_analysis.obj: imt_analysis.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object CMakeFiles/MyExecutable.dir/imt_analysis.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/imt_analysis.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/imt_analysis.c

CMakeFiles/MyExecutable.dir/imt_analysis.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/imt_analysis.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/imt_analysis.c > CMakeFiles/MyExecutable.dir/imt_analysis.i

CMakeFiles/MyExecutable.dir/imt_analysis.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/imt_analysis.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/imt_analysis.c -o CMakeFiles/MyExecutable.dir/imt_analysis.s

CMakeFiles/MyExecutable.dir/loglikelihood.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/loglikelihood.obj: loglikelihood.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/MyExecutable.dir/loglikelihood.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/loglikelihood.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/loglikelihood.c

CMakeFiles/MyExecutable.dir/loglikelihood.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/loglikelihood.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/loglikelihood.c > CMakeFiles/MyExecutable.dir/loglikelihood.i

CMakeFiles/MyExecutable.dir/loglikelihood.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/loglikelihood.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/loglikelihood.c -o CMakeFiles/MyExecutable.dir/loglikelihood.s

CMakeFiles/MyExecutable.dir/main.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/main.obj: main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/MyExecutable.dir/main.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/main.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/main.c

CMakeFiles/MyExecutable.dir/main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/main.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/main.c > CMakeFiles/MyExecutable.dir/main.i

CMakeFiles/MyExecutable.dir/main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/main.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/main.c -o CMakeFiles/MyExecutable.dir/main.s

CMakeFiles/MyExecutable.dir/mexWrap.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/mexWrap.obj: mexWrap.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building C object CMakeFiles/MyExecutable.dir/mexWrap.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/mexWrap.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/mexWrap.c

CMakeFiles/MyExecutable.dir/mexWrap.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/mexWrap.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/mexWrap.c > CMakeFiles/MyExecutable.dir/mexWrap.i

CMakeFiles/MyExecutable.dir/mexWrap.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/mexWrap.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/mexWrap.c -o CMakeFiles/MyExecutable.dir/mexWrap.s

CMakeFiles/MyExecutable.dir/nn_conv.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/nn_conv.obj: nn_conv.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building C object CMakeFiles/MyExecutable.dir/nn_conv.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/nn_conv.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/nn_conv.c

CMakeFiles/MyExecutable.dir/nn_conv.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/nn_conv.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/nn_conv.c > CMakeFiles/MyExecutable.dir/nn_conv.i

CMakeFiles/MyExecutable.dir/nn_conv.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/nn_conv.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/nn_conv.c -o CMakeFiles/MyExecutable.dir/nn_conv.s

CMakeFiles/MyExecutable.dir/onestage.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/onestage.obj: onestage.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building C object CMakeFiles/MyExecutable.dir/onestage.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/onestage.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestage.c

CMakeFiles/MyExecutable.dir/onestage.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/onestage.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestage.c > CMakeFiles/MyExecutable.dir/onestage.i

CMakeFiles/MyExecutable.dir/onestage.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/onestage.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestage.c -o CMakeFiles/MyExecutable.dir/onestage.s

CMakeFiles/MyExecutable.dir/onestagepdf_lag.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/onestagepdf_lag.obj: onestagepdf_lag.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building C object CMakeFiles/MyExecutable.dir/onestagepdf_lag.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/onestagepdf_lag.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestagepdf_lag.c

CMakeFiles/MyExecutable.dir/onestagepdf_lag.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/onestagepdf_lag.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestagepdf_lag.c > CMakeFiles/MyExecutable.dir/onestagepdf_lag.i

CMakeFiles/MyExecutable.dir/onestagepdf_lag.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/onestagepdf_lag.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestagepdf_lag.c -o CMakeFiles/MyExecutable.dir/onestagepdf_lag.s

CMakeFiles/MyExecutable.dir/onestagepdf_prime.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/onestagepdf_prime.obj: onestagepdf_prime.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building C object CMakeFiles/MyExecutable.dir/onestagepdf_prime.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/onestagepdf_prime.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestagepdf_prime.c

CMakeFiles/MyExecutable.dir/onestagepdf_prime.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/onestagepdf_prime.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestagepdf_prime.c > CMakeFiles/MyExecutable.dir/onestagepdf_prime.i

CMakeFiles/MyExecutable.dir/onestagepdf_prime.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/onestagepdf_prime.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/onestagepdf_prime.c -o CMakeFiles/MyExecutable.dir/onestagepdf_prime.s

CMakeFiles/MyExecutable.dir/tailmass.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/tailmass.obj: tailmass.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building C object CMakeFiles/MyExecutable.dir/tailmass.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/tailmass.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/tailmass.c

CMakeFiles/MyExecutable.dir/tailmass.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/tailmass.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/tailmass.c > CMakeFiles/MyExecutable.dir/tailmass.i

CMakeFiles/MyExecutable.dir/tailmass.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/tailmass.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/tailmass.c -o CMakeFiles/MyExecutable.dir/tailmass.s

CMakeFiles/MyExecutable.dir/threestage.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/threestage.obj: threestage.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building C object CMakeFiles/MyExecutable.dir/threestage.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/threestage.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/threestage.c

CMakeFiles/MyExecutable.dir/threestage.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/threestage.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/threestage.c > CMakeFiles/MyExecutable.dir/threestage.i

CMakeFiles/MyExecutable.dir/threestage.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/threestage.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/threestage.c -o CMakeFiles/MyExecutable.dir/threestage.s

CMakeFiles/MyExecutable.dir/twostage.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/twostage.obj: twostage.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building C object CMakeFiles/MyExecutable.dir/twostage.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/twostage.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/twostage.c

CMakeFiles/MyExecutable.dir/twostage.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/twostage.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/twostage.c > CMakeFiles/MyExecutable.dir/twostage.i

CMakeFiles/MyExecutable.dir/twostage.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/twostage.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/twostage.c -o CMakeFiles/MyExecutable.dir/twostage.s

CMakeFiles/MyExecutable.dir/utility.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/utility.obj: utility.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building C object CMakeFiles/MyExecutable.dir/utility.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/utility.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/utility.c

CMakeFiles/MyExecutable.dir/utility.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/utility.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/utility.c > CMakeFiles/MyExecutable.dir/utility.i

CMakeFiles/MyExecutable.dir/utility.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/utility.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/utility.c -o CMakeFiles/MyExecutable.dir/utility.s

CMakeFiles/MyExecutable.dir/window_conv.obj: CMakeFiles/MyExecutable.dir/flags.make
CMakeFiles/MyExecutable.dir/window_conv.obj: window_conv.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Building C object CMakeFiles/MyExecutable.dir/window_conv.obj"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/MyExecutable.dir/window_conv.obj   -c C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/window_conv.c

CMakeFiles/MyExecutable.dir/window_conv.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/MyExecutable.dir/window_conv.i"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/window_conv.c > CMakeFiles/MyExecutable.dir/window_conv.i

CMakeFiles/MyExecutable.dir/window_conv.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/MyExecutable.dir/window_conv.s"
	C:/MinGW/bin/gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/window_conv.c -o CMakeFiles/MyExecutable.dir/window_conv.s

# Object files for target MyExecutable
MyExecutable_OBJECTS = \
"CMakeFiles/MyExecutable.dir/binned_conv.obj" \
"CMakeFiles/MyExecutable.dir/conv.obj" \
"CMakeFiles/MyExecutable.dir/emgpdf.obj" \
"CMakeFiles/MyExecutable.dir/gp_max.obj" \
"CMakeFiles/MyExecutable.dir/imt_analysis.obj" \
"CMakeFiles/MyExecutable.dir/loglikelihood.obj" \
"CMakeFiles/MyExecutable.dir/main.obj" \
"CMakeFiles/MyExecutable.dir/mexWrap.obj" \
"CMakeFiles/MyExecutable.dir/nn_conv.obj" \
"CMakeFiles/MyExecutable.dir/onestage.obj" \
"CMakeFiles/MyExecutable.dir/onestagepdf_lag.obj" \
"CMakeFiles/MyExecutable.dir/onestagepdf_prime.obj" \
"CMakeFiles/MyExecutable.dir/tailmass.obj" \
"CMakeFiles/MyExecutable.dir/threestage.obj" \
"CMakeFiles/MyExecutable.dir/twostage.obj" \
"CMakeFiles/MyExecutable.dir/utility.obj" \
"CMakeFiles/MyExecutable.dir/window_conv.obj"

# External object files for target MyExecutable
MyExecutable_EXTERNAL_OBJECTS =

MyExecutable.exe: CMakeFiles/MyExecutable.dir/binned_conv.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/conv.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/emgpdf.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/gp_max.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/imt_analysis.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/loglikelihood.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/main.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/mexWrap.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/nn_conv.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/onestage.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/onestagepdf_lag.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/onestagepdf_prime.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/tailmass.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/threestage.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/twostage.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/utility.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/window_conv.obj
MyExecutable.exe: CMakeFiles/MyExecutable.dir/build.make
MyExecutable.exe: CMakeFiles/MyExecutable.dir/linklibs.rsp
MyExecutable.exe: CMakeFiles/MyExecutable.dir/objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles --progress-num=$(CMAKE_PROGRESS_18) "Linking C executable MyExecutable.exe"
	"C:/Program Files/CMake/bin/cmake.exe" -E remove -f CMakeFiles/MyExecutable.dir/objects.a
	C:/MinGW/bin/ar.exe cr CMakeFiles/MyExecutable.dir/objects.a @CMakeFiles/MyExecutable.dir/objects1.rsp
	C:/MinGW/bin/gcc.exe -g   -Wl,--whole-archive CMakeFiles/MyExecutable.dir/objects.a -Wl,--no-whole-archive  -o MyExecutable.exe -Wl,--out-implib,libMyExecutable.dll.a -Wl,--major-image-version,0,--minor-image-version,0 @CMakeFiles/MyExecutable.dir/linklibs.rsp

# Rule to build all files generated by this target.
CMakeFiles/MyExecutable.dir/build: MyExecutable.exe

.PHONY : CMakeFiles/MyExecutable.dir/build

CMakeFiles/MyExecutable.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MyExecutable.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MyExecutable.dir/clean

CMakeFiles/MyExecutable.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" C:/Users/remof/source/repos/IMT_FAST/IMT_FAST C:/Users/remof/source/repos/IMT_FAST/IMT_FAST C:/Users/remof/source/repos/IMT_FAST/IMT_FAST C:/Users/remof/source/repos/IMT_FAST/IMT_FAST C:/Users/remof/source/repos/IMT_FAST/IMT_FAST/CMakeFiles/MyExecutable.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MyExecutable.dir/depend

