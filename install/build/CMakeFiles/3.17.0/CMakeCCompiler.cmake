set(CMAKE_C_COMPILER "/software/gcc/11.2.0/b1/bin/gcc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "11.2.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "/software/gcc/11.2.0/b1/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/software/gcc/11.2.0/b1/bin/gcc-ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/software/gcc/11.2.0/b1/lib/gcc/x86_64-pc-linux-gnu/11.2.0/include;/software/openmpi/4.0.4/b1/include;/software/openblas/0.2.20/b1/include;/software/json11/1.0.0/include/x86_64-redhat-linux;/software/hdf5/1.8.17/b1/include;/software/armadillo/10.7.1/include;/software/libgd/2.2.5/include;/software/slurm/current/include;/gpfs/fs1/sfw2/gcc/11.2.0/b1/lib/gcc/x86_64-pc-linux-gnu/11.2.0/include-fixed;/usr/local/include;/gpfs/fs1/sfw2/gcc/11.2.0/b1/include;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/gpfs/fs1/sfw2/gcc/11.2.0/b1/lib/gcc/x86_64-pc-linux-gnu/11.2.0;/gpfs/fs1/sfw2/gcc/11.2.0/b1/lib/gcc;/software/gcc/11.2.0/b1/lib64;/software/armadillo/10.7.1/lib64;/software/slurm/current/lib64;/gpfs/fs1/sfw2/gcc/11.2.0/b1/lib64;/lib64;/usr/lib64;/software/openmpi/4.0.4/b1/lib;/software/openblas/0.2.20/b1/lib;/software/lapack/3.7.1/b1;/software/json11/1.0.0/lib/x86_64-redhat-linux;/software/hdf5/1.8.17/b1/lib;/software/libgd/2.2.5/lib;/software/slurm/current/lib;/gpfs/fs1/sfw2/gcc/11.2.0/b1/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
