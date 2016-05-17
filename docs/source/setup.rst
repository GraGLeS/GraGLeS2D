How to install the GraGLeS2D application
=============================


   
1. Check your Linux installation for cmake and make
2. Pull the project from the git repository https://github.com/GraGLeS/GraGLeS2D
3. Check your compiler. Both the Intel compiler (v14.0 or newer) and the GNU compiler can be utilized
4. Open a console and go to the src directory 
5. create a build folder: *mkdir build*
6. go to the build directory and run cmake "cmake .."
7. download the library: https://github.com/jemalloc
8. ./configure the library and copy ''libjemalloc.so'' and ''libjemalloc.so.2'' (from lib directory) to your source directory

9. Compile the program with *make*
10. Find happily the binary in *build*
11. copy one of the example paramter files ( see params folder ) into the build folder
12. run the apllication with "GraGLeS_2D parameters.xml" 

