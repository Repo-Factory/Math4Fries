# Math4Fries
A collection functions to accelerate spline fitting, regressions, digital filtering, and more.

## Building & Linking ##
There are a few convenience scripts that should help make the (cancerous) CMake build process easy
After cloning the project use this command to build the shared library for Matrix Utilities onto your system

    cd Matrix_Utilities_4_C && ../scripts/build_project -i && cd ..

Now you can use this library from anywhere by using find_library(Matrix_Utilities) in your CMake
Instead of writing this yourself, you can take advantage of the provided new_package script which will
populate a simple project onto your computer with an AutoGen'd CMake file that is already linked to this library
Use a command like this to do so 

    ./scripts/new_package my_package

You will of course have to modify this with your actual desired name for the package and make sure you add on to the relative
path to correctly call the script from the desired directory on your machine.

