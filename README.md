# CS591 C1 Twist Free Rod Implementation
Daniel Kehr
Paul Menexas

A project that, given a spline or helix function, produces a "zero-twist" visualization of the curve.
The libigl example project was used for library linking and cmake initialization.

-- LIBIGL EXAMPLE PROJECT INSTRUCTION GIVEN BELOW --
## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git


## Compile
The included 'ImGuiMenu.cpp' file should replace libigl's original version in ../libigl/include/igl/opengl/glfw/imgui/


The project can be compiled using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

Which should find and build the dependencies and create a `example_bin` binary.
If there is trouble locating libigl using CMAKE, try calling cmake as:

cmake -DLIBIGL_INCLUDE_DIR="your/path/for/libigl/include" -DCMAKE_PREFIX_PATH="the/same/path/for/libigl/include" ../

and run make as usual

## Run

From within the `build` directory call:

    ./example

Our project in the libigl viewer should appear.
