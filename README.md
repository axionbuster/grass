# Gravity (*n*-body) simulation

![Screenshot](demo/galaxies.png)

My goal is to have a performant real-time n-body simulator that works on all major platforms and that which is
reasonably, statistically accurate.

The demo shown is in the `/demo` directory.
Go to the `demo` folder documentation for build instructions.

You may need a C\+\+20-capable compiler (but C\+\+20 modules are not in use).

## Compilers

Again, C\+\+20 support is necessary (sans modules).

* MSVC
* GCC
* LLVM Clang
* (Maybe Intel ICC)

## CPU

* x86_64
* Apple ARM
* (Maybe others. I only use portable code.)

## IDE

* Visual Studio 2022 with official CMake extension
* CLion
* None

## Dependencies

* The CMake script is instructed to download and cache Raylib automatically.
* However, Raylib 5's build dependencies must be satisfied.
  On Windows and macOS, I didn't need to install anything (I believe).
  On Linux, on the other hand, I had to install a few things. See Raylib's docs for more.
