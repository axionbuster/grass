# The main demo (including "galaxies" mode)

![Galaxies Mode](galaxies.png)
Screenshot: Galaxies Mode (there's actually 2,500 particles)

Build using CMake.

Enable "galaxies" mode by running the binary with the environment variable
`GRASS_GALAXIES` set (with any value).

For example, supposing I'm in the build directory (bash):
```bash
# 1. Build (assuming Make in use; this is the usual case on Linux or macOS)
make grass

# 2. Run with optional GRASS_GALAXIES environment variable (any value;
# value is ignored as long as the key is set.)
GRASS_GALAXIES=111111 demo/grass
```

On most setups, the build directory is either called `build` or `out`, but other
choices are possible.

## What is the "galaxies" mode?

In this mode, at program startup, a random assortment of clumps of particles are generated.