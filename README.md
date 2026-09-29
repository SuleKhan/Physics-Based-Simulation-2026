# Physically-based Simulation in Computer Graphics HS2026 - Course Exercises

## Installation

The exercises build on Linux, macOS and Windows. If you would rather not install
a toolchain locally, the repository ships a dev container (see
[Dev container](#dev-container) below).

### Dependencies

**Ubuntu / Debian**

```
sudo apt-get install build-essential cmake pkg-config git gdb ca-certificates \
    libgl1-mesa-dev libglu1-mesa-dev mesa-common-dev libblas-dev \
    libx11-dev libxcursor-dev libxi-dev libxinerama-dev libxmu-dev libxrandr-dev \
    libwayland-dev libxkbcommon-dev wayland-protocols libdecor-0-dev
```

The last line covers Wayland, which is the default session on current Ubuntu.
The X11 packages above it are still needed: the viewer is built with both
backends and picks one at runtime, so it keeps working under X11 and XWayland.

On other Linux distributions the package names differ, but the set is the same:
a C++ compiler, CMake, OpenGL headers, and the X11 and Wayland client libraries.

**macOS**

```
xcode-select --install
brew install cmake
brew install libomp     # optional: enables OpenMP parallelism
```

**Windows**

Install [Visual Studio](https://visualstudio.microsoft.com/) with the
"Desktop development with C++" workload, which includes CMake.

### The Exercise Repository

This repository will be updated each week with the new assignment, along with the
solution to the previous week's assignment. We recommend creating a private
[fork](https://docs.gitlab.com/ee/gitlab-basics/fork-project.html) of this
repository such that you can store your code using Git. For this you need to have
an active [gitlab@ETH](https://gitlab.ethz.ch/) account.

Once you have forked the repository, clone the code to your machine.

```
git clone https://gitlab.ethz.ch/'{your_username}'/pbs26.git
```

Clone submodules (polyscope and dependencies) using:

```
git submodule update --init --recursive
```

### Building

From the repository root:

```
cmake --preset release
cmake --build --preset release
```

The first command downloads libigl and GLFW, so it needs a network connection and
takes a few minutes.

Three presets are available, each building into its own directory so they do not
overwrite one another:

| Preset | Directory | Use for |
| --- | --- | --- |
| `relwithdebinfo` | `build-relwithdebinfo/` | Optimized, with debug information |
| `debug` | `build-debug/` | Unoptimized, simulations run slowly |
| `release` | `build/` | Optimized, no debug information |

Both commands work identically on all three platforms. On Windows the second one
drives MSBuild, and the solution can also be opened directly in Visual Studio.

Run an exercise:

```
./build/1_mass_spring/1_mass_spring
```

Each exercise builds as `<name>` from the sources in its `template/`
directory, which is where your implementation goes. Once the solution to an
exercise has been released it also builds as `<name>_solution`, so you can
compare your implementation against it.

Remember to re-run `cmake --preset release` after pulling updates.

### Dev container

The repository includes a [dev container](https://containers.dev/) based on
Ubuntu 26.04 with every dependency preinstalled. Open the folder in VS Code with
the Dev Containers extension and choose "Reopen in Container".

The container forwards both the X11 and Wayland sockets from the host, so the
viewer window opens on either kind of session. For hardware-accelerated
rendering, uncomment the `runArgs` line in `.devcontainer/devcontainer.json`;
it is off by default because it prevents the container from starting on hosts
with no DRI device.

### Update Your Forked Repository

To update your forked repository, check this page:
[how-do-i-update-a-github-forked-repository](https://stackoverflow.com/questions/7244321/how-do-i-update-a-github-forked-repository)

Add our repository as a remote to your own:

```
git remote add upstream https://gitlab.ethz.ch/crl-pbs/pbs26.git
```

Then, fetch updates from it:

```
git fetch upstream
```

Lastly, move to your `master` branch and merge updates into yours:

```
git checkout master
git merge upstream/master
```

Remember to run `cmake` again after updating!
