FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive

# Both windowing backends are installed: GLFW compiles X11 and Wayland support
# together and selects between them at runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        pkg-config \
        gdb \
        git \
        libblas-dev \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        mesa-common-dev \
        libx11-dev \
        libxcursor-dev \
        libxi-dev \
        libxinerama-dev \
        libxmu-dev \
        libxrandr-dev \
        libdecor-0-dev \
        libwayland-dev \
        libxkbcommon-dev \
        wayland-protocols \
    && rm -rf /var/lib/apt/lists/*
