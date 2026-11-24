#!/usr/bin/env bash
set -euo pipefail

# Explicit source build; Julia finds the default installation without an override.
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-$root/deps/build}"
prefix="${PREFIX:-$root/deps/usr}"
cmake_args=()
if [[ "${OSTYPE:-}" == darwin* && -n "${CMAKE_OSX_DEPLOYMENT_TARGET:-}" ]]; then
  cmake_args+=("-DCMAKE_OSX_DEPLOYMENT_TARGET=$CMAKE_OSX_DEPLOYMENT_TARGET")
fi

cmake -S "$root/cxxBindings" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  ${cmake_args[@]+"${cmake_args[@]}"}
cmake --build "$build_dir" --config Release --parallel
cmake --install "$build_dir" --config Release --component RuntimeLibraries
printf 'Installed native library under %s\n' "$prefix"
