# stripack-sys

Raw Rust FFI bindings to STRIPACK, a Fortran library for Delaunay triangulation and Voronoi diagrams on the unit sphere.

## Overview

This crate provides low-level bindings to STRIPACK.

## Building

Requires a Fortran compiler:

```bash
# Arch
pacman -S gcc-fortran

# Ubuntu/Debian
apt install gfortran

# macOS
brew install gcc
```

## License

The Rust bindings are licensed under MIT OR Apache-2.0 at your option.

The bundled Fortran code is distributed under the GNU LGPL. See [fortran/LICENSE](fortran/LICENSE).

## Attribution

**STRIPACK** — Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere

- Original author: Robert J. Renka, University of North Texas
- Reference: R. J. Renka, "Algorithm 772: STRIPACK: Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere", *ACM Transactions on Mathematical Software*, Vol. 23, No. 3, September 1997, pp. 416-434.
- DOI: [10.1145/275323.275329](https://doi.org/10.1145/275323.275329)

The Fortran 90 version used in this crate was prepared by [John Burkardt](https://people.sc.fsu.edu/~jburkardt/f_src/stripack/stripack.html) and is distributed under the GNU LGPL.
