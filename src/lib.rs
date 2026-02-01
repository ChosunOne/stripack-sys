//! Raw Rust FFI bindings to STRIPACK.
//!
//! This crate provides low-level unsafe bindings to STRIPACK, a Fortran library for
//! constructing Delaunay triangulations and Voronoi diagrams on the surface of the unit sphere.
//!
//! # Overview
//!
//! STRIPACK (Algorithm 772) computes Delaunay triangulations and Voronoi diagrams of a set of
//! nodes on the surface of the unit sphere. The algorithms are generalizations of Robert Renka's
//! TRIPACK software for the plane.
//!
//! This `-sys` crate provides direct FFI bindings to the underlying Fortran routines.
//!
//! # Building
//!
//! This crate requires a Fortran compiler to build the bundled STRIPACK library:
//!
//! - **Arch Linux**: `pacman -S gcc-fortran`
//! - **Ubuntu/Debian**: `apt install gfortran`
//! - **macOS**: `brew install gcc`
//!
//! # Safety
//!
//! All functions in this crate are `unsafe` because they:
//! - Work with raw pointers
//! - Call into Fortran code with different calling conventions
//! - Require the caller to ensure proper memory allocation and alignment
//! - May have undefined behavior if preconditions are not met
//!
//! Users must ensure:
//! - Arrays are properly sized (e.g., `list[6*(n-2)]`, `lptr[6*(n-2)]`, `lend[n]`)
//! - Pointers are valid and properly aligned
//! - Nodes are unit vectors (x² + y² + z² = 1)
//! - Indices are 1-based (as per Fortran convention)
//!
//! # Data Structure
//!
//! The triangulation is represented using a linked list data structure:
//! - `list`: Adjacency lists
//! - `lptr`: Pointers within `list`
//! - `lend`: Pointers to ends of adjacency lists
//! - `lnew`: Pointer to first empty location
//!
//! Refer to the individual function documentation for details on array sizes and conventions.
//!
//! # Example
//!
//! ```rust
//! use stripack_sys::ffi::*;
//!
//! // Create a simple tetrahedral triangulation with 4 nodes
//! let n = 4;
//! let x = vec![1.0, 0.0, 0.0, 0.0];
//! let y = vec![0.0, 1.0, 0.0, 0.0];
//! let z = vec![0.0, 0.0, 1.0, 1.0];
//!
//! let mut list = vec![0; 6 * (n - 2) as usize];
//! let mut lptr = vec![0; 6 * (n - 2) as usize];
//! let mut lend = vec![0; n as usize];
//! let mut lnew = 0;
//! let mut near = vec![0; n as usize];
//! let mut next = vec![0; n as usize];
//! let mut dist = vec![0.0; n as usize];
//! let mut ier = 0;
//!
//! unsafe {
//!     trmesh(
//!         &raw const n,
//!         x.as_ptr(),
//!         y.as_ptr(),
//!         z.as_ptr(),
//!         list.as_mut_ptr(),
//!         lptr.as_mut_ptr(),
//!         lend.as_mut_ptr(),
//!         &raw mut lnew,
//!         near.as_mut_ptr(),
//!         next.as_mut_ptr(),
//!         dist.as_mut_ptr(),
//!         &raw mut ier,
//!     );
//! }
//!
//! assert_eq!(ier, 0);
//! ```
//! # License
//!
//! The Rust bindings are licensed under MIT OR Apache-2.0 at your option.
//!
//! The bundled Fortran code is distributed under the GNU LGPL.
//!
//! # Attribution
//!
//! STRIPACK — Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere
//!
//! - Original author: Robert J. Renka, University of North Texas
//! - Reference: R. J. Renka, "Algorithm 772: STRIPACK: Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere",
//!   ACM Transactions on Mathematical Software, Vol. 23, No. 3, September 1997, pp. 416-434.
//! - DOI: <https://doi.org/10.1145/275323.275329>
//!
//! The Fortran 90 version used in this crate was prepared by John Burkardt and is distributed under the GNU LGPL.

pub mod ffi;
