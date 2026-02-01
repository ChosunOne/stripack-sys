use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;
use which::which;

fn main() {
    let target = env::var("TARGET").unwrap();
    let host = env::var("HOST").unwrap();
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    let lib_dir = manifest_dir.join("lib").join(&target);
    let lib_file = if target.contains("windows-msvc") {
        lib_dir.join("stripack.lib")
    } else {
        lib_dir.join("libstripack.a")
    };

    if lib_file.exists() {
        println!("cargo:warning=Using pre-built library for {target}");
        println!("cargo:rustc-link-search=native={}", lib_dir.display());
        println!("cargo:rustc-link-lib=static=stripack");
    } else if target == host {
        println!("cargo:warning=Compiling from source for {target}");
        compile_from_source(&out_dir);
        create_static_lib(&out_dir, &target);
    } else {
        panic!("No pre-built library for target: {target}\nHost: {host}");
    }

    if target.contains("windows") {
        println!("cargo:rustc-link-search=native=C:/msys64/mingw64/lib");
        println!("cargo:rustc-link-search=native=C:/msys64/mingw64/lib/gcc/x86_64-w64-mingw32/14");
        println!("cargo:rustc-link-lib=dylib=gfortran");
        println!("cargo:rustc-link-lib=dylib=quadmath");
    } else if target.contains("darwin") {
        println!("cargo:rustc-link-search=native=/opt/homebrew/lib");
        println!("cargo:rustc-link-search=native=/opt/homebrew/lib/gcc/14");
        println!("cargo:rustc-link-search=native=/usr/local/lib");
        println!("cargo:rustc-link-lib=gfortran");
    } else {
        println!("cargo:rustc-link-lib=gfortran");
    }
}

fn compile_from_source(out_dir: &Path) {
    let fortran_src = "fortran/stripack.f90";

    let obj_path = out_dir.join("stripack.o");

    let Some(compiler) = detect_fortran_compiler() else {
        panic!("No fortran compiler found.");
    };

    let status = Command::new(compiler)
        .args(["-c", "-O2", "-fPIC", fortran_src, "-o"])
        .arg(&obj_path)
        .status()
        .expect("Failed to compile Fortran source");

    assert!(status.success(), "gfortran compilation failed");
}

fn create_static_lib(out_dir: &Path, target: &str) {
    let fortran_src = "fortran/stripack.f90";
    let Some(archiver) = detect_archiver() else {
        panic!("No archiver found.");
    };
    let obj_path = out_dir.join("stripack.o");
    let lib_path = if target.contains("windows-msvc") {
        out_dir.join("stripack.lib")
    } else {
        out_dir.join("libstripack.a")
    };

    let args = match archiver.as_str() {
        "lib" => vec![],
        _ => vec!["rcs"],
    };

    let lib_path = match archiver.as_str() {
        "lib" => format!("/OUT:{}", lib_path.to_string_lossy()),
        _ => format!("{}", lib_path.to_string_lossy()),
    };

    let status = std::process::Command::new(&archiver)
        .args(args)
        .arg(&lib_path)
        .arg(&obj_path)
        .status()
        .expect("Failed to create static library");

    assert!(status.success(), "{archiver} failed");

    println!(
        "cargo:rustc-link-search=native={}",
        out_dir.to_string_lossy()
    );
    println!("cargo:rustc-link-lib=static=stripack");
    println!("cargo:rerun-if-changed={fortran_src}");
}

fn detect_fortran_compiler() -> Option<String> {
    if let Ok(fc) = env::var("FC") {
        return Some(fc);
    }

    let candidates = vec![
        "gfortran",
        "gfortran-14",
        "gfortran-13",
        "gfortran-12",
        "ifort",
        "ifx",
        "flang",
    ];

    for name in candidates {
        if let Ok(path) = which(name) {
            return Some(path.to_string_lossy().to_string());
        }
    }

    None
}

fn detect_archiver() -> Option<String> {
    if let Ok(ar) = env::var("AR") {
        return Some(ar);
    }

    let candidates = vec!["ar", "x86_64-w64-mingw32-ar", "lib"];
    for name in candidates {
        if let Ok(path) = which(name) {
            return Some(path.to_string_lossy().to_string());
        }
    }

    None
}
