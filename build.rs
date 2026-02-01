fn main() {
    let fortran_src = "fortran/stripack.f90";

    let out_dir = std::env::var("OUT_DIR").unwrap();
    let obj_path = format!("{out_dir}/stripack.o");
    let lib_path = format!("{out_dir}/libstripack.a");

    let status = std::process::Command::new("gfortran")
        .args(["-c", "-O2", "-fPIC", fortran_src, "-o", &obj_path])
        .status()
        .expect("Failed to compile Fortran source");

    assert!(status.success(), "gfortran compilation failed");

    let status = std::process::Command::new("ar")
        .args(["rcs", &lib_path, &obj_path])
        .status()
        .expect("Failed to create static library");

    assert!(status.success(), "ar failed");

    println!("cargo:rustc-link-search=native={out_dir}");
    println!("cargo:rustc-link-lib=static=stripack");
    println!("cargo:rustc-link-lib=gfortran");
    println!("cargo:rerun-if-changed={fortran_src}");
}
