use std::ffi::{c_double, c_int};

unsafe extern "C" {
    /// Transform spherical coordinates into Cartesian coordinates
    /// on the unit sphere for input to `trmesh`. Storage for X and Y
    /// may coincide with storage for `rlat` and `rlon` if the latter
    /// need not be saved.
    ///
    /// # Arguments
    /// * `n` - Input. The number of nodes (points on the unit sphere) whose coordinates are to be transformed
    /// * `rlat` - Input. The latitudes of the nodes in radians.
    /// * `rlon` - Input. The longitudes of the nodes in radians.
    /// * `x`, `y`, and `z` - Output. The coordinates in the range `-1` to `1`. `x[i]**2 + y[i]**2 + z[i]**2 = 1` for `i` = 1 to `n`
    /// # Safety
    /// - All pointers must be valid and properly aligned
    /// - Arrays `rlat`, `rlon`, `x`, `y`, `z` must have length >= `*n`
    /// - Arrays `x`, `y`, `z` must not overlap with `rlat`, `rlon`
    /// - `n` must be > 0
    #[link_name = "trans_"]
    pub fn trans(
        n: *const c_int,
        rlat: *const c_double,
        rlon: *const c_double,
        x: *mut c_double,
        y: *mut c_double,
        z: *mut c_double,
    );

    /// Computes the area of a spherical triangle on the unit sphere.
    ///
    /// # Arguments
    /// * `v1[3]`, `v2[3]`, `v3[3]` - Input. The Cartesian coordinates of unit vectors (the three triangle
    /// vertices in any order). These vectors, if nonzero, are implicitly scaled to have length `1`.
    ///
    /// # Returns
    /// The area of the spherical triangle defined by `v1`, `v2`, and `v3`, in the range `0` to
    /// `2*PI` (the area of a hemisphere). `0` if and only if `v1`, `v2`, and `v3` lie in (or close to)
    /// a plane containing the origin.
    ///
    /// # Safety
    /// - All pointers must be valid and properly aligned
    /// - Arrays `v1`, `v2`, `v3` must have length == `3`
    #[link_name = "areas_"]
    pub fn areas(v1: *const c_double, v2: *const c_double, v3: *const c_double) -> c_double;
}

#[cfg(test)]
mod test {
    use rstest::rstest;

    use super::*;
    use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4, PI};

    #[rstest]
    #[case(1, &[FRAC_PI_2], &[0.0], &[0.0], &[0.0], &[1.0])]
    #[case(1, &[-FRAC_PI_2], &[0.0], &[0.0], &[0.0], &[-1.0])]
    #[case(1, &[0.0], &[0.0], &[1.0], &[0.0], &[0.0])]
    #[case(1, &[0.0], &[FRAC_PI_2], &[0.0], &[1.0], &[0.0])]
    #[case(1, &[0.0], &[PI], &[-1.0], &[0.0], &[0.0])]
    #[case(1, &[0.0], &[-FRAC_PI_2], &[0.0], &[-1.0], &[0.0])]
    #[case(1, &[std::f64::consts::FRAC_PI_4], &[0.0], &[FRAC_1_SQRT_2], &[0.0], &[FRAC_1_SQRT_2])]
    #[case(1, &[FRAC_PI_4], &[FRAC_PI_4], &[0.5], &[0.5], &[FRAC_1_SQRT_2])]
    #[case(3, &[FRAC_PI_2, -FRAC_PI_2, 0.0], &[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], &[0.0, 0.0, 0.0], &[1.0, -1.0, 0.0])]
    fn test_trans(
        #[case] n: i32,
        #[case] rlat: &[f64],
        #[case] rlon: &[f64],
        #[case] expected_x: &[f64],
        #[case] expected_y: &[f64],
        #[case] expected_z: &[f64],
    ) {
        let mut x = vec![0.0; n as usize];
        let mut y = vec![0.0; n as usize];
        let mut z = vec![0.0; n as usize];

        unsafe {
            trans(
                &raw const n,
                rlat.as_ptr(),
                rlon.as_ptr(),
                x.as_mut_ptr(),
                y.as_mut_ptr(),
                z.as_mut_ptr(),
            );
        }

        for i in 0..n as usize {
            assert!((x[i] - expected_x[i]).abs() < f64::EPSILON);
            assert!((y[i] - expected_y[i]).abs() < f64::EPSILON);
            assert!((z[i] - expected_z[i]).abs() < f64::EPSILON);
        }
    }

    #[rstest]
    #[case(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], &[0.0, 0.0, 1.0], FRAC_PI_2)]
    #[case(&[1.0, 0.0, 0.0], &[0.0, -1.0, 0.0], &[0.0, 0.0, 1.0], FRAC_PI_2)]
    #[case(&[-1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], &[0.0, 0.0, -1.0], FRAC_PI_2)]
    #[case(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], &[0.0, 0.0, -1.0], FRAC_PI_2)]
    #[case(&[1.0, 0.0, 0.0], &[1.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 0.0)]
    #[case(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0], &[-1.0, 0.0, 0.0], 0.0)]
    #[case(&[FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0], &[FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0], &[0.0, 0.0, 1.0], FRAC_PI_2)]
    fn test_areas(
        #[case] v1: &[f64; 3],
        #[case] v2: &[f64; 3],
        #[case] v3: &[f64; 3],
        #[case] expected_area: f64,
    ) {
        let area = unsafe { areas(v1.as_ptr(), v2.as_ptr(), v3.as_ptr()) };
        assert!(
            (area - expected_area).abs() < f64::EPSILON,
            "expected {expected_area} got {area}"
        );
    }
}
