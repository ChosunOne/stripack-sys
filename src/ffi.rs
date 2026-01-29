use std::ffi::{c_double, c_int};

unsafe extern "C" {
    /// Computes the area of a spherical triangle on the unit sphere.
    ///
    /// # Arguments
    /// * `v1[3]`, `v2[3]`, `v3[3]` - Input. The Cartesian coordinates of unit vectors (the three triangle
    ///   vertices in any order). These vectors, if nonzero, are implicitly scaled to have length `1`.
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

    /// Returns the boundary nodes of a triangulation. Given a triangulation of `n` nodes on the
    /// unit sphere created by `trmesh`, this subroutine returns an array containing the indexes (if any) of the counterclockwise sequence of boundary nodes, that is, the nodes on the boundary of the convex hull of the set of nodes. The boundary is empty if the nodes do not lie in a single hemisphere. The numbers of boundary nodes, arcs, and triangles are also returned.
    ///
    /// # Arguments
    ///
    /// * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    /// * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input. The data structure
    ///   defining the triangulation, created by `trmesh`.
    /// * `nodes` - Output. The ordered sequence of `nb` boundary node indexes in the range `1` to `n`.
    ///   For safety, the dimension of `nodes` should be `n`.
    /// * `nb` - Output. The number of boundary nodes.
    ///   `na`, `nt` - Output. The number of arcs and triangles, repectively, in the triangulation.
    #[link_name = "bnodes_"]
    pub fn bnodes(
        n: *const c_int,
        list: *const c_int,
        lptr: *const c_int,
        lend: *const c_int,
        nodes: *mut c_int,
        nb: *mut c_int,
        na: *mut c_int,
        nt: *mut c_int,
    );

    /// Converts from Cartesian to spherical coordinates (latitude, longitude, radius).
    ///
    /// # Arguments
    /// * `px`, `py`, `pz` - Input. The coordinates of `p`
    /// * `plat` - Output. The latitude of `p` in the range `-PI/2` to `PI/2`, or `0` if `pnrm = 0`.
    /// * `plon` - Output. The longitude of `p` in the range `-PI` to `PI`, or `0` if `p` lies on the
    ///   Z-axis.
    /// * `pnrm` - Output. The magnitude (Euclidean norm) of `p`.
    ///
    /// # Safety
    /// - All pointers must be valid and properly aligned
    #[link_name = "scoord_"]
    pub fn scoord(
        px: *const c_double,
        py: *const c_double,
        pz: *const c_double,
        plat: *mut c_double,
        plon: *mut c_double,
        pnrm: *mut c_double,
    );

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

    /// Convert a triangulation data structure into a triangle list.
    ///
    /// # Arguments
    /// * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    /// * `list`, `lptr`, `lend` - Input. Linked list data structure defining the triangulation.
    /// * `nrow` - Input. The number of rows (entries per triangle) reserved for the triangle list `ltri`. The
    ///   value must be 6 if only the vertex indexes and neighboring triangle indexes are to be stored, or
    ///   if arc indexes are also to be assigned and stored. Refer to `ltri`.
    /// * `nt` - Output. The number of triangles in the triangulation unless `ier /= 0`, which
    ///   in case `nt = 0`. `nt = 2n - nb - 2` if `nb >= 3` or `2n-4` if `nb = 0`, where `nb` is the
    ///   number of boundary nodes.
    /// * `ltri` - Output. The second dimension of `ltri` must be at least `nt`, where `nt` will be at
    ///   most `2*n - 4`. The `j`-th column contains the vertex nodal indexes (first three rows),
    ///   neighboring triangle indxes (second three rows), and, if `nrow = 9`, arc indexes (last three
    ///   rows) associated with triangle `j` for `j = 1, ..., nt`. The vertices are ordered counterclockwise
    ///   with the first vertex taken to be the one with smallest index. Thus `ltri[2][j]` and `ltri[3][j]` are larger than `ltri[1, j]` and index adjacent neighbors of node `ltri[1][j]`. For `i = 1, 2, 3`, `ltri[i + 3][j]` and `ltri[i + 6][j]` index the triangle and arc, respectively, which are opposite (not shared by) node `ltri[i][j]`, with `ltri[i + 3, j] = 0` if `ltri[i + 6][j]` indexes a boundary arc. Vertex indexes range from `1` to `n`, triangle indexes from `0` to `nt`, and, if included, arc indexes from `1` to `na`, where `na = 3n - nb - 3` if `nb >= 3` or `3n - 6` if `nb = 0`. The triangles are ordered on first (smallest) vertex indexes.
    /// * ier - Output. Error indicator.
    ///   0, if no errors were encountered.
    ///   1, if `n` or `nrow` is outside its valid range on input.
    ///   2, if the triangulation data structure (`list`, `lptr`, `lend`) is invalid. Note, however, that
    ///   these arrays are not completely tested for validity.
    #[link_name = "trlist_"]
    pub fn trlist(
        n: *const c_int,
        list: *const c_int,
        lptr: *const c_int,
        lend: *const c_int,
        nrow: *const c_int,
        nt: *mut c_int,
        ltri: *mut c_int,
        ier: *mut c_int,
    );

    /// Creates a Delaunay triangulation on the unit sphere.
    ///
    /// The Delaunay triangulation is defined as a set of (spherical) triangles with the following five
    /// properties:
    /// 1. The triangle vertices are nodes.
    /// 2. No triangle contains a node other than its vertices.
    /// 3. The interiors of the triangles are pairwise disjoint.
    /// 4. The union of triangles is the convex hull of the set of nodes (the smallest convex set that
    ///    contains the nodes). If the nodes are not contained in a single hemisphere, their convex hull
    ///    is the entire sphere and there are no boundary nodes. Otherwise, there are at least three
    ///    boundary nodes.
    /// 5. The interior of the circumcircle of each triangle contains no node.
    ///
    /// The first four properties define a triangulation, and the last property results in a
    /// triangulation which is as close as possible to equiangular in a certain sense and which is
    /// uniquely defined unless four or more nodes lie in a common plane. This property makes the
    /// triangulation well-suited for solving closest-point problems and for triangle based
    /// interpolation.
    ///
    /// Provided the nodes are randomly ordered, the algorithm has expected time complexity O(N*log(N))
    /// for most nodal distributions. Note, however, that the complexity may be as high as O(N**2) if,
    /// for example, the nodes are ordered on increasing latitude.
    ///
    /// Spherical coordinates (latitude and longitude) may be converted to Cartesian coordinates by
    /// `trans`.
    ///
    /// The following is a list of the software package modules which a user may wish to call directly:
    /// `addnod` - Updates the triangulation by appending a new node.
    /// `areas` - Returns the area of a spherical triangle
    /// `bnodes` - Returns an array containing the index of the boundary nodes (if any) in
    /// counterclockwise order. Counts of boundary nodes, triangles, and arcs are also returned.
    /// `circum` - Returns the circumcenter of a spherical triangle.
    /// `crlist` - Returns the set of triangle circumcenters (Voronoi vertices) and circumradii
    /// associated with a triangulation.
    /// `delarc` - Deletes a boundary arc from a triangulation.
    /// `edge` - Forces an arbitrary pair of nodes to be connected by an arc in the triangulation.
    /// `getnp` - Determines the ordered sequence of L closest nodes to a given node, along with the
    /// associated distances.
    /// `inside` - Locates a point relative to a polygon on the surface of the sphere.
    /// `intrsc` - Returns the point of intersection between a pair of great circle arcs.
    /// `jrand` - Generates a uniformly distributed pseudo-random integer.
    /// `left` - Locates a point relative to a great circle
    /// `nearnd` - Returns the index of the nearest node to an arbitrary point, along with its squared
    /// distance.
    /// `scoord` - Converts a point from Cartesian coordinates to spherical coordinates.
    /// `store` - Forces a value to be stored in main memory so that the precision of floating point
    /// numbers in memory locations rather than registers is computed
    /// `trans` - Transforms spherical coordinates into Cartesian coordinates on the unit sphere for
    /// input to `trmesh`
    /// `trlist` - Converts the triangulation data structure to a triangle list more suitable for use in
    /// a finite element code.
    /// `trlprt` - Creates a Delaunay triangulation of a set of nodes.
    /// `trplot` - Creates a level-2 Encapsulated Postscript (EPS) file containing a triangulation plot.
    /// `trprnt` - Prints the triangulation data structure and, optionally, the nodal coordinates.
    /// `vrplot` - Createsa level-2 Encapsulated Postscript (EPS) file containing a Voronoi diagram
    /// plot.
    ///
    /// # Arguments
    /// * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    /// * `x[n]`, `y[n]`, `z[n]` - Input. The coordinates of distinct nodes. `(x[k], y[k], z[k])` is referred to as node `k`, and `k` is referred to as a nodal index. It is required that `x[k]**2 + y[k]**2 + z[k]**2 = 1` for all `k`. The first three nodes must not be collinear (lie on a common great circle).
    /// * `list` - Output. `6 * (n - 2)` nodal indexes which, along with `lptr`, `lend`, and `lnew`
    ///   define the triangulation as a set of `n` adjacency lists; counterclockwise-ordered sequences of
    ///   neighboring nodes such that the first and last neighbors of a boundary node are boundary nodes
    ///   (the first neighbor of an interior node is arbitrary). In order to distinguish between interior
    ///   and boundary nodes, the last neighbor of each boundary node is represented by the negative of
    ///   its index.
    /// * `lptr` - Output. Set of pointers (`list` indexes) in one-to-one correspondence with the
    ///   elements of `list`. `list[lptr[i]]` indexes the node which follows `list[i]` in cyclical
    ///   counterclockwise order (the first neighbor follows the last neighbor).
    /// * `lend` - Output. `n` pointers to adjacency lists. `lend[k]` points to the last neighbor of
    ///   node `k`. `list[lend[k]] < 0` if and only if `k` is a boundary node.
    /// * `lnew` - Output. Pointer to the first empty location in `list` and `lptr` (list length plus
    ///   one). `list`, `lptr`, `lend` and `lnew` are not altered if `ier < 0`, and are incomplete if `0 <
    /// ier`.
    /// * `near` - Workspace. An array of `n` integers used to efficiently determine the nearest
    ///   triangulation node to each unprocessed node for use by `addnod`.
    /// * `next` - Workspace. An array of `n` integers used to efficiently determine the nearest triangulation node
    ///   to each unprocessed node for use by `addnod`.
    /// * `dist` - Workspace. An array of `n` floats used to efficiently determine the neareast
    ///   triangulation node to each unprocessed node for use by `addnod`.
    /// * `ier` - Output. An integer error indicator:
    ///   0, if no errors were ecountered.
    ///   -1, if `n < 3` on input.
    ///   -2, if the first three nodes are collinear.
    ///   L, if nodes L and M coincide for some L < M. The data structure represents a triangulation
    ///   of nodes 1 to M-1 in this case.
    #[link_name = "trmesh_"]
    pub fn trmesh(
        n: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
        near: *mut c_int,
        next: *mut c_int,
        dist: *mut c_double,
        ier: *mut c_int,
    );

}

#[cfg(test)]
mod test {
    use rstest::{fixture, rstest};

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

    #[rstest]
    #[case(1.0, 0.0, 0.0, 0.0, 0.0, 1.0)]
    #[case(0.0, 1.0, 0.0, 0.0, FRAC_PI_2, 1.0)]
    #[case(-1.0, 0.0, 0.0, 0.0, PI, 1.0)]
    #[case(0.0, -1.0, 0.0, 0.0, -FRAC_PI_2, 1.0)]
    #[case(0.0, 0.0, 1.0, FRAC_PI_2, 0.0, 1.0)]
    #[case(0.0, 0.0, -1.0, -FRAC_PI_2, 0.0, 1.0)]
    #[case(FRAC_1_SQRT_2, FRAC_1_SQRT_2, 0.0, 0.0, FRAC_PI_4, 1.0)]
    #[case(0.5, 0.5, FRAC_1_SQRT_2, FRAC_PI_4, FRAC_PI_4, 1.0)]
    #[case(2.0, 0.0, 0.0, 0.0, 0.0, 2.0)]
    #[case(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]
    fn test_scoord(
        #[case] px: f64,
        #[case] py: f64,
        #[case] pz: f64,
        #[case] expected_plat: f64,
        #[case] expected_plon: f64,
        #[case] expected_pnrm: f64,
    ) {
        let mut plat = 0.0;
        let mut plon = 0.0;
        let mut pnrm = 0.0;

        unsafe {
            scoord(
                &raw const px,
                &raw const py,
                &raw const pz,
                &raw mut plat,
                &raw mut plon,
                &raw mut pnrm,
            );
        };

        assert!((expected_plat - plat).abs() < f64::EPSILON);
        assert!((expected_plon - plon).abs() < f64::EPSILON);
        assert!((expected_pnrm - pnrm).abs() < f64::EPSILON);

        let n = 1;
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        unsafe {
            trans(
                &raw const n,
                &raw const plat,
                &raw const plon,
                &raw mut x,
                &raw mut y,
                &raw mut z,
            );
        };
        if pnrm == 0.0 {
            return;
        }

        let norm = (x * x + y * y + z * z).sqrt();

        assert!(
            (x / norm - px / pnrm).abs() < f64::EPSILON,
            "expected: {px} got: {x}"
        );
        assert!(
            (y / norm - py / pnrm).abs() < f64::EPSILON,
            "expected: {py} got: {y}"
        );
        assert!(
            (z / norm - pz / pnrm).abs() < f64::EPSILON,
            "expected: {pz} got: {z}"
        );
    }

    #[rstest]
    #[case(2, &[1.0, 0.0], &[0.0, 1.0], &[0.0, 0.0], -1)]
    fn test_trmesh(
        #[case] n: i32,
        #[case] x_in: &[f64],
        #[case] y_in: &[f64],
        #[case] z_in: &[f64],
        #[case] expected_ier: i32,
    ) {
        let mut x = Vec::with_capacity(n as usize);
        let mut y = Vec::with_capacity(n as usize);
        let mut z = Vec::with_capacity(n as usize);

        for i in 0..n as usize {
            let norm = (x_in[i].powi(2) + y_in[i].powi(2) + z_in[i].powi(2)).sqrt();
            x.push(x_in[i] / norm);
            y.push(y_in[i] / norm);
            z.push(z_in[i] / norm);
        }

        let list_size = if n >= 3 { 6 * (n - 2) } else { 0 } as usize;

        let mut list = vec![0i32; list_size];
        let mut lptr = vec![0i32; list_size];
        let mut lend = vec![0i32; n as usize];
        let mut lnew = 0i32;

        let mut near = vec![0i32; n as usize];
        let mut next = vec![0i32; n as usize];
        let mut dist = vec![0.0f64; n as usize];
        let mut ier = 0i32;

        unsafe {
            trmesh(
                &raw const n,
                x.as_ptr(),
                y.as_ptr(),
                z.as_ptr(),
                list.as_mut_ptr(),
                lptr.as_mut_ptr(),
                lend.as_mut_ptr(),
                &raw mut lnew,
                near.as_mut_ptr(),
                next.as_mut_ptr(),
                dist.as_mut_ptr(),
                &raw mut ier,
            );
        }

        assert_eq!(ier, expected_ier);

        if expected_ier != 0 {
            return;
        }

        assert!(lnew > 0, "lnew should be positive");
        assert!(lnew as usize <= list_size, "lnew should be within bounds");

        for i in 0..n as usize {
            let lend_val = lend[i];
            assert!(lend_val > 0, "lend[{i}] should be positive");
        }
    }

    #[fixture]
    fn tetrahedron() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let tx = vec![
            0.816_496_580_9,
            -0.408_248_290_50,
            -0.408_248_290_590_5,
            0.0,
        ];
        let ty = vec![0.0, FRAC_1_SQRT_2, -FRAC_1_SQRT_2, 0.0];
        let tz = vec![
            -0.577_350_269_2,
            0.577_350_269_2,
            0.577_350_269_2,
            -0.577_350_269_2,
        ];
        (tx, ty, tz)
    }

    #[fixture]
    fn octahedron() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let ox = vec![1.0, -1.0, 0.0, 0.0, 0.0, 0.0];
        let oy = vec![0.0, 0.0, 1.0, -1.0, 0.0, 0.0];
        let oz = vec![0.0, 0.0, 0.0, 0.0, 1.0, -1.0];

        (ox, oy, oz)
    }

    #[rstest]
    #[case(4, 6, 4)]
    #[case(4, 9, 4)]
    #[case(6, 6, 8)]
    #[case(6, 9, 8)]
    #[case(5, 6, 0)]
    fn test_trlist(
        #[case] n: i32,
        #[case] nrow: i32,
        #[case] expected_nt: i32,
        tetrahedron: (Vec<f64>, Vec<f64>, Vec<f64>),
        octahedron: (Vec<f64>, Vec<f64>, Vec<f64>),
    ) {
        let x = if n == 4 { tetrahedron.0 } else { octahedron.0 };
        let y = if n == 4 { tetrahedron.1 } else { octahedron.1 };
        let z = if n == 4 { tetrahedron.2 } else { octahedron.2 };

        let list_size = (6 * (n - 2)) as usize;
        let mut list = vec![0i32; list_size];
        let mut lptr = vec![0i32; list_size];
        let mut lend = vec![0i32; n as usize];
        let mut lnew = 0i32;

        let mut near = vec![0i32; n as usize];
        let mut next = vec![0i32; n as usize];
        let mut dist = vec![0.0f64; n as usize];
        let mut ier = 0i32;

        unsafe {
            trmesh(
                &raw const n,
                x.as_ptr(),
                y.as_ptr(),
                z.as_ptr(),
                list.as_mut_ptr(),
                lptr.as_mut_ptr(),
                lend.as_mut_ptr(),
                &raw mut lnew,
                near.as_mut_ptr(),
                next.as_mut_ptr(),
                dist.as_mut_ptr(),
                &raw mut ier,
            );
        };

        assert_eq!(ier, 0, "trmesh failed");

        let max_triangles = (2 * n - 4) as usize;
        let mut ltri = vec![0i32; (nrow as usize) * max_triangles];
        let mut nt = 0i32;
        let mut ier2 = 0i32;

        unsafe {
            trlist(
                &raw const n,
                list.as_ptr(),
                lptr.as_ptr(),
                lend.as_ptr(),
                &raw const nrow,
                &raw mut nt,
                ltri.as_mut_ptr(),
                &raw mut ier2,
            );
        }

        if nrow == 6 || nrow == 9 {
            assert_eq!(ier2, 0, "trlist should succeed");
            assert_eq!(nt, expected_nt, "triangle count mismatch");
            let v1 = ltri[0];
            assert!(v1 >= 1 && v1 <= n, "invalid vertex index");
            return;
        }

        assert_eq!(ier2, 1, "should fail with inalid nrow");
        assert_eq!(nt, 0, "nt should be 0 on error");
    }

    fn fibonacci_sphere(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let mut x_points = Vec::with_capacity(n);
        let mut y_points = Vec::with_capacity(n);
        let mut z_points = Vec::with_capacity(n);

        for i in 0..n {
            let y = 1.0 - (i as f64 / (n as f64 - 1.0)) * 2.0;
            let radius = (1.0 - y * y).sqrt();
            let theta = 2.0 * std::f64::consts::PI * i as f64 * phi;
            let x = radius * theta.cos();
            let z = radius * theta.sin();

            x_points.push(x);
            y_points.push(y);
            z_points.push(z);
        }

        (x_points, y_points, z_points)
    }

    #[fixture]
    fn hemisphere() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let x = vec![1.0, 0.0, -1.0, 0.0, 0.0];
        let y = vec![0.0, 1.0, 0.0, -1.0, 0.0];
        let z = vec![0.0, 0.0, 0.0, 0.0, 1.0];

        (x, y, z)
    }

    #[rstest]
    #[case(4, &[], 0, 6, 4)]
    #[case(6, &[], 0, 12, 8)]
    #[case(5, &[1, 2, 3, 4], 4, 9, 4)]
    fn test_bnodes(
        #[case] n: i32,
        #[case] expected_boundary: &[i32],
        #[case] expected_nb: i32,
        #[case] expected_na: i32,
        #[case] expected_nt: i32,
        tetrahedron: (Vec<f64>, Vec<f64>, Vec<f64>),
        hemisphere: (Vec<f64>, Vec<f64>, Vec<f64>),
    ) {
        let (x, y, z) = if n == 4 {
            tetrahedron
        } else if n == 5 {
            hemisphere
        } else {
            fibonacci_sphere(n as usize)
        };

        let list_size = (6 * (n - 2)) as usize;
        let mut list = vec![0i32; list_size];
        let mut lptr = vec![0i32; list_size];
        let mut lend = vec![0i32; n as usize];
        let mut lnew = 0i32;

        let mut near = vec![0i32; n as usize];
        let mut next = vec![0i32; n as usize];
        let mut dist = vec![0.0f64; n as usize];
        let mut ier = 0i32;

        unsafe {
            trmesh(
                &raw const n,
                x.as_ptr(),
                y.as_ptr(),
                z.as_ptr(),
                list.as_mut_ptr(),
                lptr.as_mut_ptr(),
                lend.as_mut_ptr(),
                &raw mut lnew,
                near.as_mut_ptr(),
                next.as_mut_ptr(),
                dist.as_mut_ptr(),
                &raw mut ier,
            );
        };

        assert_eq!(ier, 0, "trmesh failed");

        let mut nodes = vec![0i32; n as usize];
        let mut nb = 0i32;
        let mut na = 0i32;
        let mut nt = 0i32;

        unsafe {
            bnodes(
                &raw const n,
                list.as_ptr(),
                lptr.as_ptr(),
                lend.as_ptr(),
                nodes.as_mut_ptr(),
                &raw mut nb,
                &raw mut na,
                &raw mut nt,
            );
        };

        assert_eq!(nb, expected_nb, "boundary node count mismatch");
        assert_eq!(na, expected_na, "arc count mismatch");
        assert_eq!(nt, expected_nt, "triangle count mismatch");

        if nb > 0 {
            for i in 0..nb as usize {
                assert_eq!(nodes[i], expected_boundary[i], "boundary node {i} mismatch");
            }
        }
    }
}
