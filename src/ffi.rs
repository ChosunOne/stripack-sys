use std::ffi::{c_double, c_int};

unsafe extern "C" {
    /**
    Adds a node to a triangulation of the convex hull of nodes `1, ..., k-1`, producing a
    triangulation of the convex hull of nodes `1, ..., k`.

    The algorithm consists of the following steps: node `k` is located relative to the
    triangulation ([trfind]), its index is added to the data structure ([intadd] or [bdyadd]), and a sequence of swaps ([swptst] and [swap]) are applied to the arcs opposite `k` so that all arcs incident on node `k` and opposite node `k` are locally optimal (statisfy the circumcircle test).

    Thus, if a Delaunay triangulation of nodes `1` through `k-1` is input, a Delaunay
    triangulation of nodes `1` through `k` will be output.

    # Arguments

    * `nst` - Input. The index of a node at which [trfind] begins its search. Search time depends on
      the proximity of this node to `k`. If `nst < 1`, the search is begun at node `k-1`.

    * `k` - Input. The nodal index (index for `x`, `y`, and `z`, and `lend`) of the new node
      to be added. `4 <= k`.

    * `x[k]`, `y[k]`, `z[k]` - Input. The coordinates of the nodes.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[k]`, `lnew` - Input/Output. On input, the data
      structure associated with the triangulation of nodes `1` to `k-1`. On output, the data has been
      updated to include node `k`. The array lengths are assumed to be large enough to add node `k`.
      Refer to [trmesh].
    * `ier` - Output. Error indicator:
      * `0` if no errors were encountered.
      * `-1` if `k` is outside its valid range on input.
      * `-2` if all nodes (including `k`) are collinear (lie on a common geodesic).
      * `l` if nodes `l` and `k` coincide for some `l < k`
    */
    #[link_name = "addnod_"]
    pub fn addnod(
        nst: *const c_int,
        k: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
        ier: *mut c_int,
    );

    /**
    Computes the arc cosine function, with argument truncation.

    If you call your system `acos` routine with an input argument that is outside the range `[-1.0,
    1.0]` you may get an unpleasant surprise. This routine truncates arguments outside the range.
    # Arguments

    * `c` - Input. The argument

    # Returns
    * An angle whose cosine is `c`
    */
    #[link_name = "arc_cosine_"]
    pub fn arc_cosine(c: *const c_double) -> c_double;

    /**
    Computes the area of a spherical triangle on the unit sphere.

    # Arguments
    * `v1[3]`, `v2[3]`, `v3[3]` - Input. The Cartesian coordinates of unit vectors (the three triangle
      vertices in any order). These vectors, if nonzero, are implicitly scaled to have length `1`.

    # Returns
    The area of the spherical triangle defined by `v1`, `v2`, and `v3`, in the range `0` to
    `2*PI` (the area of a hemisphere). `0` if and only if `v1`, `v2`, and `v3` lie in (or close to)
    a plane containing the origin.

    # Safety
    - All pointers must be valid and properly aligned
    - Arrays `v1`, `v2`, `v3` must have length == `3`
    */
    #[link_name = "areas_"]
    pub fn areas(v1: *const c_double, v2: *const c_double, v3: *const c_double) -> c_double;

    /**
    Adds a boundary node to a triangulation.

    This subroutine adds a boundary node to a triangulation of a set
    of `kk - 1` points on the unit sphere. The data structure is
    updated with the insertion of node `kk`, but no optimizaiton is
    performed.

    This routine is identical to the similarly named routine in
    TRIPACK.

    # Arguments

    * `kk` - Input. The index of a node to be connected to the
      sequence of all visible boundary nodes. 1 <= `kk` and
      `kk` must not be equal to `i1` or `i2`.
    * `i1` - Input. The first (rightmost as viewed from `kk`)
      boundary node in the triangulation that is visible from
      node `kk` (the line segment `kk-i1` intersects no arcs).
    * `i2` - Input. The last (leftmost) boundary node that is visible
      from node `kk`. `i1` and `i2` may be determined by [trfind].
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]`, `lnew` - Input/output. The
      triangulation data structure created by [trmesh]. Nodes `i1` and `i2` must be included in the triangulation. On output, the data structure is updated with the addition of node `kk`. Node `kk` is connected to `i1`, `i2`, and all boundary nodes in between.
     */
    #[link_name = "bdyadd_"]
    pub fn bdyadd(
        kk: *const c_int,
        i1: *const c_int,
        i2: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
    );

    /**
    Returns the boundary nodes of a triangulation.

    Given a triangulation of `n` nodes on the
    unit sphere created by [trmesh], this subroutine returns an array containing the indexes (if any) of the counterclockwise sequence of boundary nodes, that is, the nodes on the boundary of the convex hull of the set of nodes. The boundary is empty if the nodes do not lie in a single hemisphere. The numbers of boundary nodes, arcs, and triangles are also returned.

    # Arguments

    * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input. The data structure
      defining the triangulation, created by [trmesh].
    * `nodes` - Output. The ordered sequence of `nb` boundary node indexes in the range `1` to `n`.
      For safety, the dimension of `nodes` should be `n`.
    * `nb` - Output. The number of boundary nodes.
      `na`, `nt` - Output. The number of arcs and triangles, repectively, in the triangulation.
    */
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

    /**
    Returns the circumcenter of a spherical triangle.

    Returns the circumcenter of a spherical triangle on the unit sphere: the point on the
    sphere surface that is equally distant from the three triangle vertices and lies in the same
    hemisphere, where distance is taken to be arc-length on the sphere surface.

    # Arguments
    * `v1[3]`, `v2[3]`, `v3[3]` - Input. The coordinates of the three triangle vertices (unit
      vectors) in counterclockwise order.
    * `c[3]` - Output. The coordinates of the circumcenter unless `0 < ier`, which in case `c` is
      not defined. `c = (v2 - v1) X (v3 - v1)` normalized to a unit vector.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered.
      * `1`, if `v1`, `v2`, and `v3` lie on a common line: `(v2 - v1) X (v3 - v1) = 0`.
    */
    #[link_name = "circum_"]
    pub fn circum(
        v1: *const c_double,
        v2: *const c_double,
        v3: *const c_double,
        c: *mut c_double,
        ier: *mut c_int,
    );

    /**
    Connects an exterior node to boundary nodes, covering the sphere.

    This subroutine connects an exterior node `kk` to all boundary nodes of a triangulation of `kk - 1` points on the unit sphere, producing a triangulation that covers the sphere. The data structure is updated with the addition of node `kk`, but no optimization is performed. All boundary nodes must be visible from node `kk`.

    # Arguments

    * `kk` - Input. Index of the node to be connected to the set of all boundary nodes. `4 <= kk`.
    * `n0` - Input. Index of a boundary node (in the range `1` to `kk - 1`). `n0` may be
      determined by [trfind].
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]`, `lnew` - Input/output. The
      triangulation data structure created by [trmesh]. Node `n0` must be included in the triangulation. On output, updated with the addition of node `kk` as the last entry. The updated triangulation contains no boundary nodes.
    */
    #[link_name = "covsph_"]
    pub fn covsph(
        kk: *const c_int,
        n0: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
    );

    /**
    Returns triangle circumcenters and other information.

    Given a Delaunay triangulation of nodes on the surface of the unit sphere, this subroutine returns the set of triangle circumcenters corresponding to Voronoi vertices, along with the circumradii and a list of triangle indexes `listc` stored in one-to-one correspondence with `list/lptr` entries.

    A triangle circumcenter is the point (unit vector) lying at the same angular distance from the three vertices and contained in the same hemisphere as the vertices. (Note that the negative of a circumcenter is also equidistant from the vertices.) If the triangulation covers the surface, the Voronoi vertices are the circumenters of the triangles in the Delaunay triangulation. `lptr`, `lend`, and `lnew` are not altered in this case.

    On the other hand, if the nodes are contained in a single hemisphere, the triangulation is implicitly extended to the entire surface by adding pseudo-arcs (of length greater than 180 degrees) between boundary nodes forming pseudo-triangles whose 'circumcenters' are included in the list. This extension to the triangulation actually consists of a triangulation of the set of boundary nodes in which the swap test is reversed (a non-empty circumcircle test). The negative circumcenters are stored as the pseudo-triangle 'circumcenters'. `listc`, `lptr`, `lend`, and `lnew` contain a data structure corresponding to the extended triangulation (Voronoi diagram), but `list` is not altered in this case. Thus, if it is necessary to retain the original (unextended) triangulation data structure, copies of `lptr` and `lnew` must be saved before calling this routine.

    # Arguments

    * `n` - Input. The number of nodes in the triangulation `3 <= n`. Note that, if `n = 3`,
      there are only two Voronoi vertices separated by 180 degrees, and the Voronoi regions are not well defined.
    * `ncol` - Input. The number of columns reserved for `ltri`. This must be at least `nb - 2`,
      where `nb` is the number of boundary nodes.
    * `x[n]`, `y[n]`, `z[n]` - Input. The coordinates of the nodes (unit vectors).
    * `list[6 * (n - 2)]` - Input. The set of adjacency lists. Refer to [trmesh].
    * `lend[n]` - Input. The set of pointers to ends of adjacency lists. Refer to [trmesh].
    * `lptr[6 * (n - 2)]` - Input/output. Pointers associated with `list`. Refer to [trmesh]. On
      output, pointers associated with `listc`. Updated for the addition of pseudo-triangles if the original triangulation contains boundary nodes (`0 < nb`).
    * `lnew` - Input/output. On input, a pointer to the first empty location in `list` and `lptr` (list length plus one). On output, pointer to the first empty location in `listc` and `lptr` (list length plus one). `lnew` is not altered if `nb = 0`.
    * `ltri[6][ncol]` - Output. Triangle list whose first `nb - 2` columns contain the indexes
      of a clockwise-ordered sequence of vertices (first three rows) followed by the `ltri` column indexes of the triangles opposite the vertices (or 0 denoting the exterior region) in the last three rows. This array is not generally of any further use outside this routine.
    * `listc[3 * nt]` - Output. Where `nt = 2 * n - 4` is the number of triangles in the
      triangulation (after extending it to cover the entire surface if necessary). Contains the triangle indexes (indexes to `xc`, `yc`, `zc`, and `rc`) stored in 1-1 correspondence with `list`/`lptr` entries (or entries that would be stored in `list` for the extended triangulation): the index of triangle (`n1`, `n2`, `n3`) is stored in `listc[k]`, `listc[l]`, and `listc[m]`, where `list[k]`, `list[l]`, and `list[m]` are the indexes of `n2` as a neighbor of `n1`, `n3` as a neighbor of `n2`, and `n1` as a neighbor of `n3`. The Voronoi region associated with a node is defined by the CCW-ordered sequence of circumcenters in one-to-one correspondence with its adjacency list (in the extended triangulation).
    * `nb` - Output. The number of boundary nodes unless `ier = 1`.
    * `xc[2 * n - 4]`, `yc[2 * n - 4]`, `zc[2 * n - 4]`, the coordinates of the triangle
      circumcenters (Voronoi vertices). `xc[i]**2 + yc[i]**2 + zc[i]**2 = 1`. The first `nb - 2`
      entries correspond to pseudo-triangles if `0 < nb`.
    * `rc[2 * n - 4]` - Output. The circumradii (the arc lengths or angles between the
      circumcenters and associated triangle vertices) in 1-1 correspondence with circumcenters.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered.
      * `1`, if `n < 3`.
      * `2`, if `ncol < nb - 2`.
      * `3`, if a triangle is degenerate (has vertices lying on a common geodesic).
    */
    #[link_name = "crlist_"]
    pub fn crlist(
        n: *const c_int,
        ncol: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        list: *const c_int,
        lend: *const c_int,
        lptr: *mut c_int,
        lnew: *mut c_int,
        ltri: *mut c_int,
        listc: *mut c_int,
        nb: *mut c_int,
        xc: *mut c_double,
        yc: *mut c_double,
        zc: *mut c_double,
        rc: *mut c_double,
        ier: *mut c_int,
    );

    /**
    Deletes a boundary arc from a triangulation.

    This subroutine deletes a boundary arc from a triangulation. It may be used to remove a null triangle from the convex hull boundary. Note, however, that if the union of triangles is rendered nonconvex, subroutines [delnod], [edge], and [trfind] (and hence [addnod]) may fail. Also, function [nearnd] should not be called following an arc deletion.

    This routine is identical to the similarly named routine in TRIPACK.

    # Arguments

    * `n` - Input. The number of nodes in the triangulation. `4 <= n`.
    * `io1`, `io2` - Input. Indexes (in the range of `1` to `n`) of a pair of adjacent boundary
      nodes defining the arc to be removed.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]`, `lnew` - Input/output. The triangulation
      data structure created by [trmesh]. On output, updated with the removal of arc `io1-io2` unless `0 < ier`.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered.
      * `1`, if `n`, `io1`, or `io2` is outside its valid range, or `io1 = io2`.
      * `2`, if `io1-io2` is not a boundary arc.
      * `3`, if the node opposite `io1-io2` is already a boundary node, and thus `io1` or `io2` has only two neighbors or a deletion would result in two triangulations sharing a single node.
      * `4`, if one of the nodes is a neighbor of the other, but not vice versa, implying an invalid triangulation data structure.
     */
    #[link_name = "delarc_"]
    pub fn delarc(
        n: *const c_int,
        io1: *const c_int,
        io2: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
        ier: *mut c_int,
    );

    /**
    Deletes a neighbor from the adjacency list.

    This subroutine deletes a neighbor `nb` from the adjacency list of node `n0` (but `n0` is not deleted from the adjacency list of `nb`) and, if `nb` is a boundary node, makes `n0` a boundary node.

    For pointer (`list` index) `lph` to `nb` as a neighbor of `n0`, the empty `list`, `lptr` location `lph` is filled in with the values at `lnew - 1`, pointer `lnew - 1` (in `lptr` and possibly in `lend`) is changed to `lph`, and `lnew` is decremented.

    This requires a search of `lend` and `lptr` entailing an expected operation count of O(n).

    This routine is identical to the similarly named routine in TRIPACK.

    # Arguments

    * `n0`, `nb` - Input. Indexes, in the range `1` to `n`, of a pair of nodes such that `nb` is
      a neighbor of `n0`. (`n0` need not be a neighbor of `nb`.)
    * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]`, `lnew` - Input/Output. The data structure defining
      the triangulation. On output, updated with the removal of `nb` from the adjacency list of
      `n0` unless `lph < 0`.
    * `lph` - Input/output. List pointer to the hole (`nb` as a neighbor of `n0`) filled in by
      the values at `lnew - 1` or error indicator:
      `> 0`, if no errors were encountered.
      `-1`, if `n0`, `nb`, or `n` is outside its valid range.
      `-2`, if `nb` is not a neighbor of `n0`.
    */
    #[link_name = "delnb_"]
    pub fn delnb(
        n0: *const c_int,
        nb: *const c_int,
        n: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
        lph: *mut c_int,
    );

    /**
    Deletes a node from a triangulation.

    This subroutine deletes node `k` (along with all arcs incident on node `k`) from a triangulation of `n` nodes on the unit sphere, and inserts arcs as necessary to produce a triangulation of the remaining `n - 1` nodes. If a Delaunay triangulation is input, a Delaunay triangulation will result, and thus, [delnod] reverses the effect of a call to [addnod].

    Note that the deletion may result in all remaining nodes being collinear. This situation is not flagged.

    # Arguments

    * `k` - Input. The index (for `x`, `y`, and `z`) of the node to be deleted. `1 <= k <= n`.
    * `n` - Input/output. The number of nodes in the triangulation. `4 <= n`. Note that `n` will
      be decremented following the deletion.
    * `x[n]`, `y[n]`, `z[n]` - Input/output. The coordinates of the nodes in the triangulation.
      On output, updated with elements `k + 1, ..., n + 1` shifted up one position, thus overwriting element `k`, unless `1 <= ier <= 4`.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]`, `lnew` - Input/output. The data
      structure defining the triangulation, created by [trmesh]. On output, updated to reflect the deletion unless `1 <= ier <= 4`. Note that the data structure may have been altered if `3 < ier`.
    * `lwk` - Input/output. The number of columns reserved for `iwk`. `lkw` must be at least
      `nnb-3`, where `nnb` is the number of neighbors of node `k`, including an extra pseudo-node if `k` is a boundary node. On output, the number of `iwk` columns required unless `ier = 1` or `ier = 3`.
    * `iwk[2][lwk]` - Output. The indexes of the endpoints of the new arcs added unless `lwk = 0` or `1 <= ier <= 4`. (Arcs are associated with columns.)
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered
      * `1`, if `k` or `n` is outside its valid range or `lwk < 0` on input.
      * `2`, if more space is required in `iwk`. Refer to `lwk`.
      * `3`, if the triangulation data structure is invalid on input.
      * `4`, if `k` indexes an interior node with four or more neighbors, none of which can be swapped out due to collinearity, and `k` cannot therefore be deleted.
      * `5`, if an error flag (other than `ier = 1`) was returned by [optim]. An error message is written to the standard output unit in this case.
      * `6`, if error flag `1` was returned by [optim]. This is not necessarily an error, but the arcs may not be optimal.
    */
    #[link_name = "delnod_"]
    pub fn delnod(
        k: *const c_int,
        n: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
        lwk: *mut c_int,
        iwk: *mut c_int,
        ier: *mut c_int,
    );

    /**
    Swaps arcs to force two nodes to be adjacent.

    Given a triangulation of `n` nodes and a pair of nodal indexes `in1` and `in2`, this routine
    swaps arcs as necessary to force `in1` and `in2` to be adjacent. Only arcs which intersect `in1
    -in2` are swapped out. If a Delaunay triangulation is input, the resulting triangulation is as
    close as possible to a Delaunay triangulation in the sense that all arcs other than `in1-in2`
    are locally optimal.

    A sequence of calls to [edge] may be used to force the presence of a set of edges defining the
    boundary of a non-convex and/or multiply connected region, or to introduce barriers into the
    triangulation. Note that [getnp] will not necessarily return closest nodes if the triangulation
    has been constrained by a call to [edge]. However, this is appropriate in some applications,
    such as triangle-based interpolation on a nonconvex domain.

    # Arguments

    * `in1`, `in2` - Input. The indexes (of `x`, `y`, and `z`) in the range `1` to `n` defining a pair of
      nodes to be connected by an arc.
    * `x[n]`, `y[n]`, `z[n]` - Input. The coordinates of the nodes.
    * `lwk` - Input/output. On input, the number of columns reserved for `iwk`. This must be at
      least `ni`, the number of arcs that intersect `in1-in2`. (`ni` is bounded by `n - 3`.) On
      output, the number of arcs which intersect `in1-in2` (but not more than the input value of `lwk`) unless `ier = 1` or `ier = 3`. `lwk = 0` if and only if `in1` and `in2` were adjacent (or `lwk = 0`) on input.
    * `iwk[2 * lwk]` - Output. The indexes of the endpoints of the new arcs other than `in1-in2`
      unless `0 < ier` or `lwk = 0`. New arcs to the left of `in1->in2` are stored in the first `k-1`
      columns (left portion of `iwk`), column `k` contains zeros, and new arcs to the right of
      `in1->in2` occupy columns `k + 1, ..., lwk`. (`k` can be determined by searching `iwk` for the
      zeros.)
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input/output. The data structure
      defining the triangulation, created by [trmesh]. On output, updated if necessary to refelct the
      presence of an arc connecting `in1` and `in2` unless `0 < ier`. The data structure has been
      altered if `4 <= ier`.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered
      * `1`, if `in1 < 1`, `in2 < 1`, `in1 = in2`, or `lwk < 0` on input.
      * `2`, if more space is required in `iwk`. Refer to `lwk`.
      * `3`, if `in1` and `in2` could not be connected due to either an invalid data structure or
        collinear nodes (and floating point error).
      * `4`, if an error flag other than `ier = 1` was returned by [optim]
      * `5`, if error flag `1` was returned by [optim]. This is not necessarily an error, but the arcs
        other than `in1-in2` may not be optimal
    */
    #[link_name = "edge_"]
    pub fn edge(
        in1: *const c_int,
        in2: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        lwk: *mut c_int,
        iwk: *mut c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        ier: *mut c_int,
    );

    /**
    Gets the next nearest node to a given node.

    Given a Delaunay triangulation of `n` nodes on the unit sphere and an array `npts` containing the indexes of `l-1` nodes ordered by angular distance from `npts[1]`, this routine sets `npts[l]` to the index of the next node in the sequence -- the node, other than `npts[1], ..., npts[l - 1]`, that is closest to `npts[1]`. Thus, the ordered sequence of `k` closest nodes to `n1` (including `n1`) may be determined by `k - 1` calls to `getnp` with `npts[1] = n1` and `l = 2,3,..., k` for `k >= 2`.

    The algorithm uses the property of a Delaunay triangulation that the `k`-th closest node to `n1` is a neighbor of one of the `k-1` closest nodes to `n1`.

    # Arguments

    * `x[n]`, `y[n]`, `z[n]` - Input. The coordinates of the nodes.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input. The triangulation data
      structure, created by [trmesh].
    * `l` - Input. The number of nodes in the sequence on output. `2 <= l <= n`.
    * `npts[l]` - Input/output. On input, the indexes of the `l - 1` closest nodes to `npts[1]`
      in the first `l - 1` locations. On output, updated with the index of the `l`-th closest
      node to `npts[1]` in position `l` unless `ier = 1`.
    * `df` - Output. The value of an increasing function (negative cosine) of the angular
      distance between `npts[1]` and `npts[l]` unless `ier = 1`.
    * `ier` - Output. Error indicator:
      `0`, if no errors were encountered.
      `1`, if `l < 2`.
    */
    #[link_name = "getnp_"]
    pub fn getnp(
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        list: *const c_int,
        lptr: *const c_int,
        lend: *const c_int,
        l: *const c_int,
        npts: *mut c_int,
        df: *mut c_double,
        ier: *mut c_int,
    );

    /**
    Inserts `k` as a neighbor of `n1`.

    This subroutine inserts `k` as a neighbor of `n1` following `n2`,
    where `lp` is the `list` pointer of `n2` as a neighbor of `n1`.
    Note that, if `n2` is the last neighbor of `n1`, K will become
    the first neighbor (even if `n1` is a boundary node).

    This routine is identical to the similarly named routine in
    TRIPACK.

    # Arguments

    * `k` - Input. The index of the node to be inserted.
    * `lp` - Input. The `list` pointer of `n2` as a neighbor of `n1`.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lnew` - Input/Output.
      The data structure defining the triangulation, created by [trmesh]. On output, updated with the addition of node `k`
    */
    #[link_name = "insert_"]
    pub fn insert(
        k: *const c_int,
        lp: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lnew: *mut c_int,
    );

    /**
    Determines if a point is inside a polygonal region.

    This function locates a point `p` relative to a polygonal region `r` on the surface of the unit sphere, returning `inside = true` if and only if `p` is contained in `r`. `r` is defined by a cyclically ordered sequence of vertices which form a positively-oriented simple closed curve. Adjacent vertices need not be distinct but the curve must not be self-intersecting. Also, while polygon edges are by definition restricted to a single hemisphere, `r` is not so restricted. Its interior is the region to the left as the vertices are traversed in order.

    The algorithm consists of selecting a point `q` in `r` and then finding all points at which the great circle defined by `p` and `q` intersects the boundary of `r`. `p` lies inside `r` if and only if there is an even number of intersection points between `q` and `p`. `q` is taken to be a point immediately to the left of a directed boundary edge -- the first one that results in no consistency-check failures.

    If `p` is close to the polygon boundary, the problem is ill-conditioned and the decision may be incorrect. Also, an incorrect decision may result from a poor choice of `q` (if, for example, a boundary edge lies on the great circle defined by `p` and `q`). A more reliable result could be obtained by a sequence of calls to [inside] with the vertices cyclically permuted before each call (to alter the choice of `q`).

    # Arguments

    * `p[3]` - Input. The coordinates of the point (unit vector) to be located.
    * `lv` - Input. The length of the arrays `xv`, `yv`, and `zv`.
    * `xv[lv]`, `yv[lv]`, `zv[lv]` - Input. The coordinates of unit vectors (points on the unit sphere).
    * `nv` - Input. The number of vertices in the polygon. `3 <= nv <= lv`.
    * `listv[nv]` - Input. The indexes (for `xv`, `yv`, and `zv`) of a cyclically-ordered (and
      CCW-ordered) sequence of vertices that define `r`. The last vertex (indexed by `listv[nv]`) is followed by the first (indexed by `listv[1]`). `listv` entries must be in the range `1` to `lv`.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered.
      * `1`, if `lv` or `nv` is outside its valid range.
      * `2`, if a `listv` entry is outside its valid range.
      * `3`, if the polygon boundary was found to be self-intersecting. This error will not necessarily be detected.
      * `4`, if every choice of `q` (one for each boundary edge) led to failure of some internal consistency check. The most likely cause of this error is invalid input: `p = (0, 0, 0)`, a null or self-intersecting polygon, etc.

    # Returns

    True if and only if `p` lies inside `r` unless `ier != 0`, in which case the value is not altered.

    */
    #[link_name = "inside_"]
    pub fn inside(
        p: *const c_double,
        lv: *const c_int,
        xv: *const c_double,
        yv: *const c_double,
        zv: *const c_double,
        nv: *const c_int,
        listv: *const c_int,
        ier: *mut c_int,
    ) -> bool;

    /**
    Adds an interior node to a triangulation.

    This subroutine adds an interior node to a triangulation of a set of points on the unit sphere. The data structure is updated with the insertion of node `kk` into the triangle whose vertices are `i1`, `i2`, and `i3`. No optimization of the triangulation is performed.

    This routine is identical to the similarly named routine in TRIPACK.

    # Arguments

    * `kk` - Input. The index of the node to be inserted. `1 <= kk` and `kk` must not be equal
      to `i1`, `i2`, or `i3`.
    * `i1`, `i2`, `i3` - Input. Indexes of the counterclockwise-orderd sequence of vertices of a
      triangle which contains node `kk`.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]`, `lnew` - Input/output. The data structure
      defining the triangulation, created by [trmesh]. Triangle (`i1`, `i2`, `i3`) must be included in the triangulation. On output, updated with the addition of node `kk`. `kk` will be connected to nodes `i1`, `i2`, and `i3`.
    */
    #[link_name = "intadd_"]
    pub fn intadd(
        kk: *const c_int,
        i1: *const c_int,
        i2: *const c_int,
        i3: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lnew: *mut c_int,
    );

    /**
    Finds the intersection of two great circles.

    Given a great circle `c` and points `p1` and `p2` defining an arc `a` on the surface of the unit sphere, where `a` is the shorter of the two portions of the great circle `c12` associated with `p1` and `p2`, this subroutine returns the point of intersection `p` between `c` and `c12` that is closer to `a`. Thus, if `p1` and `p2` lie in opposite hemispheres defined by `c`, `p` is the point of intersection of `c` with `a`.

    # Arguments

    * `p1[3]`, `p2[3]` - Input. The coordinates of unit vectors.
    * `cn[3]` - Input. The coordinates of a nonzero vector which defines `c` as the intersection
      of the plane whose normal is `cn` with the unit sphere. Thus, if `c` is to be the great
      circle defined by `p` and `q`, `cn` should be `p X q`.
    * `p` - Output. Point of intersection defined above unless `ier` is not 0, in which case `p`
      is not altered.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered.
      * `1`, if `<cn, p1> = <cn, p2>`. This occurs if and only if `p1 = p2` or `cn = 0` or there are two intersection points at the same distance from `a`.
      * `2` if `p2 = -p1` and the definition of `a` is therefore ambiguous.
    */
    #[link_name = "intrsc_"]
    pub fn intrsc(
        p1: *const c_double,
        p2: *const c_double,
        cn: *const c_double,
        p: *mut c_double,
        ier: *mut c_int,
    );

    /**
    Returns a random integer between `1` and `n`.

    This function returns a uniformly distributed pseudorandom integer in the range `1` to `n`.

    # Arguments

    * `n` - Input. The maximum value to be returned.
    * `ix`, `iy`, `iz` - Input/output. The Seeds initialized to values in the range `1` to
      `30,000` before the first call to `jrand`, and not altered between subsequent calls (unless a sequence of random numbers is to be repeated by reinitializing the seeds).

    # Returns
    A random integer in the range `1` to `n`.
    */
    #[link_name = "jrand_"]
    pub fn jrand(n: *const c_int, ix: *mut c_int, iy: *mut c_int, iz: *mut c_int) -> c_int;

    /**
    Determines whether a node is left of a plane through the origin.

    This function determines whether node `n0` is in the (closed) left hemisphere defined by the
    plane containing `n1`, `n2`, and the origin, where left is defined relative to an observer at
    `n1` facing `n2`.

    # Arguments

    * `x1`, `y1`, `z1` - Input. The coordinates of `n1`.
    * `x2`, `y2`, `z2` - Input. The coordinates of `n2`.
    * `x0`, `y0`, `z0` - Input. The coordinates of `n0`.

    # Returns
    True if and only if `n0` is in the closed left hemisphere.
    */
    #[link_name = "left_"]
    pub fn left(
        x1: *const c_double,
        y1: *const c_double,
        z1: *const c_double,
        x2: *const c_double,
        y2: *const c_double,
        z2: *const c_double,
        x0: *const c_double,
        y0: *const c_double,
        z0: *const c_double,
    ) -> bool;

    /**
    Returns the index of `nb` in the adjacency list.

    This function returns the index (`list` pointer) of `nb` in the adjacency list for `n0`, where `lpl = lend[n0]`. This function is identical to the similarly named function in TRIPACK.

    # Arguments
    * `lpl` - Input. Is `lend[n0]`
    * `nb` - Input. The index of the node whose pointer is to be returned. `nb` must be connected to
      `n0`.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]` - Input. The data structure defining the
      triangulation, created by [trmesh].

    # Returns
    Pointer `lstptr` such that `list[lstptr] = nb` or `list[lstptr] = -nb`, unless `nb`
      is not a neighbor of `n0`, in which case `lstptr = lpl`.
    */
    #[link_name = "lstptr_"]
    pub fn lstptr(
        lpl: *const c_int,
        nb: *const c_int,
        list: *const c_int,
        lptr: *const c_int,
    ) -> c_int;

    /**
    Returns the number of neighbors of a node.

    This function returns the number of neighbors of a node `n0` in a triangulation created by [trmesh].

    The number of neighbors also gives the order of the Voronoi polygon containing the point. Thus, a neighbor count of `6` means the node is contained in a `6`-sided Voronoi region.

    This function is identical to the similarly named function in TRIPACK.

    # Arguments

    * `lpl` - Input. Pointer to the last neighbor of `n0`.
    * `lptr[6 * (n - 2)]` - Input. Pointers associated with `list`.

    # Returns

    The number of neighbors of `n0`.
    */
    #[link_name = "nbcnt_"]
    pub fn nbcnt(lpl: *const c_int, lptr: *const c_int) -> c_int;

    /**
    Returns the nearest node to a given point.

    Given a point `p` on the surface of the unit sphere and a Delaunay triangulation created by [trmesh], this function returns the index of the nearest triangulation node to `p`.

    The algorithm consists of implicitly adding `p` to the triangulation, finding the nearest neighbor to `p`, and implicitly deleting `p` from the triangulation. Thus, it is based on the fact that, if `p` is a node in a Delaunay triangulation, the nearest node to `p` is a neighbor of `p`.

    For large values of `n`, this procedure will be faster than the naive approach of computing the distance from `p` to every node.

    Note that the number of candidates for [nearnd] (neighbors of `p`) is limited to `lmax` defined in the arguments below.

    # Arguments

    * `p[3]` - Input. The Cartesian coordinates of the point `p` to be located relative to the
    * triangulation. It is assumed that `p[0]**2 + p[1]**2 + p[2]**2 = 1`, that is, that the
      point lies on the unit sphere.
    * `ist` - Input. The index of the node at which the search is to begin. The search time
    * depends on the proximity of this node to `p`. If no good candidate is known, any value
      between `1` and `n` will do.
    * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    * `x[n]`, `y[n]`, `z[n]` - Input. The Cartesian coordinates of the nodes.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input. The data structure defining
      the triangulation, created by [trmesh].
    * `al` - Output. The arc length between `p` and node [nearnd]. Because both points are on
      the unit sphere, this is also the angular separation in radians.

    # Returns
    The index of the nearest node to `p`. Will be `0` if `n < 3` or the triangulation data structure is invalid.
    */
    #[link_name = "nearnd_"]
    pub fn nearnd(
        p: *const c_double,
        ist: *const c_int,
        n: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        list: *const c_int,
        lptr: *const c_int,
        lend: *const c_int,
        al: *mut c_double,
    ) -> c_int;

    /**
    Optimizes the quadrilateral portion of a triangulation.

    Given a set of `na` triangulation arcs, this subroutine optimizes the portion of the
    triangulation consisting of the quadrilaterals (pairs of adjacent triangles) which have the arcs
    as diagonals by applying the circumcircle test and appropriate swaps to the arcs.

    An iteration consists of applying the swap test and swaps to all `na` arcs in the order in which
    they are stored. The iteration is repeated until no swap occurs or `nit` iterations have
    been performed. The bound on the number of iterations may be necessary to prevent an infinite
    loop caused by cycling (reversing the effect of a previous swap) due to floating point
    inaccuracy when four or more nodes are nearly cocircular.

    # Arguments

    * `x[*]`, `y[*]`, `z[*]` - Input. The nodal coordinates.
    * `na` - Input. The number of arcs in the set. `na >= 0`.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input/output. The data structure
      defining the triangulation, created by [trmesh]. On output, updated to reflect the swaps.
    * `nit` - Input/output. On input, the maximum number of iterations to be performed. `nit = 4 *
    na` should be sufficient. `nit >= 1`. On output, the number of iterations performed.
    * `iwk[2][na]` - Input/output. The nodal indexes of the arc endpoints (paris of endpoints are
      stored in columns). On output, endpoint indexes of the new set of arcs reflecting the swaps.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered.
      * `1`, if a swap occurred on the last of `maxit` iterations, where `maxit` is the value of
        `nit` on input. The new set of arcs is not necessarily optimal in this case.
      * `2`, if `na < 0` or `nit < 1` on input
      * `3`, if `iwk[2][i]` is not a neighbor of `iwk[1][i]` for some `i` in the range `1` to `na`. A
        swap may have occurred in this case.
      * `4`, if a zero pointer was returned by subroutine [swap].
    */
    #[link_name = "optim_"]
    pub fn optim(
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        na: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        nit: *mut c_int,
        iwk: *mut c_int,
        ier: *mut c_int,
    );

    /**
    Normalizes each R83 in an R83Vec to have unit norm.

    # Arguments

    * `n` - Input. The number of nodes in the triangulation.
    * `x[n]`, `y[n]`, `z[n]` - Input/output. On input, the unnormalized coordinates of the
      triangulation. On output, the normalized coordinates of the triangulation (unit vectors).
     */
    #[link_name = "r83vec_normalize_"]
    pub fn r83vec_normalize(n: *const c_int, x: *mut c_double, y: *mut c_double, z: *mut c_double);

    /**
    Converts from Cartesian to spherical coordinates (latitude, longitude, radius).

    # Arguments
    * `px`, `py`, `pz` - Input. The coordinates of `p`
    * `plat` - Output. The latitude of `p` in the range `-PI/2` to `PI/2`, or `0` if `pnrm = 0`.
    * `plon` - Output. The longitude of `p` in the range `-PI` to `PI`, or `0` if `p` lies on the
      Z-axis.
    * `pnrm` - Output. The magnitude (Euclidean norm) of `p`.

    # Safety
    - All pointers must be valid and properly aligned
    */
    #[link_name = "scoord_"]
    pub fn scoord(
        px: *const c_double,
        py: *const c_double,
        pz: *const c_double,
        plat: *mut c_double,
        plon: *mut c_double,
        pnrm: *mut c_double,
    );

    /**
    Replaces the diagonal arc of a quadrilateral with the other diagonal.

    Given a triangulation of a set of points on the unit sphere, this subroutine replaces a diagonal
    arc in a strictly convex quadrilateral (defined by a pair of adjacent triangles) with the other
    diagonal. Equivalently, a pair of adjacent triangles is replaced by another pair having the same
    union.

    # Arguments

    * `in1`, `in2`, `io1`, `io2` - Input. Nodal indexes of the vertices of the quadrilateral. `io1 -
    io2` is replaced by `in1 - in2`. (`io1`, `io2`, `in1`) and (`io2`, `io1`, `in2`) must be
      triangles on input.

    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input/output. The data structure
      defining the triangulation, created by [trmesh]. On output, updated with the swap; triangles
      (`io1`, `io2`, `in1`) and (`io2`, `io1`, `in2`) are replaced by (`in1`, `in2`, `io2`) and (`in2`, `in1`, `io1`) unless `lp21 = 0`.
    * `lp21` - Output. Index of `in1` as a neighbor of `in2` after the swap is performed unless `in1` and `in2` are adjacent on input, in which case `lp21 = 0`.
    */
    #[link_name = "swap_"]
    pub fn swap(
        in1: *const c_int,
        in2: *const c_int,
        io1: *const c_int,
        io2: *const c_int,
        list: *mut c_int,
        lptr: *mut c_int,
        lend: *mut c_int,
        lp21: *mut c_int,
    );

    /**
    Decides whether to replace a diagonal arc by the other in a quadrilateral.

    The decision will be to swap (`swptst = true`) if and only if `n4` lies above the plane (in the half-space
    not containing the origin) defined by (`n1`, `n2`, `n3`), or equivalently, if the projection of
    `n4` onto this plane is interior to the circumcircle of (`n1`, `n2`, `n3`). The decision will be
    for no swap if the quadrilateral is not strictly convex.

    # Arguments

    * `n1`, `n2`, `n3`, `n4` - Input. The indexes of the four nodes defining the quadrilateral with
      `n1` adjacent to `n2`, and (`n1`, `n2`, `n3`) in counterclockwise order. The arc connecting `n1`
      to `n2` should be replaced by an arc connection `n3` to `n4` if `swptst = true`. Refer to
      subroutine [swap].
    * `x[n]`, `y[n]`, `z[n]` - Input. The coordinates of the nodes.

    # Returns

    True if and only if the arc connecting `n1` and `n2` should be swapped for an arc connecting `n3` and `n4`.
    */
    #[link_name = "swptst_"]
    pub fn swptst(
        n1: *const c_int,
        n2: *const c_int,
        n3: *const c_int,
        n4: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
    ) -> bool;

    /**
    Transform spherical coordinates into Cartesian coordinates
    on the unit sphere for input to [trmesh].

    Storage for X and Y may coincide with storage for `rlat` and `rlon` if the latter need not be saved.

    # Arguments
    * `n` - Input. The number of nodes (points on the unit sphere) whose coordinates are to be transformed
    * `rlat` - Input. The latitudes of the nodes in radians.
    * `rlon` - Input. The longitudes of the nodes in radians.
    * `x`, `y`, and `z` - Output. The coordinates in the range `-1` to `1`. `x[i]**2 + y[i]**2 + z[i]**2 = 1` for `i` = 1 to `n`
    # Safety
    - All pointers must be valid and properly aligned
    - Arrays `rlat`, `rlon`, `x`, `y`, `z` must have length >= `*n`
    - Arrays `x`, `y`, `z` must not overlap with `rlat`, `rlon`
    - `n` must be > 0
    */
    #[link_name = "trans_"]
    pub fn trans(
        n: *const c_int,
        rlat: *const c_double,
        rlon: *const c_double,
        x: *mut c_double,
        y: *mut c_double,
        z: *mut c_double,
    );

    /**
    Locates a point relative to a triangulation.

    This subroutine locates a point `p` relative to a triangulation created by [trmesh]. If `p` is contained
    in a triangle, the three vertex indexes and barycentric coordinates are returned. Otherwise, the
    indexes of the visible boundary nodes are returned.

    # Arguments
    * `nst` - Input. The index of a node at which [trfind] begins its search. Search time depends on
      the proximity of this node to `p`.
    * `p[3]` - Input. The x, y, and z coordinates (in that order) of the point `p` to be located.
    * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    * `x[n]`, `y[n]`, `z[n]`, the coordinates of the triangulation nodes (unit vectors).
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input. The data structure defining the
      triangulation, created by [trmesh].
      `b1`, `b2`, `b3` - Output. The unnormalized barycentric coordinates of the central projection of
      `p` onto the underlying planar triangle if `p` is in the convex hull of the nodes. These
      parameters are not altered if `i1 = 0`.
    * `i1`, `i2`, `i3` - Output. The counterclockwise-ordered vertex indexes of a triangle
      containing `p` if `p` is contained in a triangle. If `p` is not in the convex hull of
      the nodes, `i1` and `i2` are the rightmost and leftmost (boundary) nodes that are visible from
      `p`, and `i3 = 0`. (If all boundary nodes are visible from `p`, then `i1` and `i2` coincide.)
      `i1 = i2 = i3 = 0` if `p` and all of the nodes are coplanar (lie on a common great circle).
    */
    #[link_name = "trfind_"]
    pub fn trfind(
        nst: *const c_int,
        p: *const c_double,
        n: *const c_int,
        x: *const c_double,
        y: *const c_double,
        z: *const c_double,
        list: *const c_int,
        lptr: *const c_int,
        lend: *const c_int,
        b1: *mut c_double,
        b2: *mut c_double,
        b3: *mut c_double,
        i1: *mut c_int,
        i2: *mut c_int,
        i3: *mut c_int,
    );

    /**
    Convert a triangulation data structure into a triangle list.

    # Arguments
    * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
    * `list`, `lptr`, `lend` - Input. Linked list data structure defining the triangulation.
    * `nrow` - Input. The number of rows (entries per triangle) reserved for the triangle list `ltri`. The
      value must be 6 if only the vertex indexes and neighboring triangle indexes are to be stored, or
      if arc indexes are also to be assigned and stored. Refer to `ltri`.
    * `nt` - Output. The number of triangles in the triangulation unless `ier != 0`, which
      in case `nt = 0`. `nt = 2n - nb - 2` if `nb >= 3` or `2n-4` if `nb = 0`, where `nb` is the
      number of boundary nodes.
    * `ltri` - Output. The second dimension of `ltri` must be at least `nt`, where `nt` will be at
      most `2*n - 4`. The `j`-th column contains the vertex nodal indexes (first three rows),
      neighboring triangle indxes (second three rows), and, if `nrow = 9`, arc indexes (last three
      rows) associated with triangle `j` for `j = 1, ..., nt`. The vertices are ordered counterclockwise
      with the first vertex taken to be the one with smallest index. Thus `ltri[2][j]` and `ltri[3][j]` are larger than `ltri[1, j]` and index adjacent neighbors of node `ltri[1][j]`. For `i = 1, 2, 3`, `ltri[i + 3][j]` and `ltri[i + 6][j]` index the triangle and arc, respectively, which are opposite (not shared by) node `ltri[i][j]`, with `ltri[i + 3, j] = 0` if `ltri[i + 6][j]` indexes a boundary arc. Vertex indexes range from `1` to `n`, triangle indexes from `0` to `nt`, and, if included, arc indexes from `1` to `na`, where `na = 3n - nb - 3` if `nb >= 3` or `3n - 6` if `nb = 0`. The triangles are ordered on first (smallest) vertex indexes.
    * ier - Output. Error indicator.
      * `0`, if no errors were encountered.
      * `1`, if `n` or `nrow` is outside its valid range on input.
      * `2`, if the triangulation data structure (`list`, `lptr`, `lend`) is invalid. Note, however, that
        these arrays are not completely tested for validity.
    */
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

    /**
    Converts a triangulation data structure to a triangle list.

    This subroutine converts a triangulation data structure from the linked list created by [trmesh] to a triangle list.

    It is a version of [trlist] for the special case where the triangle list should only include the nodes that define each triangle.

    # Arguments

    * `n` - Input. The number of nodes in the triangulation.
    * `list[6 * (n - 2)]`, `lptr[6 * (n - 2)]`, `lend[n]` - Input. The linked list data
      structure defining the triangulation. Refer to [trmesh].
    * `nt` - Output. The number of triangles in the triangulation unless `ier != 0`, in which
      case `nt = 0`. `nt = 2n - nb - 2` if `nb >= 3` or `2n - 4` if `nb = 0`, where `nb` is the number of boundary nodes.
    * `ltri[3][*]` - Output. The second dimension of `ltri` must be at least `nt`, where `nt`
      will be at most `2 * n - 4`. The `j`-th column contains the vertex nodal indexes
      associated with triangle `j` for `j = 1, ..., nt`. The vertices are ordered counterclockwise with the first vertex taken to be the one with the smallest index. Thus, `ltri[2][j]` and `ltri[3][j]` are larger than `ltri[1][j]` and index adjacent neighbors of node `ltri[1][j]`. The triangles are ordered on first (smallest) vertex indexes.
    * `ier` - Output. Error indicator:
      * `0`, if no errors were encountered.
      * `1`, if `n` is outside its valid range on input.
      * `2`, if the triangulation data structure (`list`, `lptr`, `lend`) is invalid. Note, however, that these arrays are not completely tested for validity.
    */
    #[link_name = "trlist2_"]
    pub fn trlist2(
        n: *const c_int,
        list: *const c_int,
        lptr: *const c_int,
        lend: *const c_int,
        nt: *mut c_int,
        ltri: *mut c_int,
        ier: *mut c_int,
    );

    /**
     Creates a Delaunay triangulation on the unit sphere.

     The Delaunay triangulation is defined as a set of (spherical) triangles with the following five
     properties:
     1. The triangle vertices are nodes.
     2. No triangle contains a node other than its vertices.
     3. The interiors of the triangles are pairwise disjoint.
     4. The union of triangles is the convex hull of the set of nodes (the smallest convex set that
        contains the nodes). If the nodes are not contained in a single hemisphere, their convex hull
        is the entire sphere and there are no boundary nodes. Otherwise, there are at least three
        boundary nodes.
     5. The interior of the circumcircle of each triangle contains no node.

     The first four properties define a triangulation, and the last property results in a
     triangulation which is as close as possible to equiangular in a certain sense and which is
     uniquely defined unless four or more nodes lie in a common plane. This property makes the
     triangulation well-suited for solving closest-point problems and for triangle based
     interpolation.

     Provided the nodes are randomly ordered, the algorithm has expected time complexity O(N*log(N))
     for most nodal distributions. Note, however, that the complexity may be as high as O(N**2) if,
     for example, the nodes are ordered on increasing latitude.

     Spherical coordinates (latitude and longitude) may be converted to Cartesian coordinates by
     [trans].

     The following is a list of the software package modules which a user may wish to call directly:
     * [addnod] - Updates the triangulation by appending a new node.
     * [areas] - Returns the area of a spherical triangle
     * [bnodes] - Returns an array containing the index of the boundary nodes (if any) in
       counterclockwise order. Counts of boundary nodes, triangles, and arcs are also returned.
     * [circum] - Returns the circumcenter of a spherical triangle.
     * [crlist] - Returns the set of triangle circumcenters (Voronoi vertices) and circumradii
       associated with a triangulation.
     * [delarc] - Deletes a boundary arc from a triangulation.
     * [edge] - Forces an arbitrary pair of nodes to be connected by an arc in the triangulation.
     * [getnp] - Determines the ordered sequence of L closest nodes to a given node, along with the
       associated distances.
     * [inside] - Locates a point relative to a polygon on the surface of the sphere.
     * [intrsc] - Returns the point of intersection between a pair of great circle arcs.
     * [jrand] - Generates a uniformly distributed pseudo-random integer.
     * [left] - Locates a point relative to a great circle
     * [nearnd] - Returns the index of the nearest node to an arbitrary point, along with its squared
       distance.
     * [scoord] - Converts a point from Cartesian coordinates to spherical coordinates.
     * [trans] - Transforms spherical coordinates into Cartesian coordinates on the unit sphere for
       input to [trmesh]
     * [trlist] - Converts the triangulation data structure to a triangle list more suitable for use in
       a finite element code.

     # Arguments
     * `n` - Input. The number of nodes in the triangulation. `3 <= n`.
     * `x[n]`, `y[n]`, `z[n]` - Input. The coordinates of distinct nodes. `(x[k], y[k], z[k])` is referred to as node `k`, and `k` is referred to as a nodal index. It is required that `x[k]**2 + y[k]**2 + z[k]**2 = 1` for all `k`. The first three nodes must not be collinear (lie on a common great circle).
     * `list` - Output. `6 * (n - 2)` nodal indexes which, along with `lptr`, `lend`, and `lnew`
       define the triangulation as a set of `n` adjacency lists; counterclockwise-ordered sequences of
       neighboring nodes such that the first and last neighbors of a boundary node are boundary nodes
       (the first neighbor of an interior node is arbitrary). In order to distinguish between interior
       and boundary nodes, the last neighbor of each boundary node is represented by the negative of
       its index.
    * `lptr` - Output. Set of pointers (`list` indexes) in one-to-one correspondence with the
      elements of `list`. `list[lptr[i]]` indexes the node which follows `list[i]` in cyclical
      counterclockwise order (the first neighbor follows the last neighbor).
    * `lend` - Output. `n` pointers to adjacency lists. `lend[k]` points to the last neighbor of
      node `k`. `list[lend[k]] < 0` if and only if `k` is a boundary node.
    * `lnew` - Output. Pointer to the first empty location in `list` and `lptr` (list length plus
      one). `list`, `lptr`, `lend` and `lnew` are not altered if `ier < 0`, and are incomplete if `0 <
     ier`.
    * `near` - Workspace. An array of `n` integers used to efficiently determine the nearest
      triangulation node to each unprocessed node for use by [addnod].
    * `next` - Workspace. An array of `n` integers used to efficiently determine the nearest triangulation node
      to each unprocessed node for use by [addnod].
    * `dist` - Workspace. An array of `n` floats used to efficiently determine the neareast
      triangulation node to each unprocessed node for use by [addnod].
    * `ier` - Output. An integer error indicator:
       * `0`, if no errors were encountered.
       * `-1`, if `n < 3` on input.
       * `-2`, if the first three nodes are collinear.
       * `L`, if nodes L and M coincide for some L < M. The data structure represents a triangulation
         of nodes 1 to M-1 in this case.
     */
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
    use proptest::prelude::*;
    use rstest::{fixture, rstest};

    use super::*;
    use std::{
        collections::HashSet,
        f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4, PI, TAU},
    };

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
    fn hemisphere_fixed() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let x = vec![1.0, 0.0, -1.0, 0.0, 0.0];
        let y = vec![0.0, 1.0, 0.0, -1.0, 0.0];
        let z = vec![0.0, 0.0, 0.0, 0.0, 1.0];

        (x, y, z)
    }

    fn hemisphere(n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let (x, y, z) = fibonacci_sphere(n * 2);

        let mut hx = Vec::with_capacity(n);
        let mut hy = Vec::with_capacity(n);
        let mut hz = Vec::with_capacity(n);

        for i in 0..x.len() {
            if z[i] >= 0.0 && hx.len() < n {
                hx.push(x[i]);
                hy.push(y[i]);
                hz.push(z[i]);
            }
        }

        (hx, hy, hz)
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
        hemisphere_fixed: (Vec<f64>, Vec<f64>, Vec<f64>),
    ) {
        let (x, y, z) = if n == 4 {
            tetrahedron
        } else if n == 5 {
            hemisphere_fixed
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

    fn create_triangulation(
        n: i32,
        x: &[f64],
        y: &[f64],
        z: &[f64],
    ) -> (Vec<i32>, Vec<i32>, Vec<i32>, i32) {
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

        (list, lptr, lend, lnew)
    }

    proptest! {
        #[test]
        fn trfind_locates_any_interior_point(px in -1.0f64..1.0, py in -1.0f64..1.0, pz in -1.0f64..1.0, n in 6..20) {
            let norm = (px * px + py * py + pz * pz).sqrt();
            if norm < f64::EPSILON {
                return Ok(());
            }
            let p = [px / norm, py / norm, pz / norm];
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (list, lptr, lend, _) = create_triangulation(n, &x, &y, &z);

            let nst = 1;
            let mut b1 = 0.0;
            let mut b2 = 0.0;
            let mut b3 = 0.0;
            let mut i1 = 0i32;
            let mut i2 = 0i32;
            let mut i3 = 0i32;

            unsafe {
                trfind(
                    &raw const nst,
                    p.as_ptr(),
                    &raw const n,
                    x.as_ptr(),
                    y.as_ptr(),
                    z.as_ptr(),
                    list.as_ptr(),
                    lptr.as_ptr(),
                    lend.as_ptr(),
                    &raw mut b1,
                    &raw mut b2,
                    &raw mut b3,
                    &raw mut i1,
                    &raw mut i2,
                    &raw mut i3,
                );
            };

            prop_assert!(i3 > 0, "i3 should be > 0");
            prop_assert!(i1 >= 1 && i1 <= n, "i1 should be a valid node index");
            prop_assert!(i2 >= 1 && i2 <= n, "i2 should be a valid node index");
            prop_assert!(i3 >= 1 && i3 <= n, "i3 should be a valid node index");
            prop_assert!(b1 > 0.0 && b2 > 0.0 && b3 > 0.0, "barycentric coords should be positive for an interior point");
        }
    }

    fn unit_vector() -> impl Strategy<Value = [f64; 3]> {
        (-1.0f64..1.0, -1.0f64..1.0, -1.0f64..1.0).prop_map(|(x, y, z)| {
            let norm = (x * x + y * y + z * z).sqrt();
            if norm < f64::EPSILON {
                [1.0, 0.0, 0.0]
            } else {
                [x / norm, y / norm, z / norm]
            }
        })
    }

    fn spherical_distance(a: &[f64], b: &[f64]) -> f64 {
        (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
            .clamp(-1.0, 1.0)
            .acos()
    }

    fn normalize(v: &[f64]) -> [f64; 3] {
        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if norm < f64::EPSILON {
            [1.0, 0.0, 0.0]
        } else {
            [v[0] / norm, v[1] / norm, v[2] / norm]
        }
    }

    proptest! {
        #[test]
        fn circum_equidistant(v1 in unit_vector(), v2 in unit_vector(), v3 in unit_vector()) {
            let mut c = [0.0; 3];
            let mut ier = 0i32;

            unsafe {
                circum(
                    v1.as_ptr(),
                    v2.as_ptr(),
                    v3.as_ptr(),
                    c.as_mut_ptr(),
                    &raw mut ier
                );
            };

            if ier != 0 {
                return Ok(());
            }

            let d1 = spherical_distance(&c, &v1);
            let d2 = spherical_distance(&c, &v2);
            let d3 = spherical_distance(&c, &v3);

            prop_assert!((d1 - d2).abs() < 1e-10, "d1={d1} d2={d2}");
            prop_assert!((d2 - d3).abs() < 1e-10, "d2={d2} d3={d3}");
            prop_assert!((d1 - d3).abs() < 1e-10, "d1={d1} d3={d3}");
        }

        #[test]
        fn circum_on_unit_sphere(v1 in unit_vector(), v2 in unit_vector(), v3 in unit_vector()) {
            let mut c = [0.0; 3];
            let mut ier = 0i32;

            unsafe {
                circum(
                    v1.as_ptr(),
                    v2.as_ptr(),
                    v3.as_ptr(),
                    c.as_mut_ptr(),
                    &raw mut ier
                );
            };

            if ier != 0 {
                return Ok(());
            }

            let norm = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
            prop_assert!((norm - 1.0).abs() < 1e-10, "norm: {norm}");
        }

        #[test]
        fn circum_collinear_error(value in 0f64..TAU) {
            let v1 = [value, 0.0, 0.0];
            let v2 = [value * 2.0, 0.0, 0.0];
            let v3 = [value * 3.0, 0.0, 0.0];

            let mut c = [0.0; 3];
            let mut ier = 0i32;

            unsafe {
                circum(
                    v1.as_ptr(),
                    v2.as_ptr(),
                    v3.as_ptr(),
                    c.as_mut_ptr(),
                    &raw mut ier
                );
            };

            prop_assert!(ier == 1, "expected 1 got {ier}");
        }
    }

    fn add_node_to_triangulation(
        nst: i32,
        k: i32,
        x: &[f64],
        y: &[f64],
        z: &[f64],
        list: &mut [i32],
        lptr: &mut [i32],
        lend: &mut [i32],
        lnew: &mut i32,
    ) -> i32 {
        let mut ier = 0i32;
        unsafe {
            addnod(
                &raw const nst,
                &raw const k,
                x.as_ptr(),
                y.as_ptr(),
                z.as_ptr(),
                list.as_mut_ptr(),
                lptr.as_mut_ptr(),
                lend.as_mut_ptr(),
                &raw mut *lnew,
                &raw mut ier,
            );
        };

        ier
    }

    #[rstest]
    fn test_addnod_adds_fourth_node_to_tetrahedron(tetrahedron: (Vec<f64>, Vec<f64>, Vec<f64>)) {
        let mut x = tetrahedron.0[0..3].to_vec();
        let mut y = tetrahedron.1[0..3].to_vec();
        let mut z = tetrahedron.2[0..3].to_vec();

        let n = 3i32;
        let (mut list, mut lptr, mut lend, mut lnew) = create_triangulation(n, &x, &y, &z);

        x.push(tetrahedron.0[3]);
        y.push(tetrahedron.1[3]);
        z.push(tetrahedron.2[3]);

        let new_list_size = 6 * (n + 1 - 2) as usize;
        list.resize(new_list_size, 0);
        lptr.resize(new_list_size, 0);
        lend.resize(n as usize + 1, 0);

        let k = n + 1;
        let nst = 1i32;

        let add_ier = add_node_to_triangulation(
            nst, k, &x, &y, &z, &mut list, &mut lptr, &mut lend, &mut lnew,
        );

        assert_eq!(add_ier, 0, "addnod should succeed");
        assert!(lnew > 0, "lnew should be positive after adding node");
    }

    proptest! {
        #[test]
        fn test_arc_cosine(c in -10.0f64..10.0f64) {
            let acos = unsafe { arc_cosine(&raw const c) };
            prop_assert!((acos - c.clamp(-1.0, 1.0).acos()).abs() < f64::EPSILON);
        }
    }

    fn scalar_triple_product(n0: &[f64], n1: &[f64], n2: &[f64]) -> f64 {
        let cross_x = n1[1] * n2[2] - n1[2] * n2[1];
        let cross_y = n1[2] * n2[0] - n1[0] * n2[2];
        let cross_z = n1[0] * n2[1] - n1[1] * n2[0];

        n0[0] * cross_x + n0[1] * cross_y + n0[2] * cross_z
    }

    fn is_left(n1: &[f64], n2: &[f64], n0: &[f64]) -> bool {
        unsafe {
            left(
                &raw const n1[0],
                &raw const n1[1],
                &raw const n1[2],
                &raw const n2[0],
                &raw const n2[1],
                &raw const n2[2],
                &raw const n0[0],
                &raw const n0[1],
                &raw const n0[2],
            )
        }
    }

    proptest! {
        #[test]
        fn left_agrees_with_scalar_triple_product(n1 in unit_vector(), n2 in unit_vector(), n0 in unit_vector()) {
            let stp = scalar_triple_product(&n0, &n1, &n2);
            let result = is_left(&n1, &n2, &n0);

            prop_assert_eq!(result, stp >= -f64::EPSILON);
        }

        #[test]
        fn left_boundary_points_return_true(n1 in unit_vector(), n2 in unit_vector(), coeff1 in -1.0f64..1.0, coeff2 in -1.0f64..1.0) {
            let n0 = [
                coeff1 * n1[0] + coeff2 * n2[0],
                coeff1 * n1[1] + coeff2 * n2[1],
                coeff1 * n1[2] + coeff2 * n2[2],
            ];

            let stp = scalar_triple_product(&n0, &n1, &n2);
            let result = is_left(&n1, &n2, &n0);
            prop_assert_eq!(result, stp >= 0.0);
        }
    }

    fn should_swap(n1: i32, n2: i32, n3: i32, n4: i32, x: &[f64], y: &[f64], z: &[f64]) -> bool {
        unsafe {
            swptst(
                &raw const n1,
                &raw const n2,
                &raw const n3,
                &raw const n4,
                x.as_ptr(),
                y.as_ptr(),
                z.as_ptr(),
            )
        }
    }

    fn swap_determinant(n1: &[f64], n2: &[f64], n3: &[f64], n4: &[f64]) -> f64 {
        let dx1 = n1[0] - n4[0];
        let dy1 = n1[1] - n4[1];
        let dz1 = n1[2] - n4[2];

        let dx2 = n2[0] - n4[0];
        let dy2 = n2[1] - n4[1];
        let dz2 = n2[2] - n4[2];

        let dx3 = n3[0] - n4[0];
        let dy3 = n3[1] - n4[1];
        let dz3 = n3[2] - n4[2];

        dx3 * (dy2 * dz1 - dy1 * dz2) - dy3 * (dx2 * dz1 - dx1 * dz2)
            + dz3 * (dx2 * dy1 - dx1 * dy2)
    }

    proptest! {
        #[test]
        fn test_swptst(
            n1 in unit_vector(),
            n2 in unit_vector(),
            n3 in unit_vector(),
            n4 in unit_vector()
        ) {
            let x = vec![n1[0], n2[0], n3[0], n4[0]];
            let y = vec![n1[1], n2[1], n3[1], n4[1]];
            let z = vec![n1[2], n2[2], n3[2], n4[2]];

            let det = swap_determinant(&n1, &n2, &n3, &n4);
            let result = should_swap(1, 2, 3, 4, &x, &y, &z);

            prop_assert_eq!(result, det > 0.0);
        }
    }

    fn find_node_pointer(lpl: i32, nb: i32, list: &[i32], lptr: &[i32]) -> i32 {
        unsafe { lstptr(&raw const lpl, &raw const nb, list.as_ptr(), lptr.as_ptr()) }
    }

    fn check_triangulation(n: i32, list: &[i32], lend: &[i32], lptr: &[i32]) {
        for node in 1..=n {
            let lpl = lend[(node - 1) as usize];
            assert!(lpl > 0, "Node {node} should have valid lend after edge");

            let mut current = lpl;
            let mut count = 0;
            loop {
                let neighbor = list[(current - 1) as usize].abs();
                assert!(
                    neighbor >= 1 && neighbor <= n,
                    "Node {node} neighbor {neighbor} out of range"
                );
                count += 1;
                current = lptr[(current - 1) as usize];
                if current == lpl {
                    break;
                }
                assert!(count <= 6 * (n - 2), "Infinite loop in adjacency list");
            }
        }
    }

    proptest! {
        fn lsptr_finds_all_neighbors_in_triangulation(n in 6..50i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (list, lptr, lend, _) = create_triangulation(n, &x, &y, &z);

            for node_idx in 1..=n {
                let lpl = lend[(node_idx - 1) as usize];

                let mut current = lpl;
                let mut neighbor_count = 0;

                loop {
                    let neighbor = list[(current - 1) as usize].abs();

                    let found_ptr = find_node_pointer(lpl, neighbor, &list, &lptr);

                    prop_assert!(found_ptr > 0 && found_ptr <= 6 * (n - 2), "Node {node_idx}: lstptr should return valid poiter for neighbor {neighbor}");
                    prop_assert_eq!(list[(found_ptr - 1) as usize].abs(), neighbor, "Node {}: lstptr should find correct neighbor {} at position {}", node_idx, neighbor, found_ptr);

                    let mut traversal_ptr = lpl;
                    let mut found_in_traversal = false;
                    for _ in 0..(6 * (n - 2)) {
                        if traversal_ptr == found_ptr {
                            found_in_traversal = true;
                            break;
                        }
                        traversal_ptr = lptr[(traversal_ptr - 1) as usize];
                        if traversal_ptr == lpl {
                            break;
                        }
                    }

                    prop_assert!(found_in_traversal, "Node {node_idx}: lstptr result {found_ptr} should be reachable from lpl {lpl}");

                    neighbor_count += 1;
                    current = lptr[(current - 1) as usize];

                    if current == lpl {
                        break;
                    }
                }

                prop_assert!(neighbor_count >= 2, "Node {node_idx} should have at least 2 neighbors, found {neighbor_count}");
            }
        }
    }

    proptest! {
         #[test]
         fn test_swap(n in 6..20i32) {
             let (x, y, z) = fibonacci_sphere(n as usize);
             let (mut list, mut lptr, mut lend, _) = create_triangulation(n, &x, &y, &z);

             let nrow = 6i32;
             let max_triangles = (2 * n - 4) as usize;
             let mut ltri = vec![0i32; (nrow as usize) * max_triangles];
             let mut nt = 0i32;
    let mut ier = 0i32;

             unsafe {
                 trlist(
                     &raw const n,
                     list.as_ptr(),
                     lptr.as_ptr(),
                     lend.as_ptr(),
                     &raw const nrow,
                     &raw mut nt,
                     ltri.as_mut_ptr(),
                     &raw mut ier,
                 );
             }

             prop_assert_eq!(ier, 0, "trlist failed");
             prop_assert!(nt >= 2, "Need at least 2 triangles for a swap");


             'outer: for t1 in 0..(nt as usize) {
                 for t2 in (t1 + 1)..(nt as usize) {
                     let v1_0 = ltri[t1 * 6 + 0];
                     let v1_1 = ltri[t1 * 6 + 1];
                     let v1_2 = ltri[t1 * 6 + 2];

                     let v2_0 = ltri[t2 * 6 + 0];
                     let v2_1 = ltri[t2 * 6 + 1];
                     let v2_2 = ltri[t2 * 6 + 2];

                     let verts1 = [v1_0, v1_1, v1_2];
                     let verts2 = [v2_0, v2_1, v2_2];

                     let mut shared = vec![];
                     let mut unique1 = None;
                     let mut unique2 = None;

                     for &v in &verts1 {
                         if verts2.contains(&v) {
                             shared.push(v);
                         } else {
                             unique1 = Some(v);
                         }
                     }

                     for &v in &verts2 {
                         if !verts1.contains(&v) {
                             unique2 = Some(v);
                         }
                     }

                     if !(shared.len() == 2 && unique1.is_some() && unique2.is_some()) {
                         continue;
                     }

                     let io1 = shared[0];
                     let io2 = shared[1];
                     let in1 = unique1.unwrap();
                     let in2 = unique2.unwrap();

                     let lpl_io1 = lend[(io1 - 1) as usize];
                     let ptr_io2 = find_node_pointer(lpl_io1, io2, &list, &lptr);

                     if !(ptr_io2 > 0 && list[(ptr_io2 - 1) as usize].abs() == io2) {
                         continue;
                     }

                     let lpl_in1 = lend[(in1 - 1) as usize];
                     let ptr_in2_check = find_node_pointer(lpl_in1, in2, &list, &lptr);
                     if list[(ptr_in2_check - 1) as usize].abs() == in2 {
                         continue;
                     }

                     if !should_swap(in1, in2, io1, io2, &x, &y, &z) {
                         continue;
                     }

                     let mut lp21 = 0i32;
                     unsafe {
                         swap(
                             &raw const in1,
                             &raw const in2,
                             &raw const io1,
                             &raw const io2,
                             list.as_mut_ptr(),
                             lptr.as_mut_ptr(),
                             lend.as_mut_ptr(),
                             &raw mut lp21,
                         );
                     };

                     prop_assert!(lp21 > 0, "Swap should succeed with lp21 > 0");

                     let lpl_in2 = lend[(in2 - 1) as usize];
                     let ptr_in1 = find_node_pointer(lpl_in2, in1, &list, &lptr);
                     prop_assert!(ptr_in1 > 0, "IN1 should be neighbor of IN2 after swap");
                     prop_assert_eq!(list[(ptr_in1 - 1) as usize].abs(), in1);

                     break 'outer;
                 }
             }
         }
     }

    proptest! {
        #[test]
        fn test_optim(n in 8..30i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (mut list, mut lptr, mut lend, _) = create_triangulation(n, &x, &y, &z);

            let mut na = 0i32;
            let max_arcs = (3 * n) as usize;
            let mut iwk = vec![0i32; 2 * max_arcs];

            for node in 1..n {
                let lpl = lend[(node - 1) as usize];
                let mut current = lpl;

                loop {
                    let neighbor = list[(current - 1) as usize].abs();
                    if node < neighbor {
                        iwk[na as usize * 2] = node;
                        iwk[na as usize * 2 + 1] = neighbor;
                        na += 1;
                        if na >= 20 {
                            break;
                        }
                    }

                    current = lptr[(current - 1) as usize];
                    if current == lpl {
                        break;
                    }
                }

                if na >= 20 {
                    break;
                }
            }

            prop_assert!(na > 0, "Should have at least one arc to optimize");

            let mut nit = 100i32;
            let mut ier = 0i32;

            unsafe {
                optim(
                    x.as_ptr(),
                    y.as_ptr(),
                    z.as_ptr(),
                    &raw const na,
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    lend.as_mut_ptr(),
                    &raw mut nit,
                    iwk.as_mut_ptr(),
                    &raw mut ier,
                );
            }

            prop_assert!(ier == 0 || ier == 1, "optim failed");

            check_triangulation(n, &list, &lend, &lptr);

            if ier == 0 {
                prop_assert!(nit < 100, "Converged but used all iterations");
            }
        }
    }

    proptest! {
        #[test]
        fn test_edge(n in 8..25i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (mut list, mut lptr, mut lend, _) = create_triangulation(n, &x, &y, &z);

            let mut in1 = 1i32;
            let mut in2 = 2i32;
            let mut found_non_adjacent = false;

            'search: for i in 1..=n {
                let lpl = lend[(i - 1) as usize];
                let mut current = lpl;
                let mut neighbors = vec![];

                loop {
                    let neighbor = list[(current - 1) as usize].abs();
                    neighbors.push(neighbor);
                    current = lptr[(current - 1) as usize];
                    if current == lpl {
                        break;
                    }
                }

                for j in (i + 1)..=n {
                    if !neighbors.contains(&j) && i != j {
                        in1 = i;
                        in2 = j;
                        found_non_adjacent = true;
                        break 'search;
                    }
                }
            }

            if !found_non_adjacent {
                return Ok(());
            }

            let max_lwk = (n - 3) as usize;
            let mut lwk = max_lwk as i32;
            let mut iwk = vec![0i32; 2 * max_lwk];
            let mut ier = 0i32;

            unsafe {
                edge(
                    &raw const in1,
                    &raw const in2,
                    x.as_ptr(),
                    y.as_ptr(),
                    z.as_ptr(),
                    &raw mut lwk,
                    iwk.as_mut_ptr(),
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    lend.as_mut_ptr(),
                    &raw mut ier
                );
            }

            prop_assert!(ier == 0 || ier == 5, "edge failed");

            let lpl_in1 = lend[(in1 - 1) as usize];
            let ptr_in2 = find_node_pointer(lpl_in1, in2, &list, &lptr);
            let is_adjacent = list[(ptr_in2 - 1) as usize].abs() == in2;
            prop_assert!(is_adjacent, "Nodes {in1} and {in2} should be adjacent after edge, lwk={lwk}");

            check_triangulation(n, &list, &lend, &lptr);

            prop_assert!(lwk >= 0 && lwk <= max_lwk as i32, "lwk should be in valid range [0, {max_lwk}], got {lwk}");
        }
    }

    proptest! {
        #[test]
        fn test_insert(n in 6..15i32, insertion_point in 1..3usize) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (mut list, mut lptr, lend, mut lnew) = create_triangulation(n, &x, &y, &z);

            let lpl = lend[0];

            let mut lp = lpl;
            for _ in 0..insertion_point {
                lp = lptr[(lp - 1) as usize];
            }

            let k = n;

            let mut current = lpl;
            let mut already_neighbor = false;
            loop {
                if list[(current - 1) as usize].abs() == k as i32 {
                    already_neighbor = true;
                    break;
                }

                current = lptr[(current - 1) as usize];
                if current == lpl {
                    break;
                }
            }

            if already_neighbor {
                return Ok(());
            }

            let lnew_before = lnew;

            let current_size = list.len();
            list.resize(current_size + 1, 0);
            lptr.resize(current_size + 1, 0);

            unsafe {
                insert(
                    &raw const k,
                    &raw const lp,
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    &raw mut lnew
                );
            }

            prop_assert_eq!(lnew, lnew_before + 1, "lnew should increment by 1");
            prop_assert_eq!(list[(lnew_before - 1) as usize], k, "k should be at position {}", lnew_before);
            let next_after_insert = lptr[(lp - 1) as usize];
            prop_assert_eq!(next_after_insert, lnew_before, "lp should now point ot new entry");

            let next_of_new = lptr[(lnew_before - 1) as usize];
            prop_assert!(next_of_new > 0 && next_of_new <= list.len() as i32, "New entry should point to valid position");
        }
    }

    #[test]
    fn test_bdyadd() {
        let (x, y, z) = hemisphere_fixed();
        let n_base = 5i32;

        let (mut list, mut lptr, mut lend, mut lnew) = create_triangulation(n_base, &x, &y, &z);

        let mut nodes = vec![0i32; n_base as usize];
        let mut nb = 0i32;
        let mut na = 0i32;
        let mut nt = 0i32;

        unsafe {
            bnodes(
                &raw const n_base,
                list.as_ptr(),
                lptr.as_ptr(),
                lend.as_ptr(),
                nodes.as_mut_ptr(),
                &raw mut nb,
                &raw mut na,
                &raw mut nt,
            );
        }

        if nb < 2 || lnew as usize >= list.len() {
            return;
        }

        let i1 = nodes[0];
        let i2 = nodes[1];
        let kk = n_base + 1;

        let new_size = list.len() + 10;
        list.resize(new_size, 0);
        lptr.resize(new_size, 0);
        lend.resize(n_base as usize + 1, 0);

        let lnew_before = lnew;

        unsafe {
            bdyadd(
                &raw const kk,
                &raw const i1,
                &raw const i2,
                list.as_mut_ptr(),
                lptr.as_mut_ptr(),
                lend.as_mut_ptr(),
                &raw mut lnew,
            );
        }

        assert!(lnew > lnew_before, "bdyadd should increment lnew");
        assert!(
            lend[(kk - 1) as usize] > 0,
            "New node should have valid lend"
        );

        let lpl_kk = lend[(kk - 1) as usize];
        let ptr_i1 = find_node_pointer(lpl_kk, i1, &list, &lptr);
        assert!(ptr_i1 > 0, "i1 should be a neighbor of kk");
        assert_eq!(list[(ptr_i1 - 1) as usize].abs(), i1);

        let ptr_i2 = find_node_pointer(lpl_kk, i2, &list, &lptr);
        assert!(ptr_i2 > 0, "i2 should be a neighbor of kk");
        assert_eq!(list[(ptr_i2 - 1) as usize].abs(), i2);

        check_triangulation(kk, &list, &lend, &lptr);
    }

    proptest! {
        #[test]
        fn test_intadd(n in 6..20i32) {
            let (mut x, mut y, mut z) = fibonacci_sphere(n as usize);
            let (mut list, mut lptr, mut lend, mut lnew) = create_triangulation(n, &x, &y, &z);

            let nrow = 6i32;
            let max_triangles = (2 * n - 4) as usize;
            let mut ltri = vec![0i32; (nrow as usize) * max_triangles];
            let mut nt = 0i32;
            let mut ier = 0i32;

            unsafe {
                trlist(
                    &raw const n,
                    list.as_ptr(),
                    lptr.as_ptr(),
                    lend.as_ptr(),
                    &raw const nrow,
                    &raw mut nt,
                    ltri.as_mut_ptr(),
                    &raw mut ier,
                );
            }

            if ier != 0 || nt < 1 {
                return Ok(());
            }

            let i1 = ltri[0];
            let i2 = ltri[1];
            let i3 = ltri[2];

            let x1 = x[(i1 - 1) as usize];
            let y1 = y[(i1 - 1) as usize];
            let z1 = z[(i1 - 1) as usize];
            let x2 = x[(i2 - 1) as usize];
            let y2 = y[(i2 - 1) as usize];
            let z2 = z[(i2 - 1) as usize];
            let x3 = x[(i3 - 1) as usize];
            let y3 = y[(i3 - 1) as usize];
            let z3 = z[(i3 - 1) as usize];

            let px = 0.25 * x1 + 0.25 * x2 + 0.5 * x3;
            let py = 0.25 * y1 + 0.25 * y2 + 0.5 * y3;
            let pz = 0.25 * z1 + 0.25 * z2 + 0.5 * z3;

            let [px, py, pz] = normalize(&[px, py, pz]);
            let p = [px, py, pz];
            let mut b1 = 0.0;
            let mut b2 = 0.0;
            let mut b3 = 0.0;
            let mut tf_i1 = 0i32;
            let mut tf_i2 = 0i32;
            let mut tf_i3 = 0i32;
            let nst = 1i32;

            unsafe {
                trfind(
                    &raw const nst,
                    p.as_ptr(),
                    &raw const n,
                    x.as_ptr(),
                    y.as_ptr(),
                    z.as_ptr(),
                    list.as_ptr(),
                    lptr.as_ptr(),
                    lend.as_ptr(),
                    &raw mut b1,
                    &raw mut b2,
                    &raw mut b3,
                    &raw mut tf_i1,
                    &raw mut tf_i2,
                    &raw mut tf_i3,
                );
            }

            if tf_i3 == 0 {
                return Ok(());
            }

            let kk = n + 1;
            x.push(px);
            y.push(py);
            z.push(pz);

            let new_list_size = list.len() + 6;
            list.resize(new_list_size, 0);
            lptr.resize(new_list_size, 0);
            lend.resize(kk as usize, 0);

            let lnew_before = lnew;

            unsafe {
                intadd(
                    &raw const kk,
                    &raw const tf_i1,
                    &raw const tf_i2,
                    &raw const tf_i3,
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    lend.as_mut_ptr(),
                    &raw mut lnew,
                );
            }

            prop_assert_eq!(lnew, lnew_before + 6, "intadd should increment lnew by 3");
            prop_assert!(lend[(kk - 1) as usize] > 0, "New node should have valid lend");

            for &vertex in &[tf_i1, tf_i2, tf_i3] {
                let lpl_kk = lend[(kk - 1) as usize];
                let ptr = find_node_pointer(lpl_kk, vertex, &list, &lptr);
                prop_assert!(ptr > 0, "Vertex {vertex} should be a neighbor of kk");
                prop_assert_eq!(list[(ptr - 1) as usize].abs(), vertex, "kk should have vertex {} as a neighbor", vertex);
            }

            for &vertex in &[tf_i1, tf_i2, tf_i3] {
                let lpl_vertex = lend[(vertex - 1) as usize];
                let ptr = find_node_pointer(lpl_vertex, kk, &list, &lptr);
                prop_assert!(ptr > 0, "kk should be a neighbor of vertex {}", vertex);
                prop_assert_eq!(list[(ptr - 1) as usize].abs(), kk, "Vertex {} should have kk as a neighbor", vertex);
            }

            check_triangulation(kk, &list, &lend, &lptr);
        }
    }

    proptest! {
        #[test]
        fn test_covsph(n in 5..25i32) {
            let (x, y, z) = hemisphere(n as usize);
            let n_hemi = x.len() as i32;

            if n_hemi < 4 {
                return Ok(());
            }

            let (mut list, mut lptr, mut lend, mut lnew) = create_triangulation(n_hemi, &x, &y, &z);

            let mut nodes = vec![0i32; n_hemi as usize];
            let mut nb = 0i32;
            let mut na = 0i32;
            let mut nt = 0i32;

            unsafe {
                bnodes(
                    &raw const n_hemi,
                    list.as_ptr(),
                    lptr.as_ptr(),
                    lend.as_ptr(),
                    nodes.as_mut_ptr(),
                    &raw mut nb,
                    &raw mut na,
                    &raw mut nt,
                );
            }

            if nb < 3 {
                return Ok(());
            }

            let new_size = list.len() + (2 * nb as usize);
            list.resize(new_size, 0);
            lptr.resize(new_size, 0);
            lend.resize(n_hemi as usize + 1, 0);

            let n0 = nodes[0];
            let kk = n_hemi + 1;

            unsafe {
                covsph(
                    &raw const kk,
                    &raw const n0,
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    lend.as_mut_ptr(),
                    &raw mut lnew,
                );
            }

            for i in 0..nb as usize {
                let lpl = lend[(nodes[i] - 1) as usize];
                prop_assert!(list[(lpl - 1) as usize] > 0);
            }

            let lpl_kk = lend[(kk - 1) as usize];
            for i in 0..nb as usize {
                let ptr = find_node_pointer(lpl_kk, nodes[i], &list, &lptr);
                prop_assert!(ptr > 0);
            }

            check_triangulation(kk, &list, &lend, &lptr);
        }
    }

    proptest! {
        #[test]
        fn test_delnod(mut n in 5..25i32) {
            let (mut x, mut y, mut z) = fibonacci_sphere(n as usize);
            let (mut list, mut lptr, mut lend, mut lnew) = create_triangulation(n, &x, &y, &z);

            let k = n / 2;
            let starting_n = n;

            let lwk_max = n;
            let mut lwk = lwk_max;
            let mut iwk = vec![0i32; 2 * lwk_max as usize];
            let mut ier = 0i32;

            unsafe {
                delnod(
                    &raw const k,
                    &raw mut n,
                    x.as_mut_ptr(),
                    y.as_mut_ptr(),
                    z.as_mut_ptr(),
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    lend.as_mut_ptr(),
                    &raw mut lnew,
                    &raw mut lwk,
                    iwk.as_mut_ptr(),
                    &raw mut ier,
                );
            }

            prop_assert!(ier == 0 || ier == 6, "delnod failed");
            prop_assert_eq!(n, starting_n - 1, "n should be decremented by 1");

            check_triangulation(n, &list, &lend, &lptr);
        }
    }

    proptest! {
        #[test]
        fn test_jrand(n in 100..150, mut ix in 1..30000, mut iy in 1..30000, mut iz in 1..30000) {
            let x = unsafe { jrand(&raw const n, &raw mut ix, &raw mut iy, &raw mut iz) };
            prop_assert!(x <= n);
        }
    }

    proptest! {
        #[test]
        fn test_nearnd(n in 5..25i32, p in unit_vector()) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (list, lptr, lend, _) = create_triangulation(n, &x, &y, &z);

            let mut al = 0.0;
            let ist = 1;
            let nearest = unsafe {
                nearnd(
                    p.as_ptr(),
                    &raw const ist,
                    &raw const n,
                    x.as_ptr(),
                    y.as_ptr(),
                    z.as_ptr(),
                    list.as_ptr(),
                    lptr.as_ptr(),
                    lend.as_ptr(),
                    &raw mut al,
                )
            };

            let mut naive_nearest = 0;
            let mut naive_nearest_distance = f64::MAX;
            for i in 0..x.len() {
                let dist = ((p[0] - x[i]).powi(2) + (p[1] - y[i]).powi(2) + (p[2] - z[i]).powi(2)).sqrt();
                if dist < naive_nearest_distance {
                    naive_nearest = i;
                }
                naive_nearest_distance = naive_nearest_distance.min(dist);
            }
            prop_assert_eq!(nearest, naive_nearest as i32 + 1);
            prop_assert!((naive_nearest_distance - al).abs() < 0.10, "Incorrect al value: {al}, expected: {naive_nearest_distance}");

        }
    }

    proptest! {
        #[test]
        fn test_getnp(n in 10..25i32, l in 3..8i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (list, lptr, mut lend, _) = create_triangulation(n, &x, &y, &z);
            let mut npts = vec![0i32; l as usize];
            npts[0] = n / 2;
            let mut distances = vec![0.0f64; l as usize];
            distances[0] = -1.0;

            for i in 2..=l {
                let current_l = i;
                let mut df = 0.0f64;
                let mut ier = 0i32;

                unsafe {
                    getnp(
                        x.as_ptr(),
                        y.as_ptr(),
                        z.as_ptr(),
                        list.as_ptr(),
                        lptr.as_ptr(),
                        lend.as_mut_ptr(),
                        &raw const current_l,
                        npts.as_mut_ptr(),
                        &raw mut df,
                        &raw mut ier,
                    );
                }

                prop_assert_eq!(ier, 0, "getnp failed");
                prop_assert!(df >= -1.0 && df <= 1.0, "df: {df}");
                distances[(i - 1) as usize] = df;
                prop_assert!(df >= distances[(i - 2) as usize], "distance should be increasing: df={df}, prev={}", distances[(i - 2) as usize]);
            }

            prop_assert_eq!(npts.into_iter().collect::<HashSet<_>>().len(), l as usize, "All {} nodes in sequence should be distinct", l);
        }
    }

    proptest! {
        #[test]
        fn test_nbcnt(n in 6..30i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (_, lptr, lend, _) = create_triangulation(n, &x, &y, &z);

            for node_idx in 1..=n {
                let lpl = lend[(node_idx - 1) as usize];
                let count = unsafe { nbcnt(&raw const lpl, lptr.as_ptr()) };
                let mut manual_count = 0;
                let mut current = lpl;
                loop {
                    manual_count += 1;
                    current = lptr[(current - 1) as usize];
                    if current == lpl {
                        break;
                    }
                    prop_assert!(manual_count <= 6 * (n - 2), "infinite loop detected");
                }
                prop_assert_eq!(count, manual_count, "nbcnt shoudl match manula count for node {}", node_idx);
            }
        }
    }

    proptest! {
        #[test]
        fn test_delnb(n in 6..25i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (mut list, mut lptr, mut lend, mut lnew) = create_triangulation(n, &x, &y, &z);

            let n0 = n / 2;
            let lpl_n0 = lend[(n0 - 1) as usize];

            let first_neighbor_pos = lptr[(lpl_n0 - 1) as usize];
            let nb = list[(first_neighbor_pos - 1) as usize];

            if nb < 1 || nb > n {
                return Ok(());
            }

            let lnew_before = lnew;
            let mut lph = 0i32;

            unsafe {
                delnb(
                    &raw const n0,
                    &raw const nb,
                    &raw const n,
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    lend.as_mut_ptr(),
                    &raw mut lnew,
                    &raw mut lph,
                );
            }

            prop_assert!(lph > 0, "delnb failed");
            prop_assert!(lnew <= lnew_before, "lnew should decrease or stay the same, got {lnew} expected <= {lnew_before}");

            check_triangulation(n, &list, &lend, &lptr);
        }
    }

    proptest! {
        #[test]
        fn test_delarc(n in 5..25i32) {
            let (x, y, z) = hemisphere(n as usize);
            let (mut list, mut lptr, mut lend, mut lnew) = create_triangulation(n, &x, &y, &z);

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
            }

            if nb < 2 {
                return Ok(());
            }

            let io1 = nodes[0];
            let io2 = nodes[1];

            let lpl_io1 = lend[(io1 - 1) as usize];
            let ptr_io2 = find_node_pointer(lpl_io1, io2, &list, &lptr);
            let are_neighbors = list[(ptr_io2 - 1) as usize].abs() == io2;

            if !are_neighbors {
                return Ok(());
            }

            let mut ier = 0i32;

            unsafe {
                delarc(
                    &raw const n,
                    &raw const io1,
                    &raw const io2,
                    list.as_mut_ptr(),
                    lptr.as_mut_ptr(),
                    lend.as_mut_ptr(),
                    &raw mut lnew,
                    &raw mut ier,
                );
            }

            prop_assert_eq!(ier, 0, "delarc failed");

            let lpl_io1_new = lend[(io1 - 1) as usize];
            let ptr_io2_new = find_node_pointer(lpl_io1_new, io2, &list, &lptr);
            let still_neighbors = list[(ptr_io2_new - 1) as usize].abs() == io2;
            prop_assert!(!still_neighbors, "io1 and io2 should no longer be neighbors");

            check_triangulation(n, &list, &lend, &lptr);
        }
    }

    proptest! {
        #[test]
        fn test_inside(n in 6..25i32, p in unit_vector()) {
            let (x, y, z) = fibonacci_sphere(n as usize);

            let nv = 4i32;
            let listv = vec![1i32, 2, 3, 4];

            let mut ier = 0i32;

            let _ = unsafe {
                inside(
                    p.as_ptr(),
                    &raw const n,
                    x.as_ptr(),
                    y.as_ptr(),
                    z.as_ptr(),
                    &raw const nv,
                    listv.as_ptr(),
                    &raw mut ier,
                )
            };

            prop_assert_eq!(ier, 0, "inside failed");
        }
    }

    proptest! {
        #[test]
        fn test_crlist(n in 5i32..25i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (list, lptr, lend, lnew) = create_triangulation(n, &x, &y, &z);

            if lnew <= 0 {
                return Ok(());
            }

            let nrow = 6i32;
            let max_triangles = (2 * n - 4) as usize;
            let mut ltri = vec![0i32; (nrow as usize) * max_triangles];


            let ncol = n;
            let max_centers = (2 * n - 4) as usize;
            let mut listc = vec![0i32; 3 * max_triangles];
            let mut nb = 0i32;
            let mut xc = vec![0.0f64; max_centers];
            let mut yc = vec![0.0f64; max_centers];
            let mut zc = vec![0.0f64; max_centers];
            let mut rc = vec![0.0f64; max_centers];
            let mut ier_crlist = 0i32;
            let mut lptr_copy = lptr.clone();
            let mut lnew_copy = lnew;

            unsafe {
                crlist(
                    &raw const n,
                    &raw const ncol,
                    x.as_ptr(),
                    y.as_ptr(),
                    z.as_ptr(),
                    list.as_ptr(),
                    lend.as_ptr(),
                    lptr_copy.as_mut_ptr(),
                    &raw mut lnew_copy,
                    ltri.as_mut_ptr(),
                    listc.as_mut_ptr(),
                    &raw mut nb,
                    xc.as_mut_ptr(),
                    yc.as_mut_ptr(),
                    zc.as_mut_ptr(),
                    rc.as_mut_ptr(),
                    &raw mut ier_crlist,
                );
            };

            if nb > 0 {
                return Ok(());
            }

            prop_assert_eq!(ier_crlist, 0, "crlist failed");

            let nt = 2 * n - 4;
            for i in 0..nt as usize {
                let norm_sq = xc[i] * xc[i] + yc[i] * yc[i] + zc[i] * zc[i];
                prop_assert!((norm_sq - 1.0).abs() < 1e-10, "circumcenter {i} should be a unit vector, norm_sq={norm_sq}");

                let base = i * 6;
                let v1 = (ltri[base] - 1) as usize;
                let v2 = (ltri[base + 1] - 1) as usize;
                let v3 = (ltri[base + 2] - 1) as usize;

                if v1 >= n as usize || v2 >= n as usize || v3 >= n as usize {
                    continue;
                }

                let c = [xc[i], yc[i], zc[i]];
                let p1 = [x[v1], y[v1], z[v1]];
                let p2 = [x[v2], y[v2], z[v2]];
                let p3 = [x[v3], y[v3], z[v3]];

                let d1 = spherical_distance(&c, &p1);
                let d2 = spherical_distance(&c, &p2);
                let d3 = spherical_distance(&c, &p3);

                prop_assert!((d1 - d2).abs() < 1e-9 && (d2 - d3).abs() < 1e-9, "circumcenter {i} not equidistant: d1={d1}, d2={d2}, d3={d3}, rc={}", rc[i]);
                prop_assert!((rc[i] - d1).abs() < 1e-9, "circumradius mismatch for triangle {i}: rc={}, actual distance={d1}", rc[i]);
            }
        }
    }

    fn cross_product(v1: &[f64; 3], v2: &[f64; 3]) -> [f64; 3] {
        [
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0],
        ]
    }

    fn norm(v: &[f64; 3]) -> f64 {
        (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt()
    }

    fn dot(v1: &[f64; 3], v2: &[f64; 3]) -> f64 {
        v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    }

    proptest! {
        #[test]
        fn test_intrsc(p1 in unit_vector(), p2 in unit_vector(), q in unit_vector()) {
            let c12 = cross_product(&p1, &p2);
            let norm_c12 = norm(&c12);
            if norm_c12 < f64::EPSILON {
                return Ok(());
            }
            let c12 = normalize(&c12);

            let cn = cross_product(&p1, &q);
            let norm_cn = norm(&cn);
            if norm_cn < f64::EPSILON {
                return Ok(());
            }

            let mut p = [0f64; 3];
            let mut ier = 0i32;

            unsafe {
                intrsc(
                    p1.as_ptr(),
                    p2.as_ptr(),
                    cn.as_ptr(),
                    p.as_mut_ptr(),
                    &raw mut ier,
                );
            };

            prop_assert_eq!(ier, 0, "intrsc failed");

            let norm_p = norm(&p);
            prop_assert!((norm_p - 1.0).abs() < 1e-10, "intersection point should be on unit sphere, norm={norm_p}");

            let dot_p_cn = dot(&p, &cn);
            prop_assert!(dot_p_cn.abs() < 1e-10, "intersection should be perpendicular to cn, dot={dot_p_cn}");

            let dot_p_c12 = dot(&p, &c12);
            prop_assert!(dot_p_c12.abs() < 1e-10, "intersection should be perpendicular to c12 normal, dot={dot_p_c12}");
        }
    }

    proptest! {
        #[test]
        fn test_r83vec_normalize(x in -10.0f64..10.0, y in -10.0f64..10.0, z in -10.0f64..10.0f64) {
            let expected_vector = normalize(&[x, y, z]);
            let mut xv = [x];
            let mut yv = [y];
            let mut zv = [z];
            let n = 1i32;
            unsafe {
                r83vec_normalize(
                    &raw const n,
                    xv.as_mut_ptr(),
                    yv.as_mut_ptr(),
                    zv.as_mut_ptr()
                );
            }

            prop_assert!((xv[0] - expected_vector[0]).abs() < 1e-10, "normalize failed. Expected x={} got x={}", expected_vector[0], xv[0]);
            prop_assert!((yv[0] - expected_vector[1]).abs() < 1e-10, "normalize failed. Expected y={} got y={}", expected_vector[1], yv[0]);
            prop_assert!((zv[0] - expected_vector[2]).abs() < 1e-10, "normalize failed. Expected z={} got z={}", expected_vector[2], zv[0]);
        }
    }

    proptest! {
        #[test]
        fn test_trlist2(n in 5i32..25i32) {
            let (x, y, z) = fibonacci_sphere(n as usize);
            let (list, lptr, lend, lnew) = create_triangulation(n, &x, &y, &z);

            if lnew <= 0 {
                return Ok(());
            }

            let max_triangles = (2 * n - 4) as usize;
            let mut ltri = vec![0i32; 3 * max_triangles];
            let mut nt = 0i32;
            let mut ier = 0i32;

            unsafe {
                trlist2(
                    &raw const n,
                    list.as_ptr(),
                    lptr.as_ptr(),
                    lend.as_ptr(),
                    &raw mut nt,
                    ltri.as_mut_ptr(),
                    &raw mut ier,
                );
            }

            prop_assert_eq!(ier, 0, "trlist2 failed");
            prop_assert!(nt > 0, "should have at least one triangle");
            prop_assert!(nt as usize <= max_triangles, "nt exceeds maximum");

            let expected_nt = 2 * n - 4;
            prop_assert_eq!(nt, expected_nt, "triangle count mismatch: expected {} for n={}, got {}", expected_nt, n, nt);

            for t in 0..nt as usize {
                let base = t * 3;
                let v1 = ltri[base];
                let v2 = ltri[base + 1];
                let v3 = ltri[base + 2];

                prop_assert!(v1 >= 1 && v1 <= n, "invalid v1: {v1}");
                prop_assert!(v2 >= 1 && v2 <= n, "invalid v2: {v2}");
                prop_assert!(v3 >= 1 && v3 <= n, "invalid v3: {v3}");

                prop_assert_ne!(v1, v2, "duplicate vertices in triangle");
                prop_assert_ne!(v2, v3, "duplicate vertices in triangle");
                prop_assert_ne!(v1, v3, "duplicate vertices in triangle");

                prop_assert!(v1 < v2, "v1 should be smallest: v1={v1}, v2={v2}");
                prop_assert!(v1 < v3, "v1 should be smallest: v1={v1}, v3={v3}");
            }

            for t in 1..nt as usize {
                let prev_base = (t - 1) * 3;
                let curr_base = t * 3;
                let prev_v1 = ltri[prev_base];
                let curr_v1 = ltri[curr_base];

                prop_assert!(curr_v1 >= prev_v1, "triangles not ordered: {prev_v1} before {curr_v1}");
            }
        }
    }
}
