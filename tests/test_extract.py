try:
    from . import common as c
except BaseException:
    import common as c


class extractTest(c.unittest.TestCase):
    def test_extract_edges(self):
        """Test all possible scenarios in the extract edges routine"""
        # 1D Parametric (single)
        bezier_line = c.splinepy.Bezier(
            degrees=[2], control_points=[[0, 1], [1, 1], [2, 1]]
        )
        edges = bezier_line.extract.edges(resolution=20)
        self.assertTrue(
            c.np.allclose(edges.vertices[:, 0], c.np.linspace(0, 2, 20))
        )

        # 1D Parametric (MP)
        bezier_line_2 = bezier_line.copy()
        bezier_line_2.cps += [2, 0]
        mp_line = c.splinepy.Multipatch(splines=[bezier_line, bezier_line_2])
        edges = mp_line.extract.edges(resolution=20)
        self.assertTrue(
            c.np.allclose(
                edges.vertices[:, 0],
                c.np.hstack(
                    (c.np.linspace(0, 2, 20), c.np.linspace(2, 4, 20))
                ),
            )
        )

        # Extract patch-borders (single patch)
        bezier_cube = c.splinepy.Bezier(
            [1, 1, 1],
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [0, 1, 1],
                [1, 1, 1],
            ],
        ).bspline
        bezier_cube.insert_knots(0, [i * 0.25 for i in range(1, 4)])
        bezier_cube.insert_knots(1, [i * 0.33333 for i in range(1, 3)])
        bezier_cube.insert_knots(2, [i * 0.2 for i in range(1, 5)])
        edges = bezier_cube.extract.edges(resolution=2, all_knots=False)
        expected = c.np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
            ]
        )
        self.assertTrue(c.np.allclose(edges.vertices, expected))

        # Extract Patch border (MP)
        bezier_cube_2 = bezier_cube.copy()
        bezier_cube_2.cps += [0, 0, 1]
        mp_cubes = c.splinepy.Multipatch(splines=[bezier_cube, bezier_cube_2])
        edges = mp_cubes.extract.edges(resolution=2, all_knots=False)
        self.assertTrue(
            c.np.allclose(
                edges.vertices, c.np.vstack([expected, expected + [0, 0, 1]])
            )
        )

    def test_extract_faces(self):
        """Test extraction of faces"""
        # Two dimensional single patch
        bspline = c.splinepy.BSpline(
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
            degrees=[1, 1],
            control_points=[[0, 0], [2, 0], [1, 1], [3, 1]],
        )
        # Modify to make more interesting
        bspline.elevate_degrees([0, 0, 1, 1])
        faces = bspline.extract.faces(resolution=6)
        self.assertTrue(c.np.allclose(faces.vertices, bspline.sample(6)))
        # Check a random face center
        face_id = c.np.random.randint(0, 5 * 5)
        face_col = face_id // 5
        face_row = face_id % 5
        d_x = 0.2
        query_points = c.np.array([[0, 0], [d_x, 0], [0, d_x], [d_x, d_x]]) + [
            face_row * d_x,
            face_col * d_x,
        ]
        face_center = c.np.sum(bspline.evaluate(query_points), axis=0) / 4
        self.assertTrue(
            c.np.allclose(face_center, faces.centers()[face_id, :])
        )

        # Two dimensional single patch
        bspline_2 = bspline.copy()
        bspline_2.cps += [1, 1]
        mp_2d = c.splinepy.Multipatch(splines=[bspline, bspline_2])
        mp_faces = mp_2d.extract.faces(resolution=6, watertight=False)
        self.assertTrue(c.np.allclose(mp_faces.vertices, mp_2d.sample(6)))

        # Check faces
        self.assertTrue(
            c.np.allclose(face_center, faces.centers()[face_id, :])
        )
        self.assertTrue(
            c.np.allclose(
                face_center + [1, 1],
                mp_faces.centers()[face_id + faces.elements.shape[0], :],
            )
        )

        # Watertighten faces
        mp_faces_wt = mp_2d.extract.faces(resolution=6, watertight=True)

        # Check faces (remain true)
        face_id = c.np.random.randint(0, 5 * 5)
        self.assertTrue(
            c.np.allclose(
                mp_faces_wt.centers()[face_id, :], faces.centers()[face_id, :]
            )
        )
        self.assertTrue(
            mp_faces.vertices.shape[0] - mp_faces_wt.vertices.shape[0] == 6
        )

        # p3D
        nurbs = c.splinepy.NURBS(
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]],
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [2, 0, 0],
                [1, 1, 0],
                [3, 1, 0],
                [0, 0, 3],
                [2, 0, 3],
                [1, 1, 3],
                [3, 1, 3],
            ],
            weights=c.np.ones(8),
        )
        faces = nurbs.extract.faces(resolution=8, watertight=False)
        self.assertTrue(c.np.allclose(faces.vertices.shape[0], 8 * 8 * 6))

        # Check watertighten
        faces = nurbs.extract.faces(resolution=8, watertight=True)
        self.assertTrue(faces.vertices.shape[0] == 8 * 8 * 6 - 6 * 12 - 2 * 8)

        # Check multipatch
        nurbs_2 = nurbs.copy()
        nurbs_2.cps += [2, 0, 0]
        mp_3d = c.splinepy.Multipatch(splines=[nurbs, nurbs_2])
        faces = mp_3d.extract.faces(resolution=8, watertight=False)
        self.assertTrue(c.np.allclose(faces.vertices.shape[0], 2 * 8 * 8 * 5))

        # Check watertighten
        mp_3d = c.splinepy.Multipatch(splines=[nurbs, nurbs_2])
        faces = mp_3d.extract.faces(resolution=8, watertight=True)
        self.assertTrue(
            faces.vertices.shape[0] == 2 * 8 * 8 * 5 - 6 * 20 - 2 * 8 - 3 * 4
        )

        # Check if the boundary centers correspond to the expected faces
        boundaries = (
            nurbs.extract.boundaries() + nurbs_2.extract.boundaries()[1:]
        )
        boundaries.pop(1)  # First id of first spline is duplicate
        expected_faces = c.gus.Faces.concat(
            [b.extract.faces(resolution=8) for b in boundaries]
        )
        self.assertTrue(
            c.np.allclose(
                c.np.sort(expected_faces.centers(), axis=0),
                c.np.sort(faces.centers(), axis=0),
            ),
            "Wrong",
        )

    def test_extract_volumes(self):
        """Check the extraction of volumetric meshes from splines"""
        error_spline = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [0, 0], [0, 0], [1, 1]]
        )
        self.assertRaisesRegex(
            ValueError,
            "Volume extraction from a spline is only valid for para_dim: 3 dim"
            ": 3 splines.",
            error_spline.extract.volumes,
            resolution=5,
        )
        # Single Patch
        nurbs = c.splinepy.NURBS(
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1], [0, 0, 1, 1]],
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [2, 0, 0],
                [1, 1, 0],
                [3, 1, 0],
                [0, 0, 3],
                [2, 0, 3],
                [1, 1, 3],
                [3, 1, 3],
            ],
            weights=c.np.ones(8),
        )
        volumes_sp = nurbs.extract.volumes(resolution=6)

        # Use linearity of spline to check against evaluations
        queries = c.splinepy.utils.data.cartesian_product(
            [c.np.linspace(0.1, 0.9, 5) for _ in range(3)]
        )

        # Check length and compare centers
        self.assertTrue(volumes_sp.elements.shape[0] == 5**3)
        self.assertTrue(
            c.np.allclose(nurbs.evaluate(queries), volumes_sp.centers())
        )

        # Watertight flag does not have an influence on single patch geometries
        self.assertTrue(
            c.np.allclose(
                nurbs.extract.volumes(resolution=6, watertight=True).vertices,
                volumes_sp.vertices,
            )
        )

        # Check multipatch geometries
        nurbs_2 = nurbs.copy()
        nurbs_2.cps += [0, 0, 3]
        mp_nurbs = c.splinepy.Multipatch(splines=[nurbs, nurbs_2])

        # Retrieve volumes from multipatch
        volumes_mp = mp_nurbs.extract.volumes(resolution=6, watertight=False)

        # Check length and compare centers
        self.assertTrue(volumes_mp.elements.shape[0] == 2 * 5**3)
        self.assertTrue(
            c.np.allclose(mp_nurbs.evaluate(queries), volumes_mp.centers())
        )

        # Copy mp_nurbs because it currently has issues #198
        mp_nurbs = c.splinepy.Multipatch(splines=[nurbs, nurbs_2])
        volumes_mp_wt = mp_nurbs.extract.volumes(resolution=6, watertight=True)

        # Watertight flag does not have an influence on single patch geometries
        self.assertTrue(volumes_mp_wt.vertices.shape[0] == 2 * 6**3 - 6**2)
        self.assertTrue(
            c.np.allclose(volumes_mp.centers(), volumes_mp_wt.centers())
        )

    def test_extract_control_points(self):
        """Extract control points of all splines"""
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
        )
        bez_el1 = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[1, 0], [2, 0], [1, 1], [2, 1]]
        )

        # Extract ctps from single patches
        ctps0 = bez_el0.extract.control_points()
        ctps1 = bez_el1.extract.control_points()

        # Extract from mp geometry
        mp = c.splinepy.Multipatch(splines=[bez_el0, bez_el1])
        ctps_mp = mp.extract.control_points()

        self.assertTrue(
            c.np.allclose(
                c.np.vstack([ctps0.vertices, ctps1.vertices]), ctps_mp.vertices
            )
        )

    def test_extract_boundaries_2D(self):
        """Create a Spline, extract boundaries and check validity"""
        # uniform
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1], control_points=[[0, 0], [1, 0], [0, 1], [1, 1]]
        )
        boundary_0 = c.splinepy.Bezier(
            degrees=[1], control_points=[[0, 0], [0, 1]]
        )
        boundary_3 = c.splinepy.Bezier(
            degrees=[1], control_points=[[0, 1], [1, 1]]
        )
        list_of_boundaries = bez_el0.extract.boundaries()
        self.assertTrue(c.are_splines_equal(list_of_boundaries[0], boundary_0))
        self.assertTrue(c.are_splines_equal(list_of_boundaries[3], boundary_3))

        # non-uniform
        bez_el0 = c.splinepy.NURBS(
            knot_vectors=[[0, 0, 0.5, 1, 1], [0, 0, 1, 1]],
            weights=[1.0, 0.5, 1, 1, 0.5, 1],
            degrees=[1, 1],
            control_points=[
                [0, 0],
                [0.5, 0.5],
                [1, 0],
                [0, 1],
                [0.5, 1.5],
                [1, 1],
            ],
        )
        boundary_2 = c.splinepy.NURBS(
            knot_vectors=[[0, 0, 0.5, 1, 1]],
            weights=[1.0, 0.5, 1],
            degrees=[1],
            control_points=[[0, 0], [0.5, 0.5], [1, 0]],
        )
        list_of_boundaries = bez_el0.extract.boundaries()
        self.assertTrue(c.are_splines_equal(list_of_boundaries[2], boundary_2))

    def test_extract_boundaries_3D(self):
        """Create a Spline, extract boundaries and check validity"""
        # uniform
        bez_el0 = c.splinepy.Bezier(
            degrees=[1, 1, 1],
            control_points=[
                [0, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [0, 1, 1],
                [1, 1, 1],
            ],
        )
        boundary_0 = c.splinepy.Bezier(
            degrees=[1, 1],
            control_points=[
                [0, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [0, 1, 1],
            ],
        )
        list_of_boundaries = bez_el0.extract.boundaries()[0]
        self.assertTrue(c.are_splines_equal(list_of_boundaries, boundary_0))

    def test_extract_arrow_data(self):
        """
        see if arrow_data is created as we wished
        """
        spline = c.splinepy.helpme.create.box(10, 20, 30)
        # try to avoid zero division below
        spline.cps -= 0.0001

        # set self data and resolutions
        spline.spline_data["me"] = spline
        res = [11, 13, 17]
        spline.show_options["resolutions"] = res

        # 1. check number of edges
        me_edges = spline.extract.arrow_data("me")
        assert isinstance(me_edges, c.gus.Edges)
        # 3d should give you arrows on the faces
        total = 0
        for i in range(spline.para_dim):
            r = res.copy()
            r.pop(i)
            total += c.np.prod(r) * 2

        assert len(me_edges.edges) == total

        # 2. check if they point to the direction they need to point
        def check_values(a, b):
            # normalize - this may do some zero division
            a /= c.np.linalg.norm(a, axis=1).reshape(-1, 1)
            b /= c.np.linalg.norm(b, axis=1).reshape(-1, 1)
            a[c.np.isnan(a)] = 0
            b[c.np.isnan(b)] = 0

            assert c.np.allclose(a, b)

        direction = c.np.diff(
            me_edges.vertices[me_edges.edges], axis=1
        ).reshape(-1, me_edges.vertices.shape[1])
        sample = spline.extract.faces(res, watertight=False).vertices

        check_values(direction, sample)

        # 3. do the same with locations
        locations = c.np.random.random((10, 3))
        spline.spline_data["me_at"] = c.splinepy.SplineDataAdaptor(
            spline, locations=locations
        )
        me_at_edges = spline.extract.arrow_data("me_at")

        # 3.1 len test
        assert isinstance(me_at_edges, c.gus.Edges)
        assert len(me_at_edges.edges) == len(locations)

        # 3.2 arrow validity
        direction = c.np.diff(
            me_at_edges.vertices[me_at_edges.edges], axis=1
        ).reshape(-1, me_at_edges.vertices.shape[1])
        sample = spline.evaluate(locations)

        check_values(direction, sample)

        # 4. at last, scaling
        spline.show_options["arrow_data_scale"] = 200
        me_at_edges = spline.extract.arrow_data("me_at")

        # 3.1 len test
        assert isinstance(me_at_edges, c.gus.Edges)
        assert len(me_at_edges.edges) == len(locations)

        # 3.2 arrow validity
        direction = c.np.diff(
            me_at_edges.vertices[me_at_edges.edges], axis=1
        ).reshape(-1, me_at_edges.vertices.shape[1])
        sample = spline.evaluate(locations)
        assert c.np.allclose(
            direction, sample * spline.show_options["arrow_data_scale"]
        )


if __name__ == "__main__":
    c.unittest.main()
