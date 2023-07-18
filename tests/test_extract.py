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
        list_of_boundaries = bez_el0.extract_boundaries([0, 3])
        self.assertTrue(c.are_splines_equal(list_of_boundaries[0], boundary_0))
        self.assertTrue(c.are_splines_equal(list_of_boundaries[1], boundary_3))

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
        list_of_boundaries = bez_el0.extract_boundaries([2])
        self.assertTrue(c.are_splines_equal(list_of_boundaries[0], boundary_2))

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
        list_of_boundaries = bez_el0.extract_boundaries([0])
        self.assertTrue(c.are_splines_equal(list_of_boundaries[0], boundary_0))


if __name__ == "__main__":
    c.unittest.main()
