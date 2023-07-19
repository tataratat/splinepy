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
                [3, 1, 3]
            ],
            weights=c.np.ones((8))
        )
        faces = nurbs.extract.faces(resolution=8, watertight=False)
        self.assertTrue(c.np.allclose(faces.vertices.shape[0], 8*8*6))
        
        # Check watertighten
        faces = nurbs.extract.faces(resolution=8, watertight=True)
        self.assertTrue(c.np.allclose(faces.vertices.shape[0], 8*8*6 - 6*12 - 2*8))    

        # Check multipatch
        nurbs_2 = nurbs.copy()    
        nurbs_2.cps += [2 ,0 ,0]
        mp_3d = c.splinepy.Multipatch(splines=[nurbs, nurbs_2])
        faces = mp_3d.extract.faces(resolution=8, watertight=False)
        self.assertTrue(c.np.allclose(faces.vertices.shape[0], 2*8*8*5))
        
        # Check watertighten       
        mp_3d = c.splinepy.Multipatch(splines=[nurbs, nurbs_2])
        faces = mp_3d.extract.faces(resolution=8, watertight=True)
        self.assertTrue(c.np.allclose(faces.vertices.shape[0], 2*8*8*5 - 6*20 - 2*8 - 3*4)) 


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
