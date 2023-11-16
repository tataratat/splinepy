try:
    from . import common as c
except BaseException:
    import common as c


def cps_are_synced(spl):
    assert c.np.allclose(
        spl.cps, spl.current_core_properties()["control_points"]
    )


class InplaceModificationTest(c.SplineBasedTestCase):
    def test_inplace_change_degrees(self):
        """inplace change of degrees should not be allowed if core spline is
        initialized"""
        z = c.dict_bezier_2p2d()
        r = c.dict_rational_bezier_2p2d()
        b = c.dict_bspline_2p2d()
        n = c.dict_nurbs_2p2d()

        Z, R, B, N = (
            c.splinepy.Bezier,
            c.splinepy.RationalBezier,
            c.splinepy.BSpline,
            c.splinepy.NURBS,
        )

        for props, SClass in zip((z, r, b, n), (Z, R, B, N)):
            # test no core spline
            # this should be fine
            s = SClass()
            s.degrees = props["degrees"]
            s.degrees += 1

            # this shoundn't be fine
            s._new_core(**props)
            with self.assertRaises(ValueError):
                s.degrees += 1

    def test_inplace_change_knot_vectors(self):
        """test inplace change of knot_vectors"""
        # let's test 3D splines
        dim = 3
        box_data = c.nd_box(dim)
        n = c.splinepy.NURBS(**box_data)
        box_data.pop("weights")
        b = c.splinepy.BSpline(**box_data)

        # elevate degrees and insert some knots
        knot_insert_dims = [1, 2]
        for s in (n, b):
            s.elevate_degrees([0, 0, 1, 2])
            for kid in knot_insert_dims:
                s.insert_knots(kid, c.np.linspace(0.1, 0.9, 9))

        qres = 5
        raster_query = c.raster([[0] * dim, [1] * dim], [qres] * dim)

        for s in (n, b):
            # is this valid spline?
            assert c.np.allclose(raster_query, s.evaluate(raster_query))

            # modify knots and corresponding queries, so that
            # evaluated points are same as raster_query
            modified_query = raster_query.copy()
            factor = 2
            for kid in knot_insert_dims:
                s.knot_vectors[kid][:] = c.np.multiply(
                    s.knot_vectors[kid], factor
                ).tolist()
                modified_query[:, kid] *= factor

            # evaluation check
            assert c.np.allclose(raster_query, s.evaluate(modified_query))

    def test_inplace_change_control_points(self):
        """test inplace changes of control points"""
        dim = 3
        res = [3] * dim
        box_data = c.nd_box(dim)

        n = c.splinepy.NURBS(**box_data)
        weights = box_data.pop("weights")
        b = c.splinepy.BSpline(**box_data)
        box_data.pop("knot_vectors")
        z = c.splinepy.Bezier(**box_data)
        r = c.splinepy.RationalBezier(**box_data, weights=weights)

        for s in (z, r, b, n):
            # init
            orig = s.copy()

            # check if we are good to start
            assert c.np.allclose(orig.sample(res), s.sample(res))

            # hold reference to control points
            s_cps = s.control_points
            # modify cps
            s.control_points /= 2
            assert c.np.allclose(orig.sample(res) / 2, s.sample(res))
            # check if they are still the same array
            assert s_cps is s.control_points

            # modify ws
            if s.is_rational:
                s_ws = s.weights
                s.weights *= 0.5
                assert s_ws is s.weights

    def test_inplace_change_weights(self):
        """test inplace change of weights by comparing quarter circle"""
        n_q_circle = c.nurbs_2p2d_quarter_circle()
        res = [4] * 2
        # init rational
        n = c.splinepy.NURBS(**n_q_circle)
        kvs = n_q_circle.pop("knot_vectors")
        r = c.splinepy.RationalBezier(**n_q_circle)
        # init non-rational
        n_q_circle.pop("weights")
        b = c.splinepy.BSpline(**n_q_circle, knot_vectors=kvs)
        z = c.splinepy.Bezier(**n_q_circle)

        # make sure they aren't the same
        assert not c.np.allclose(b.sample(res), n.sample(res))
        assert not c.np.allclose(z.sample(res), r.sample(res))

        # modify weights of rational splines
        for rs in (n, r):
            # set them all to 1
            rs.weights[abs(rs.weights - 1.0) > 1e-10] = 1.0

        # now, they should be the same
        assert c.np.allclose(b.sample(res), n.sample(res))
        assert c.np.allclose(z.sample(res), r.sample(res))
        assert c.np.allclose(b.sample(res), r.sample(res))
        assert c.np.allclose(z.sample(res), n.sample(res))

    def test_physical_space_array(self):
        """tests if physical space array is syncing correctly."""
        n = c.splinepy.NURBS(**c.nurbs_2p2d_quarter_circle())

        # all setitem cases
        # 1.
        n.cps[0] = 100.0
        cps_are_synced(n)

        n.cps[c.np.array(1)] = [0.99, 0.11]  # here, weight is not 1
        cps_are_synced(n)

        # 2.
        n.cps[c.np.array([1, 4, 3])] = 123
        cps_are_synced(n)

        b_mask = c.np.ones(len(n.cps), dtype=bool)
        b_mask[[0, 4]] = False
        n.cps[b_mask] = -1235.5
        cps_are_synced(n)

        # 3.
        n.cps[[2, 0, 5]] = 299354.0
        cps_are_synced(n)

        n.cps[b_mask.tolist()] = 975.5
        cps_are_synced(n)

        # 4.
        n.cps[3, :2] = [-10, -20]
        cps_are_synced(n)

        n.cps[5, [1, 0]] = [1.722, 24.29]
        cps_are_synced(n)

        n.cps[c.np.array(2), 1] = 7
        cps_are_synced(n)

        # 5.
        n.cps[[0, 5, 2], 0] = 92314.111234
        cps_are_synced(n)
        n.cps[b_mask, :2] = n.cps[[2, 5, 0, 1], -2:]
        cps_are_synced(n)
        n.cps[b_mask.tolist(), 0] = n.cps[b_mask[::-1], 1]
        cps_are_synced(n)

        # 6.
        n.cps[:1, :2] = 160.4
        cps_are_synced(n)
        n.cps[2:-1, :2] = 78883.2
        n.cps[-1:-4:-1, :1] = -1687
        cps_are_synced(n)

        # 7.
        n.cps[...] = 1237792.4634527
        n.cps[..., :1] = -62.5
        cps_are_synced(n)

        # 8. if you have ideas for use case not listed above, please add!

        # now, child arrays
        carr = n.cps[0]
        assert carr.base is n.cps

        carr += [1324235.1, -234.5]
        assert c.np.allclose(carr, n.cps[0])
        cps_are_synced(n)

        # copy should return a normal np.ndarray
        ndarr = n.cps.copy()
        assert isinstance(ndarr, c.np.ndarray)
        assert not isinstance(ndarr, type(n.cps))
        assert not isinstance(ndarr, c.splinepy.utils.data.PhysicalSpaceArray)

    def test_new_control_point_pointers_creation(self):
        """Test cp pointers are newly created if they change.
        _sync_source_ptr compares len of cp pointers and PhysicalSpaceArray
        before syncing.
        """
        for s in self.all_2p2d_splines():
            # initial check
            s.cps._sync_source_ptr()
            cps_are_synced(s)

            # check after degree elevation
            s.elevate_degrees([0, 1])
            s.cps._sync_source_ptr()
            cps_are_synced(s)

            # and after knot vectors and reduce degree
            if s.has_knot_vectors:
                s.insert_knots(0, [0.2, 0.5, 0.7])
                s.cps._sync_source_ptr()
                cps_are_synced(s)

                s.insert_knots(1, [0.2, 0.5, 0.7])
                s.cps._sync_source_ptr()
                cps_are_synced(s)

                s.remove_knots(0, [0.7])
                s.cps._sync_source_ptr()
                cps_are_synced(s)

                # reduce degree only available for bsplines
                s.reduce_degrees([0])
                s.cps._sync_source_ptr()
                cps_are_synced(s)


if __name__ == "__main__":
    c.unittest.main()
