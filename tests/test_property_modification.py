try:
    from . import common as c
except BaseException:
    import common as c


class InplaceModificationTest(c.unittest.TestCase):
    def test_inplace_change_degrees(self):
        """inplace change of degrees should not be allowed if core spline is
        initialized"""
        z = c.z2p2d()
        r = c.r2p2d()
        b = c.b2p2d()
        n = c.n2p2d()

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
            s.new_core(**props)
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
                s.knot_vectors[kid] *= factor
                modified_query[:, kid] *= factor

                # check modified flag
                assert s.knot_vectors[kid]._modified
                # full modified check
                assert c.splinepy.spline.is_modified(s)

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

            # modify cps
            s.control_points /= 2
            assert c.np.allclose(orig.sample(res) / 2, s.sample(res))

    def test_inplace_change_weights(self):
        """test inplace change of weights by compareing quarter circle"""
        n_q_circle = c.n2p2d_quarter_circle()
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


class CoordinateReferencesModificationTest(c.unittest.TestCase):
    def test_coordinate_references(self):
        dim = 3
        res = [3] * dim
        box_data = c.nd_box(dim)

        n = c.splinepy.NURBS(**box_data)
        weights = box_data.pop("weights")
        b = c.splinepy.BSpline(**box_data)
        box_data.pop("knot_vectors")
        z = c.splinepy.Bezier(**box_data)
        r = c.splinepy.RationalBezier(**box_data, weights=weights)

        query = c.raster([[0] * dim, [1] * dim], res)

        for s in (z, r, b, n):
            orig = s.copy()

            # assign all to one value - tests broadcast_scalar
            # should be okay for rational splines,
            # since weights are all one
            s.coordinate_references[:] = 1.0
            evaluated = s.evaluate(query)
            assert c.np.allclose(evaluated, c.np.ones_like(evaluated))
            assert not c.splinepy.spline.is_modified(s)

            s = orig.copy()
            dim_to_modify = int(dim - 1)
            factor = 2.3
            # inplace modification is not supported
            s.coordinate_references[:, dim_to_modify] = (
                s.coordinate_references.numpy()[:, dim_to_modify] * factor
            )
            ref_eval = query.copy()
            ref_eval[:, dim_to_modify] *= factor
            assert c.np.allclose(s.evaluate(query), ref_eval)

            # let's do exactly the same, but see if it also updates
            s = orig.copy()
            s.coordinate_references[:, dim_to_modify] = (
                s.control_points[:, dim_to_modify] * factor
            )
            assert c.np.allclose(
                s.coordinate_references.numpy(), s.control_points
            )
            assert c.np.allclose(s.evaluate(query), ref_eval)

    def test_coordinate_references_weight_handling(self):
        """test if coordiante references reflect weights
        and handles weight as expected is apply_weight is True."""

        n_q_circle = c.n2p2d_quarter_circle()
        res = [4] * 2
        # init rational
        n = c.splinepy.NURBS(**n_q_circle)
        n_q_circle.pop("knot_vectors")
        r = c.splinepy.RationalBezier(**n_q_circle)

        for rs in (n, r):
            orig_cps = c.np.array(n_q_circle["control_points"])
            orig_ws = c.np.array(n_q_circle["weights"])

            orig_rs = rs.copy()
            # check if references are actually weighted
            assert c.np.allclose(
                rs.coordinate_references.numpy(), orig_cps * orig_ws
            )

            # explicitly set apply_weight to False
            rs.coordinate_references.apply_weight = False
            rs.coordinate_references[:] = orig_cps
            # this will result in unweighted reference
            assert c.np.allclose(rs.coordinate_references.numpy(), orig_cps)
            # but this means it won't be halfhalf circle
            assert not c.np.allclose(rs.sample(res), orig_rs.sample(res))
            # explicitly set apply_weight to True
            rs = orig_rs.copy()
            rs.coordinate_references.apply_weight = True
            weighted = rs.coordinate_references.numpy()
            rs.coordinate_references[:] = orig_cps  # this should apply weight

            assert c.np.allclose(rs.coordinate_references.numpy(), weighted)
            assert c.np.allclose(rs.control_points, orig_cps)

            # this will result in same eval
            assert c.np.allclose(rs.sample(res), orig_rs.sample(res))


if __name__ == "__main__":
    c.unittest.main()
