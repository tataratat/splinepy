try:
    from . import common as c
except BaseException:
    import common as c


class SplineListTest(c.unittest.TestCase):
    def test_spline_list_evaluate(self):
        """list evaluation"""
        # get list of splines
        splines = c.all2p2d()

        # get reference solution
        ref_solutions = c.np.vstack([s.evaluate(c.q2D) for s in splines])

        # call eval, don't check dim
        list_solutions = c.splinepy.splinepy_core.lists.evaluate(
            splines, c.q2D, nthreads=2, check_dims=False
        )
        assert c.np.allclose(ref_solutions, list_solutions)

        # call eval, check dim
        list_solutions = c.splinepy.splinepy_core.lists.evaluate(
            splines, c.q2D, nthreads=2, check_dims=True
        )
        assert c.np.allclose(ref_solutions, list_solutions)

        # call eval, single thread
        list_solutions = c.splinepy.splinepy_core.lists.evaluate(
            splines, c.q2D, nthreads=1, check_dims=True
        )
        assert c.np.allclose(ref_solutions, list_solutions)

    def test_spline_list_sample(self):
        """spline sample"""
        splines = c.all2p2d()

        # set resolution. this is int,
        # because it is implemented for equal resolution only
        resolution = 3

        ref_solutions = c.np.vstack(
            [s.sample([resolution] * splines[0].para_dim) for s in splines]
        )

        # check everything and compute
        list_samples = c.splinepy.splinepy_core.lists.sample(
            splines,
            resolution,
            nthreads=2,
            same_parametric_bounds=True,
            check_dims=True,
        )
        assert c.np.allclose(ref_solutions, list_samples)

        # don't check p bounds
        list_samples = c.splinepy.splinepy_core.lists.sample(
            splines,
            resolution,
            nthreads=2,
            same_parametric_bounds=True,
            check_dims=True,
        )
        assert c.np.allclose(ref_solutions, list_samples)

        # update p bounds for bsplines
        bspline_id = 2
        nurbs_id = 3
        for i in [bspline_id, nurbs_id]:
            new_kvs = [kv * 2 for kv in splines[i].kvs]
            splines[i].kvs = new_kvs

        # check p bounds,
        list_samples = c.splinepy.splinepy_core.lists.sample(
            splines,
            resolution,
            nthreads=2,
            same_parametric_bounds=False,
            check_dims=True,
        )
        assert c.np.allclose(ref_solutions, list_samples)

    def test_extract_boundaries(self):
        """extract boundaries. allowed for splines with unmatching dims"""
        # prepare 2d splines
        splines2d = c.all2p2d()

        # prepare mixed dim splines
        splines2d3d = [*splines2d, *c.all3p3d()]

        # prepare test func
        def _test(pure_list, same_p_dims):
            """actual test routine"""
            ref_boundaries = []
            for s in pure_list:
                ref_boundaries.extend(
                    c.splinepy.splinepy_core.extract_boundaries(s, [])
                )

            list_boundaries = (
                c.splinepy.splinepy_core.lists.extract_boundaries(
                    pure_list, nthreads=2, same_para_dims=False
                )
            )

            # for same p_dims,
            # try same_para_dims=True by extending test splines
            if same_p_dims:
                ref_boundaries.extend(ref_boundaries)
                list_boundaries.extend(
                    c.splinepy.splinepy_core.lists.extract_boundaries(
                        pure_list, nthreads=2, same_para_dims=True
                    )
                )

            # should have same size
            assert len(ref_boundaries) == len(list_boundaries)

            # same spline?
            for rb, lb in zip(ref_boundaries, list_boundaries):
                # currently, core function does not have getter
                # for individual attr, so cast
                assert c.are_splines_equal(c.to_derived(rb), c.to_derived(lb))

        _test(splines2d, True)
        _test(splines2d3d, False)

    def test_boundary_centers(self):
        """boundary centers"""
        splines = c.all2p2d()

        ref_centers = c.np.vstack(
            [c.splinepy.splinepy_core.boundary_centers(s) for s in splines]
        )

        # test same bounds
        list_centers = c.splinepy.splinepy_core.lists.boundary_centers(
            splines, nthreads=2, same_parametric_bounds=True, check_dims=False
        )
        assert c.np.allclose(ref_centers, list_centers)

        # different bounds
        # update p bounds for bsplines
        bspline_id = 2
        nurbs_id = 3
        for i in [bspline_id, nurbs_id]:
            new_kvs = [kv * 2 for kv in splines[i].kvs]
            splines[i].kvs = new_kvs

        list_centers = c.splinepy.splinepy_core.lists.boundary_centers(
            splines, nthreads=2, same_parametric_bounds=False, check_dims=False
        )
        assert c.np.allclose(ref_centers, list_centers)

    def test_raise_dim_mismatch(self):
        """see if function raises dim mismatch correctly"""
        slist2d = c.all2p2d()

        # prepare mixed
        slist2d3d = [*c.all2p2d(), *c.all3p3d()]

        # nothing happens
        c.splinepy.splinepy_core.lists.raise_dim_mismatch(slist2d, nthreads=2)

        with self.assertRaises(RuntimeError):
            c.splinepy.splinepy_core.lists.raise_dim_mismatch(
                slist2d3d, nthreads=2
            )

    def _bezier_noisy_boxes_and_test_shapes(self):
        # prepare boxes with some noise
        box2d = c.nd_box(2)
        box2d.pop("knot_vectors")
        rbox2d = c.splinepy.RationalBezier(**box2d)
        box2d.pop("weights")
        zbox2d = c.splinepy.Bezier(**box2d)
        box3d = c.nd_box(3)
        box3d.pop("knot_vectors")
        rbox3d = c.splinepy.RationalBezier(**box3d)
        box3d.pop("weights")
        zbox3d = c.splinepy.Bezier(**box2d)

        outer = [zbox2d, rbox2d, zbox3d, rbox3d]
        # add some noise
        for o in outer:
            o.cps = o.cps + c.np.random.normal(0, 0.025, o.cps.shape)

        # test shapes
        z2 = c.splinepy.Bezier(**c.z2p2d())
        z3 = c.splinepy.Bezier(**c.z3p3d())
        r2 = c.splinepy.RationalBezier(**c.r2p2d())
        r3 = c.splinepy.RationalBezier(**c.r3p3d())

        return outer[:2], outer[2:], [z2, r2], [z3, r3]

    def test_list_compose(self):
        """check if list compose yield same spline as spline compose"""
        pass

    def _bezier_noisy_boxes_and_test_shapes(self):
        # prepare boxes with some noise
        box2d = c.nd_box(2)
        box2d.pop("knot_vectors")
        rbox2d = c.splinepy.RationalBezier(**box2d)
        box2d.pop("weights")
        zbox2d = c.splinepy.Bezier(**box2d)
        box3d = c.nd_box(3)
        box3d.pop("knot_vectors")
        rbox3d = c.splinepy.RationalBezier(**box3d)
        box3d.pop("weights")
        zbox3d = c.splinepy.Bezier(**box2d)

        outer = [zbox2d, rbox2d, zbox3d, rbox3d]
        # add some noise
        for o in outer:
            o.cps = o.cps + c.np.random.normal(0, 0.025, o.cps.shape)

        # test shapes
        z2 = c.splinepy.Bezier(**c.z2p2d())
        z3 = c.splinepy.Bezier(**c.z3p3d())
        r2 = c.splinepy.RationalBezier(**c.r2p2d())
        r3 = c.splinepy.RationalBezier(**c.r3p3d())

        return outer[:2], outer[2:], [z2, r2], [z3, r3]

    def test_list_compose(self):
        """check if list compose yield same spline as spline compose"""
        # prepare boxes with some noise
        outer_2d, outer_3d, inner_2d, inner_3d = [
            c.splinepy.splinepy_core.SplineList(beziers)
            for beziers in self._bezier_noisy_boxes_and_test_shapes()
        ]

        t = []
        for beziers in self._bezier_noisy_boxes_and_test_shapes():
            sl = c.splinepy.splinepy_core.SplineList(beziers)
            print(sl[0])
            t.append(sl)


        ref_composed = [o.compose(i) for o, i in zip(outer_2d, inner_2d)]

        # list compose
        list_composed = outer_2d.compose(
            inner_2d, cartesian_product=False, nthreads=2
        )

        for r, l in zip(ref_composed, list_composed):
            assert c.are_splines_equal(r, l)


if __name__ == "__main__":
    c.unittest.main()
