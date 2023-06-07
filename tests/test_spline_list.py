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
        zbox3d = c.splinepy.Bezier(**box3d)

        outer = [zbox2d, rbox2d, zbox3d, rbox3d]
        # add some noise
        for o in outer:
            o.cps = o.cps + c.np.random.normal(0, 0.025, o.cps.shape)

        # test shapes
        z2 = c.splinepy.Bezier(**c.z2p2d())
        z3 = c.splinepy.Bezier(**c.z3p3d())
        r2 = c.splinepy.RationalBezier(**c.r2p2d())
        r3 = c.splinepy.RationalBezier(**c.r3p3d())

        # fit into unit box
        for bez in [z2, z3, r2, r3]:
            # offset
            bez.cps -= bez.cps.min(axis=0)
            # scale
            bez.cps = bez.cps * (1 / bez.cps.max(axis=0))

        return outer[:2], outer[2:], [z2, r2], [z3, r3]

    def test_list_compose(self):
        """check if list compose yield same spline as spline compose"""
        # prepare boxes with some noise
        (
            outer_2d,
            outer_3d,
            inner_2d,
            inner_3d,
        ) = self._bezier_noisy_boxes_and_test_shapes()

        ref_composed = [o.compose(i) for o, i in zip(outer_2d, inner_2d)]

        # list compose
        list_composed = c.splinepy.splinepy_core.lists.compose(
            outer_2d, inner_2d, cartesian_product=False, nthreads=2
        )

        # add 3d
        ref_composed.extend([o.compose(i) for o, i in zip(outer_3d, inner_3d)])
        list_composed.extend(
            c.splinepy.splinepy_core.lists.compose(
                outer_3d, inner_3d, cartesian_product=False, nthreads=2
            )
        )

        for r, l in zip(ref_composed, list_composed):
            assert c.are_splines_equal(c.to_derived(r), c.to_derived(l))

        # test cartesian
        ref_composed = []
        for o in outer_2d:
            for i in inner_2d:
                ref_composed.append(o.compose(i))
        for o in outer_3d:
            for i in inner_3d:
                ref_composed.append(o.compose(i))

        list_composed = c.splinepy.splinepy_core.lists.compose(
            outer_2d, inner_2d, cartesian_product=True, nthreads=2
        )
        list_composed.extend(
            c.splinepy.splinepy_core.lists.compose(
                outer_3d, inner_3d, cartesian_product=True, nthreads=2
            )
        )

        for r, l in zip(ref_composed, list_composed):
            assert c.are_splines_equal(c.to_derived(r), c.to_derived(l))

    def test_list_compose_derivatives(self):
        # prepare boxes with some noise
        (
            outer_2d,
            outer_3d,
            inner_2d,
            inner_3d,
        ) = self._bezier_noisy_boxes_and_test_shapes()

        # create 2d ref
        ref_composed = [
            o.composition_derivative(i, i) for o, i in zip(outer_2d, inner_2d)
        ]

        # list compose
        list_composed = c.splinepy.splinepy_core.lists.composition_derivative(
            outer_2d, inner_2d, inner_2d, cartesian_product=False, nthreads=2
        )

        # add 3d
        ref_composed.extend(
            [
                o.composition_derivative(i, i)
                for o, i in zip(outer_3d, inner_3d)
            ]
        )
        list_composed.extend(
            c.splinepy.splinepy_core.lists.composition_derivative(
                outer_3d,
                inner_3d,
                inner_3d,
                cartesian_product=False,
                nthreads=2,
            )
        )

        for r, l in zip(ref_composed, list_composed):
            assert c.are_splines_equal(c.to_derived(r), c.to_derived(l))

        # test cartesian
        #
        # rational_bz.composition_derivative(poly, poly_der) not supported
        # so, remove
        for inn in [inner_2d, inner_3d]:
            # pop poly
            inn.pop(0)

            # clone rat
            rat = inn[-1].copy()
            rat.cps[:, 0] = 0.5
            rat.cps[:, 0] += c.np.random.uniform(-0.4, 0.4, len(rat.cps))
            # force sync - list exe doesn't sync before exe
            rat.new_core(**rat._data["properties"])
            inn.append(rat)

        ref_composed = []
        for o in outer_2d:
            for i in inner_2d:
                ref_composed.append(o.composition_derivative(i, i))
        for o in outer_3d:
            for i in inner_3d:
                ref_composed.append(o.composition_derivative(i, i))

        list_composed = c.splinepy.splinepy_core.lists.composition_derivative(
            outer_2d, inner_2d, inner_2d, cartesian_product=True, nthreads=2
        )

        list_composed.extend(
            c.splinepy.splinepy_core.lists.composition_derivative(
                outer_3d,
                inner_3d,
                inner_3d,
                cartesian_product=True,
                nthreads=2,
            )
        )

        for r, l in zip(ref_composed, list_composed):
            assert c.are_splines_equal(c.to_derived(r), c.to_derived(l))


if __name__ == "__main__":
    c.unittest.main()
