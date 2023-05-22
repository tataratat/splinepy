try:
    from . import common as c
except BaseException:
    import common as c


class SplineListTest(c.unittest.TestCase):
    def test_spline_list(self):
        """
        test SplineList and some cpp vector manipulation calls.
        std::vector wrapping is provided by pybind,
        so, minimal test to make sure it works in our context
        """
        # Define some splines
        splines = c.all2p2d()
        b_id = splines.index(splines[2])  # this should be b_spline

        # short cut
        SplineList = c.splinepy.splinepy_core.SplineList

        # test if creating a spline is at the end the same thing
        # 0. append
        splist0 = SplineList()
        splist0.append(splines[0])
        splist0.append(splines[1])
        splist0.append(splines[2])
        splist0.append(splines[3])

        # 1. extend
        splist1 = SplineList()
        splist1.extend(splines)

        # 2. iterable init - list
        splist2 = SplineList(splines)

        # 3. iterable init - tuple
        splist3 = SplineList(tuple(splines))

        # are they the same?
        assert splist0 == splist1 == splist2 == splist3

        # are they really ref?
        # edit cps
        splines[b_id].cps[0] += 0.5
        for splist in [splist0, splist1, splist2, splist3]:
            # if they are ref, this should be modified
            assert splist[b_id].cps._modified

        # test resize
        splist4 = SplineList()

        # should fill with nullptr
        splist4.resize(len(splines))
        assert len(splist4) == len(splines)

        for s4 in splist4:
            # thankfully nullptr is casted to None
            assert s4 is None

        for i, s in enumerate(splines):
            splist4[i] = s
        assert splist0 == splist4

    def test_spline_list_evaluate(self):
        """list evaluation"""
        # get list of splines
        splines = c.all2p2d()

        # create SplineList
        slist = c.splinepy.splinepy_core.SplineList(splines)

        # get reference solution
        ref_solutions = c.np.vstack([s.evaluate(c.q2D) for s in splines])

        # call eval, don't check dim
        list_solutions = slist.evaluate(c.q2D, nthreads=2, check_dims=False)
        assert c.np.allclose(ref_solutions, list_solutions)

        # call eval, check dim
        list_solutions = slist.evaluate(c.q2D, nthreads=2, check_dims=True)
        assert c.np.allclose(ref_solutions, list_solutions)

        # call eval, single thread
        list_solutions = slist.evaluate(c.q2D, nthreads=1, check_dims=True)
        assert c.np.allclose(ref_solutions, list_solutions)

    def test_spline_list_sample(self):
        """spline sample"""
        splines = c.all2p2d()
        slist = c.splinepy.splinepy_core.SplineList(splines)

        # set resolution. this is int,
        # because it is implemented for equal resolution only
        resolution = 3

        ref_solutions = c.np.vstack(
            [s.sample([resolution] * splines[0].para_dim) for s in splines]
        )

        # check everything and compute
        list_samples = slist.sample(
            resolution,
            nthreads=2,
            same_parametric_bounds=True,
            check_dims=True,
        )
        assert c.np.allclose(ref_solutions, list_samples)

        # don't check p bounds
        list_samples = slist.sample(
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
            new_kvs = [kv * 2 for kv in slist[i].kvs]
            slist[i].kvs = new_kvs

        # check p bounds,
        list_samples = slist.sample(
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
        slist2d = c.splinepy.splinepy_core.SplineList(splines2d)

        # prepare mixed dim splines
        splines2d3d = [*splines2d, *c.all3p3d()]
        slist2d3d = c.splinepy.splinepy_core.SplineList(splines2d3d)

        # prepare test func
        def _test(pure_list, spline_list, same_p_dims):
            """actual test routine"""
            ref_boundaries = []
            for s in pure_list:
                ref_boundaries.extend(
                    c.splinepy.splinepy_core.extract_boundaries(s, [])
                )

            list_boundaries = spline_list.extract_boundaries(
                nthreads=2, same_para_dims=False
            )

            # for same p_dims,
            # try same_para_dims=True by extending test splines
            if same_p_dims:
                ref_boundaries.extend(ref_boundaries)
                list_boundaries.extend(
                    spline_list.extract_boundaries(
                        nthreads=2, same_para_dims=True
                    )
                )

            # should have same size
            assert len(ref_boundaries) == len(list_boundaries)

            # same spline?
            for rb, lb in zip(ref_boundaries, list_boundaries):
                # currently, core function does not have getter
                # for individual attr, so cast
                assert c.are_splines_equal(c.to_derived(rb), c.to_derived(lb))

        _test(splines2d, slist2d, True)
        _test(splines2d3d, slist2d3d, False)

    def test_boundary_centers(self):
        """boundary centers"""
        splines = c.all2p2d()
        slist = c.splinepy.splinepy_core.SplineList(splines)

        ref_centers = c.np.vstack(
            [c.splinepy.splinepy_core.boundary_centers(s) for s in splines]
        )

        # test same bounds
        list_centers = slist.boundary_centers(
            nthreads=2, same_parametric_bounds=True, check_dims=False
        )
        assert c.np.allclose(ref_centers, list_centers)

        # different bounds
        # update p bounds for bsplines
        bspline_id = 2
        nurbs_id = 3
        for i in [bspline_id, nurbs_id]:
            new_kvs = [kv * 2 for kv in slist[i].kvs]
            slist[i].kvs = new_kvs

        list_centers = slist.boundary_centers(
            nthreads=2, same_parametric_bounds=False, check_dims=False
        )
        assert c.np.allclose(ref_centers, list_centers)

    def test_raise_dim_mismatch(self):
        """see if function raises dim mismatch correctly"""
        slist2d = c.splinepy.splinepy_core.SplineList(c.all2p2d())

        # prepare mixed
        slist2d3d = c.splinepy.splinepy_core.SplineList(
            [*c.all2p2d(), *c.all3p3d()]
        )

        # nothing happens
        slist2d.raise_dim_mismatch(nthreads=2)

        with self.assertRaises(RuntimeError):
            slist2d3d.raise_dim_mismatch(nthreads=2)


if __name__ == "__main__":
    c.unittest.main()
