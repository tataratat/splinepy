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
        z = c.splinepy.Bezier(**c.z2p2d())
        r = c.splinepy.RationalBezier(**c.r2p2d())
        b = c.splinepy.BSpline(**c.b2p2d())
        n = c.splinepy.NURBS(**c.n2p2d())
        splines = [z, r, b, n]
        b_id = splines.index(b)

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


if __name__ == "__main__":
    c.unittest.main()
