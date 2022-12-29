try:
    from . import common as c
except BaseException:
    import common as c


class ContiguousArrayInputTest(c.unittest.TestCase):
    def test_c_contiguous_array_input(self):
        """
        Test whether cpp sides receives contiguous array
        """
        spl_ts = (
            c.splinepy.Bezier,
            c.splinepy.RationalBezier,
            c.splinepy.BSpline,
            c.splinepy.NURBS,
        )
        spl_props = (c.z2P2D, c.r2P2D, c.b2P2D, c.n2P2D)

        for spline_t, props in zip(spl_ts, spl_props):
            p = props.copy()  # create shallow copy of the dict

            # make f contiguous
            f_contig_orig = c.np.asarray(p["control_points"], order="F")
            assert not f_contig_orig.flags["C_CONTIGUOUS"]
            assert f_contig_orig.flags["F_CONTIGUOUS"]

            # set non contiguous array as prop
            p["control_points"] = f_contig_orig

            # initialize spl
            spl = spline_t(**props)
            # check round_trip cps
            round_trip_cps = spl.cps

            # they should be the same
            assert c.np.allclose(round_trip_cps, props["control_points"])

            # test setter
            # make sure this array is still not c contiguous
            assert not f_contig_orig.flags["C_CONTIGUOUS"]

            spl.cps = f_contig_orig  # setter

            assert spl.cps.flags["C_CONTIGUOUS"]
            assert c.np.allclose(spl.cps, props["control_points"])


if __name__ == "__main__":
    c.unittest.main()
