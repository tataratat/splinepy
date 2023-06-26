try:
    from . import common as c
except BaseException:
    import common as c


class ContiguousArrayInputTest(c.unittest.TestCase):
    def test_c_contiguous_array_input(self):
        """
        Test whether cpp sides receives contiguous array
        """

        for spl, props in zip(
            c.get_all_spline_typs_as_list(),
            c.get_all_splines_as_dict_as_list(),
        ):
            # make f contiguous
            f_contig_orig = c.np.asarray(props["control_points"], order="F")
            assert not f_contig_orig.flags["C_CONTIGUOUS"]
            assert f_contig_orig.flags["F_CONTIGUOUS"]

            # set non contiguous array as prop
            props["control_points"] = f_contig_orig

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
