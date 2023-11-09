import platform
import tempfile

try:
    from . import common as c
except BaseException:
    import common as c


class IGESExportTest(c.unittest.TestCase):
    def test_iges_round_trip_bsplines(self):
        """
        Test iges export-import routine
        """
        # temporarily disable windows - debug test (help wanted!)
        if platform.system().startswith(
            "Win"
        ) and c.splinepy.splinepy_core.build_type().startswith("debug"):
            c.splinepy.utils.log.warning(
                "Skipping io.iges test.",
                "Try release mode, if you need to use io.iges",
            )
            return None

        bsp_el2 = c.splinepy.BSpline(
            degrees=[1, 1],
            control_points=[[0, 1], [1, 1], [0, 2], [1, 2]],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )
        nur_el3 = c.splinepy.NURBS(
            degrees=[1, 1],
            control_points=[[1, 1], [2, 1], [1, 2], [2, 2]],
            weights=[1, 1, 1, 1],
            knot_vectors=[[0, 0, 1, 1], [0, 0, 1, 1]],
        )

        # Make it more tricky
        bsp_el2.elevate_degrees(0)
        bsp_el2.elevate_degrees(1)
        nur_el3.elevate_degrees(0)
        nur_el3.elevate_degrees(1)
        bsp_el2.insert_knots(1, [0.5])
        nur_el3.insert_knots(1, [0.5])

        # Test output against input
        list_of_splines = [bsp_el2, nur_el3]
        with tempfile.TemporaryDirectory() as tmpd:
            tmpf = c.to_tmpf(tmpd)
            c.splinepy.io.iges.export(tmpf, list_of_splines)
            list_of_splines_loaded = c.splinepy.io.iges.load(tmpf)

            # to compare we need 3d "bumped" splines of originals
            splines_in_3d = []
            for s in list_of_splines:
                tmp_s = s.copy()
                tmp_s.cps = c.np.hstack(
                    [tmp_s.cps, c.np.zeros((len(tmp_s.cps), 1))]
                )
                splines_in_3d.append(tmp_s)

            assert all(
                c.are_splines_equal(a, b)
                for a, b in zip(splines_in_3d, list_of_splines_loaded)
            )


if __name__ == "__main__":
    c.unittest.main()
