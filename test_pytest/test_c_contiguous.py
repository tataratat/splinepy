import numpy as np

def test_c_contiguous_array_input(request):
        """ Test whether cpp sides receives contiguous array """
        
        splinelist = request.getfixturevalue("spline_types_as_list")
        splinedict = request.getfixturevalue("get_spline_dictionaries")
        
        for spl, properties in zip(splinelist,splinedict,):
            # make f contiguous
            f_contig_orig = np.asarray(
                properties["control_points"], order="F"
            )
            assert not f_contig_orig.flags["C_CONTIGUOUS"]
            assert f_contig_orig.flags["F_CONTIGUOUS"]

            # set non contiguous array as prop
            properties["control_points"] = f_contig_orig

            # check round_trip cps
            round_trip_cps = spl.cps

            # they should be the same
            assert np.allclose(round_trip_cps, properties["control_points"])

            # test setter
            # make sure this array is still not c contiguous
            assert not f_contig_orig.flags["C_CONTIGUOUS"]

            spl.cps = f_contig_orig  # setter

            assert spl.cps.flags["C_CONTIGUOUS"]
            assert np.allclose(spl.cps, properties["control_points"])
