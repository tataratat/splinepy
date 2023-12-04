try:
    from . import common as c
except BaseException:
    import common as c


class KnotVectorsTest(c.SplineBasedTestCase):
    def test_knot_vectors(self):
        for spline in (
            self.bspline_2p2d(),
            self.bspline_3p3d(),
            self.nurbs_2p2d(),
            self.nurbs_3p3d(),
        ):
            # test unique_knots
            copy_knot_vectors = spline.knot_vectors[:]
            unique_knots = [c.np.unique(ckvs) for ckvs in copy_knot_vectors]

            for uk, uk_fct in zip(unique_knots, spline.unique_knots):
                assert c.np.allclose(uk, uk_fct)

            # test knot_multiplicities
            multiplicity = [
                c.np.unique(ckvs, return_counts=True)
                for ckvs in copy_knot_vectors
            ]

            for m_u, m_fct in zip(multiplicity, spline.knot_multiplicities):
                assert c.np.allclose(m_u[1], m_fct)

            # test knot_vector creation
            for u_kv, kn_m, spl_kv in zip(
                spline.unique_knots,
                spline.knot_multiplicities,
                spline.knot_vectors,
            ):
                assert c.np.allclose(
                    c.np.array(spl_kv), c.np.repeat(u_kv, kn_m)
                )


if __name__ == "__main__":
    c.unittest.main()
