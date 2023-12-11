try:
    from . import common as c
except BaseException:
    import common as c


class ArrayTest(c.unittest.TestCase):
    def test_uniform_query(self):
        """
        Test uniform query against gustaf's vertice creation
        """
        for i in range(1, 11):
            l_bound = c.np.random.random(i)
            bounds = [l_bound, l_bound + 1]
            res = c.np.random.randint([3] * i) + 2  # minimum 2
            gus_q = c.gus.create.vertices.raster(bounds, res).vertices
            spp_q = c.splinepy.utils.data.uniform_query(bounds, res)

            assert c.np.allclose(gus_q, spp_q)


if __name__ == "__main__":
    c.unittest.main()
