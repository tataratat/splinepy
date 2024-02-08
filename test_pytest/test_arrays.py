import gustaf as gus
import numpy as np

import splinepy


def test_uniform_query():
    """
    Test uniform query against gustaf's vertice creation
    """
    for i in range(1, 11):
        l_bound = np.random.random(i)
        bounds = [l_bound, l_bound + 1]
        res = np.random.randint([3] * i) + 2  # minimum 2
        gus_q = gus.create.vertices.raster(bounds, res).vertices
        spp_q = splinepy.utils.data.uniform_query(bounds, res)

        assert np.allclose(gus_q, spp_q)
