import os
import re

import gustaf as gus
import numpy as np
import pytest

import splinepy

__all__ = [
    "pytest",
    "np",
    "splinepy",
    "gus",
]


@pytest.fixture
def np_rng():
    return np.random.default_rng()


@pytest.fixture
def get_2d_control_points_b_spline():
    return [
        [0, 0],
        [0, 1],
        [1, 1.5],
        [3, 1.5],
        [-1, 0],
        [-1, 2],
        [1, 4],
        [3, 4],
        [-2, 0],
        [-2, 2],
        [1, 5],
        [3, 5],
    ]


@pytest.fixture
def get_2d_control_points_nurbs():
    return [
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ]


@pytest.fixture
def get_2d_control_points_bezier():
    return [
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ]


@pytest.fixture
def get_3d_control_points():
    return [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, -1.0, 1.0],
        [1.0, 0.0, 1.0],
        [-1.0, 1.0, 2.0],
        [2.0, 2.0, 2.0],
    ]


@pytest.fixture
def get_knotvectors_2():
    return [
        [0, 0, 0, 0.5, 1, 1, 1],
        [0, 0, 0, 1, 1, 1],
    ]


@pytest.fixture
def get_knotvectors_3():
    return [
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
    ]


# query points
@pytest.fixture
def get_queries_2D():
    return [
        [0.01, 0.01],
        [0.01, 0.5],
        [0.9, 0.1],
        [0.8, 0.7],
        [0.4, 0.99],
    ]


@pytest.fixture
def get_queries_3D():
    return [
        [0.1, 0.1, 0.1],
        [0.734, 0.525, 0.143],
        [0.9666, 0.991, 0.003],
        [0.5623, 0.0089, 0.99],
        [0.0431, 0.2, 0.523],
    ]


@pytest.fixture
def dict_bspline_2p2d(get_knotvectors_2, get_2d_control_points_b_spline):
    return {
        "degrees": [2, 2],
        "knot_vectors": get_knotvectors_2,
        "control_points": get_2d_control_points_b_spline,
    }


@pytest.fixture
def dict_nurbs_2p2d(get_2d_control_points_nurbs):
    return {
        "degrees": [2, 1],
        "knot_vectors": [
            [0, 0, 0, 1, 1, 1],
            [0, 0, 1, 1],
        ],
        "control_points": get_2d_control_points_nurbs,
        "weights": [
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    }


@pytest.fixture
def nurbs_2p2d_quarter_circle(dict_nurbs_2p2d):
    """explicit function for quarter circle
    in case n2p2d changes in the future..."""
    return dict_nurbs_2p2d


#
@pytest.fixture
def dict_bezier_2p2d(get_2d_control_points_nurbs):
    return {
        "degrees": [2, 1],
        "control_points": get_2d_control_points_nurbs,
    }


#
@pytest.fixture
def dict_rational_bezier_2p2d(get_2d_control_points_nurbs):
    return {
        "degrees": [2, 1],
        "control_points": get_2d_control_points_nurbs,
        "weights": [
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    }


# 3D
@pytest.fixture
def dict_bezier_3p3d(get_3d_control_points):
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points,
    }


@pytest.fixture
def dict_rational_bezier_3p3d(get_3d_control_points):
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points,
        "weights": [1.0] * len(get_3d_control_points),
    }


@pytest.fixture
def dict_bspline_3p3d(get_3d_control_points, get_knotvectors_3):
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points,
        "knot_vectors": get_knotvectors_3,
    }


@pytest.fixture
def dict_nurbs_3p3d(get_knotvectors_3, get_3d_control_points):
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points,
        "weights": [1.0] * len(get_3d_control_points),
        "knot_vectors": get_knotvectors_3,
    }


# initializing a spline should be a test itself, so provide `dict_spline`
# this is "iga-book"'s fig 2.15.
@pytest.fixture
def bspline_2p2d(get_knotvectors_2, get_2d_control_points_b_spline):
    return splinepy.BSpline(
        degrees=[2, 2],
        knot_vectors=get_knotvectors_2,
        control_points=get_2d_control_points_b_spline,
    )


# half-half circle.
@pytest.fixture
def nurbs_2p2d(get_2d_control_points_nurbs):
    return splinepy.NURBS(
        degrees=[2, 1],
        knot_vectors=[
            [0, 0, 0, 1, 1, 1],
            [0, 0, 1, 1],
        ],
        control_points=get_2d_control_points_nurbs,
        weights=[
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    )


@pytest.fixture
def bezier_2p2d(get_2d_control_points_nurbs):
    return splinepy.Bezier(
        degrees=[2, 1],
        control_points=get_2d_control_points_nurbs,
    )


@pytest.fixture
def rational_bezier_2p2d(get_2d_control_points_bezier):
    return splinepy.RationalBezier(
        degrees=[2, 1],
        control_points=get_2d_control_points_bezier,
        weights=[
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    )


@pytest.fixture
def bezier_3p3d(get_3d_control_points):
    return splinepy.Bezier(
        degrees=[1, 1, 1], control_points=get_3d_control_points
    )


@pytest.fixture
def rational_bezier_3p3d(get_3d_control_points):
    return splinepy.RationalBezier(
        degrees=[1, 1, 1],
        control_points=get_3d_control_points,
        weights=[1.0] * len(get_3d_control_points),
    )


@pytest.fixture
def bspline_3p3d(get_3d_control_points, get_knotvectors_3):
    return splinepy.BSpline(
        degrees=[1, 1, 1],
        control_points=get_3d_control_points,
        knot_vectors=get_knotvectors_3,
    )


@pytest.fixture
def nurbs_3p3d(get_3d_control_points, get_knotvectors_3):
    return splinepy.NURBS(
        degrees=[1, 1, 1],
        control_points=get_3d_control_points,
        weights=[1.0] * len(get_3d_control_points),
        knot_vectors=get_knotvectors_3,
    )


@pytest.fixture
def raster():
    def _raster(bounds, resolutions):
        """prepares raster points using np.meshgrid"""
        l_bounds, u_bounds = bounds[0], bounds[1]
        pts = np.meshgrid(
            *[
                np.linspace(lo, up, re)
                for lo, up, re in zip(l_bounds, u_bounds, resolutions)
            ],
            indexing="ij",
        )
        # return pts
        return np.hstack([p.reshape(-1, 1) for p in pts[::-1]])

    return _raster


@pytest.fixture
def nd_box(raster):
    def _nd_box(dim):
        """creates simple box in n-d"""
        ds = [1 for _ in range(dim)]
        cps = raster(
            [[0 for _ in range(dim)], [1 for _ in range(dim)]],
            [2 for _ in range(dim)],
        )
        kvs = [[0, 0, 1, 1] for _ in range(dim)]
        ws = np.ones((len(cps), 1))
        return {
            "degrees": ds,
            "control_points": cps,
            "knot_vectors": kvs,
            "weights": ws,
        }

    return _nd_box


@pytest.fixture
def to_tmpf():
    def _to_tmpf(tmpd):
        """given tmpd, returns tmpf"""
        return os.path.join(tmpd, "nqv248p90")

    return _to_tmpf


@pytest.fixture
def are_splines_equal():
    def _are_splines_equal(a, b, print_=False):
        """returns True if Splines are equivalent"""
        if a.whatami != b.whatami:
            return False
        for req_prop in a.required_properties:
            if req_prop == "knot_vectors":
                for aa, bb in zip(a.knot_vectors, b.knot_vectors):
                    if not np.allclose(aa.numpy(), bb.numpy()):
                        if print_:
                            print("a.kvs", a.kvs)
                            print("b.kvs", b.kvs)
                        return False
            elif not np.allclose(getattr(a, req_prop), getattr(b, req_prop)):
                if print_:
                    print(f"a.{req_prop}", getattr(a, req_prop))
                    print(f"b.{req_prop}", getattr(b, req_prop))
                return False
        return True

    return _are_splines_equal


@pytest.fixture
def are_items_close():
    def _are_items_close(a, b):
        """returns True if items in a and b are close"""
        all_close = True

        for i, (aa, bb) in enumerate(zip(a, b)):
            this_is_close = all(np.isclose(aa, bb))
            if not this_is_close:
                # print to inform
                print(f"elements in index-{i} are not close")
                print(f"  from first: {aa}")
                print(f"  from second: {bb}")

                all_close = False

        return all_close

    return _are_items_close


@pytest.fixture
def are_items_same():
    def _are_items_same(a, b):
        """returns True if items in a and b are same"""
        all_same = True

        for i, (aa, bb) in enumerate(zip(a, b)):
            this_is_same = aa == bb
            if not this_is_same:
                # print to inform
                print(f"element in index-{i} are not same")
                print(f"  from first: {aa}")
                print(f"  from second: {bb}")

                all_same = False

        return all_same

    return _are_items_same


@pytest.fixture
def are_stripped_lines_same():
    def _are_stripped_lines_same(a, b, ignore_order=False):
        """returns True if items in a and b same, preceding and tailing whitespaces
        are ignored and strings are joined"""
        all_same = True

        for i, (line_a, line_b) in enumerate(zip(a, b)):
            # check stripped string
            stripped_a, stripped_b = line_a.strip(), line_b.strip()
            this_is_same = stripped_a == stripped_b

            # print general info
            if not this_is_same:
                print(f"stripped line at index-{i} are not the same")
                print(f"  from first: {line_a}")
                print(f"  from second: {line_b}")

            # give one more chance if ignore_order
            if not this_is_same and ignore_order:
                print("  checking again, while ignoring word order:")

                # This is meant for attributes
                delimiters = r" |\>|\<|\t"
                splitted_a = list(
                    filter(None, re.split(delimiters, stripped_a))
                )
                splitted_b = list(
                    filter(None, re.split(delimiters, stripped_b))
                )
                # first, len check
                len_a, len_b = len(splitted_a), len(splitted_b)
                if len(splitted_a) != len(splitted_b):
                    print(f"    different word counts: a-{len_a}, b-{len_b}")
                    all_same = False
                else:
                    # word order
                    a_to_b = []
                    for word_a in splitted_a:
                        try:
                            a_to_b.append(splitted_b.index(word_a))
                        except BaseException:
                            print(f"    second does not contain ({word_a})")
                            all_same = False

        return all_same

    return _are_stripped_lines_same


@pytest.fixture
def spline_types_as_list(
    bspline_2p2d, nurbs_2p2d, bezier_2p2d, rational_bezier_2p2d
):
    return [bspline_2p2d, nurbs_2p2d, bezier_2p2d, rational_bezier_2p2d]


@pytest.fixture
def get_all_spline_types_empty_as_list():
    return [
        splinepy.BSpline(),
        splinepy.NURBS(),
        splinepy.Bezier(),
        splinepy.RationalBezier(),
    ]


@pytest.fixture
def get_spline_dictionaries(
    dict_bspline_2p2d,
    dict_nurbs_2p2d,
    dict_bezier_2p2d,
    dict_rational_bezier_2p2d,
):
    return [
        dict_bspline_2p2d,
        dict_nurbs_2p2d,
        dict_bezier_2p2d,
        dict_rational_bezier_2p2d,
    ]


# @pytest.fixture
# def all_2p2d_splines(self):
#     return (
#         self.bezier_2p2d(),
#         self.rational_bezier_2p2d(),
#         self.bspline_2p2d(),
#         self.nurbs_2p2d(),
#     )


# @pytest.fixture
# def all_3p3d_splines(self):
#     return (
#         self.bezier_3p3d(),
#         self.rational_bezier_3p3d(),
#         self.bspline_3p3d(),
#         self.nurbs_3p3d(),
#     )
