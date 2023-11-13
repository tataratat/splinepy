import os
import re
import unittest

import gustaf as gus
import numpy as np

import splinepy

__all__ = [
    "unittest",
    "np",
    "splinepy",
    "gus",
]


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


def get_2d_control_points_nurbs():
    return [
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ]


def get_2d_control_points_bezier():
    return [
        [-1.0, 0.0],
        [-1.0, 1.0],
        [0.0, 1.0],
        [-2.0, 0.0],
        [-2.0, 2.0],
        [0.0, 2.0],
    ]


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


def get_knotvectors_2():
    return [
        [0, 0, 0, 0.5, 1, 1, 1],
        [0, 0, 0, 1, 1, 1],
    ]


def get_knotvectors_3():
    return [
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
        [0.0, 0.0, 1.0, 1.0],
    ]


# query points
def get_queries_2D():
    return [
        [0.01, 0.01],
        [0.01, 0.5],
        [0.9, 0.1],
        [0.8, 0.7],
        [0.4, 0.99],
    ]


def get_queries_3D():
    return [
        [0.1, 0.1, 0.1],
        [0.734, 0.525, 0.143],
        [0.9666, 0.991, 0.003],
        [0.5623, 0.0089, 0.99],
        [0.0431, 0.2, 0.523],
    ]


def dict_bspline_2p2d():
    return {
        "degrees": [2, 2],
        "knot_vectors": get_knotvectors_2(),
        "control_points": get_2d_control_points_b_spline(),
    }


def dict_nurbs_2p2d():
    return {
        "degrees": [2, 1],
        "knot_vectors": [
            [0, 0, 0, 1, 1, 1],
            [0, 0, 1, 1],
        ],
        "control_points": get_2d_control_points_nurbs(),
        "weights": [
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    }


#
def dict_bezier_2p2d():
    return {"degrees": [2, 1], "control_points": get_2d_control_points_nurbs()}


#
def dict_rational_bezier_2p2d():
    return {
        "degrees": [2, 1],
        "control_points": get_2d_control_points_nurbs(),
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
def dict_bezier_3p3d():
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points(),
    }


def dict_rational_bezier_3p3d():
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points(),
        "weights": [1.0] * len(get_3d_control_points()),
    }


def dict_bspline_3p3d():
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points(),
        "knot_vectors": get_knotvectors_3(),
    }


def dict_nurbs_3p3d():
    return {
        "degrees": [1, 1, 1],
        "control_points": get_3d_control_points(),
        "weights": [1.0] * len(get_3d_control_points()),
        "knot_vectors": get_knotvectors_2(),
    }


# initializing a spline should be a test itself, so provide `dict_spline`
# this is "iga-book"'s fig 2.15.
def bspline_2p2d():
    return splinepy.BSpline(
        degrees=[2, 2],
        knot_vectors=get_knotvectors_2(),
        control_points=get_2d_control_points_b_spline(),
    )


# half-half circle.
def nurbs_2p2d():
    return splinepy.NURBS(
        degrees=[2, 1],
        knot_vectors=[
            [0, 0, 0, 1, 1, 1],
            [0, 0, 1, 1],
        ],
        control_points=get_2d_control_points_nurbs(),
        weights=[
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    )


def nurbs_2p2d_quarter_circle():
    """explicit function for quarter circle
    in case n2p2d changes in the future..."""
    return dict_nurbs_2p2d()


def bezier_2p2d():
    return splinepy.Bezier(
        degrees=[2, 1],
        control_points=get_2d_control_points_nurbs(),
    )


def rational_bezier_2p2d():
    return splinepy.RationalBezier(
        degrees=[2, 1],
        control_points=get_2d_control_points_bezier(),
        weights=[
            [1.0],
            [2**-0.5],
            [1.0],
            [1.0],
            [2**-0.5],
            [1.0],
        ],
    )


def bezier_3p3d():
    return splinepy.Bezier(
        degrees=[1, 1, 1], control_points=get_3d_control_points()
    )


def rational_bezier_3p3d():
    return splinepy.RationalBezier(
        degrees=[1, 1, 1],
        control_points=get_3d_control_points(),
        weights=[1.0] * len(get_3d_control_points()),
    )


def bspline_3p3d():
    return splinepy.BSpline(
        degrees=[1, 1, 1],
        control_points=get_3d_control_points(),
        knot_vectors=get_knotvectors_3(),
    )


def nurbs_3p3d():
    return splinepy.NURBS(
        degrees=[1, 1, 1],
        control_points=get_3d_control_points(),
        weights=[1.0] * len(get_3d_control_points()),
        knot_vectors=get_knotvectors_3(),
    )


def raster(bounds, resolutions):
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


def nd_box(dim):
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


def to_tmpf(tmpd):
    """given tmpd, returns tmpf"""
    return os.path.join(tmpd, "nqv248p90")


def are_splines_equal(a, b, print_=False):
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


def are_items_close(a, b):
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


def are_items_same(a, b):
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


def are_stripped_lines_same(a, b, ignore_order=False):
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
            splitted_a = list(filter(None, re.split(delimiters, stripped_a)))
            splitted_b = list(filter(None, re.split(delimiters, stripped_b)))
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


def spline_types_as_list():
    return [
        bspline_2p2d(),
        nurbs_2p2d(),
        bezier_2p2d(),
        rational_bezier_2p2d(),
    ]


def get_all_spline_types_empty_as_list():
    return [
        splinepy.BSpline(),
        splinepy.NURBS(),
        splinepy.Bezier(),
        splinepy.RationalBezier(),
    ]


def get_spline_dictionaries():
    return [
        dict_bspline_2p2d(),
        dict_nurbs_2p2d(),
        dict_bezier_2p2d(),
        dict_rational_bezier_2p2d(),
    ]


class SplineBasedTestCase(unittest.TestCase):
    def bezier_2p2d(self):
        return bezier_2p2d()

    def rational_bezier_2p2d(self):
        return rational_bezier_2p2d()

    def bspline_2p2d(self):
        return bspline_2p2d()

    def nurbs_2p2d(self):
        return nurbs_2p2d()

    def all_2p2d_splines(self):
        return (
            self.bezier_2p2d(),
            self.rational_bezier_2p2d(),
            self.bspline_2p2d(),
            self.nurbs_2p2d(),
        )

    def bezier_3p3d(self):
        return bezier_3p3d()

    def rational_bezier_3p3d(self):
        return rational_bezier_3p3d()

    def bspline_3p3d(self):
        return bspline_3p3d()

    def nurbs_3p3d(self):
        return nurbs_3p3d()

    def all_3p3d_splines(self):
        return (
            self.bezier_3p3d(),
            self.rational_bezier_3p3d(),
            self.bspline_3p3d(),
            self.nurbs_3p3d(),
        )
