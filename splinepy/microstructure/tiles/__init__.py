"""splinepy/microstructure/tiles/__init__.py.

Interface for tools and generators creating simple microstructures.
"""

# required in TileLib
from splinepy.microstructure.tiles import (
    armadillo,
    chi,
    cross_2d,
    cross_3d_linear,
    cube_void,
    double_lattice,
    ellips_v_oid,
    hollow_cube,
    hollow_octagon,
    hollow_octagon_extrude,
    inverse_cross_3d,
    snappy,
    tile_base,
)
from splinepy.microstructure.tiles.armadillo import Armadillo
from splinepy.microstructure.tiles.chi import Chi
from splinepy.microstructure.tiles.cross_2d import Cross2D
from splinepy.microstructure.tiles.cross_3d import Cross3D
from splinepy.microstructure.tiles.cross_3d_linear import Cross3DLinear
from splinepy.microstructure.tiles.cube_void import CubeVoid
from splinepy.microstructure.tiles.double_lattice import DoubleLattice
from splinepy.microstructure.tiles.ellips_v_oid import EllipsVoid
from splinepy.microstructure.tiles.hollow_cube import HollowCube
from splinepy.microstructure.tiles.hollow_octagon import HollowOctagon
from splinepy.microstructure.tiles.hollow_octagon_extrude import (
    HollowOctagonExtrude,
)
from splinepy.microstructure.tiles.inverse_cross_3d import InverseCross3D
from splinepy.microstructure.tiles.snappy import Snappy
from splinepy.microstructure.tiles.tile_base import TileBase
from splinepy.utils import log as _log

__all__ = [
    "armadillo",
    "chi",
    "cross_2d",
    "cube_void",
    "cross_3d_linear",
    "double_lattice",
    "ellips_v_oid",
    "hollow_cube",
    "hollow_octagon",
    "hollow_octagon_extrude",
    "inverse_cross_3d",
    "snappy",
    "tile_base",
    "Armadillo",
    "Chi",
    "Cross2D",
    "Cross3D",
    "Cross3DLinear",
    "CubeVoid",
    "DoubleLattice",
    "EllipsVoid",
    "HollowCube",
    "HollowOctagon",
    "HollowOctagonExtrude",
    "InverseCross3D",
    "Snappy",
    "TileBase",
    "by_dim",
    "everything",
    "show",
    "get",
]


def _summarize_tiles():
    """
    Get and save tile types
    """
    tile_types = {}
    d1 = {}
    d2 = {}
    d3 = {}

    # loop subclasses and classify
    for SubClass in TileBase.__subclasses__():
        key = SubClass.__qualname__
        # save types and sort with direction
        tile_types[key] = SubClass
        dim = SubClass.dim
        if dim == 1:
            d1[key] = SubClass
        elif dim == 2:
            d2[key] = SubClass
        elif dim == 3:
            d3[key] = SubClass

    return tile_types, d1, d2, d3


_tile_types, _1d, _2d, _3d = _summarize_tiles()


def by_dim(para_dim=None, dim=None):
    """
    Returns names of tiles that satisfies dimension inputs.
    Per default, it returns list of all the name of available tiles.

    Parameters
    ----------
    para_dim: int
    dim: int

    Returns
    -------
    tiles: dict<str, type>
        dict of tile names and types that fulfills dimension inputs.
    """
    # initialize pool of tiles
    pool = _tile_types.copy() if dim is None else eval(f"_{int(dim)}d").copy()

    # filter
    if para_dim is not None:
        para_dim = int(para_dim)
        filtered = {}
        for key, value in pool.items():
            if value.para_dim == para_dim:
                filtered[key] = value

        # overwrite pool with filtered
        pool = filtered

    if len(pool) == 0:
        _log.error(
            "Tiles does not exist with given dimension combination - ",
            f"(para_dim {para_dim}, dim {dim})",
        )

    return pool


def everything():
    """
    Returns all predefined tiles.

    Parameters
    ----------
    None

    Returns
    -------
    tiles: dict<str, type>
    """
    return _tile_types.copy()


def show(**kwargs):
    """
    Shows name and default tile (with default parameter values).

    Parameters
    ----------
    **kwargs: kwargs
      show options that splinepy.show() accepts

    Returns
    -------
    show: Any
      show() return based on kwargs. Look gus.show() documentations.
    """
    from splinepy import Multipatch, show

    to_show = []
    for key, value in _tile_types.items():
        to_show.append([key, Multipatch(value().create_tile()[0])])

    # turn off control points if kwargs doesn't have it
    if not any(k.startswith("control") for k in kwargs):
        kwargs["control_points"] = False

    return show(*to_show, **kwargs)


def get(key):
    """
    Returns initialized tile object based that corresponds to given str key.

    Parameters
    ----------
    key: str
      Name of a tile in str

    Returns
    -------
    tile: TileBase
      Initialized, derived tile.
    """
    if key in _tile_types:
        # return initialized as they don't expect any variables
        return _tile_types[key]()

    raise KeyError(f"{key}-tile does not exist.")
