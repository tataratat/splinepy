"""splinepy/microstructure/tiles/__init__.py.

Interface for tools and generators creating simple microstructures.
"""

# required in TileLib
from splinepy._base import SplinepyBase as _SplinepyBase
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
]


class TileLib(_SplinepyBase):
    """Class that provides easier access available tiles.
    All the member functions are classmethods, meaning you can access them
    directly through the class.
    """

    __slots__ = ()

    _tile_types = None

    # physical dimensions
    _1d = None
    _2d = None
    _3d = None

    @classmethod
    def _summarize(cls):
        """
        Get and save tile types
        """
        if cls._tile_types is not None:
            return None

        import inspect

        cls._tile_types = {}
        cls._1d = {}
        cls._2d = {}
        cls._3d = {}
        for a in __all__:
            # eval to turn str into object
            obj = eval(a)

            # if not class, not interested
            if not inspect.isclass(obj):
                continue

            # if TileBase, also not interested
            if obj is TileBase:
                continue

            # is class, and subclass of TileBase, but not TileBase
            if issubclass(obj, TileBase):
                # save types and sort with direction
                cls._tile_types[a] = obj
                if obj.dim == 1:
                    cls._1d[a] = obj
                elif obj.dim == 2:
                    cls._2d[a] = obj
                elif obj.dim == 3:
                    cls._3d[a] = obj

    @classmethod
    def by_dim(cls, para_dim=None, dim=None):
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
        cls._summarize()

        # initialize pool of tiles
        pool = None
        if dim is None:
            pool = cls._tile_types.copy()
        else:
            # cast int and copy
            pool = getattr(cls, f"_{int(dim)}d").copy()

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
            cls._loge(
                "Tiles does not exist with given dimension combination - ",
                f"(para_dim {para_dim}, dim {dim})",
            )

        return pool

    @classmethod
    def everything(cls):
        """
        Returns all available
        """
        cls._summarize()
        return cls._tile_types.copy()

    @classmethod
    def show(cls, **kwargs):
        """
        Shows name and default tile (with default parameter values).
        """
        from splinepy import Multipatch, show

        cls._summarize()

        to_show = []
        for key, value in cls._tile_types.items():
            to_show.append([key, Multipatch(value().create_tile()[0])])

        # turn off control points if kwargs doesn't have it
        if not any(k.startswith("control") for k in kwargs):
            kwargs["control_points"] = False

        show(*to_show, **kwargs)

    @classmethod
    def __getitem__(cls, key):
        cls._summarize()

        if key in cls._tile_types:
            # return initialized as they don't expect any variables
            return cls._tile_types[key]()

        raise KeyError(f"{key}-tile does not exist.")
