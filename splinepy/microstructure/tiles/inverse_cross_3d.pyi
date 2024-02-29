from _typeshed import Incomplete

from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase

class InverseCross3D(_TileBase):
    def __init__(self) -> None: ...
    def create_tile(
        self,
        parameters: Incomplete | None = None,
        parameter_sensitivities: Incomplete | None = None,
        seperator_distance: Incomplete | None = None,
        center_expansion: float = 1.0,
        closure: Incomplete | None = None,
        **kwargs,
    ): ...
