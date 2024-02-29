from _typeshed import Incomplete

from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase

class Cross3D(_TileBase):
    def __init__(self) -> None: ...
    def create_tile(
        self,
        parameters: Incomplete | None = None,
        parameter_sensitivities: Incomplete | None = None,
        center_expansion: float = 1.0,
        closure: Incomplete | None = None,
        **kwargs,
    ): ...
