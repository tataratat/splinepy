from _typeshed import Incomplete

from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase

class Snappy(_TileBase):
    def __init__(self) -> None: ...
    def create_tile(
        self,
        parameters: Incomplete | None = None,
        parameter_sensitivities: Incomplete | None = None,
        contact_length: float = 0.1,
        a: float = 0.1,
        b: float = 0.2,
        c: float = 0.3,
        r: float = 0.15,
        closure: Incomplete | None = None,
        **kwargs,
    ): ...
