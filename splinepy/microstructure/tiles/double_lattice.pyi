from _typeshed import Incomplete

from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase

class DoubleLattice(_TileBase):
    def __init__(self) -> None: ...
    def create_tile(
        self,
        parameters: Incomplete | None = None,
        parameter_sensitivities: Incomplete | None = None,
        contact_length: float = 0.5,
        **kwargs,
    ): ...
