from _typeshed import Incomplete

from splinepy.microstructure.tiles.tile_base import TileBase as _TileBase

class Armadillo(_TileBase):
    def __init__(self) -> None: ...
    def closing_tile(
        self,
        parameters: Incomplete | None = None,
        parameter_sensitivities: Incomplete | None = None,
        contact_length: float = 0.3,
        closure: Incomplete | None = None,
        **kwargs,
    ): ...
    def create_tile(
        self,
        parameters: Incomplete | None = None,
        parameter_sensitivities: Incomplete | None = None,
        contact_length: float = 0.3,
        closure: Incomplete | None = None,
        **kwargs,
    ): ...
