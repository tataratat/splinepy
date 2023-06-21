from splinepy._base import SplinepyBase


class TileBase(SplinepyBase):
    """
    Base class for tile objects
    """

    def __init__(self):
        """
        Init Values to None
        """
        self._dim = None
        self._evaluation_points = None
        self._parameter_space_dimension = None

    @property
    def parameter_space_dimension(self):
        """Number of parameters per evaluation point."""
        if self._parameter_space_dimension is None:
            raise TypeError(
                "Inherited Tile-types need to provide "
                "_parameter_space_dimension, see documentation."
            )
        return self._parameter_space_dimension

    @property
    def evaluation_points(self):
        """Positions in the parametrization function to be evaluated when tile
        " "is constructed prior to composition.

        Parameters
        ----------
        None

        Returns
        -------
        evaluation_points : np.ndarray(6,3)
        """
        if self._evaluation_points is None:
            raise TypeError(
                "Inherited Tile-types need to provide _evaluation_points, see "
                "documentation."
            )
        return self._evaluation_points

    @property
    def dim(self):
        """Returns dimensionality in physical space of the Microtile.

        Parameters
        ----------
        None

        Returns
        -------
        dim : int
        """
        if self._dim is None:
            raise TypeError(
                "Inherited Tile-types need to provide _dim, see documentation."
            )
        return self._dim
