import numpy as _np

from splinepy._base import SplinepyBase as _SplinepyBase


class TileBase(_SplinepyBase):
    """
    Base class for tile objects
    """

    _dim = None
    _para_dim = None

    def __init__(self):
        """
        Init Values to None
        """
        self._evaluation_points = None
        self._n_info_per_eval_point = None

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

    @classmethod
    @property
    def dim(cls):
        """Returns dimensionality in physical space of the Microtile.

        Parameters
        ----------
        None

        Returns
        -------
        dim : int
        """
        if cls._dim is None:
            raise TypeError(
                "Inherited Tile-types need to provide _dim, see documentation."
            )
        return cls._dim

    @classmethod
    @property
    def para_dim(cls):
        """Returns dimensionality in parametric space of the Microtile.

        Parameters
        ----------
        None

        Returns
        -------
        para_dim : int
        """
        if cls._para_dim is None:
            raise TypeError(
                "Inherited Tile-types need to provide _para_dim, see documentation."
            )
        return cls._para_dim

    def check_params(self, params):
        """Checks if the parameters have the correct format and shape

        Parameters
        ----------
        params: np.ndarray
            the parameters for the tile

        Returns
        -------
        True: Boolean
        """
        # check if tuple

        if not (isinstance(params, _np.ndarray) and params.ndim == 2):
            raise TypeError("parameters must be two-dimensional np array")

        if not (
            (self._evaluation_points.shape[0] == params.shape[0])
            and (self._n_info_per_eval_point == params.shape[1])
        ):
            raise TypeError(
                f"Mismatch in parameter size, expected "
                f"{self._evaluation_points.shape[0]} x "
                f"{self._n_info_per_eval_point}"
            )

        return True

    def check_param_derivatives(self, derivatives):
        """Checks if all derivatives have the correct format and shape

        Parameters
        ----------
        derivatives: np.ndarray
            all derivatives as 3D array

        Returns
        -------
        True: Boolean
        """
        if derivatives is None:
            return False

        if not (
            isinstance(derivatives, _np.ndarray) and derivatives.ndim == 3
        ):
            raise TypeError("parameters must be three-dimensional np array")

        if not (
            (self._evaluation_points.shape[0] == derivatives.shape[0])
            and (self._n_info_per_eval_point == derivatives.shape[1])
        ):
            raise TypeError(
                f"Mismatch in parameter size, expected "
                f"{self._evaluation_points.shape[0]} x "
                f"{self._n_info_per_eval_point} x n_sensitivities"
            )

        return True

    def create_tile(self, **kwargs):
        """Tile creation interface for derived classes.
        All types should implement this.

        Parameters
        ----------
        kwargs: **kwargs

        Returns
        -------
        tiles: list
          List of splines that create a tile
        derivatives: list<list<Spline>> or None
          List of list of splines that represents parameter sensitivities.
          If it is not implemented, returns None.
        """
        raise NotImplementedError(
            f"create_tile() not implemented for {type(self)}"
        )
