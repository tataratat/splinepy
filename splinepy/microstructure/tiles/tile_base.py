import numpy as _np

from splinepy._base import SplinepyBase as _SplinepyBase


class TileBase(_SplinepyBase):
    """
    Base class for tile objects

    Attributes
    ---------------
    _dim: int
        Dimension in physical space
    _para_dim: int
        Dimension in parametric space
    _evaluation_points: np.ndarray (2D)
        Points in parametric space where tile parameters are evaluated. Each parameter
        is attributed one evaluation point. This is used for the parametrization spline.
    _n_info_per_eval_point: int
        Number of tile parameters per evaluation point
    _sensitivities_implemented: bool
        Whether sensitivities w.r.t. tile parameters are implemented
    _closure_directions: list<str>
        List of directions in which the closure has been implemented
    _parameter_bounds: list<list<float>>
        List of bounds for the tile parameters
    _parameters_shape: tuple<int>
        Shape of parameters array
    _default_parameter_value: float/np.ndarray
        Default values for all tile parameters
    """

    _dim = None
    _para_dim = None
    _evaluation_points = None
    _n_info_per_eval_point = None
    _sensitivities_implemented = None
    _closure_directions = None
    _parameter_bounds = None
    _parameters_shape = None
    _default_parameter_value = None

    def __init__(self):
        if type(self) is TileBase:
            raise TypeError(
                "Can't initialize a TileBase class! Please use one of the "
                f"derived classes: {TileBase.__subclasses__()}"
            )

    @classmethod
    def _raise_if_not_set_else_return(cls, attr_name):
        attr = getattr(cls, attr_name, None)
        if attr is None and cls is not TileBase:
            raise NotImplementedError(
                f"Inherited Tile-types need to provide {attr_name}, see "
                "documentation."
            )
        return attr

    @property
    def evaluation_points(cls):
        """Positions in the parametrization function to be evaluated when tile
        is constructed prior to composition.

        Parameters
        ----------
        None

        Returns
        -------
        evaluation_points : np.ndarray(6,3)
        """
        return cls._raise_if_not_set_else_return("_evaluation_points")

    @property
    def n_info_per_eval_point(cls):
        """Number of parameters per evaluation point

        Parameters
        ----------
        None

        Returns
        -------
        n_info : int
        """
        return cls._raise_if_not_set_else_return("_n_info_per_eval_point")

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
        return cls._raise_if_not_set_else_return("_dim")

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
        return cls._raise_if_not_set_else_return("_para_dim")

    @property
    def sensitivities_implemented(cls):
        """Returns whether sensitivities are implemented for the microtile

        Parameters
        ----------
        None

        Returns
        -------
        is_implemented: bool
        """
        return cls._raise_if_not_set_else_return("_sensitivities_implemented")

    @property
    def closure_directions(self):
        """Returns the available closure directions of the microtile

        Parameters
        ----------
        None

        Returns
        -------
        directions: None/list<str>
        """
        return self._closure_directions

    @property
    def parameter_bounds(cls):
        """Returns the bounds for the microtiles' parameters

        Parameters
        ----------
        None

        Returns
        -------
        bounds: list<list<float>>
        """
        return cls._raise_if_not_set_else_return("_parameter_bounds")

    @property
    def parameters_shape(cls):
        """Returns the shape of the microtile's parameters array

        Parameters
        ----------
        None

        Returns
        -------
        shape: tuple<int>
        """
        return cls._raise_if_not_set_else_return("_parameters_shape")

    @property
    def default_parameter_value(cls):
        """Returns the default value of the microtile's parameters

        Parameters
        ----------
        None

        Returns
        -------
        default_value: float/None
        """
        return cls._default_parameter_value

    def check_params(self, parameters):
        """Checks if the parameters have the correct format and shape and are within
        defined bounds

        Parameters
        ----------
        params: np.ndarray
            the parameters for the tile

        Returns
        -------
        True: Boolean
        """
        # Check correct format
        if not (isinstance(parameters, _np.ndarray) and parameters.ndim == 2):
            raise TypeError("parameters must be two-dimensional np array")

        # Check correct shape
        if not (
            (self._evaluation_points.shape[0] == parameters.shape[0])
            and (self._n_info_per_eval_point == parameters.shape[1])
        ):
            raise TypeError(
                f"Mismatch in parameter size, expected "
                f"{self._evaluation_points.shape[0]} x "
                f"{self._n_info_per_eval_point}"
            )

        # Check if all parameters are within bounds
        if self._parameter_bounds is not None:
            bounds = _np.array(self._parameter_bounds)
            lower_bounds = bounds[:, 0]
            upper_bounds = bounds[:, 1]
            within_bounds = (parameters.ravel() > lower_bounds) & (
                parameters.ravel() < upper_bounds
            )
            if not _np.all(within_bounds):
                raise ValueError(
                    f"The parameters {parameters} must be within the following bounds: "
                    + f"lower: {lower_bounds} and upper: {upper_bounds}"
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
        **kwargs: kwargs

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

    def _process_input(self, parameters, parameter_sensitivities):
        """Processing input for create_tile and _closing_tile

        Parameters
        -------------
        parameters: np.ndarray
            Tile parameters
        parameter_sensitivities: np.ndarray
            Tile parameter sensitivities to be calculated

        Returns
        ---------
        parameters: np.ndarray
            Tile parameters
        n_derivatives: np.ndarray
            Number of different derivatives to compute
        derivatives: list/None
            Initialized list of derivatives
        """
        # Set parameters to default values if not user-given
        if parameters is None:
            default_value = self.default_parameter_value
            self._logd(
                f"Setting parameters to default values ({default_value})"
            )
            if isinstance(default_value, float):
                parameters = _np.full(
                    (len(self.evaluation_points), self.n_info_per_eval_point),
                    default_value,
                )
            elif isinstance(default_value, _np.ndarray):
                parameters = default_value

        # Initialize list of derivatives
        if parameter_sensitivities is not None:
            n_derivatives = parameter_sensitivities.shape[2]
            derivatives = []
        else:
            n_derivatives = 0
            derivatives = None

        # Validity check of parameters and their sensitivities
        self.check_params(parameters)
        self.check_param_derivatives(parameter_sensitivities)

        return parameters, n_derivatives, derivatives

    def _check_custom_parameter(self, value, param_name, bounds):
        """Check if a custom tile parameter (e.g. contact length) is within bounds

        Parameters
        ------------------
        value: float
            Value of the custom parameter
        param_name: str
            Name of custom parameter
        bounds: list<int>
            List of min. and max. bound
        """
        assert isinstance(bounds, list), "Bounds has to be a list"
        assert (
            len(bounds) == 2
        ), "Bounds must consist of a min. and a max. value"
        min_bound, max_bound = bounds

        if not isinstance(value, float):
            raise ValueError(f"Invalid type for {param_name}")

        if not ((value > min_bound) and (value < max_bound)):
            raise ValueError(
                f"{param_name} must be in ({min_bound}, {max_bound})"
            )
