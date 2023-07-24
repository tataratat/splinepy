import gustaf as gus
import numpy as np

from splinepy._base import SplinepyBase
from splinepy.bezier import Bezier
from splinepy.multipatch import Multipatch
from splinepy.splinepy_core import PySpline
from splinepy.utils.data import cartesian_product


class Microstructure(SplinepyBase):
    """Helper class to facilitatae the construction of microstructures."""

    def __init__(
        self,
        deformation_function=None,
        tiling=None,
        microtile=None,
        parametrization_function=None,
    ):
        """Helper class to facilitatae the construction of microstructures.

        Parameters
        ----------
        deformation_function : spline
          Outer function that describes the contour of the microstructured
          geometry
        tiling : list of integers
          microtiles per parametric dimension
        microtile : spline or list of splines
          Representation of the building block defined in the unit cube
        parametrization_function : Callable (optional)
          Function to describe spline parameters
        """
        if deformation_function is not None:
            self.deformation_function = deformation_function

        if tiling is not None:
            self.tiling = tiling

        if microtile is not None:
            self.microtile = microtile

        if parametrization_function is not None:
            self.parametrization_function = parametrization_function

    @property
    def deformation_function(self):
        """Deformation function defining the outer geometry (contour) of the
        microstructure.

        Parameters
        ----------
        None

        Returns
        -------
        deformation_function : spline
        """
        if hasattr(self, "_deformation_function"):
            return self._deformation_function
        else:
            return None

    @deformation_function.setter
    def deformation_function(self, deformation_function):
        """Deformation function setter defining the outer geometry of the
        microstructure. Must be spline type and as such inherit from
        splinepy.Spline.

        Parameters
        ----------
        deformation_function : spline

        Returns
        -------
        None
        """
        if not isinstance(deformation_function, PySpline):
            raise ValueError(
                "Deformation function must be splinepy-Spline."
                " e.g. splinepy.NURBS"
            )
        self._deformation_function = deformation_function
        self._sanity_check()

    @property
    def tiling(self):
        """Number of microtiles per parametric dimension.

        Parameters
        ----------
        None

        Returns
        -------
        tiling : list<int>
        """
        if hasattr(self, "_tiling"):
            return self._tiling
        else:
            return None

    @tiling.setter
    def tiling(self, tiling):
        """Setter for the tiling attribute, defining the number of microtiles
        per parametric dimension.

        Parameters
        ----------
        tiling : int / list<int>
          Number of tiles for each dimension respectively
        Returns
        -------
        None
        """
        if not isinstance(tiling, list):
            if not isinstance(tiling, int):
                raise ValueError(
                    "Tiling mus be either list of integers of integer " "value"
                )
        self._tiling = tiling
        # Is defaulted to False using function arguments
        self._sanity_check()
        self._logd(f"Successfully set tiling to : {self.tiling}")

    @property
    def microtile(self):
        """Microtile that is either a spline, a list of splines, or a class
        that provides a `create_tile` function."""
        if hasattr(self, "_microtile"):
            return self._microtile
        else:
            return None

    @microtile.setter
    def microtile(self, microtile):
        """Setter for microtile.

        Microtile must be either a spline, a list of splines, or a class that
        provides (at least) a `create_tile` function and a `dim` member.

        Parameters
        ----------
        microtile : spline / list<splines> / user-object
          arbitrary long list of splines that define the microtile

        Returns
        -------
        None
        """
        # place single tiles into a list to provide common interface
        if isinstance(microtile, list) or isinstance(microtile, PySpline):
            microtile = self._make_microtilable(microtile)
        # Assign Microtile object to member variable
        self._microtile = microtile

        self._sanity_check()

    @property
    def parametrization_function(self):
        """Function, that - if required - parametrizes the microtiles.

        In order to use said function, the Microtile needs to provide a couple
        of attributes:

         - evaluation_points - a list of points defined in the unit cube
           that will be evaluated in the parametrization function to provide
           the required set of data points
         - parameter_space_dimension - dimensionality of the parametrization
           function and number of design variables for said microtile

        Parameters
        ----------
        None

        Returns
        -------
        parametrization_function : Callable
          Function that descibes the local tile parameters
        """
        if hasattr(self, "_parametrization_function"):
            return self._parametrization_function
        else:
            return None

    @parametrization_function.setter
    def parametrization_function(self, parametrization_function):
        if not callable(parametrization_function):
            raise TypeError("parametrization_function must be callable")
        self._parametrization_function = parametrization_function
        self._sanity_check()

    @property
    def parameter_sensitivity_function(self):
        """Function, that - if required - determines the parameter sensitivity
        of a set of microtiles

        In order to use said function, the Microtile needs to provide a couple
        of attributes:

         - evaluation_points - a list of points defined in the unit cube
           that will be evaluated in the parametrization function to provide
           the required set of data points
         - parameter_space_dimension - dimensionality of the parametrization
           function and number of design variables for said microtile
         - parametrization_function - a function that calculates the microtile
           parameters based on the position of the tile within the deformation
           functions parametric space

        Parameters
        ----------
        None

        Returns
        -------
        parameter_sensitivity_function : Callable
          Function that descibes the local tile parameters
        """
        if hasattr(self, "_parameter_sensitivity_function"):
            return self._parameter_sensitivity_function
        else:
            return None

    @parameter_sensitivity_function.setter
    def parameter_sensitivity_function(self, parameter_sensitivity_function):
        if parameter_sensitivity_function is None:
            self._parameter_sensitivity_function = None
            self._sanity_check()
            return
        if not callable(parameter_sensitivity_function):
            raise TypeError("parametrization_function must be callable")
        self._parameter_sensitivity_function = parameter_sensitivity_function
        self._sanity_check()

    def create(self, closing_face=None, knot_span_wise=None, **kwargs):
        """Create a Microstructure.

        Parameters
        ----------
        closing_face : string
          If not None, Microtile must provide a function `closing_tile`
          Represents coordinate to be a closed surface {"x", "y", "z"}
        knot_span_wise : bool
          Insertion per knotspan vs. total number per paradim
        **kwargs
          will be passed to `create_tile` function

        Returns
        -------
        Microstructure : list<spline>
          finished microstructure based on object requirements
        """
        # Check if all information is gathered
        if not self._sanity_check():
            raise ValueError("Not enough information provided, abort")

        # Set default values
        if knot_span_wise is None:
            knot_span_wise = True

        # check if user wants closed structure
        if (closing_face is not None) and (
            not hasattr(self.microtile, "closing_tile")
        ):
            raise ValueError(
                "Microtile does not provide closing tile definition"
            )
        closing_face_dim = {"x": 0, "y": 1, "z": 2}.get(closing_face)
        is_closed = closing_face_dim is not None
        if not is_closed and (closing_face is not None):
            raise ValueError(
                "Invalid format for closing_face argument, (handed: "
                f"{closing_face}), must be one of"
                "{'x', 'y', 'z'}"
            )

        if is_closed:
            is_closed = True
            if closing_face_dim >= self._deformation_function.para_dim:
                raise ValueError(
                    "closing face must be smaller than the deformation "
                    "function's parametric dimension"
                )
            if self._parametrization_function is None:
                raise ValueError(
                    "Faceclosure is currently only implemented for "
                    "parametrized microstructures"
                )

        # Prepare the deformation function
        # Transform into a non-uniform splinetype and make sure to work on copy
        if "BSpline" in self._deformation_function.whatami:
            deformation_function_copy = self._deformation_function.copy()
        elif hasattr(self._deformation_function, "bspline"):
            deformation_function_copy = self._deformation_function.bspline
        else:
            deformation_function_copy = self._deformation_function.nurbs

        # Create Spline that will be used to iterate over parametric space
        ukvs = deformation_function_copy.unique_knots
        if knot_span_wise:
            for i_pd in range(deformation_function_copy.para_dim):
                if self.tiling[i_pd] == 1:
                    continue
                inv_t = 1 / self.tiling[i_pd]
                new_knots = [
                    ukvs[i_pd][i - 1]
                    + j * inv_t * (ukvs[i_pd][i] - ukvs[i_pd][i - 1])
                    for i in range(1, len(ukvs[i_pd]))
                    for j in range(1, self.tiling[i_pd])
                ]
                # insert knots in both the deformation function
                deformation_function_copy.insert_knots(i_pd, new_knots)
        else:
            self._logd(
                "New knots will be inserted one by one with the objective"
                " to evenly distribute tiles within the parametric domain"
            )
            for i_pd in range(deformation_function_copy.para_dim):
                n_current_spans = len(ukvs[i_pd]) - 1
                if self.tiling[i_pd] == n_current_spans:
                    continue
                elif self.tiling[i_pd] < n_current_spans:
                    self._logw(
                        f"The requested tiling can not be provided, as "
                        f"there are too many knotspans in the deformation"
                        f" function. The tiling in parametric dimension "
                        f"{i_pd} will be set to {n_current_spans}"
                    )
                    self.tiling[i_pd] = n_current_spans
                else:
                    # Determine new knots
                    n_k_span = np.zeros(n_current_spans, dtype=int)
                    span_measure = np.diff(ukvs[i_pd])
                    for _ in range(self.tiling[i_pd] - n_current_spans):
                        add_knot = np.argmax(span_measure)
                        n_k_span[add_knot] += 1
                        span_measure[add_knot] *= n_k_span[add_knot] / (
                            n_k_span[add_knot] + 1
                        )

                    new_knots = []
                    for i, nks in enumerate(n_k_span):
                        new_knots.extend(
                            np.linspace(
                                ukvs[i_pd][i], ukvs[i_pd][i + 1], nks + 2
                            )[1:-1]
                        )
                    deformation_function_copy.insert_knots(i_pd, new_knots)

        # Bezier Extraction for composition
        def_fun_patches = deformation_function_copy.extract.beziers()

        # Calculate parametric space representation for parametrized
        # microstructures
        is_parametrized = self.parametrization_function is not None
        if is_parametrized:
            para_space_dimensions = [np.array([u[0], u[-1]]) for u in ukvs]
            def_fun_para_space = Bezier(
                degrees=[1] * deformation_function_copy.para_dim,
                control_points=cartesian_product(para_space_dimensions),
            ).bspline
            for i_pd in range(deformation_function_copy.para_dim):
                if len(deformation_function_copy.unique_knots[i_pd][1:-1]) > 0:
                    def_fun_para_space.insert_knots(
                        i_pd,
                        deformation_function_copy.unique_knots[i_pd][1:-1],
                    )
            def_fun_para_space = def_fun_para_space.extract.beziers()

        # Determine element resolution
        element_resolutions = [
            len(c) - 1 for c in deformation_function_copy.unique_knots
        ]

        # Initialize return values
        n_sensitivities = 0
        if self.parameter_sensitivity_function is not None:
            n_sensitivities = self.parameter_sensitivity_function(
                self._microtile.evaluation_points
            ).shape[2]
            self._microstructure_derivatives = [
                [] for i in range(n_sensitivities)
            ]
        self._microstructure = []

        # Start actual composition
        if is_parametrized:
            for i, (def_fun, def_fun_para) in enumerate(
                zip(def_fun_patches, def_fun_para_space)
            ):
                # Evaluate tile parameters
                positions = def_fun_para.evaluate(
                    self._microtile.evaluation_points
                )
                tile_parameters = self.parametrization_function(positions)

                if self.parameter_sensitivity_function is not None:
                    tile_sensitivities = self.parameter_sensitivity_function(
                        positions
                    )
                    support = np.where(
                        np.any(
                            np.any(tile_sensitivities != 0.0, axis=0), axis=0
                        )
                    )[0]
                    anti_support = np.where(
                        np.any(
                            np.all(tile_sensitivities == 0.0, axis=0), axis=0
                        )
                    )[0]
                    tile_sens_on_support = tile_sensitivities[:, :, support]
                else:
                    tile_sens_on_support = None

                # Check if center or closing tile
                if is_closed:
                    # check index
                    index = i
                    for ipd in range(closing_face_dim):
                        index -= index % element_resolutions[ipd]
                        index = int(index / element_resolutions[ipd])
                    index = index % element_resolutions[closing_face_dim]
                    if index == 0:
                        # Closure at minimum id
                        tile = self._microtile.closing_tile(
                            parameters=tile_parameters,
                            parameter_sensitivities=tile_sens_on_support,
                            closure=closing_face + "_min",
                            **kwargs,
                        )
                    elif (index + 1) == element_resolutions[closing_face_dim]:
                        # Closure at maximum
                        tile = self._microtile.closing_tile(
                            parameters=tile_parameters,
                            parameter_sensitivities=tile_sens_on_support,
                            closure=closing_face + "_max",
                            **kwargs,
                        )
                    else:
                        tile = self._microtile.create_tile(
                            parameters=tile_parameters,
                            parameter_sensitivities=tile_sens_on_support,
                            **kwargs,
                        )
                else:
                    tile = self._microtile.create_tile(
                        parameters=tile_parameters,
                        parameter_sensitivities=tile_sens_on_support,
                        **kwargs,
                    )

                if isinstance(tile, tuple):
                    # Returned tile and derivatives
                    (splines, derivatives) = tile
                    for tile_patch in splines:
                        self._microstructure.append(
                            def_fun.compose(tile_patch)
                        )
                    n_splines_per_tile = len(splines)
                    empty_splines = [None] * n_splines_per_tile
                    for i in anti_support:
                        self._microstructure_derivatives[i].extend(
                            empty_splines
                        )
                    for i, deris in enumerate(derivatives):
                        for j, (tile_v, tile_deriv) in enumerate(
                            zip(splines, deris)
                        ):
                            self._microstructure_derivatives[
                                support[i]
                            ].append(
                                def_fun.composition_derivative(
                                    tile_v, tile_deriv
                                )
                            )
                else:
                    # Perform composition
                    for tile_patch in tile:
                        self._microstructure.append(
                            def_fun.compose(tile_patch)
                        )

        # Not parametrized
        else:
            # Tile can be computed once (prevent to many evaluations)
            tile = self._microtile.create_tile(**kwargs)
            for def_fun in def_fun_patches:
                for t in tile:
                    self._microstructure.append(def_fun.compose(t))

        # return copy of precomputed member
        if n_sensitivities == 0:
            return self._microstructure.copy()
        else:
            return (
                self._microstructure.copy(),
                self._microstructure_derivatives.copy(),
            )

    def show(self, use_saved=False, return_gustaf=False, **kwargs):
        """
        Shows microstructure. Consists of deformation_function, microtile, and
        microstructure. Supported only by vedo.

        Parameters
        ----------
        use_saved: bool
        return_gustaf: bool
        **kwargs: kwargs
          Will be passed to show function

        Returns
        -------
        gustaf_obj: dict
          keys are deformation_function, microtile, and microstructure.
          Iff return_gustaf is True.
        plt: vedo.Plotter
        """
        if use_saved:
            if hasattr(self, "_microstructure"):
                microstructure = Multipatch(self._microstructure)
            else:
                raise ValueError("No previous microstructure saved")
        else:
            # Create on the fly
            microstructure = Multipatch(splines=self.create(**kwargs))

        # Precompute splines
        microtile = self.microtile.create_tile(**kwargs)
        deformation_function = self.deformation_function

        # Show in vedo
        return gus.show(
            ["Deformation Function", deformation_function],
            ["Microtile", microtile],
            ["Composed Microstructure", microstructure],
            **kwargs,
        )

    def _sanity_check(self):
        """Check all members and consistency of user data.

        Parameters
        ----------
        updated_properties : bool
          Sets the updated_properties variable to value, which indicates,
          wheither the microstructure needs to be rebuilt

        Returns
        -------
        passes: bool
        """
        if (
            (self.deformation_function is None)
            or (self.microtile is None)
            or (self.tiling is None)
        ):
            self._logd(
                "Current information not sufficient,"
                " awaiting further assignments"
            )
            return False
        # Check if microtile object fulfils requirements
        if not hasattr(self._microtile, "create_tile"):
            raise ValueError(
                "Microtile class does not provide the necessary "
                "attribute `create_tile`, that is required for "
                "microstructure construction"
            )
        if not hasattr(self._microtile, "dim"):
            raise ValueError(
                "Microtile class does not provide the necessary "
                "attribute `dim`, defining the dimensionality of "
                "the created tile"
            )

        # Check if parametric dimensions are consistent
        if not self.deformation_function.para_dim == self._microtile.dim:
            raise ValueError(
                "Microtile dimension must match parametric dimension of "
                "deformation function to enable composition"
            )

        # Check if tiling is consistent
        if isinstance(self.tiling, int):
            self.tiling = [self.tiling] * self.deformation_function.para_dim
        if len(self.tiling) != self.deformation_function.para_dim:
            raise ValueError(
                "Tiling list must have one entry per parametric dimension"
                " of the deformation function"
            )

        if self.parametrization_function is None:
            # All checks have passed
            return True

        self._logd("Checking compatibility of parametrization function")
        if not hasattr(self._microtile, "evaluation_points"):
            raise ValueError(
                "Microtile class does not provide the necessary "
                "attribute `evaluation_points`, that is required for"
                " a parametrized microstructure construction"
            )
        result = self._parametrization_function(
            self._microtile.evaluation_points
        )
        if not isinstance(result, np.ndarray):
            raise ValueError(
                "Function outline of parametrization function must be "
                "`f(np.ndarray)->np.ndarray`"
            )
        if self.parameter_sensitivity_function is None:
            # All checks have passed
            return True

        # Check sensitivity function
        result = self.parameter_sensitivity_function(
            self._microtile.evaluation_points
        )
        if (not isinstance(result, np.ndarray)) or (not result.ndim == 3):
            raise ValueError(
                "Function outline of parameter sensitivity function must "
                "be `f(np.ndarray)->np.ndarray`, with dimension 3 and shape:"
                "n_evaluation_points x n_info_per_evaluation_point x "
                "n_sensitivities"
            )

        return True

    def _make_microtilable(self, microtile):
        """Creates a Microtile object on the fly if user only provides (a list
        of) splines. Internal use only.

        Parameters
        ----------
        microtile : spline / list<splines>
          Microtile definition of a spline
        """
        return _UserTile(microtile)


class _UserTile:
    def __init__(self, microtile):
        """
        On the fly created class of a user tile
        Parameters
        ----------
        microtile : spline , list<splines>
        """
        # Assign microtiles
        self._user_tile = []

        if not isinstance(microtile, list):
            microtile = [microtile]

        for m in microtile:
            if not isinstance(m, PySpline):
                raise ValueError(
                    "Microtiles must be (list of) "
                    "splinepy-Splines. e.g. splinepy.NURBS"
                )
            # Extract beziers for every non Bezier patch else this just
            # returns itself
            self._user_tile.extend(m.extract.beziers())
        self._dim = microtile[0].dim
        for m in microtile:
            if m.dim != self._dim:
                raise ValueError("Dimensions of spline lists inconsistent")

    @property
    def dim(self):
        return self._dim

    def create_tile(self, **kwargs):
        """Create a tile on the fly."""
        return self._user_tile.copy()
