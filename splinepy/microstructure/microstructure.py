import gustaf as _gus
import numpy as _np

from splinepy._base import SplinepyBase as _SplinepyBase
from splinepy.bspline import BSpline as _BSpline
from splinepy.multipatch import Multipatch as _Multipatch
from splinepy.splinepy_core import PySpline as _PySpline
from splinepy.utils.data import cartesian_product as _cartesian_product


class Microstructure(_SplinepyBase):
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

        if not isinstance(deformation_function, _PySpline):
            raise ValueError(
                "Deformation function must be splinepy-Spline."
                " e.g. splinepy.NURBS"
            )
        self._deformation_function = deformation_function

        if self._deformation_function.has_knot_vectors:
            self._deformation_function.normalize_knot_vectors()

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
        if not isinstance(tiling, list) and not isinstance(tiling, int):
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
        if isinstance(microtile, (_PySpline, list)):
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
          Function that describes the local tile parameters
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
          Function that describes the local tile parameters
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

    def _additional_knots(self, knot_span_wise):
        """Determine all required additional knots to be inserted into the
        original spline representation of the deformation function

        Parameters
        ----------
        knot_span_wise: bool
          Decicive, whether tiling is applied to all knot-spans respectively,
          or globally

        Returns
        -------
        additional_knots : list<list>
          returns a list of all required additional knots as a list of list to
          the respective parametric dimension. If no additional knot is to be
          inserted, the list is empty
        """
        additional_knots = []
        # Create Spline that will be used to iterate over parametric space
        ukvs = self.deformation_function.unique_knots
        if knot_span_wise:
            for tt, ukv in zip(self.tiling, ukvs):
                inv_t = 1 / tt
                new_knots = [
                    ukv[i - 1] + j * inv_t * (ukv[i] - ukv[i - 1])
                    for i in range(1, len(ukv))
                    for j in range(1, tt)
                ]
                additional_knots.append(new_knots)
        else:
            self._logd(
                "New knots will be inserted one by one with the objective"
                " to evenly distribute tiles within the parametric domain"
            )
            for i_pd, (tt, ukv) in enumerate(zip(self.tiling, ukvs)):
                n_current_spans = len(ukv) - 1
                if tt == n_current_spans:
                    continue
                elif tt < n_current_spans:
                    self._logw(
                        f"The requested tiling can not be provided, as "
                        f"there are too many knotspans in the deformation"
                        f" function. The tiling in parametric dimension "
                        f"{i_pd} will be set to {n_current_spans}"
                    )
                    self.tiling[i_pd] = n_current_spans
                else:
                    # Determine new knots
                    n_k_span = _np.zeros(n_current_spans, dtype=int)
                    span_measure = _np.diff(ukv)
                    for _ in range(tt - n_current_spans):
                        add_knot = _np.argmax(span_measure)
                        n_k_span[add_knot] += 1
                        span_measure[add_knot] *= n_k_span[add_knot] / (
                            n_k_span[add_knot] + 1
                        )

                    new_knots = []
                    for i, nks in enumerate(n_k_span):
                        new_knots.extend(
                            _np.linspace(ukv[i], ukv[i + 1], nks + 2)[1:-1]
                        )
                    additional_knots.append(new_knots)

        # Return all knots
        return additional_knots

    def _compute_tiling_prerequisites(
        self, knot_span_wise, is_parametrized, macro_sensitivity
    ):
        """
        Prepare all required information to start the construction process. All
        additional properties are returned None, to avoid argument dependent
        return values

        Parameters
        ----------
        knot_span_wise: bool
          Decides if the number of tiles in tiling is supposed to be
          interpreted as the total amount or per knot span
        is_parametrized : bool
          Decides if the microtiles are adapted before insertion into the
          microstructure
        macro_sensitivity :bool
          Provides information required for the calculation of derivatives with
          respect to the control points of the deformation function

        Returns
        -------
        def_fun_patches: list<BezierBasis>
          Extracted Bezier patches, that define the position of the individual
          microtiles
        bezier_matrices: list<numpy/scipy.sparse>
          (None if `macro_sensitivity` is `False`) Matrices that can be used to
          map from the original control points to the control points of the
          extracted bezier patches
        para_space_patches: list<BezierBasis>
          (None if `is_parametrized` is `False`) Linear representation of the
          parametric space that can be used to evaluate the parametrization
          function on the individual elements
        element_resolutions: list<int>
          Number of elements in the respective parametric dimension
        """

        # Prepare the deformation function
        # Transform into a non-uniform splinetype and make sure to work on copy
        if self._deformation_function.is_rational:
            deformation_function_copy = self._deformation_function.nurbs
        else:
            deformation_function_copy = self._deformation_function.bspline

        additional_knots = self._additional_knots(knot_span_wise)

        # First step : insert all required knots into the deformation function
        # spline
        if macro_sensitivity:
            knot_insertion_matrix = (
                deformation_function_copy.knot_insertion_matrix(
                    0, additional_knots[0]
                )
            )
        deformation_function_copy.insert_knots(0, additional_knots[0])
        for i_pd, akv in enumerate(additional_knots[1:], start=1):
            if macro_sensitivity:
                knot_insertion_matrix = (
                    deformation_function_copy.knot_insertion_matrix(i_pd, akv)
                    @ knot_insertion_matrix
                )
            deformation_function_copy.insert_knots(i_pd, akv)

        # Second step (if MS is parametrized)
        if is_parametrized:
            para_bounds = self.deformation_function.parametric_bounds.T
            def_fun_para_space = _BSpline(
                degrees=[1] * self.deformation_function.para_dim,
                knot_vectors=para_bounds.repeat(2, 1),
                control_points=_cartesian_product(para_bounds),
            )
            for i_pd in range(deformation_function_copy.para_dim):
                if len(deformation_function_copy.unique_knots[i_pd][1:-1]) > 0:
                    def_fun_para_space.insert_knots(
                        i_pd,
                        deformation_function_copy.unique_knots[i_pd][1:-1],
                    )

        # Extract the bezier patches
        def_fun_patches = deformation_function_copy.extract.beziers()

        # Precompute the extraction matrices
        if macro_sensitivity:
            bezier_matrices = deformation_function_copy.knot_insertion_matrix(
                beziers=True
            )

            # However, these refer to the macro spline that has knots already
            # inserted into it, so we need to apply the matrix constructed
            # earlier
            for i_matrix in range(len(bezier_matrices)):
                bezier_matrices[i_matrix] @= knot_insertion_matrix
        else:
            bezier_matrices = None

        # Do the same for the parametric space representation
        if is_parametrized:
            para_space_patches = def_fun_para_space.extract.beziers()
        else:
            para_space_patches = None

        # Determine element resolution
        element_resolutions = [
            len(c) - 1 for c in deformation_function_copy.unique_knots
        ]

        return (
            def_fun_patches,
            bezier_matrices,
            para_space_patches,
            element_resolutions,
        )

    def create(
        self,
        closing_face=None,
        knot_span_wise=None,
        macro_sensitivities=None,
        **kwargs,
    ):
        """Create a Microstructure.

        Parameters
        ----------
        closing_face : string
          If not None, Microtile must provide a function `closing_tile`
          Represents coordinate to be a closed surface {"x", "y", "z"}
        knot_span_wise : bool
          Insertion per knotspan vs. total number per paradim
        macro_sensitivities: bool
          Calculate the derivatives of the structure with respect to the outer
          control point variables
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
        if macro_sensitivities is None:
            macro_sensitivities = False

        # check if user wants closed structure
        if closing_face is not None:
            if not hasattr(self.microtile, "_closing_tile"):
                raise ValueError(
                    "Microtile does not provide closing tile definition"
                )
            closing_face_dim = {"x": 0, "y": 1, "z": 2}.get(closing_face)
            if closing_face is None:
                raise ValueError(
                    "Invalid format for closing_face argument, (handed: "
                    f"{closing_face}), must be one of"
                    "{'x', 'y', 'z'}"
                )
            if closing_face_dim >= self._deformation_function.para_dim:
                raise ValueError(
                    "closing face must be smaller than the deformation "
                    "function's parametric dimension"
                )
            is_closed = True
        else:
            is_closed = False

        # Check if parametrized
        is_parametrized = self.parametrization_function is not None
        parameter_sensitivities = (
            self.parameter_sensitivity_function is not None
        )

        if parameter_sensitivities and not is_parametrized:
            raise ValueError(
                "Parameter sensitivities can not be interpreted without"
                "parametrization function"
            )

        # Using the additional knots and the bezier extraction operator, we can
        # provide all necessary information on the construction process
        (
            bezier_patches,
            knot_insertion_matrices,
            para_space_patches,
            element_resolutions,
        ) = self._compute_tiling_prerequisites(
            knot_span_wise=knot_span_wise,
            is_parametrized=is_parametrized,
            macro_sensitivity=macro_sensitivities,
        )

        # Use element_resolutions to get global indices of patches
        local_indices = _cartesian_product(
            [_np.arange(er) for er in element_resolutions]
        )

        # Determine the number of geometric derivatives with respect to the
        # design variables
        if parameter_sensitivities:
            n_parameter_sensitivities = self.parameter_sensitivity_function(
                self._microtile.evaluation_points
            ).shape[2]
        else:
            n_parameter_sensitivities = 0

        # Determine the number of geometric derivatives with respect to the
        # deformation function's control points
        if macro_sensitivities:
            n_macro_sensitivities = self.deformation_function.cps.size
        else:
            n_macro_sensitivities = 0

        # Prepare field for derivatives
        if macro_sensitivities or parameter_sensitivities:
            spline_list_derivs = [
                []
                for i in range(
                    n_parameter_sensitivities + n_macro_sensitivities
                )
            ]

        spline_list_ms = []

        # Start actual composition
        for i_patch, def_fun in enumerate(bezier_patches):
            # If the structure is parametrized, provide parameters for the tile
            # construction
            if is_parametrized:
                # Retrieve position from the parametric space patches
                positions = para_space_patches[i_patch].evaluate(
                    self._microtile.evaluation_points
                )
                # Evaluate parametrization function
                tile_parameters = self.parametrization_function(positions)
            else:
                tile_parameters = None

            # If the sensitivities are requested, evaluate the sensitivity
            # function, which must be provided by the user
            if parameter_sensitivities:
                tile_sensitivities = self.parameter_sensitivity_function(
                    positions
                )
                # To avoid unnecessary constructions (if a parameter
                # sensitivity evaluates to zero), perform only those that are
                # in the support of a specific design variable
                support = _np.where(
                    _np.any(_np.any(tile_sensitivities != 0.0, axis=0), axis=0)
                )[0]
                anti_support = _np.where(
                    _np.any(_np.all(tile_sensitivities == 0.0, axis=0), axis=0)
                )[0]
                tile_sens_on_support = tile_sensitivities[:, :, support]
            else:
                tile_sens_on_support = None

            # Check if center or closing tile
            if is_closed:
                # check index
                index = local_indices[i_patch, closing_face_dim]
                face_closure = (
                    closing_face + "_min"
                    if (index == 0)
                    else (
                        closing_face + "_max"
                        if (index + 1) == element_resolutions[closing_face_dim]
                        else None
                    )
                )
            else:
                face_closure = None

            # Determine the microtile (prior to insertion) derivative is set to
            # None if the parameter sensitivity is not requested
            tile, derivatives = self._microtile.create_tile(
                parameters=tile_parameters,
                parameter_sensitivities=tile_sens_on_support,
                closure=face_closure,
                **kwargs,
            )

            # Perform the composition
            if macro_sensitivities:
                basis_function_compositions = []
                for tile_patch in tile:
                    composed, basis_function_composition = def_fun.compose(
                        tile_patch, compute_sensitivities=True
                    )
                    spline_list_ms.append(composed)
                    basis_function_compositions.append(
                        basis_function_composition
                    )
            else:
                for tile_patch in tile:
                    spline_list_ms.append(def_fun.compose(tile_patch))

            # Determine the geometric sensitivities using the chain rule
            if parameter_sensitivities:
                n_splines_per_tile = len(tile)
                empty_splines = [None] * n_splines_per_tile
                for j in anti_support:
                    spline_list_derivs[j].extend(empty_splines)
                for j, deris in enumerate(derivatives):
                    for tile_v, tile_deriv in zip(tile, deris):
                        spline_list_derivs[support[j]].append(
                            def_fun.composition_derivative(tile_v, tile_deriv)
                        )

            # Compute the derivatives with respect to the deformation
            # function's ctps
            if macro_sensitivities:
                for patch_info in basis_function_compositions:
                    # patch_info is a list of splines all representing the
                    # derivatives of the composed spline with respect to a
                    # component of a control point within the bezier patch
                    # of the deformation function.
                    cps = _np.hstack([p.cps for p in patch_info])
                    # We use the matrices to map the contributions of the
                    # bezier ctps to the deformation functions
                    mapped_cps = cps @ knot_insertion_matrices[i_patch]

                    # Last step: create beziers to be added to the derivative
                    # fields
                    for j_cc in range(n_macro_sensitivities):
                        control_points = _np.zeros(
                            (cps.shape[0], self._deformation_function.dim)
                        )
                        ii_ctps, jj_dim = divmod(
                            j_cc, self._deformation_function.dim
                        )
                        control_points[:, jj_dim] = mapped_cps[:, ii_ctps]
                        spline_list_derivs[
                            n_parameter_sensitivities + j_cc
                        ].append(
                            type(patch_info[0])(
                                patch_info[0].degrees, control_points
                            )
                        )

        # Use a multipatch object to bundle all information
        self._microstructure = _Multipatch(splines=spline_list_ms)

        # Add fields if requested
        if macro_sensitivities or parameter_sensitivities:
            self._microstructure.add_fields(
                spline_list_derivs, self.deformation_function.dim
            )

        return self._microstructure

    def show(self, use_saved=False, **kwargs):
        """
        Shows microstructure. Consists of deformation_function, microtile, and
        microstructure. Supported only by vedo.

        Parameters
        ----------
        use_saved: bool
        **kwargs: kwargs
          Will be passed to show function

        Returns
        -------
        plt: vedo.Plotter
        """
        if use_saved:
            if hasattr(self, "_microstructure"):
                microstructure = self._microstructure
            else:
                raise ValueError("No previous microstructure saved")
        else:
            # Create on the fly
            microstructure = self.create(**kwargs)

        # Precompute splines
        microtile = self.microtile.create_tile(**kwargs)
        deformation_function = self.deformation_function

        # Show in vedo
        return _gus.show(
            ["Deformation Function", deformation_function],
            ["Microtile", microtile[0]],
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
        if not isinstance(result, _np.ndarray):
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
        if (not isinstance(result, _np.ndarray)) or (not result.ndim == 3):
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


class _UserTile(_SplinepyBase):
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
            if not isinstance(m, _PySpline):
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

    def create_tile(self, **kwargs):  # noqa ARG002
        """Create a tile on the fly.
        Derivative of a constant tile can not be requested, so derivative must
        be None
        Parameters
        ----------
        None

        Returns
        -------
        Microtile: list<Bezier>
        Derivative : None
        """
        for k, v in kwargs.items():
            self._logd(f"Additional argument {k} : {v} will be ignored")
        return (self._user_tile.copy(), None)
