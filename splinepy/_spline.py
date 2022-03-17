import abc
import logging
import copy
import os
import itertools

import numpy as np

from splinepy import utils
from splinepy import io

class InputDimensionError(Exception):
    """
    Raised when input dimension does not match spline's internal dimension.
    """
    pass


class Spline(abc.ABC):
    """
    Abstract Spline Class.
    """

    __slots__ = [
        "_degrees",
        "_knot_vectors",
        "_control_points",
        "_weights",
        "_c_spline",
        "_para_dim",
        "_dim",
    ]

    def __init__(
            self,
            degrees=None,
            knot_vectors=None,
            control_points=None,
            weights=None,
    ):
        """
        Spline.

        Parameters
        -----------
        degrees: (para_dim,) list-like
        knot_vectors: (para_dim, n) list
        control_points: (m, dim) list-like
        weights: (m,) list-like

        Attributes
        -----------
        whatami: str
        para_dim: int
        dim: int
        degrees: np.ndarray
        knot_vectors: list
        control_points: np.ndarray
        knot_vector_bounds: np.ndarray
        control_point_bounds: np.ndarray
        skip_update: bool

        Returns
        --------
        None
        """
        self.degrees = degrees
        self.knot_vectors = knot_vectors
        self.control_points = control_points

        # Only for NURBS
        if weights is not None:
            self.weights = weights

    def clear(self):
        """
        Clear all properties.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        for s in self.__slots__:
            if hasattr(self, s):
                delattr(self, s)

        logging.debug("Spline - All attributes are cleared!")

    @property
    def whatami(self):
        """
        Answers a deep philosophical question of "what am i?"

        Parameters
        -----------
        None

        Returns
        --------
        whatami: str
        """
        if hasattr(self, "_c_spline"):
            return self._c_spline.whatami

        else:
            return "Nothing"

    @whatami.setter
    def whatami(self, how_dare_you):
        """
        Sassily refuse to tell me what i am.
        Totally unnecessary.

        Parameters
        -----------
        how_dare_you: seriously

        Returns
        --------
        None
        """
        logging.warning("Spline - Excuse me, you cannot tell me what I am.")

    @property
    def skip_update(self):
        """
        Returns the status of skip_update switch

        Parameters
        -----------
        None

        Returns
        --------
        skip_update: bool
        """
        if hasattr(self, "_c_spline"):
            return self._c_spline.skip_update

        else:
            return None

    @skip_update.setter
    def skip_update(self, skip_update):
        """
        Default is that cpp syncs itself with python each time before it
        executes any function. This is to deal with in-place changes of knot
        vectors and control_points.
        If you know that you aren't doing this (in-place changes),
        set this to `True`, to avoid update/syncing.

        Parameters
        -----------
        skip_update: bool
        """
        if hasattr(self, "_c_spline"):
            self._c_spline.skip_update = skip_update

        else:
            return None

    @property
    def para_dim(self):
        """
        Returns parametric dimension.

        Parameters
        -----------
        None

        Returns
        --------
        para_dim: int
        """
        if hasattr(self, "_para_dim"):
            return self._para_dim

        else:
            return None

    @property
    def dim(self):
        """
        Returns physical dimension.

        Parameters
        -----------
        None

        Returns
        --------
        dim: int
        """
        if hasattr(self, "_dim"):
            return self._dim

        else:
            return None

    @property
    def degrees(self):
        """
        Returns Degrees.

        Parameters
        -----------
        None

        Returns
        --------
        degrees: (para_dim,) np.ndarray
        """
        if hasattr(self, "_degrees"):
            return self._degrees

        else:
            return None

    @degrees.setter
    def degrees(self, degrees):
        """
        Sets Degrees.

        Parameters
        -----------
        degrees: list-like

        Returns
        --------
        None
        """
        if degrees is None:
            if hasattr(self, "_degrees"):
                delattr(self, "_degrees")
                if hasattr(self, "_para_dim") and self.knot_vectors is None:
                    delattr(self, "_para_dim")
            return None

        # flatten copies
        degrees = utils.make_c_contiguous(degrees, "int32").flatten()

        if self.para_dim is None:
            self._para_dim = degrees.size

        else:
            if self.para_dim != degrees.size:
                raise InputDimensionError(
                    "Input dimension does not match spline's para_dim "
                )

        self._degrees = degrees

        logging.debug(f"Spline - Degrees set: {self.degrees}")

        self._update_c()

    @property
    def knot_vectors(self):
        """
        Returns knot vectors.

        Parameters
        -----------
        None

        Returns
        --------
        knot_vectors: list
        """
        if hasattr(self, "_knot_vectors"):
            return self._knot_vectors

        else:
            return None

    @knot_vectors.setter
    def knot_vectors(self, knot_vectors):
        """
        Set knot vectors.

        Parameters
        -----------
        knot_vectors: list

        Returns
        --------
        None
        """
        if knot_vectors is None:
            if hasattr(self, "_knot_vectors"):
                delattr(self, "_knot_vectors")
                if hasattr(self, "_para_dim") and self.degrees is None:
                    delattr(self, "_para_dim")
            return None

        if not isinstance(knot_vectors, list):
            raise TypeError("knot_vectors should be a `list`!")
            
        if self.para_dim is None:
            self._para_dim = len(knot_vectors)

        else:
            if self.para_dim != len(knot_vectors):
                raise InputDimensionError(
                    "Input dimension does not match spline's para_dim "
                )

        self._knot_vectors = copy.deepcopy(knot_vectors)

        logging.debug("Spline - Knot vectors set:")
        for i, kv in enumerate(self.knot_vectors):
            logging.debug(
                f"Spline -   {i}"
                ". knot vector length: "
                f"{len(kv)}"
            )

        self._update_c()

    @property
    def unique_knots(self,):
        """
        Returns unique knots.
        Computed using `np.unique`.
        Does not store results.

        Parameters
        -----------
        None

        Returns
        --------
        unique_knots: list
        """
        logging.debug("Spline - Computing unique knots using `np.unique`.")
        unique_knots = []
        for k in self.knot_vectors:
            unique_knots.append(np.unique(k).tolist())

        return unique_knots

    @property
    def knot_vector_bounds(self,):
        """
        Returns bounds of knot vectors.

        Parameters
        -----------
        None

        Returns
        --------
        knot_vector_bounds: (2, para_dim) np.ndarray
        """
        if self.knot_vectors is None:
            return None

        logging.debug("Spline - Computing knot_vectors_bounds")
        lower_bounds = []
        upper_bounds = []
        kvs = self.knot_vectors
        for i in range(self.para_dim):
            lower_bounds.append(min(kvs[i]))
            upper_bounds.append(max(kvs[i]))

        return np.vstack((lower_bounds, upper_bounds))

    @property
    def control_points(self,):
        """
        Returns control points.

        Parameters
        -----------
        None

        Returns
        --------
        control_points_: (n, dim) np.ndarray
        """
        if hasattr(self, "_control_points"):
            return self._control_points

        else:
            return None

    @control_points.setter
    def control_points(self, control_points):
        """
        Set control points.

        Parameters
        -----------
        control_points: (n, dim) list-like

        Returns
        --------
        None
        """
        if control_points is None:
            if hasattr(self, "_control_points"):
                delattr(self, "_control_points")
                delattr(self, "_dim")
            return None

        control_points = utils.make_c_contiguous(control_points, "float64").copy()

        if self.dim is None:
            self._dim = control_points.shape[1]

        else:
            if self.dim != control_points.shape[1]:
                raise InputDimensionError(
                    "Input dimension does not match spline's dim "
                )

        self._control_points = control_points
        logging.debug(
            f"Spline - {self.control_points.shape[0]} Control points set."
        )

        self._update_c()

    def _lexsort_control_points(self, order):
        """
        Lexsorts control points.
        Not always the solution.
        Leading underscore to mark development/internal use.

        Parameters
        -----------
        order: list

        Returns
        --------
        None
        """
        ind = np.lexsort([self.control_points[:,i] for i in order])
        logging.debug(f"Spline - `lexsort` control points ({order})")
        self.control_points = self.control_points[ind]

    @property
    def control_point_bounds(self,):
        """
        Returns bounds of control points.

        Parameters
        -----------
        None

        Returns
        --------
        control_point_bounds: (2, dim) np.ndarray
        """
        logging.debug("Spline - Computing control_point_bounds")
        cps = self.control_points

        return np.vstack(
            (
                cps.min(axis=0),
                cps.max(axis=0),
            )
        )

    @property
    def control_mesh_resolutions(self,):
        """
        Returns control mesh resolutions.

        Parameters
        -----------
        None

        Returns
        --------
        control_mesh_resolutions: list
        """
        cmr = []
        #for i in range(self.para_dim):
        #    cmr.append(
        #        len(self.knot_vectors[i]) - self.degrees[i] - 1
        #    )
        for kv, d in zip(self.knot_vectors, self.degrees):
            cmr.append(len(kv) - d - 1)

        return cmr

    @abc.abstractmethod
    def _update_c(self,):
        """
        Updates/Init cpp spline, if it is ready to be updated.
        Checks if all the entries are filled before updating.

        This needs to be separately implemented, as BSpline and NURBS have 
        different update process. 

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        pass

    @abc.abstractmethod
    def _update_p(self,):
        """
        Reads cpp spline and writes it here.
        Probably get an error if cpp isn't ready for this.

        This needs to be separately implemented, as BSpline and NURBS have 
        different update process. 

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        pass

    def evaluate(self, queries, n_threads=1):
        """
        Evaluates spline.

        Parameters
        -----------
        queries: (n, para_dim) list-like
        n_threads: int
          Default is 1. If n_thread > 1, it calls multithreading version of
          `evaluate`. Parallel loop implemented in cpp.

        Returns
        --------
        results: (n, dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        queries = utils.make_c_contiguous(queries, dtype="float64")

        if queries.shape[1] != self.para_dim:
            raise InputDimensionError(
                "`queries` does not match current pametric dimension."
            )

        logging.debug("Spline - Evaluating spline...")

        if int(n_threads) > 1:
            return self._c_spline.p_evaluate(
                queries=queries,
                n_threads=n_threads
            )

        else:
            return self._c_spline.evaluate(queries=queries)

    def derivative(self, queries, orders, n_threads=1):
        """
        Evaluates derivatives of spline.

        Parameters
        -----------
        queries: (n, para_dim) list-like
        orders: (para_dim,) list-like
        n_threads: int
          Default is 1. If n_thread > 1, it calls multithreading version of
          `derivate`. Parallel loop implemented in cpp.

        Returns
        --------
        results: (n, dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        queries = utils.make_c_contiguous(queries, dtype="float64")
        orders = utils.make_c_contiguous(orders, dtype="float64")

        if queries.shape[1] != self.para_dim:
            raise InputDimensionError(
                "`queries` does not match current pametric dimension."
            )
        if orders.size != self.para_dim:
            raise InputDimensionError(
                "`orders` does not match current pametric dimension."
            )

        logging.debug("Spline - Evaluating derivatives of the spline...")


        if int(n_threads) > 1:
            return self._c_spline.p_derivative(
                queries=queries,
                orders=orders,
                n_threads=n_threads
            )

        else:
            return self._c_spline.derivative(
                queries=queries,
                orders=orders
            )

    def basis_functions(self, queries, n_threads=1):
        """
        Returns basis function values and their support ids of given queries.

        Parameters
        -----------
        queries: (n, para_dim) list-like
        n_threads: int
          Default is 1. Higher number is currently not supported. #TODO: DOIT

        Returns
        --------
        results: tuple
          tuple of two elements.
          first: (prod(degrees .+ 1)) array of basis function values.
          second: support ids.
        """
        if self.whatami == "Nothing":
          return None

        queries = utils.make_c_contiguous(queries, dtype="float64")

        if queries.shape[1] != self.para_dim:
            raise InputDimensionError(
                "`queries` does not match current pametric dimension."
            )

        logging.debug("Spline - evaluating basis functions")

        return self._c_spline.basis_functions(queries)

    def insert_knots(self, parametric_dimension, knots):
        """
        Inserts knots. 

        Parameters
        -----------
        parametric_dimension: int
        knots: list or float

        Returns
        --------
        None
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to insert knots."
            )

        if isinstance(knots, float):
            knots = [knots]

        elif isinstance(knots, np.ndarray):
            knots = knots.tolist()

        if not isinstance(knots, list):
            raise TypeError(
                "We couldn't convert input to `list`. Please give us `list`."
            )

        if max(knots) > max(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too big)"
            )

        if min(knots) < min(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too small)"
            )

        self._c_spline.insert_knots(
            int(parametric_dimension),
            knots
        )

        logging.debug(f"Spline - Inserted {len(knots)} knot(s).")

        self._update_p()

    def remove_knots(self, parametric_dimension, knots, tolerance=1e-8):
        """
        Tries to removes knots. If you've compiled `splinepy` in `Debug`
        and your removal request is not "accepted", you will get an error.
        See the comments for `Nurbs::remove_knots` @ 
        `splinepy/src/nurbs.hpp` for more info.

        Parameters
        -----------
        parametric_dimension: int
        knots: list or float
        tolerance: float

        Returns
        --------
        None
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to remove knots."
            )

        if isinstance(knots, float):
            knots = [knots]

        elif isinstance(knots, np.ndarray):
            knots = knots.tolist()

        if not isinstance(knots, list):
            raise TypeError(
                "We couldn't convert input to `list`. Please give us `list`."
            )

        if max(knots) > max(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too big)"
            )

        if min(knots) < min(self.knot_vectors[parametric_dimension]):
            raise ValueError(
                "One of the query knots not in valid knot range. (Too small)"
            )

        total_knots_before = len(self.knot_vectors[int(parametric_dimension)])
        self._c_spline.remove_knots(
            int(parametric_dimension),
            knots,
            tolerance,
        )

        self._update_p()

        logging.debug(
            f"Spline - Tried to remove {len(knots)} knot(s)."
        )
        logging.debug(
            "Spline - Actually removed {nk} knot(s).".format(
                nk=(
                    total_knots_before
                    - len(self.knot_vectors[int(parametric_dimension)])
                )
            )
        )

    def elevate_degree(self, parametric_dimension):
        """
        Elevate degree.

        Parameters
        -----------
        parametric_dimension: int

        Returns
        --------
        None
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to elevate_degree."
            )

        self._c_spline.elevate_degree(parametric_dimension)
        logging.debug(
            f"Spline - Elevated {parametric_dimension}.-dim. "
            "degree of the spline."
        )

        self._update_p()

    def reduce_degree(self, parametric_dimension, tolerance=1e-8):
        """
        Tries to reduce degree.

        Parameters
        -----------
        parametric_dimension: int
        tolerance: float

        Returns
        --------
        reduced: bool
        """
        if self.whatami == "Nothing":
            return None

        if parametric_dimension >= self.para_dim:
            raise ValueError(
                "Invalid parametric dimension to reduce degree."
            )

        reduced = self._c_spline.reduce_degree(
            parametric_dimension,
            tolerance
        )

        logging.debug(
            f"Spline - Tried to reduce {parametric_dimension}.-dim. "
            "degree of the spline."
        )

        if reduced:
            logging.debug(
                f"Spline - Successfully reduced {parametric_dimension}.-dim. "
                "degree"
            )
            self._update_p()

        else:
            logging.debug(
                f"Spline - Could not reduce {parametric_dimension}.-dim. "
                "degree"
            )

    def sample(self, query_resolutions):
        """
        Uniformly sample along each parametric dimensions from spline.

        Parameters
        -----------
        query_resolutions: (n,) list-like

        Returns
        --------
        results: (math.product(query_resolutions), dim) np.ndarray
        """
        if self.whatami == "Nothing":
            return None

        query_resolutions = utils.make_c_contiguous(
            query_resolutions,
            dtype="int32",
        ).flatten()

        if query_resolutions.size != self.para_dim:
            raise InputDimensionError(
                "Query resolutions don't match current parametric dimension."
            )

        is_one_or_less = [int(qr) <= 1 for qr in query_resolutions]
        if any(is_one_or_less):
            logging.debug(
                "Spline - You cannot sample less than 2 points per each "
                "parametric dimension."
            )
            logging.debug("Spline - Applying minimum sampling resolution 2.")

            query_resolutions[is_one_or_less] = int(2)

        logging.debug(
            f"Spline - Sampling {np.product(query_resolutions)} "
            "points from spline."
        )

        return self._c_spline.sample(query_resolutions)

    def export(self, fname):
        """
        Export spline. Please be aware of the limits of `.iges`

        Parameters
        -----------
        fname: str

        Returns
        --------
        None
        """
        self._update_c()

        if self.whatami == "Nothing":
            return None

        fname = str(fname)

        dirname = os.path.dirname(fname)
        if not os.path.isdir(dirname) and dirname != "":
            os.makedirs(dirname)

        ext = os.path.splitext(fname)[1]
    
        if ext == ".iges":
            self._c_spline.write_iges(fname)

        elif ext == ".xml":
            self._c_spline.write_xml(fname)

        elif ext == ".itd":
            self._c_spline.write_irit(fname)

        elif ext == ".npz":
            io.npz.write_npz(fname, self)

        elif ext == ".mesh":
            # mfem nurbs mesh.
            io.mfem.write_mfem(fname, self, precision=12)

        else:
            raise ValueError(
                "We can only export "
                + "< .iges | .xml | .itd | .npz | .mesh > extentions"
            )

        logging.info(f"Spline - Exported current spline as {fname}.")

    @abc.abstractmethod
    def copy(self,):
        """
        Returns freshly initialized Spline of self.

        Needs to be implemented separately

        Parameters
        -----------
        None

        Returns
        --------
        new_spline: `Spline`
        """
        pass

    @abc.abstractmethod
    def todict(self,):
        """
        Return current spline as dict

        Parameters
        -----------
        None

        Returns
        --------
        spline_as_dict: `Spline`
        """
        pass

    # member function aliases for urgent situations
    # `weights` is defined in `NURBS`
    ds = degrees
    kvs = knot_vectors
    cps = control_points
