import abc
import logging
import copy
import os
import itertools

import numpy as np

from splinelibpy import utils
from splinelibpy import io

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

        logging.debug("Spline - Degrees set: {d}".format(d=self.degrees))

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
        kvs = self.knot_vectors
        for i in range(len(kvs)):
            logging.debug(
                "Spline -   "
                + str(i)
                + ". knot vector length: "
                + "{kv}".format(kv=len(kvs[i]))
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
                "Spline - {n_cps} Control points set.".format(
                    n_cps=self.control_points.shape[0]
                )
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
        logging.debug(
            "Spline - `lexsort` control points ({l})".format(l=order)
        )
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
    def control_net_resolutions(self,):
        """
        Returns control new resolutions.

        Parameters
        -----------
        None

        Returns
        --------
        control_net_resolutions: list
        """
        cnr = []
        for i in range(self.para_dim):
            cnr.append(
                len(self.knot_vectors[i]) - self.degrees[i] - 1
            )

        return cnr

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

    def evaluate(self, queries):
        """
        Evaluates spline.

        Parameters
        -----------
        queries: (n, para_dim) list-like

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

        return self._c_spline.evaluate(queries=queries)

    def derivative(self, queries, orders):
        """
        Evaluates derivatives of spline.

        Parameters
        -----------
        queries: (n, para_dim) list-like
        orders: (para_dim,) list-like

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

        return self._c_spline.derivative(
            queries=queries,
            orders=orders
        )

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

        logging.debug("Spline - Inserted {nk} knot(s).".format(nk=len(knots)))

        self._update_p()

    def remove_knots(self, parametric_dimension, knots, tolerance=1e-8):
        """
        Tries to removes knots. If you've compiled `splinelibpy` in `Debug`
        and your removal request is not "accepted", you will get an error.
        See the comments for `Nurbs::remove_knots` @ 
        `splinelibpy/src/nurbs.hpp` for more info.

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
            "Spline - Tried to remove {nk} knot(s).".format(nk=len(knots))
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
            "Spline - Elevated {p}.-dim. degree of the spline.".format(
                p=parametric_dimension
            )
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
            "Spline - Tried to reduce {p}.-dim. degree of the spline.".format(
                p=parametric_dimension
            )
        )

        if reduced:
            logging.debug(
                "Spline - Successfully reduced {p}.-dim. degree".format(
                    p=parametric_dimension
                )
            )
            self._update_p()

        else:
            logging.debug(
                "Spline - Could not reduce {p}.-dim. degree".format(
                    p=parametric_dimension
                )
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
                + "parametric dimension."
            )
            logging.debug("Spline - Applying minimum sampling resolution 2.")

            query_resolutions[is_one_or_less] = int(2)

        logging.debug(
            "Spline - Sampling {t} points from spline.".format(
                t=np.product(query_resolutions)
            )
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
            # Save knot vectors as size=1 string array
            # --> to avoid any floating point issues
            # whatami is also size=1 string array
            property_dicts = dict(
                degrees=self.degrees,
                knot_vectors=np.array([str(self.knot_vectors)]),
                control_points=self.control_points,
                whatami=np.array([self.whatami]),
            )
            if self.whatami.startswith("NURBS"):
                property_dicts.update(weights=self.weights)

            np.savez(
                fname,
                **property_dicts,
            )

        elif ext == ".mesh":
            # mfem nurbs mesh.
            io.mfem.write_mfem(self, fname, precision=12)

        else:
            raise ValueError(
                "We can only export "
                + "< .iges | .xml | .itd | .npz | .mesh > extentions"
            )

        logging.info("Spline - Exported current spline as {f}.".format(f=fname))

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

    # member function aliases for urgent situations
    # `weights` is defined in `NURBS`
    ds = degrees
    kvs = knot_vectors
    cps = control_points
