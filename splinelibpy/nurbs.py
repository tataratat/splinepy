import logging
import copy

from splinelibpy import utils
from splinelibpy._spline import Spline

class NURBS(Spline):

    def __init__(self,
            degrees=None,
            knot_vectors=None,
            control_points=None,
            weights=None,
    ):
        """
        NURBS.

        Parameters
        -----------
        degrees: (para_dim,) list-like
        knot_vectors: (para_dim, n) list
        control_points: (m, dim) list-like
        weights: (m,) list-like

        Returns
        --------
        None
        """
        super().__init__(
            degrees=degrees,
            knot_vectors=knot_vectors,
            control_points=control_points,
            weights=weights,
        )

    @property
    def weights(self,):
        """
        Returns weights.

        Parameters
        -----------
        None

        Returns
        --------
        self._weights: (n, 1) list-like
        """
        if not self._is_property("weight"):
            return None

        return self._properties["weight"]

    @weights.setter
    def weights(self, weights):
        """
        Set weights.

        Parameters
        -----------
        weights: (n,) list-like

        Returns
        --------
        None
        """
        if weights is None:
            return None

        weights = utils.make_c_contiguous(
            weights,
            dtype=np.double
        ).reshape(-1,1)

        if self.control_points is not None:
            if self.control_points.shape[0] != weights.shape[0]:
                raise ValueError(
                    "Number of control points and number of weights does not "
                    + "match."
                )

        self._properties["weights"] = weights
        
        logging.debug(
            "Spline - {nws} Weights set.".format(nws=self.weights.shape[0])
        )

        self._update_c()

    def _update_c(self,):
        """
        Updates/Init cpp spline, if it is ready to be updated.
        Checks if all the entries are filled before updating.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        if (
            self.dim is None
            or self.para_dim is None
            or self.degrees is None
            or self.knot_vectors is None
            or self.control_points is None
            or self.weights is None
        ):
            logging.debug(
                "Spline - Not enough information to update cpp spline. "
                + "Skipping update."
            )

            return None

        c_spline_class = f"NURBS{self.para_dim}P{self.dim}D()"
        c_spline = eval(c_spline_class)
        c_spline.knot_vectors = self.knot_vectors.view()
        c_spline.degrees = self.degrees.view()
        c_spline.control_points = self.control_points.view()
        c_spline.weights = self.weights.view()
        self._properties["c_spline"] = c_spline
        self._properties["c_bspline"].update_c()

        logging.debug("Spline - Your spline is {w}.".format(w=self._whatami))

    def _update_p(self,):
        """
        Reads cpp spline and writes it here.
        Probably get an error if cpp isn't ready for this.

        Parameters
        -----------
        None

        Returns
        --------
        None
        """
        self.degrees = self._properties["c_spline"].degrees
        self.knot_vectors = self._properties["c_spline"].knot_vectors
        self.control_points = self._properties["c_spline"].control_points
        self.weights = self._properties["c_spline"].weights
        logging.debug(
            "Spline - Updated python spline. CPP spline and python spline are "
            + "now identical."
        )

    def copy(self,):
        """
        Returns freshly initialized Nurbs of self.

        Parameters
        -----------
        None

        Returns
        --------
        new_nurbs: `NURBS`
        """
        new_nurbs = NURBS()
        new_nurbs._properties = copy.deepcopy(self._properties)

        return new_nurbs

