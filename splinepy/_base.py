"""splinepy/splinepy/_base.py
Useful base class for splinepy.
Plus some useful decorators.
"""

import abc

from splinepy.log import debug, info, warning

class SplinepyBase(abc.ABC):
    """
    Base class for splinepy, where logging is nicely wrapped, and some useful
    methods are defined as classmethods..
    Since attributes are predefined with __slots__, we can pre define
    all the properties that could have been saved.
    Adding `get_` in front of such properties are names for functions that
    freshly compute and save the properties and their byproducts.
    Other more complex operations will be a separate function.
    TODO: maybe add explicit `use_saved` switch to avoid recomputing
    """

    __slots__ = []

    def _logd(self, *log):
        """
        Debug logger wrapper for Mesh.
        Parameters
        -----------
        *log: *str
        Returns
        --------
        None
        """
        debug(type(self).__qualname__, "-", *log)

    def _logi(self, *log):
        """
        Info logger wrapper for Mesh.
        Parameters
        -----------
        *log: *str
        Returns
        --------
        None
        """
        info(type(self).__qualname__, "-", *log)

    def _logw(self, *log):
        """
        Warning logger wrapper for Mesh.
        Parameters
        -----------
        *log: *str
        Returns
        --------
        None
        """
        warning(type(self).__qualname__, "-", *log)