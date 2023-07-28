API References
==============

Overview
--------
- :class:`splinepy.splinepy_core` is a module that includes all the c++ bindings. Most importantly, :class:`.PySpline` and :class:`.PyMultipatch` are implemented here, which serve bases for all the derived spline types with in the library. It only contains minimal documentation, as the classes and methods are derived / integrated throughout the library with proper documentations on python side.
- :class:`.Spline` extends / documents :class:`.PySpline`. Common features for all splines are implemented here. We recommend using :class:`.Bezier`, :class:`.RationalBezier`, :class:`.BSpline`, or :class:`.NURBS` in front-end applications.
- :class:`.BezierBase` extends :class:`.Spline` with features available for Bezier-families. It becomes a basis for :class:`.Bezier` and :class:`.RationalBezier`, where the derived classes implement appropriate initializations.
- :class:`.BSplineBase` extends :class:`.Spline` with features available for BSpline-famlilies. It becomes a basis for :class:`.BSpline` and :class:`.NURBS`, where derived classes implements appropriate initializations.
- :class:`.Multipatch` extends / documents :class:`.PyMultipatch`.
- :class:`splinepy.helpme` module contains useful application-driven tools / utility functions.


Python API
----------
.. toctree::
   python_api
