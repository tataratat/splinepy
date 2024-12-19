---
title: 'splinepy : Library for prototyping spline geometries of arbitrary dimensions and degrees and IGA'
tags:
  - Python
  - C++
  - Spline
  - NURBS
  - (rational) Bezier
  - IsoGeometric Analysis
authors:
  - name: Jaewook Lee
    orcid: 0000-0002-9321-4176
    equal-contrib: true
    corresponding: true 
    affiliation: 1
  - name: Jacques Zwar
    orcid: 0009-0000-1285-2364
    equal-contrib: true
    affiliation: 1
  - name: Clemens Fricke
    affiliation: 1
  - name: Stefanie Elgeti
    affiliation: 1
affiliations:
 - name: Institute for Lightweight Design and Structural Biomechanics, TU Wien, Vienna, Austria
   index: 1
date: 15 December 2023
bibliography: paper.bib
---

# Summary


Some random reference [@piegl1996nurbs].

# Statement of need

The field of computational geometry and isogeometric analysis (IGA) has experienced significant advancements in recent years. As a conequence, there is a growing demand for versatile libraries that facilitate the prototyping and analysis of spline-based geometries and IGA-equations in arbitrary dimensions and degrees. We introduce 'splinepy,' an open-source library that efficiently handles tensor product splines, specifically NURBS, B-splines, and (rational) Bezier splines [@piegl2012nurbs,@cohen2001geometric,@rogers2001introduction]. Beyond that, generating and manipulating spline geometries is a crucial task in various fields in engineering and applied sciences, from data interpolation to CAD-design. 

Although many tools for spline creation and manipulation exist to date, there remains a need for a library capable of producing splines in higher parametric dimensions, that provide fast evaluation in multi-query scenarios and are capable of handling large, computational demanding geometries. To this end, splinepy leverages a combination of Python for ease of use and a high-performance C++ backend for computational efficiency, particularly crucial for handling high polynomial degrees and high dimensions. Furthermore, splinepy is tailored towards IGA-applications, thus providing functionality for multipatch geometries - consisting of a high number of splines - like interface information and orientations between splines to facilitate coupling and imposition of boundary conditions.

While other python-based libraries provide similar functionality, for instance `PyNurbs` _citation_ or `NURSB-python` _citation_, they do not provide the required performance to evaluate a high number of spline elements, e.g. for matrix assembly for IGA-prototyping. There are highly specialized implementations for 1D- and 2D-splines of up to degree 3 *TO BE VERIFIED* in `SciPy` _citation_, that can outperform `splinepy`, however, there is no multivariate implementation for higher orders availabe.Other packages like `ndsplines` _citation_ that are implemented with a precompiled backend lack functionality beyond evaluation and fitting. Furthermore, none of these packages seem to be under active development. To the authers' best knowledge, there exists no library today, with n-dimensional splines of arbitrary degrees with a focus on geometry modeling for numerical analysis.

`splinepy` provides an intuitive syntax, making it accessible to a wide audience, including novice programmers and those focused on teaching. At the same time, its computational efficiency and versatility address both research and practical applications, providing valuable functionality for geometry manipulation and preprocessing for analysis such as finite elements or IsoGeometric Analysis [@hughes2005isogeometric, @cottrell2009isogeometric]. `splinepy` is designed for rapid prototyping of complex geometries and offers visualization options to analyze the new geometry and incorporate changes in short time.

# Feature Overview
```Ich dachte wir könnten an dieser Stelle vielleicht so vorgehen, dass wir jeweils das Feature kurz präsentieren und dann ein Bild ein streuseln```

`splinepy` provides four different types of splines, that is (rational) Bézier splines, BSplines and NURBS, all of which are provided up to 10 parametric dimensions in the default build (which is provided via pip) but this can be extended using custom builds. The Bézier splines build on the `bezman` _citation_ C++ library, where as the BSpline and NURBS backend builds on `BSplineLib` _citation_, both of which are developed and maintained by the same group. While some of the features (like arithmetic operations) are only available for Bézier type splines, most of the operations are available for independent of the spline architecture.

### Spline refinement
### Derivatives
### Spline creation, extraction and modification
### Fitting and Interpolation
### Mappings
### Arithmetic operations
### Microstructure construction
### Visualization
### Multipatch geometries and outputs
### Miscellaneous


### 
# State of the field
# Acknowledgement


This belongs into a seperate file:
```
@book{piegl2012nurbs,
  title={The NURBS book},
  author={Piegl, Les and Tiller, Wayne},
  year={2012},
  publisher={Springer Science \& Business Media}
},
@book{cohen2001geometric,
  title={Geometric modeling with splines: an introduction},
  author={Cohen, Elaine and Riesenfeld, Richard F and Elber, Gershon},
  year={2001},
  publisher={CRC Press}
},
@book{rogers2001introduction,
  title={An introduction to NURBS: with historical perspective},
  author={Rogers, David F},
  year={2001},
  publisher={Morgan Kaufmann}
},
@article{hughes2005isogeometric,
  title={Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and mesh refinement},
  author={Hughes, Thomas JR and Cottrell, John A and Bazilevs, Yuri},
  journal={Computer methods in applied mechanics and engineering},
  volume={194},
  number={39-41},
  pages={4135--4195},
  year={2005},
  publisher={Elsevier}
},
@book{cottrell2009isogeometric,
  title={Isogeometric analysis: toward integration of CAD and FEA},
  author={Cottrell, J Austin and Hughes, Thomas JR and Bazilevs, Yuri},
  year={2009},
  publisher={John Wiley \& Sons}
}
```
