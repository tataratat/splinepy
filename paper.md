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

The field of computational geometry and isogeometric analysis (IGA) has experienced significant advancements. There is a growing demand for versatile libraries that facilitate the prototyping and analysis of spline geometries in arbitrary dimensions and degrees. We introduce 'splinepy,' an open-source library that efficiently handles tensor product splines, specifically NURBS, B-splines, and (rational) Bezier splines [@piegl2012nurbs,@cohen2001geometric,@rogers2001introduction]. Generating and manipulating spline geometries is a crucial task in various fields in engineering and applied sciences. While many tools in this field exist already, there remains a need for a library that providing fast and versatile evaluation of spline geometries. To this end, splinepy leverages a combination of Python for ease of use and a high-performance C++ backend for computational efficiency, particularly crucial for handling high polynomial degrees.

Splinepy provides an intuitive syntax, making it accessible to a wide audience, including novice programmers and those focused on teaching. At the same time, its computational efficiency and versatility address both research and practical applications, providing valuable functionality for geometry manipulation and preprocessing for analysis such as finite elements or IsoGeometric Analysis [@hughes2005isogeometric, @cottrell2009isogeometric]. splinepy is designed for rapid prototyping of complex geometries and offers visualization options to analyze the new geometry and incorporate changes in short time.

# Applications (Some fancy pictures)
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
