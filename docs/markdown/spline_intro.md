# A general introduction to splines
A spline is a very flexible and in particular smooth way of representing geometry. As such, splines have enjoyed great success in particular in design, be it engineering design or architecture. The following description is inspired by the books of Rogers as well as Piegl and Tiller (see section further reading for the exact reference).

Splines come in many colors and facets, but there are a few properties that the most common ones (Bézier splines, B-splines, and NURBS) share:

1. Splines are parametric geometry descriptions. This means that each spline coordinate is expressed as a function of a local parameter. Assume that we are dealing with a spline curve, which is embedded in a two-dimensional space. Then the coordinates :math:`(x,y)` are computed as :math:`(x(u),y(u))`, where :math:`u` is the local parameter. An example could be :math:`x(u)=\cos(u),y(u)=\sin(u)`, which represents the unit circle. Usually, the local parameter is normalized, meaning that it takes on the value of :math:`0` at one end of the curve and the value of :math:`1` at the other end of the curve. Notice that parameterizations are not unique and one can represent the same geometry using different parameterizations.
2. The dimension of the geometric entity (i.e., 1D for a curve, 2D for a surface, 3D for a volume, etc.) determines the dimension of the parameter space. This means that a curve always has one local parameter, a curve two, etc.
3. In many cases, one is faced with so-called immersions, i.e., a spline of lower dimension is embedded into a higher-dimensional space, such as a curve into a 2D space. For this reason, the dimension of the parameter space (defined in :math:`\pmb{u}`-coordinates) does not necessarily match the dimension of embedding space, which we will call the global space (defined in :math:`\mathbf{x}`-coordinates).
4. The general equation for a spline curve :math:`\mathbf{S}` is

.. math::
    \mathbf{S} = \sum_{i=1}^n N_i(u) \mathbf{P}

Here, :math:`N` denotes the basis functions, :math:`\mathbf{P}` the coordinates of the control points, and :math:`n` the number of basis functions / control points. Geometrically, this means that one selects a certain number of points in the global space (i.e., the control points) and these points are "connected" using the basis functions. Notice that for splines, the basis functions can take on values within the interval :math:`[0,1]`. However, they will usually never reach the value of one. As a consequence, the control points are not interpolated (meaning that the spline does not go through them), but they only guide the spline.
1. Higher-dimensional splines are derived from 1D-splines in a tensor-product fashion. This means that the control points are then arranged in structured grids and the basis functions are computed as a product of the univariate basis functions :math:`N_i(u)` and :math:`M_j(\eta)`. Consider the example of a splines surface:

.. math::
  \mathbf{S}(u, \eta) = \sum_{i=1}^n \sum_{j=1}^m N_i(u) M_j(\eta) \mathbf{P}\_{i,j}

Notice that the structure of :math:`N_i` and :math:`M_j` is exactly the same. They are simply defined in different coordinate directions.

## In which situations are splines helpful?

As previously mentioned, splines, or in particular non-uniform rational B-spline (NURBS), constitute the basis of computer aided design. However, this is not where the stories of splines ends. They are also very helpful in geometry representation for shape optimization, fitting of data points, but most of all, have become an important means also in numerical analysis. Here, we refer to the idea of isogeometric analysis (IGA), proposed by Hughes, Cottrell and Bazilevs in 2005.

# An introduction to different spline types

As described in the introduction, the principal structure of most spline types is the same. The only difference between the considered spline types are the definition of the basis functions :math:`N_i`. In the following, we will introduce the basis functions for four different spline types and discuss the specific properties they yield. We will restrict ourselves to the univariate basis functions and ask the reader to keep in mind that from there, one can derive the higher-dimensional basis functions as tensor products.

## Bézier splines

Named after Pierre Bézier (1910-1999), a French engineer at Renault and one of the founders of CAD/CAM systems, Bézier splines were one of the earliest spline representations. They use Bernstein polynomials as their basis function, which are defined as follows:

.. math::
  N_{i,p}(u)=\frac{(p)!}{(i-1)!(p+1-i)!}u^{i-1}(1-u)^{p+1-i} ,

with

.. math::
  u \in \[0;1\]; 0^{0}\equiv 1; 0!=1

and the fixed relation between degree of the polynomial :math:`p` and number of control points :math:`n`

.. math::
  n=p+1 .

This means that the number of control points cannot be chosen freely, but is fixed the moment the degree has been chosen. Actually, this is a good moment to point to a terminology conflict that becomes very apparent, when geometers and analysts try to communicate: While a geometer uses the terms `degree` and `order` differently, assuming that order= degree+1, an analyst will tend to use those two terms interchangeably, meaning that order = degree.

After this basic definition, we will continue with a few examples.

#### Examples

**Example 1**: Bézier curve with :math:`n=2` control points and degree :math:`p=1`
In this case, there are two basis functions, both of degree :math:`1`. We denote the basis functions as :math:`N_{i,p}`, where :math:`i` indicates the number of the basis function and :math:`p` the degree of the basis function.

.. math::

    N_{1,1}(u)&=\frac{1!}{0!1!}u^{0}(1-u)^{1}=1-u ,

    N_{2,1}(u)&=\frac{1!}{1!0!}u^{1}(1-u)^{0}=u .

The resulting spline curve can then be computed as follows:

.. math::

  \mathbf{S}(u)=(1-u)\mathbf{P}\_{1}+u \mathbf{P}\_{2} .

**Example 2**: Bézier curve with :math:`n=3` control points and degree :math:`p=2`

\noindent In this case, there are three basis functions, all of degree :math:`2`.

.. math::

    N_{1,2}(u)&=\frac{2!}{0!(2!)}u^{0}(1-u)^{2}=(1-u)^{2}

    N_{2,2}(u)&=\frac{2!}{1!(1!)}u^{1}(1-u)^{1}=2u (1-u)

    N_{3,2}(u)&=\frac{2!}{2!(0!)}u^{1}(1-u)^{0}=u^{2}

The resulting spline curve can then be computed as follows:

.. math::

   \mathbf{S}(u)=(1-u)^{2}\mathbf{P}\_{1}+2u (1-u)\mathbf{P}\_{2}+u^{2}\mathbf{P}\_{3}

#### Properties of Bézier curves

Bézier curves have the following properties:

-  The curve is contained in the convex hull of the control polygon (i.e. the largest convex polygon defined by the CPs).
-  The basis functions are real and (non-zero) throughout the entire definition space.
-  :math:`\mathbf{P}_{1}=\mathbf{S}(0)` and :math:`\mathbf{P}_{n}=\mathbf{S}(1)`. This follows from:

.. math::

  N_{1,p}(0)&=\underbrace{\frac{p!}{0!p!}}\_{=1}\underbrace{0^{0}}\_{=1}\underbrace{(1-0)^{p}}\_{=1}=1

  N_{i,p}(0)&=\frac{p!}{(i-1)!p!}0^{i-1}(1-0)^{p+1-i}=0\qquad\text{if }i\neq 1

  C(0)&=1\cdot P_{1}

A similar argument can be made for :math:`\mathbf{P}_n`.

-  Bernstein polynomials form a partition  of unity.

#### Drawbacks of Bézier curves

In modern CAD systems, Bézier curves no longer play a role. This is due to two drawbacks:

1. Important types of curves cannot be represented using polynomials and thus also not with a Bézier spline, which has a polynomial basis. These include conic sections like circles, cylinders, cones, spheres, etc..
2. The control over the curve geometry is not sufficiently local (single segment curve). A high polynomial degree is required to satisfy a large number of constraints, i.e., complex shapes.


#### Rational Bézier splines

One of the drawbacks of Bézier spline listed above is that it cannot be used to represent certain conic sections exactly. To this end, rational functions are needed. In this case, we will use rational functions defined via the ratio of two polynomials. The new basis function becomes:

.. math::

  R_{i,p}=\frac{N_{i,p}(u)\omega_{i}}{W(u)}=\frac{N_{i,p}(u)\omega_{i}}{\sum_{\hat{i}=1}^{n}N_{\hat{i},p}(u)\omega_{\hat{i}}} .



Here, :math:`N_{i,p}` are the Bernstein polynomials discussed above. :math:`\omega_i` is a so-called weight. :math:`\omega_i` is greater than zero and provides a measure to scale the influence of a certain control point. If the weight tends to zero, the control point is technically excluded from the spline computation. If the weight becomes large, the spline moves very close to the respective control point. If all weights are equal, one retrieves a the Bézier spline (the weights cancel and the denominator than becomes zero due to the partition of unity property.

## B-splines

Just like the Bézier basis functions, B-spline basis functions are polynomials. The key difference is that they are non-zero in only a portion of the parameter space. In other words, the basis functions have local support and this is what, in contrast to Bézier splines, gives us local control over the spline. In order to indicate where the support of each basis function begins and ends, we define a series of break points. These break points are called knots. The collection of knots forms the knot vector :math:`u`:

.. math::
  u=[u_{1}=0,u_{2},\dots,u_{n}=1]\qquad \underbrace{0\le u\le 1}\_{\text{ parameter space}} .



The entries of the knot vector form a non-decreasing sequence: :math:`u_{i}\le u_{i+1}`. It is however possible to repeat knot values. The interval :math:`[u_{i},u_{i+1}]` is called the :math:`i`-th knot span. The length (meaning number of entries) of the knot vector, :math:`m`, is fixed as

.. math::

  &m=n+p+1

  &n:\text{ number of basis functions/CPs}

  &p:\text{ degree} \nonumber



An example of a knot vector could be:

.. math::
  u=[0, 0, 0.5, 1, 1] .



It has a length of :math:`m=5` and contains four knot spans. The first knot span is :math:`[0,0]` with a length of :math:`0`, the second knot span is :math:`[0,0.5]` with a length of :math:`0.5`, etc..

There are different types of knot vectors:


-  open  :math:`\leftrightarrow` periodic

In an open knot vector, the first and last value appears :math:`p+1` times. A knot vector is called periodic, if it is not open. :math:`u=[0,0,0,0.5,1,1,1]` is an open knot vector if :math:`p=2`.
  -  uniform  :math:`\leftrightarrow` non-uniform \\
In a uniform knot vector, the knots are equally spaced in the parameter space, e.g., :math:`u=[0,0.2,0.4,0.6,0.8,1]`. A knot vector is called non-uniform, if this condition is not fulfilled. An exception to this rule are repeated knot values in the beginning or end of the knot vector, e.g., :math:`[0,0,0,0.5,1,1,1]`. This knot vector is then still called uniform, if all middle knots are equally spaced.

Once the knot vector has been specified, the basis functions :math:`N_{i,p}` can be defined. Again, :math:`i` indicates the number of the basis function and :math:`p` the degree.  :math:`N_{i,p}` is defined by the Cox de Boor recurrence formula. The starting point are constant basis functions (:math:`p=0`) defined as follows:

.. math::
   N_{i,0}(u)=
    \begin{cases}
      1,  & \quad u_{i}\le u < u_{i+1} \qquad (i\text{th knot span}), \\\\
      0,  & \quad \text{otherwise.}
    \end{cases}



Subsequently, the degree of the basis function can be raised bit by bit using

.. math::
  N_{i,p}(u)=\frac{u-u_{i}}{u_{i+p}-u_{i}}N_{i,p-1}(u)+\frac{u_{i+p+1}-u}{u_{i+p+1}-u_{i+1}}N_{i+1,p-1}(u) ,



until the desired degree has been reached.


####    Properties of B-spline basis functions:
-  constitute a partition of unity :math:`\sum_{i=1}^{n}N_{i,p}(u)=1`
-  non-negative
-  continuously differentiable within knot spans, :math:`p\cdots` times differentiable across knots (multiplicity of knot)
-  computation of a set of basis functions requires a knot vector :math:`u` and a degree :math:`p`
-  basis functions are non-zero in :math:`[u_{i},u_{i+p+1})`, i.e., in :math:`p+1` knot spans
-  in each knot span, :math:`p+1` basis function will be nonzero
-  if the number of basis functions is :math:`p+1`, the B-spline basis reduces to the Bézier basis




## Non-uniform rational B-splines (NURBS)

In an argument similar to rational Bézier splines, NURBS are the rational counterpart of B-splines. In fact, a NURBS entity in :math:`\mathbb{R}^{d}` is obtained by the projective transformation of a B-spline entity in :math:`\mathbb{R}^{d+1}`.

The NURBS basis functions :math:`R_{i,p}` are computed using the B-spline basis functions :math:`N_{i,p}`:

.. math::
  R_{i,p}=\frac{N_{i,p}(u)\omega_{i}}{W(u)}=\frac{N_{i,p}(u)\omega_{i}}{\sum_{\hat{i}=1}^{n}N_{\hat{i},p}(u)\omega_{\hat{i}}} .



Again, :math:`\omega_i` is a weight.

<!-- # A selection of spline operations offered by splinepy -->

# Further reading

As further reading, we suggest

1. Rogers, David F. An introduction to NURBS: with historical perspective. Morgan Kaufmann, 2001.
2. Piegl, Les, and Wayne Tiller. The NURBS book. Springer Science and Business Media, 1996.
3. Hughes, Thomas JR, John A. Cottrell, and Yuri Bazilevs. Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and mesh refinement. Computer methods in applied mechanics and engineering 194, no. 39-41 (2005): 4135-4195.
