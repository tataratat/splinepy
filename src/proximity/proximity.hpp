#pragma once

#include "arrutils.hpp"
#include "../helpers.hpp"

namespace splinepy::proximity {

// a helper class to perform proximity operations to nurbs.
template<typename SplineType>
class Proximity {
public:

  // Frequently used array alias
  using DArray_ = std::array<double, SplineType::kDim>;
  using PArray_ = std::array<double, SplineType::kParaDim>;
  using PxPMatrix_ = std::array<
      std::array<double, SplineType::kParaDim>, SplineType::kParaDim
  >;
  using PxDMatrix_ = std::array<
      std::array<double, SplineType::kDim>, SplineType::kParaDim
  >;
  using ClipArray_ = std::array<int, SplineType::kParaDim>;


  // constructor.
  Proximity(SplineType const &spline) : spline_(spline) {};

  // kdtree related variables
  GridPoints<double, int, SplineType::kParaDim> grid_points_;
  using Coordinates_ = typename std::unique_ptr<SplineType::Coordinate_[]>;
  using Cloud_ =
      typename napf::CoordinatesCloud<Coordinates_, int, SplineType::kDim>;
  using Tree_ = typename std::conditional<
    (SplineType::kDim < 4),
    napf::CoordinatesTree<double,  /* DataT */
                          double,  /* DistT */
                          int,     /* IndexT */
                          dim,     /* dim */
                          1,       /* metric (L1) */
                          Cloud_>, /* CloudT */
    napf::CoordinatesHighDimTree<double, double, int, dim, 1, Cloud_>
  >::type;

  /* Computes difference between physical query and current parametric guess.
   *
   * Parameters
   * -----------
   * query: in
   * guess: in
   * difference: out
   */
  void QueryMinusGuess(const double* query,
                       const SplineType::ParametricCoordinate_& guess,
                       DArray_& difference) {

    FirstMinusSecondEqualsThird(
        query,
        SplineType::Base_::operator()(guess),
        difference
    );

  }

protected:
  SplineType const &spline_;
  bool tree_planted_ = false;

}

} /* namespace proximity */

namespace py = pybind11;

template<int para_dim, int dim>
class NurbsExt : public splinelib::sources::splines::Nurbs<para_dim, dim> {
public:
  using Base_ = splinelib::sources::splines::Nurbs<para_dim, dim>;
  using Coordinate_ = typename Base_::Coordinate_;
  using ScalarCoordinate_ = typename Coordinate_::value_type;
  using Derivative_ = typename Base_::Derivative_;
  using Knot_ = typename Base_::Knot_;
  using ParameterSpace_ = typename Base_::ParameterSpace_;
  using ParametricCoordinate_ = typename Base_::ParametricCoordinate_;
  using ScalarParametricCoordinate_ =
      typename ParametricCoordinate_::value_type;
  using WeightedVectorSpace_ = typename Base_::WeightedVectorSpace_;
  using OutputInformation_ = splinelib::Tuple<
      typename ParameterSpace_::OutputInformation_,
      typename WeightedVectorSpace_::OutputInformation_
  >;
  using VectorSpace_ = typename WeightedVectorSpace_::Base_; // <dim + 1>

  // Some private ones. Here we make it public
  using Index_ = typename Base_::Base_::Index_;
  using IndexLength_ = typename Index_::Length_;
  using IndexValue_ = typename Index_::Value_;
  using KnotRatios_ = typename Base_::ParameterSpace_::KnotRatios_;
  using Knots_ = typename Base_::Base_::Knots_;
  using HomogeneousBSpline_ = typename Base_::HomogeneousBSpline_;

  // Constructor
  using Base_::Base_;

  // kdtree, for smart initial guess
  using VSCoords = typename std::unique_ptr<Coordinate_[]>;
  using CloudT = napf::VSCoordCloud<VSCoords, int, dim>;
  using TreeT = typename std::conditional<
      (dim < 4),
      napf::SplineLibCoordinatesTree<double,  /* DataT */
                                     double,  /* DistT */
                                     int,     /* IndexT */
                                     dim,     /* dim */
                                     1,       /* metric */
                                     CloudT>, /* CloudT */
      napf::SplineLibCoordinatesHighDimTree<double,
                                            double,
                                            int,
                                            dim,
                                            1,
                                            CloudT>
  >::type; // takes L1 metric

  // add members for parameter coordinate search
  RasterPoints<double, int, para_dim> pc_sampler_; /* sampling locations */
  VSCoords cloud_data_; /* evaluated */
  std::unique_ptr<CloudT> cloud_;
  std::unique_ptr<TreeT> tree_;
  bool tree_planted_ = false;

  /***************************************************************************/
  /******** [start] Helper functions for `FindParametricCoordinate` **********/
  /***************************************************************************/
  /* TODO: Generalize and extract them to use it also for bsplines! */

/* Computes distance between goal and query.
 * physical_goal - spline(parametric_query)
 *
 * Parameters
 * -----------
 * spline: in
 * q: in
 * goal: in
 * dist: out
 */
template<typename SplineT, typename ParaCoordT, int dim>
void distance(SplineT& spline,
              ParaCoordT& q,
              double* goal,
              std::array<double, dim>& dist) {

  auto const &physc = spline(q); // eval
  ew_minus(physcgoal, dist);
}


/* Returns parametric bounds
 *
 * Parameters
 * -----------
 * spline: in
 * pbounds: out
 */
template<typename SplineT, int para_dim>
void parametric_bounds(SplineT& spline,
                       std::array<std::array<double, para_dim>, 2>& pbounds) {

  auto const &parameter_space = *spline::Base_::Base_::parameter_space_;
  int i = 0;
  for (auto& knotvector : parameter_space.knot_vectors_) {
    // first and last knots of each kv are the bounds
    auto const &kf = knotvector->GetFront();
    auto const &kb = knotvector->GetBack();

    pbounds[0][i] = kf.Get();
    pbounds[1][i] = kb.Get();

    i++;
  }
   
}

/* Builds a kdtree based on given resolution.
 * This needs to be built before queires with `goodguess=1`
 * Uses "nthread-for", so no need to be executed by thread-child.
 *
 * Parameters
 * -----------
 * spline: in
 * resolutions: in
 * nthread: in
 * 
 */
template<typename SplineT,
         int para_dim,
         typename PointSamplerT
         typename CloudDataT,
         typename CloudT,
         typename TreeT>
void newtree(SplineT& spline,
             std::array<int, para_dim>& resolutions,
             int nthread,
             PointSamplerT& p_sampler,
             std::unique_ptr<CloudDataT>& cloud_data, // Coordinate
             std::unique_ptr<CloudT>& cloud,
             std::unique_ptr<TreeT>& tree) {

  using ParaCoordT = typename SplineT::ParametricCoordiate_;

  // get parametric bounds and build raster points 
  std::array<std::array<double, para_dim>, 2> bounds;
  parametric_bounds(spline, bounds);

  // create sampler
  p_sampler = PointSamplerT<double, int, para_dim>(bounds, resolutions);

  // cloud_data
  cloud_data = std::unique_ptr<CloudDataT>(new CloudDataT[p_sampler.size()]);
  // lambda to fill cloud_data
  auto prepare_cloud = [&] (int begin, int end) {
    ParaCoordT_ pc;
    for (int i{begin}; i < end; i++) {
      pc_sampler_.id_to_paracoord(i, pc);
      cloud_data_[i] = std::move(Base_::operator()(pc));
    }
  };
  // fill/prepare cloud data
  nthread_execution(prepare_cloud, p_sampler.size(), nthread);

  // create cloud and tree
  cloud = std::unique_ptr<CloudT>(new CloudT(cloud_data, p_sampler_.size()));
  tree = std::unique_ptr<TreeT>(new TreeT(spline.dim_, *cloud));
}

  /* Good guess including bounds guess
   * Options:
   *   0: Gets mid-point of parametric space and parametric bounds.
   *   1: KDTree query for nearest discrete point. 
   *      Currently returns parametric bounbds as bound guess.
   *      However, we could instead give kdtree hit +- tree sample step size.
   *      That'd be aggressive, but would work as long as kdtree is build dense
   *      enough.
   *
   * Parameters
   * -----------
   * goal: in
   * option: in
   * goodguess: out
   * boundguess: out
   *
   * Raises
   * -------
   * runtime_error:
   *   If there's no kdtree planted and option is 1.
   */
  void _good_guess(double* goal,
                   int option,
                   ParametricCoordinate_& goodguess,
                   std::array<std::array<double, para_dim>, 2>& boundguess) {

    /* always start with para_bounds */
    _parametric_bounds(boundguess);

    // return mid point. will serve until better guessers arrive.
    if (option == 0) {
       std::array<double, para_dim> tmpgoodguess;
       ew_mean(boundguess[0], boundguess[1], tmpgoodguess);

       for (int i{0}; i < para_dim; i++) {
         goodguess[i] = ScalarParametricCoordinate_{tmpgoodguess[i]};
       }

    } else if (option == 1) {
      if (!tree_planted_) {
        // hate to be aggresive, but here we raise error
        throw std::runtime_error(
            "to use kd-tree initial guess, please first plant the kd-tree!"
        );
      }
      int vote = 1;
      int id[vote];
      double dist[vote];
      tree_->knnSearch(goal,
                       vote /* only closest */,
                       &id[0],
                       &dist[0]);

      pc_sampler_.id_to_paracoord(id[0], goodguess);

      // this would be a good part to adjust boundguess based on
      // step size!
        
    }
  }

  /* Builds rhs array and returns "eyeders" matrix, so that it can be used in
   * d2jdu2. "eyeders" refer to collection of spline derivatives where only one
   * parametric dimension is set to one and everything else zero, resulting in
   * "eye"-like vectorized derivative query.
   *
   * Parameters
   * -----------
   * guess: in
   *   current parametric coordinate guess
   * dist: in
   *   distance to coordinate goal. Same as result of `_distance(--)`.
   * eyeders: out
   * rhs: out
   */
  void _build_djdu(
      ParametricCoordinate_& guess, /* in */
      std::array<double, dim>& dist, /* in */
      std::array<std::array<double, dim>, para_dim>& eyeders, /* out */
      std::array<double, para_dim>& rhs /* out */) {
    // rhs is djdu
    double tmp;
    int j;
    Derivative_ derq; // derivative query 
    for (int i{0}; i < para_dim; i++) { /* fill */
      // get deriv of order 1.
      derq.fill(splinelib::Derivative{0});
      derq[i] = splinelib::Derivative{1}; 
      auto const &der = Base_::operator()(guess, derq);


      j = 0;
      tmp = 0.;
      for (const auto& d : der) { /* matmul */
        tmp += dist[j] * d.Get();
        eyeders[i][j] = d.Get(); /* needed for 2nd ders. save! */
        j++;
      }
      rhs[i] = -2. * tmp; // apply minus here already!
    }
  }

  /* Builds lhs array. Parameters are same as djdu, except this one only takes
   * eyeders as `in`.
   *
   * Parameters
   * -----------
   * guess: in
   * dist: in
   * eyeders: in
   * lhs: out
   */
  void _build_d2jdu2(
      ParametricCoordinate_& guess,
      std::array<double, dim>& dist,
      std::array<std::array<double, dim>, para_dim>& eyeders,
      std::array<std::array<double, para_dim>, para_dim>& lhs) {
    
    // prepare AAt of eyeders
    std::array<std::array<double, para_dim>, para_dim> eyedersAAt;
    AAt(eyeders, eyedersAAt);

    Derivative_ derq; // derivative query
    for (int i{0}; i < para_dim; i++) {
      for (int j{0}; j < para_dim; j++) {
        // it results in symmetric matrix.
        // until we implement something special for it,
        // fill in bottom half same as upper half
        if (i > j) continue; /* skip bottom half */

        derq.fill(splinelib::Derivative{0});
        ++derq[i];
        ++derq[j];
        auto const &der = Base_::operator()(guess, derq);

        lhs[i][j] = 2 * (dot(dist, der) + eyedersAAt[i][j]);

        // fill symmetric part
        if (i != j) lhs[j][i] = lhs[i][j];
      }
    }
  }

  /***************************************************************************/
  /********** [end] Helper functions for `FindParametricCoordinate` **********/
  /***************************************************************************/

  /* Given physical query coordinate, finds closest parametric coordinate.
   * Take a look at `_good_guess(--)` for initial guess options.
   * 
   * Paramters
   * ----------
   * query: in
   * guessoption: in
   * para_coord: out
   */
  void ClosestParametricCoordinate(double* query, /* <- from physical space */
                                   int guessoption,
                                   double* para_coord) {

    // everything we need
    // here, we try nested array
    double tolerance = 1e-13;
    ParametricCoordinate_ current_guess{};
    std::array<std::array<double, para_dim>, 2> searchbounds{};
    std::array<std::array<double, para_dim>, para_dim> d2jdu2; /* lhs */
    std::array<double, para_dim> djdu; /* rhs */
    std::array<double, para_dim> dx; /* sol */
    std::array<double, dim> dist;
    std::array<std::array<double, dim>, para_dim> eyeders;
    std::array<int, para_dim> clipped;
    std::array<int, para_dim> prevclipped;
    std::array<int, para_dim> solverclip;
    bool solverclipped = false;
    clipped.fill(0);
    prevclipped.fill(0);
    solverclip.fill(0);
    
    
    // start with some sort of guess
    // TODO: better option than this?
    _good_guess(query,
                guessoption,
                current_guess,
                searchbounds);

    constexpr int max_iter = para_dim * 20; /* max newton iter */

    double prevnorm{123456789.}, curnorm;

    /* newton loops */
    for (int i{0}; i < max_iter; i++) {
      /* build system to solve*/
      _distance(current_guess, query, dist);
      _build_djdu(current_guess, dist, eyeders, djdu); /* rhs */
      _build_d2jdu2(current_guess, dist, eyeders, d2jdu2); /* lhs */

      /* debug prints */
      std::cout << "**************** Newton iteration: " << i << " **********"; 
      std::cout << "\n distance (goal - spline(para_c): ";
      for (int i{0}; i < para_dim; i++) {
        std::cout << dist[i] << " ";
      }
      std::cout << "\n";
      std::cout << "System Matrix: A | b  (d2jdu2 | djdu)\n";
      for (int i{0}; i < para_dim; i++) {
        for (int j{0}; j < para_dim; j++) {
          std::cout << d2jdu2[i][j] << " ";
        }
        std::cout << "| " << djdu[i];
        std::cout << "\n";
      }
      /* end debug prints */

      /* solve and update */
      // solver clip only if it is clipped at the same place twice.
      if(prevclipped == clipped && !solverclipped) {
        solverclip = clipped;
        solverclipped = true;
      }
      // solve. pivoting is done inplace, so don't use these afterwards.
      // dx is always reorganized.
      // TODO: why not reorganize d2jdu2 and djdu too?
      gauss_with_pivot(d2jdu2, djdu, solverclip, dx); /* solve */
      ew_iplus(dx, current_guess); /* update */


      /* debug prints */
      std::cout << "Solution (dx): ";
      for (int i{0}; i < para_dim; i++) {
        std::cout << dx[i] << " ";
      }
      std::cout << "\nUpdated: ";
      for (int i{0}; i < para_dim; i++) {
        std::cout << current_guess[i].Get() << " ";
      }
      /* end debug prints */

      /* clip */
      prevclipped = clipped;
      clip(searchbounds, current_guess, clipped); /* clip */


      /* debug prints */
      std::cout << "\nClip info: ";
      for (int i{0}; i < para_dim; i++) {
        std::cout << clipped[i] << " ";
      }
      std::cout  << "\n Current guess after clip:";
      for (int i{0}; i < para_dim; i++) {
        std::cout << current_guess[i].Get() << " ";
      }
      /* end debug prints */

      // Customer satisfaction survey
      // ============================
      // 1. Converged?
      //   -> should be zero, if query point is on the spline.
      // 2. All clipped?
      //   -> should be the case if query point lies outside the spline.
      curnorm = norm2(djdu);
      std::cout  << "\n norm2(djdu): " << curnorm << "\n";
      if (std::abs(prevnorm - curnorm) < tolerance
          || nonzeros(solverclip) == para_dim) break;
      prevnorm = curnorm;
    }

    // fill up the para_coord
    for (int i{0}; i < para_dim; i++) {
      para_coord[i] = current_guess[i].Get();
    }
  }
};
