#pragma once

#include <iostream>
#include <map>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// SplineLib
#include <Sources/Splines/nurbs.hpp>

// napf
#include <napf.hpp>

#include "arrutils.hpp"
#include "../helpers.hpp"

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
  using ScalarParametricCoordinate_ = typename ParametricCoordinate_::value_type;
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
      napf::SplineLibCoordinatesTree<double,
                                     double,
                                     int,
                                     dim,
                                     1,
                                     CloudT>,
      napf::SplineLibCoordinatesHighDimTree<double,
                                            double,
                                            int,
                                            dim,
                                            1,
                                            CloudT>
  >::type; // takes L1 metric

  // add members
  RasterPoints<double, int, para_dim> pc_sampler_; /* samples paracoord */
  VSCoords cloud_data_; /* evaluated */
  std::unique_ptr<CloudT> cloud_;
  std::unique_ptr<TreeT> tree_;
  bool tree_planted_ = false;

  // update degrees since its size never changes
  void UpdateDegrees(int* ds_buf_ptr) {

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;
    for (int i = 0; i < para_dim; i++) {
      ds_buf_ptr[i] = parameter_space.degrees_[i].Get();
    }
  }

  // update current knot vectors to python
  // since list is mutable, update works
  void UpdateKnotVectors(py::list &p_knot_vectors) {

    // start clean
    p_knot_vectors.attr("clear")();

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;
    for (auto& knotvector : parameter_space.knot_vectors_) {
      auto const &kv = *knotvector; // in
      py::list p_kv; // out
      for (int i = 0; i < kv.GetSize(); i++) {
        auto const &knot = kv[splinelib::Index{i}];
        p_kv.append(knot.Get());
      }
      p_knot_vectors.append(p_kv);
    }

  }

  int GetNCps() {
    WeightedVectorSpace_ const &vector_space = *Base_::weighted_vector_space_;
    return vector_space.GetNumberOfCoordinates();
  }

  // Update cps and weights at the same time since they belong together in
  // weightedvectorspace
  void UpdateControlPointsAndWeights(double* cps_buf_ptr,
                                     double* ws_buf_ptr) {
    WeightedVectorSpace_ const &vector_space = *Base_::weighted_vector_space_;
    int numcps = vector_space.GetNumberOfCoordinates();


    // fill it up, phil!
    for (int i = 0; i < numcps; i++) {
      auto const &coord_named_phil = vector_space[splinelib::Index{i}];
      // phil needs to first project before filling.
      auto const &projected_phil = WeightedVectorSpace_::Project(coord_named_phil);
      for (int j = 0; j < dim; j++) {
        cps_buf_ptr[i * dim + j] = projected_phil[j].Get();
      }
      ws_buf_ptr[i] = coord_named_phil[dim].Get();
    }

  }

  // Computes (degree + 1) X ...
  void BasisFunctionsAndIDs(ParametricCoordinate_ const &parametric_coordinate,
                            double* basis_function_values,
                            int* support_control_point_ids) const {

    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;

    int i = 0;
    double W = 0.0;
    for (Index_ non_zero_basis_function{parameter_space.First()};
         non_zero_basis_function != parameter_space.Behind();
         ++non_zero_basis_function)
    {

      Index_ const &basis_function = (
          parameter_space.FindFirstNonZeroBasisFunction(parametric_coordinate) +
          non_zero_basis_function.GetIndex()
      );

      // general basis fn
      const auto evaluated = parameter_space.EvaluateBasisFunction(
          basis_function,
          parametric_coordinate
      );

      // get `w` and add to `W`
      const auto& support_id = basis_function.GetIndex1d(); // not int yet

      VectorSpace_ const &vector_space = *Base_::weighted_vector_space_;
      const auto& w = vector_space[support_id][dim].Get(); // dim: last elem

      const double N_times_w = evaluated * w;
 
      W += N_times_w;
      basis_function_values[i] = N_times_w; // not yet final
      support_control_point_ids[i] = support_id.Get();
      i++;
    }

    // Loop and divide entries by W
    int end = i;
    double W_inv = 1 / W;
    for (i = 0; i < end; i++) {
      basis_function_values[i] *= W_inv;
    }
  }


  /* [start] Helper functions for `FindParametricCoordinate` */
  void _distance(ParametricCoordinate_& q,
                 double* goal,
                 std::array<double, dim>& dist) {

    auto const &physc = Base_::operator()(q);
    ew_minus(physc, goal, dist);
  }

  void _parametric_bounds(
      std::array<std::array<double, para_dim>, 2>& pbounds) {
    ParameterSpace_ const &parameter_space = *Base_::Base_::parameter_space_;
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

  /* build a kdtree based on given resolution. */
  /* this needs to be buildt before queires! */
  void _newtree(std::array<int, para_dim>& resolutions,
                int nthread) {
    // get parametric bounds and build raster points 
    std::array<std::array<double, para_dim>, 2> bounds;
    _parametric_bounds(bounds);

    pc_sampler_ = RasterPoints<double, int, para_dim>(bounds, resolutions);
    cloud_data_ = VSCoords(new Coordinate_[pc_sampler_.size()]);

    auto prepare_cloud = [&] (int begin, int end) {
      ParametricCoordinate_ pc;
      for (int i{begin}; i < end; i++) {
        pc_sampler_.id_to_paracoord(i, pc);
        cloud_data_[i] = std::move(Base_::operator()(pc));
      }
    };

    // make sure this is not called from a child !
    nthread_execution(prepare_cloud, pc_sampler_.size(), nthread);

    cloud_ = std::unique_ptr<CloudT>(new CloudT(cloud_data_, pc_sampler_.size()));
    tree_ = std::unique_ptr<TreeT>(new TreeT(dim, *cloud_));
    tree_planted_ = true;
  }

  /* goodguess including bounds guess */
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
    }
  }

  void _build_djdu(
      ParametricCoordinate_& guess, /* in */
      std::array<double, dim>& dist, /* in */
      std::array<std::array<double, dim>, para_dim>& eyeders, /* out */
      std::array<double, para_dim>& lhs /* out */) {
    // lhs is djdu
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
      lhs[i] = -2. * tmp; // apply minus here already!
    }
  }

  void _build_d2jdu2(
      ParametricCoordinate_& guess,
      std::array<double, dim>& dist,
      std::array<std::array<double, dim>, para_dim>& eyeders,
      std::array<std::array<double, para_dim>, para_dim>& rhs) {
    
    // prepare AAt of eyeders
    std::array<std::array<double, para_dim>, para_dim> eyedersAAt;
    AAt(eyeders, eyedersAAt);

    Derivative_ derq; // derivative query
    for (int i{0}; i < para_dim; i++) {
      for (int j{0}; j < para_dim; j++) {
        // it results in symmetric matrix.
        // until we implement something special for it,
        // fill in bottom half same as upper half
        if (i > j) continue; 
        derq.fill(splinelib::Derivative{0});
        ++derq[i];
        ++derq[j];
        auto const &der = Base_::operator()(guess, derq);

        rhs[i][j] = 2 * (dot(dist, der) + eyedersAAt[i][j]);

        // fill symmetric part
        if (i != j) rhs[j][i] = rhs[i][j];
      }
    }
  }

  /* [end] Helper functions for `FindParametricCoordinate` */

  void ClosestParametricCoordinate(double* query, /* <- from physical space */
                                   double* para_coord,
                                   int guessoption = 0) {

    // everything we need
    // here, we try nested array
    double tolerance = 1e-10;
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
    clipped.fill(0);
    solverclip.fill(0);
    
    
    // start with some sort of guess
    // TODO: better option than this?
    _good_guess(query,
                guessoption,
                current_guess,
                searchbounds);

    int max_iter = para_dim * 20; /* max newton iter */

    double prevnorm{123456789.}, curnorm;
    /* newton loops */
    for (int i{0}; i < max_iter; i++) {
      _distance(current_guess, query, dist);
      std::cout << "dist===============\n";
      for (int i{0}; i < para_dim; i++) {
        std::cout << dist[i] << " ";
      }
      std::cout << "===============\n";

      _build_djdu(current_guess, dist, eyeders, djdu); /* rhs */
      _build_d2jdu2(current_guess, dist, eyeders, d2jdu2); /* lhs */

      std::cout << "d2jdu2===============\n";
      for (int i{0}; i < para_dim; i++) {
        for (int j{0}; j < para_dim; j++) {
          std::cout << d2jdu2[i][j] << " ";
        }
        std::cout << "\n";
      }
      std::cout << "===============\n";
      std::cout << "djdu===============\n";
      for (int i{0}; i < para_dim; i++) {
        std::cout << djdu[i] << " ";
      }
      std::cout << "===============\n";

      if(prevclipped == clipped) solverclip = clipped;
      gauss_with_pivot(d2jdu2, djdu, solverclip, dx); /* solve */
                                                   // don't trust d2jdu2, djdu,
                                                   // they've been altered
                                                   // inplace.
      std::cout << "d2jdu2===============\n";
      for (int i{0}; i < para_dim; i++) {
        for (int j{0}; j < para_dim; j++) {
          std::cout << d2jdu2[i][j] << " ";
        }
        std::cout << "\n";
      }
      std::cout << "===============\n";
      std::cout << "djdu===============\n";
      for (int i{0}; i < para_dim; i++) {
        std::cout << djdu[i] << " ";
      }
      std::cout << "===============\n";


      ew_iplus(dx, current_guess); /* update */
      for (int i{0}; i < para_dim; i++) {
        std::cout << current_guess[i].Get() << " ";
      }
      std::cout  << "<- current\n";
      for (int i{0}; i < para_dim; i++) {
        std::cout << dx[i] << " ";
      }
      std::cout  << "<- dx\n";

      prevclipped = clipped;
      clip(searchbounds, current_guess, clipped); /* clip */
      for (int i{0}; i < para_dim; i++) {
        std::cout << clipped[i] << " ";
      }
      std::cout  << "<- clipped\n";
      for (int i{0}; i < para_dim; i++) {
        std::cout << current_guess[i].Get() << " ";
      }
      std::cout  << "<- after clip\n";


      // Customer satisfaction survey
      // ============================
      // 1. Converged?
      // 2. All clipped?
      curnorm = norm2(djdu);
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

