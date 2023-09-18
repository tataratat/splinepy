module SplinepySceneM
   use, intrinsic :: ISO_C_BINDING

   interface

      subroutine SplinepyNearestBoundaryPoint( &
         queries, n_queries, n_threads, results &
         ) bind(c, name="nearest_boundary_point")
         import
         real(c_double), dimension(*), intent(in) :: queries
         integer(c_int), intent(in) :: n_queries
         integer(c_int), intent(in) :: n_threads
         real(c_double), dimension(*), intent(out) :: results
      end subroutine SplinepyNearestBoundaryPoint

      subroutine SplinepyNearestBoundaryPointPrevious( &
         queries, n_queries, n_threads, results &
         ) bind(c, name="nearest_boundary_point_previous")
         import
         real(c_double), dimension(*), intent(in) :: queries
         integer(c_int), intent(in) :: n_queries
         integer(c_int), intent(in) :: n_threads
         real(c_double), dimension(*), intent(out) :: results
      end subroutine SplinepyNearestBoundaryPointPrevious


      subroutine SplinepyIsInSpline( &
         queries, n_queries, n_threads, true_false &
         ) bind(c, name="is_in_spline")
         import
         real(c_double), dimension(*), intent(in) :: queries
         integer(c_int), intent(in) :: n_queries
         integer(c_int), intent(in) :: n_threads
         integer(c_int), dimension(*), intent(out) :: true_false
      end subroutine SplinepyIsInSpline

   end interface

end module SplinepySceneM
