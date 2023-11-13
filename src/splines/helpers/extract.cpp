#include "splinepy/splines/helpers/extract.hpp"

namespace splinepy::splines::helpers {

std::vector<std::vector<int>>
ExtractBezierPatchIDs(const std::vector<std::vector<int>>& knot_multiplicities,
                      const int* degrees,
                      const int* n_patches_per_para_dim) {

  const int para_dim = knot_multiplicities.size();

  // Number of total patches and ctps per patch
  int n_total_patches = n_patches_per_para_dim[0];
  int n_ctps_per_patch = degrees[0] + 1;

  // Offsets for start values of individual patches
  std::vector<std::size_t> bezier_index_offsets{};
  bezier_index_offsets.reserve(para_dim);
  bezier_index_offsets.push_back(1);
  for (int i_para_dim{1}; i_para_dim < para_dim; i_para_dim++) {
    n_total_patches *= n_patches_per_para_dim[i_para_dim];
    n_ctps_per_patch *= degrees[i_para_dim] + 1;
    bezier_index_offsets.push_back(bezier_index_offsets[i_para_dim - 1]
                                   * (degrees[i_para_dim - 1] + 1));
  }

  // Init return values
  std::vector<std::size_t> patch_ctp_id_offsets(para_dim, 0);
  std::vector<std::vector<int>> list_of_id_lists(n_total_patches);
  std::vector<int> patch_coords(para_dim, 0);
  for (int i_patch{}; i_patch < n_total_patches; i_patch++) {
    // Determine internal positions in local coord system
    int ii{i_patch};
    // Determine the parameter wise ids of the patch (i.e. the
    // patch-coordinates) and calculate the required index offsets
    for (int i{}; i < para_dim; i++) {
      // ID in spline coordinate system of current patch
      const int patch_coord = static_cast<int>(ii % n_patches_per_para_dim[i]);

      // Coordinate offset of the control points indices
      if (patch_coords[i] != patch_coord) {
        patch_coords[i] = patch_coord;
        patch_ctp_id_offsets[i] =
            (patch_coord == 0)
                ? 0
                : patch_ctp_id_offsets[i] + knot_multiplicities[i][patch_coord];
      }

      ii -= patch_coord;
      ii /= n_patches_per_para_dim[i];
    }

    // Init vectors required for initialization
    std::vector<int>& ids = list_of_id_lists[i_patch];
    ids.reserve(n_ctps_per_patch);

    // Extract relevant coordinates
    for (int i_local_id{}; i_local_id < static_cast<int>(n_ctps_per_patch);
         i_local_id++) {
      int global_id{};
      int n_ctps_in_previous_layers{1};
      // Determine index of local point in global spline
      for (int i_para_dim{}; i_para_dim < para_dim; i_para_dim++) {
        // First id in local system
        const int local_id = (i_local_id / bezier_index_offsets[i_para_dim])
                             % (degrees[i_para_dim] + 1);
        // Add patch offsets
        global_id += (local_id + patch_ctp_id_offsets[i_para_dim])
                     * n_ctps_in_previous_layers;
        // Multiply to index offset
        n_ctps_in_previous_layers *=
            n_patches_per_para_dim[i_para_dim] * degrees[i_para_dim] + 1;
      }
      ids.push_back(global_id);
    }
  }
  return list_of_id_lists;
}

} // namespace splinepy::splines::helpers
