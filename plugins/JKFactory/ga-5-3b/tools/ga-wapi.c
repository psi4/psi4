
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>
#include "ga-papi.h"
#include "typesf2c.h"

static int me;
static int nproc;


static long count_pnga_abs_value = 0;
static long count_pnga_abs_value_patch = 0;
static long count_pnga_acc = 0;
static long count_pnga_access_block_grid_idx = 0;
static long count_pnga_access_block_grid_ptr = 0;
static long count_pnga_access_block_idx = 0;
static long count_pnga_access_block_ptr = 0;
static long count_pnga_access_block_segment_idx = 0;
static long count_pnga_access_block_segment_ptr = 0;
static long count_pnga_access_ghost_element = 0;
static long count_pnga_access_ghost_element_ptr = 0;
static long count_pnga_access_ghost_ptr = 0;
static long count_pnga_access_ghosts = 0;
static long count_pnga_access_idx = 0;
static long count_pnga_access_ptr = 0;
static long count_pnga_add = 0;
static long count_pnga_add_constant = 0;
static long count_pnga_add_constant_patch = 0;
static long count_pnga_add_diagonal = 0;
static long count_pnga_add_patch = 0;
static long count_pnga_allocate = 0;
static long count_pnga_bin_index = 0;
static long count_pnga_bin_sorter = 0;
static long count_pnga_brdcst = 0;
static long count_pnga_check_handle = 0;
static long count_pnga_cluster_nnodes = 0;
static long count_pnga_cluster_nodeid = 0;
static long count_pnga_cluster_nprocs = 0;
static long count_pnga_cluster_proc_nodeid = 0;
static long count_pnga_cluster_procid = 0;
static long count_pnga_comp_patch = 0;
static long count_pnga_compare_distr = 0;
static long count_pnga_copy = 0;
static long count_pnga_copy_patch = 0;
static long count_pnga_copy_patch_dp = 0;
static long count_pnga_create = 0;
static long count_pnga_create_bin_range = 0;
static long count_pnga_create_config = 0;
static long count_pnga_create_ghosts = 0;
static long count_pnga_create_ghosts_config = 0;
static long count_pnga_create_ghosts_irreg = 0;
static long count_pnga_create_ghosts_irreg_config = 0;
static long count_pnga_create_handle = 0;
static long count_pnga_create_irreg = 0;
static long count_pnga_create_irreg_config = 0;
static long count_pnga_create_mutexes = 0;
static long count_pnga_ddot_patch_dp = 0;
static long count_pnga_deregister_type = 0;
static long count_pnga_destroy = 0;
static long count_pnga_destroy_mutexes = 0;
static long count_pnga_diag = 0;
static long count_pnga_diag_reuse = 0;
static long count_pnga_diag_seq = 0;
static long count_pnga_diag_std = 0;
static long count_pnga_diag_std_seq = 0;
static long count_pnga_distribution = 0;
static long count_pnga_dot = 0;
static long count_pnga_dot_patch = 0;
static long count_pnga_duplicate = 0;
static long count_pnga_elem_divide = 0;
static long count_pnga_elem_divide_patch = 0;
static long count_pnga_elem_maximum = 0;
static long count_pnga_elem_maximum_patch = 0;
static long count_pnga_elem_minimum = 0;
static long count_pnga_elem_minimum_patch = 0;
static long count_pnga_elem_multiply = 0;
static long count_pnga_elem_multiply_patch = 0;
static long count_pnga_elem_step_divide_patch = 0;
static long count_pnga_elem_stepb_divide_patch = 0;
static long count_pnga_error = 0;
static long count_pnga_fence = 0;
static long count_pnga_fill = 0;
static long count_pnga_fill_patch = 0;
static long count_pnga_gather = 0;
static long count_pnga_gather2d = 0;
static long count_pnga_get = 0;
static long count_pnga_get_block_info = 0;
static long count_pnga_get_debug = 0;
static long count_pnga_get_diag = 0;
static long count_pnga_get_dimension = 0;
static long count_pnga_get_field = 0;
static long count_pnga_get_ghost_block = 0;
static long count_pnga_get_pgroup = 0;
static long count_pnga_get_pgroup_size = 0;
static long count_pnga_get_proc_grid = 0;
static long count_pnga_get_proc_index = 0;
static long count_pnga_ghost_barrier = 0;
static long count_pnga_gop = 0;
static long count_pnga_has_ghosts = 0;
static long count_pnga_init_fence = 0;
static long count_pnga_initialize = 0;
static long count_pnga_initialize_ltd = 0;
static long count_pnga_inquire = 0;
static long count_pnga_inquire_memory = 0;
static long count_pnga_inquire_name = 0;
static long count_pnga_inquire_type = 0;
static long count_pnga_is_mirrored = 0;
static long count_pnga_list_nodeid = 0;
static long count_pnga_llt_solve = 0;
static long count_pnga_locate = 0;
static long count_pnga_locate_nnodes = 0;
static long count_pnga_locate_num_blocks = 0;
static long count_pnga_locate_region = 0;
static long count_pnga_lock = 0;
static long count_pnga_lu_solve = 0;
static long count_pnga_lu_solve_alt = 0;
static long count_pnga_lu_solve_seq = 0;
static long count_pnga_mask_sync = 0;
static long count_pnga_matmul = 0;
static long count_pnga_matmul_mirrored = 0;
static long count_pnga_matmul_patch = 0;
static long count_pnga_median = 0;
static long count_pnga_median_patch = 0;
static long count_pnga_memory_avail = 0;
static long count_pnga_memory_avail_type = 0;
static long count_pnga_memory_limited = 0;
static long count_pnga_merge_distr_patch = 0;
static long count_pnga_merge_mirrored = 0;
static long count_pnga_msg_brdcst = 0;
static long count_pnga_msg_pgroup_sync = 0;
static long count_pnga_msg_sync = 0;
static long count_pnga_nbacc = 0;
static long count_pnga_nbget = 0;
static long count_pnga_nbget_field = 0;
static long count_pnga_nbget_ghost_dir = 0;
static long count_pnga_nblock = 0;
static long count_pnga_nbput = 0;
static long count_pnga_nbput_field = 0;
static long count_pnga_nbtest = 0;
static long count_pnga_nbwait = 0;
static long count_pnga_ndim = 0;
static long count_pnga_nnodes = 0;
static long count_pnga_nodeid = 0;
static long count_pnga_norm1 = 0;
static long count_pnga_norm_infinity = 0;
static long count_pnga_pack = 0;
static long count_pnga_patch_enum = 0;
static long count_pnga_patch_intersect = 0;
static long count_pnga_periodic = 0;
static long count_pnga_pgroup_absolute_id = 0;
static long count_pnga_pgroup_brdcst = 0;
static long count_pnga_pgroup_create = 0;
static long count_pnga_pgroup_destroy = 0;
static long count_pnga_pgroup_get_default = 0;
static long count_pnga_pgroup_get_mirror = 0;
static long count_pnga_pgroup_get_world = 0;
static long count_pnga_pgroup_gop = 0;
static long count_pnga_pgroup_nnodes = 0;
static long count_pnga_pgroup_nodeid = 0;
static long count_pnga_pgroup_set_default = 0;
static long count_pnga_pgroup_split = 0;
static long count_pnga_pgroup_split_irreg = 0;
static long count_pnga_pgroup_sync = 0;
static long count_pnga_print = 0;
static long count_pnga_print_distribution = 0;
static long count_pnga_print_file = 0;
static long count_pnga_print_patch = 0;
static long count_pnga_print_patch2d = 0;
static long count_pnga_print_patch_file = 0;
static long count_pnga_print_patch_file2d = 0;
static long count_pnga_print_stats = 0;
static long count_pnga_proc_topology = 0;
static long count_pnga_put = 0;
static long count_pnga_put_field = 0;
static long count_pnga_randomize = 0;
static long count_pnga_read_inc = 0;
static long count_pnga_recip = 0;
static long count_pnga_recip_patch = 0;
static long count_pnga_register_type = 0;
static long count_pnga_release = 0;
static long count_pnga_release_block = 0;
static long count_pnga_release_block_grid = 0;
static long count_pnga_release_block_segment = 0;
static long count_pnga_release_ghost_element = 0;
static long count_pnga_release_ghosts = 0;
static long count_pnga_release_update = 0;
static long count_pnga_release_update_block = 0;
static long count_pnga_release_update_block_grid = 0;
static long count_pnga_release_update_block_segment = 0;
static long count_pnga_release_update_ghost_element = 0;
static long count_pnga_release_update_ghosts = 0;
static long count_pnga_scale = 0;
static long count_pnga_scale_cols = 0;
static long count_pnga_scale_patch = 0;
static long count_pnga_scale_rows = 0;
static long count_pnga_scan_add = 0;
static long count_pnga_scan_copy = 0;
static long count_pnga_scatter = 0;
static long count_pnga_scatter2d = 0;
static long count_pnga_scatter_acc = 0;
static long count_pnga_scatter_acc2d = 0;
static long count_pnga_select_elem = 0;
static long count_pnga_set_array_name = 0;
static long count_pnga_set_block_cyclic = 0;
static long count_pnga_set_block_cyclic_proc_grid = 0;
static long count_pnga_set_chunk = 0;
static long count_pnga_set_data = 0;
static long count_pnga_set_debug = 0;
static long count_pnga_set_diagonal = 0;
static long count_pnga_set_ghost_corner_flag = 0;
static long count_pnga_set_ghost_info = 0;
static long count_pnga_set_ghosts = 0;
static long count_pnga_set_irreg_distr = 0;
static long count_pnga_set_irreg_flag = 0;
static long count_pnga_set_memory_limit = 0;
static long count_pnga_set_pgroup = 0;
static long count_pnga_set_restricted = 0;
static long count_pnga_set_restricted_range = 0;
static long count_pnga_set_update4_info = 0;
static long count_pnga_set_update5_info = 0;
static long count_pnga_shift_diagonal = 0;
static long count_pnga_solve = 0;
static long count_pnga_spd_invert = 0;
static long count_pnga_step_bound_info = 0;
static long count_pnga_step_bound_info_patch = 0;
static long count_pnga_step_mask_patch = 0;
static long count_pnga_step_max = 0;
static long count_pnga_step_max_patch = 0;
static long count_pnga_strided_acc = 0;
static long count_pnga_strided_get = 0;
static long count_pnga_strided_put = 0;
static long count_pnga_summarize = 0;
static long count_pnga_symmetrize = 0;
static long count_pnga_sync = 0;
static long count_pnga_terminate = 0;
static long count_pnga_timer = 0;
static long count_pnga_total_blocks = 0;
static long count_pnga_transpose = 0;
static long count_pnga_type_c2f = 0;
static long count_pnga_type_f2c = 0;
static long count_pnga_unlock = 0;
static long count_pnga_unpack = 0;
static long count_pnga_update1_ghosts = 0;
static long count_pnga_update2_ghosts = 0;
static long count_pnga_update3_ghosts = 0;
static long count_pnga_update44_ghosts = 0;
static long count_pnga_update4_ghosts = 0;
static long count_pnga_update55_ghosts = 0;
static long count_pnga_update5_ghosts = 0;
static long count_pnga_update6_ghosts = 0;
static long count_pnga_update7_ghosts = 0;
static long count_pnga_update_ghost_dir = 0;
static long count_pnga_update_ghosts = 0;
static long count_pnga_uses_ma = 0;
static long count_pnga_uses_proc_grid = 0;
static long count_pnga_valid_handle = 0;
static long count_pnga_verify_handle = 0;
static long count_pnga_wtime = 0;
static long count_pnga_zero = 0;
static long count_pnga_zero_diagonal = 0;
static long count_pnga_zero_patch = 0;

static double time_pnga_abs_value = 0;
static double time_pnga_abs_value_patch = 0;
static double time_pnga_acc = 0;
static double time_pnga_access_block_grid_idx = 0;
static double time_pnga_access_block_grid_ptr = 0;
static double time_pnga_access_block_idx = 0;
static double time_pnga_access_block_ptr = 0;
static double time_pnga_access_block_segment_idx = 0;
static double time_pnga_access_block_segment_ptr = 0;
static double time_pnga_access_ghost_element = 0;
static double time_pnga_access_ghost_element_ptr = 0;
static double time_pnga_access_ghost_ptr = 0;
static double time_pnga_access_ghosts = 0;
static double time_pnga_access_idx = 0;
static double time_pnga_access_ptr = 0;
static double time_pnga_add = 0;
static double time_pnga_add_constant = 0;
static double time_pnga_add_constant_patch = 0;
static double time_pnga_add_diagonal = 0;
static double time_pnga_add_patch = 0;
static double time_pnga_allocate = 0;
static double time_pnga_bin_index = 0;
static double time_pnga_bin_sorter = 0;
static double time_pnga_brdcst = 0;
static double time_pnga_check_handle = 0;
static double time_pnga_cluster_nnodes = 0;
static double time_pnga_cluster_nodeid = 0;
static double time_pnga_cluster_nprocs = 0;
static double time_pnga_cluster_proc_nodeid = 0;
static double time_pnga_cluster_procid = 0;
static double time_pnga_comp_patch = 0;
static double time_pnga_compare_distr = 0;
static double time_pnga_copy = 0;
static double time_pnga_copy_patch = 0;
static double time_pnga_copy_patch_dp = 0;
static double time_pnga_create = 0;
static double time_pnga_create_bin_range = 0;
static double time_pnga_create_config = 0;
static double time_pnga_create_ghosts = 0;
static double time_pnga_create_ghosts_config = 0;
static double time_pnga_create_ghosts_irreg = 0;
static double time_pnga_create_ghosts_irreg_config = 0;
static double time_pnga_create_handle = 0;
static double time_pnga_create_irreg = 0;
static double time_pnga_create_irreg_config = 0;
static double time_pnga_create_mutexes = 0;
static double time_pnga_ddot_patch_dp = 0;
static double time_pnga_deregister_type = 0;
static double time_pnga_destroy = 0;
static double time_pnga_destroy_mutexes = 0;
static double time_pnga_diag = 0;
static double time_pnga_diag_reuse = 0;
static double time_pnga_diag_seq = 0;
static double time_pnga_diag_std = 0;
static double time_pnga_diag_std_seq = 0;
static double time_pnga_distribution = 0;
static double time_pnga_dot = 0;
static double time_pnga_dot_patch = 0;
static double time_pnga_duplicate = 0;
static double time_pnga_elem_divide = 0;
static double time_pnga_elem_divide_patch = 0;
static double time_pnga_elem_maximum = 0;
static double time_pnga_elem_maximum_patch = 0;
static double time_pnga_elem_minimum = 0;
static double time_pnga_elem_minimum_patch = 0;
static double time_pnga_elem_multiply = 0;
static double time_pnga_elem_multiply_patch = 0;
static double time_pnga_elem_step_divide_patch = 0;
static double time_pnga_elem_stepb_divide_patch = 0;
static double time_pnga_error = 0;
static double time_pnga_fence = 0;
static double time_pnga_fill = 0;
static double time_pnga_fill_patch = 0;
static double time_pnga_gather = 0;
static double time_pnga_gather2d = 0;
static double time_pnga_get = 0;
static double time_pnga_get_block_info = 0;
static double time_pnga_get_debug = 0;
static double time_pnga_get_diag = 0;
static double time_pnga_get_dimension = 0;
static double time_pnga_get_field = 0;
static double time_pnga_get_ghost_block = 0;
static double time_pnga_get_pgroup = 0;
static double time_pnga_get_pgroup_size = 0;
static double time_pnga_get_proc_grid = 0;
static double time_pnga_get_proc_index = 0;
static double time_pnga_ghost_barrier = 0;
static double time_pnga_gop = 0;
static double time_pnga_has_ghosts = 0;
static double time_pnga_init_fence = 0;
static double time_pnga_initialize = 0;
static double time_pnga_initialize_ltd = 0;
static double time_pnga_inquire = 0;
static double time_pnga_inquire_memory = 0;
static double time_pnga_inquire_name = 0;
static double time_pnga_inquire_type = 0;
static double time_pnga_is_mirrored = 0;
static double time_pnga_list_nodeid = 0;
static double time_pnga_llt_solve = 0;
static double time_pnga_locate = 0;
static double time_pnga_locate_nnodes = 0;
static double time_pnga_locate_num_blocks = 0;
static double time_pnga_locate_region = 0;
static double time_pnga_lock = 0;
static double time_pnga_lu_solve = 0;
static double time_pnga_lu_solve_alt = 0;
static double time_pnga_lu_solve_seq = 0;
static double time_pnga_mask_sync = 0;
static double time_pnga_matmul = 0;
static double time_pnga_matmul_mirrored = 0;
static double time_pnga_matmul_patch = 0;
static double time_pnga_median = 0;
static double time_pnga_median_patch = 0;
static double time_pnga_memory_avail = 0;
static double time_pnga_memory_avail_type = 0;
static double time_pnga_memory_limited = 0;
static double time_pnga_merge_distr_patch = 0;
static double time_pnga_merge_mirrored = 0;
static double time_pnga_msg_brdcst = 0;
static double time_pnga_msg_pgroup_sync = 0;
static double time_pnga_msg_sync = 0;
static double time_pnga_nbacc = 0;
static double time_pnga_nbget = 0;
static double time_pnga_nbget_field = 0;
static double time_pnga_nbget_ghost_dir = 0;
static double time_pnga_nblock = 0;
static double time_pnga_nbput = 0;
static double time_pnga_nbput_field = 0;
static double time_pnga_nbtest = 0;
static double time_pnga_nbwait = 0;
static double time_pnga_ndim = 0;
static double time_pnga_nnodes = 0;
static double time_pnga_nodeid = 0;
static double time_pnga_norm1 = 0;
static double time_pnga_norm_infinity = 0;
static double time_pnga_pack = 0;
static double time_pnga_patch_enum = 0;
static double time_pnga_patch_intersect = 0;
static double time_pnga_periodic = 0;
static double time_pnga_pgroup_absolute_id = 0;
static double time_pnga_pgroup_brdcst = 0;
static double time_pnga_pgroup_create = 0;
static double time_pnga_pgroup_destroy = 0;
static double time_pnga_pgroup_get_default = 0;
static double time_pnga_pgroup_get_mirror = 0;
static double time_pnga_pgroup_get_world = 0;
static double time_pnga_pgroup_gop = 0;
static double time_pnga_pgroup_nnodes = 0;
static double time_pnga_pgroup_nodeid = 0;
static double time_pnga_pgroup_set_default = 0;
static double time_pnga_pgroup_split = 0;
static double time_pnga_pgroup_split_irreg = 0;
static double time_pnga_pgroup_sync = 0;
static double time_pnga_print = 0;
static double time_pnga_print_distribution = 0;
static double time_pnga_print_file = 0;
static double time_pnga_print_patch = 0;
static double time_pnga_print_patch2d = 0;
static double time_pnga_print_patch_file = 0;
static double time_pnga_print_patch_file2d = 0;
static double time_pnga_print_stats = 0;
static double time_pnga_proc_topology = 0;
static double time_pnga_put = 0;
static double time_pnga_put_field = 0;
static double time_pnga_randomize = 0;
static double time_pnga_read_inc = 0;
static double time_pnga_recip = 0;
static double time_pnga_recip_patch = 0;
static double time_pnga_register_type = 0;
static double time_pnga_release = 0;
static double time_pnga_release_block = 0;
static double time_pnga_release_block_grid = 0;
static double time_pnga_release_block_segment = 0;
static double time_pnga_release_ghost_element = 0;
static double time_pnga_release_ghosts = 0;
static double time_pnga_release_update = 0;
static double time_pnga_release_update_block = 0;
static double time_pnga_release_update_block_grid = 0;
static double time_pnga_release_update_block_segment = 0;
static double time_pnga_release_update_ghost_element = 0;
static double time_pnga_release_update_ghosts = 0;
static double time_pnga_scale = 0;
static double time_pnga_scale_cols = 0;
static double time_pnga_scale_patch = 0;
static double time_pnga_scale_rows = 0;
static double time_pnga_scan_add = 0;
static double time_pnga_scan_copy = 0;
static double time_pnga_scatter = 0;
static double time_pnga_scatter2d = 0;
static double time_pnga_scatter_acc = 0;
static double time_pnga_scatter_acc2d = 0;
static double time_pnga_select_elem = 0;
static double time_pnga_set_array_name = 0;
static double time_pnga_set_block_cyclic = 0;
static double time_pnga_set_block_cyclic_proc_grid = 0;
static double time_pnga_set_chunk = 0;
static double time_pnga_set_data = 0;
static double time_pnga_set_debug = 0;
static double time_pnga_set_diagonal = 0;
static double time_pnga_set_ghost_corner_flag = 0;
static double time_pnga_set_ghost_info = 0;
static double time_pnga_set_ghosts = 0;
static double time_pnga_set_irreg_distr = 0;
static double time_pnga_set_irreg_flag = 0;
static double time_pnga_set_memory_limit = 0;
static double time_pnga_set_pgroup = 0;
static double time_pnga_set_restricted = 0;
static double time_pnga_set_restricted_range = 0;
static double time_pnga_set_update4_info = 0;
static double time_pnga_set_update5_info = 0;
static double time_pnga_shift_diagonal = 0;
static double time_pnga_solve = 0;
static double time_pnga_spd_invert = 0;
static double time_pnga_step_bound_info = 0;
static double time_pnga_step_bound_info_patch = 0;
static double time_pnga_step_mask_patch = 0;
static double time_pnga_step_max = 0;
static double time_pnga_step_max_patch = 0;
static double time_pnga_strided_acc = 0;
static double time_pnga_strided_get = 0;
static double time_pnga_strided_put = 0;
static double time_pnga_summarize = 0;
static double time_pnga_symmetrize = 0;
static double time_pnga_sync = 0;
static double time_pnga_terminate = 0;
static double time_pnga_timer = 0;
static double time_pnga_total_blocks = 0;
static double time_pnga_transpose = 0;
static double time_pnga_type_c2f = 0;
static double time_pnga_type_f2c = 0;
static double time_pnga_unlock = 0;
static double time_pnga_unpack = 0;
static double time_pnga_update1_ghosts = 0;
static double time_pnga_update2_ghosts = 0;
static double time_pnga_update3_ghosts = 0;
static double time_pnga_update44_ghosts = 0;
static double time_pnga_update4_ghosts = 0;
static double time_pnga_update55_ghosts = 0;
static double time_pnga_update5_ghosts = 0;
static double time_pnga_update6_ghosts = 0;
static double time_pnga_update7_ghosts = 0;
static double time_pnga_update_ghost_dir = 0;
static double time_pnga_update_ghosts = 0;
static double time_pnga_uses_ma = 0;
static double time_pnga_uses_proc_grid = 0;
static double time_pnga_valid_handle = 0;
static double time_pnga_verify_handle = 0;
static double time_pnga_wtime = 0;
static double time_pnga_zero = 0;
static double time_pnga_zero_diagonal = 0;
static double time_pnga_zero_patch = 0;


void wnga_abs_value(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_abs_value;
    local_start = MPI_Wtime();
    pnga_abs_value(g_a);
    local_stop = MPI_Wtime();
    time_pnga_abs_value += local_stop - local_start;
}


void wnga_abs_value_patch(Integer g_a, Integer *lo, Integer *hi)
{
    double local_start, local_stop;
    ++count_pnga_abs_value_patch;
    local_start = MPI_Wtime();
    pnga_abs_value_patch(g_a, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_abs_value_patch += local_stop - local_start;
}


void wnga_acc(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_acc;
    local_start = MPI_Wtime();
    pnga_acc(g_a, lo, hi, buf, ld, alpha);
    local_stop = MPI_Wtime();
    time_pnga_acc += local_stop - local_start;
}


void wnga_access_block_grid_idx(Integer g_a, Integer *subscript, AccessIndex *index, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_access_block_grid_idx;
    local_start = MPI_Wtime();
    pnga_access_block_grid_idx(g_a, subscript, index, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_block_grid_idx += local_stop - local_start;
}


void wnga_access_block_grid_ptr(Integer g_a, Integer *index, void *ptr, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_access_block_grid_ptr;
    local_start = MPI_Wtime();
    pnga_access_block_grid_ptr(g_a, index, ptr, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_block_grid_ptr += local_stop - local_start;
}


void wnga_access_block_idx(Integer g_a, Integer idx, AccessIndex *index, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_access_block_idx;
    local_start = MPI_Wtime();
    pnga_access_block_idx(g_a, idx, index, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_block_idx += local_stop - local_start;
}


void wnga_access_block_ptr(Integer g_a, Integer idx, void *ptr, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_access_block_ptr;
    local_start = MPI_Wtime();
    pnga_access_block_ptr(g_a, idx, ptr, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_block_ptr += local_stop - local_start;
}


void wnga_access_block_segment_idx(Integer g_a, Integer proc, AccessIndex *index, Integer *len)
{
    double local_start, local_stop;
    ++count_pnga_access_block_segment_idx;
    local_start = MPI_Wtime();
    pnga_access_block_segment_idx(g_a, proc, index, len);
    local_stop = MPI_Wtime();
    time_pnga_access_block_segment_idx += local_stop - local_start;
}


void wnga_access_block_segment_ptr(Integer g_a, Integer proc, void *ptr, Integer *len)
{
    double local_start, local_stop;
    ++count_pnga_access_block_segment_ptr;
    local_start = MPI_Wtime();
    pnga_access_block_segment_ptr(g_a, proc, ptr, len);
    local_stop = MPI_Wtime();
    time_pnga_access_block_segment_ptr += local_stop - local_start;
}


void wnga_access_ghost_element(Integer g_a, AccessIndex *index, Integer subscript[], Integer ld[])
{
    double local_start, local_stop;
    ++count_pnga_access_ghost_element;
    local_start = MPI_Wtime();
    pnga_access_ghost_element(g_a, index, subscript, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_ghost_element += local_stop - local_start;
}


void wnga_access_ghost_element_ptr(Integer g_a, void *ptr, Integer subscript[], Integer ld[])
{
    double local_start, local_stop;
    ++count_pnga_access_ghost_element_ptr;
    local_start = MPI_Wtime();
    pnga_access_ghost_element_ptr(g_a, ptr, subscript, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_ghost_element_ptr += local_stop - local_start;
}


void wnga_access_ghost_ptr(Integer g_a, Integer dims[], void *ptr, Integer ld[])
{
    double local_start, local_stop;
    ++count_pnga_access_ghost_ptr;
    local_start = MPI_Wtime();
    pnga_access_ghost_ptr(g_a, dims, ptr, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_ghost_ptr += local_stop - local_start;
}


void wnga_access_ghosts(Integer g_a, Integer dims[], AccessIndex *index, Integer ld[])
{
    double local_start, local_stop;
    ++count_pnga_access_ghosts;
    local_start = MPI_Wtime();
    pnga_access_ghosts(g_a, dims, index, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_ghosts += local_stop - local_start;
}


void wnga_access_idx(Integer g_a, Integer *lo, Integer *hi, AccessIndex *index, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_access_idx;
    local_start = MPI_Wtime();
    pnga_access_idx(g_a, lo, hi, index, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_idx += local_stop - local_start;
}


void wnga_access_ptr(Integer g_a, Integer *lo, Integer *hi, void *ptr, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_access_ptr;
    local_start = MPI_Wtime();
    pnga_access_ptr(g_a, lo, hi, ptr, ld);
    local_stop = MPI_Wtime();
    time_pnga_access_ptr += local_stop - local_start;
}


void wnga_add(void *alpha, Integer g_a, void *beta, Integer g_b, Integer g_c)
{
    double local_start, local_stop;
    ++count_pnga_add;
    local_start = MPI_Wtime();
    pnga_add(alpha, g_a, beta, g_b, g_c);
    local_stop = MPI_Wtime();
    time_pnga_add += local_stop - local_start;
}


void wnga_add_constant(Integer g_a, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_add_constant;
    local_start = MPI_Wtime();
    pnga_add_constant(g_a, alpha);
    local_stop = MPI_Wtime();
    time_pnga_add_constant += local_stop - local_start;
}


void wnga_add_constant_patch(Integer g_a, Integer *lo, Integer *hi, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_add_constant_patch;
    local_start = MPI_Wtime();
    pnga_add_constant_patch(g_a, lo, hi, alpha);
    local_stop = MPI_Wtime();
    time_pnga_add_constant_patch += local_stop - local_start;
}


void wnga_add_diagonal(Integer g_a, Integer g_v)
{
    double local_start, local_stop;
    ++count_pnga_add_diagonal;
    local_start = MPI_Wtime();
    pnga_add_diagonal(g_a, g_v);
    local_stop = MPI_Wtime();
    time_pnga_add_diagonal += local_stop - local_start;
}


void wnga_add_patch(void *alpha, Integer g_a, Integer *alo, Integer *ahi, void *beta, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_add_patch;
    local_start = MPI_Wtime();
    pnga_add_patch(alpha, g_a, alo, ahi, beta, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_add_patch += local_stop - local_start;
}


logical wnga_allocate(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_allocate;
    local_start = MPI_Wtime();
    return_value = pnga_allocate(g_a);
    local_stop = MPI_Wtime();
    time_pnga_allocate += local_stop - local_start;
    return return_value;
}


void wnga_bin_index(Integer g_bin, Integer g_cnt, Integer g_off, Integer *values, Integer *subs, Integer n, Integer sortit)
{
    double local_start, local_stop;
    ++count_pnga_bin_index;
    local_start = MPI_Wtime();
    pnga_bin_index(g_bin, g_cnt, g_off, values, subs, n, sortit);
    local_stop = MPI_Wtime();
    time_pnga_bin_index += local_stop - local_start;
}


void wnga_bin_sorter(Integer g_bin, Integer g_cnt, Integer g_off)
{
    double local_start, local_stop;
    ++count_pnga_bin_sorter;
    local_start = MPI_Wtime();
    pnga_bin_sorter(g_bin, g_cnt, g_off);
    local_stop = MPI_Wtime();
    time_pnga_bin_sorter += local_stop - local_start;
}


void wnga_brdcst(Integer type, void *buf, Integer len, Integer originator)
{
    double local_start, local_stop;
    ++count_pnga_brdcst;
    local_start = MPI_Wtime();
    pnga_brdcst(type, buf, len, originator);
    local_stop = MPI_Wtime();
    time_pnga_brdcst += local_stop - local_start;
}


void wnga_check_handle(Integer g_a, char *string)
{
    double local_start, local_stop;
    ++count_pnga_check_handle;
    local_start = MPI_Wtime();
    pnga_check_handle(g_a, string);
    local_stop = MPI_Wtime();
    time_pnga_check_handle += local_stop - local_start;
}


Integer wnga_cluster_nnodes()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_cluster_nnodes;
    local_start = MPI_Wtime();
    return_value = pnga_cluster_nnodes();
    local_stop = MPI_Wtime();
    time_pnga_cluster_nnodes += local_stop - local_start;
    return return_value;
}


Integer wnga_cluster_nodeid()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_cluster_nodeid;
    local_start = MPI_Wtime();
    return_value = pnga_cluster_nodeid();
    local_stop = MPI_Wtime();
    time_pnga_cluster_nodeid += local_stop - local_start;
    return return_value;
}


Integer wnga_cluster_nprocs(Integer node)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_cluster_nprocs;
    local_start = MPI_Wtime();
    return_value = pnga_cluster_nprocs(node);
    local_stop = MPI_Wtime();
    time_pnga_cluster_nprocs += local_stop - local_start;
    return return_value;
}


Integer wnga_cluster_proc_nodeid(Integer proc)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_cluster_proc_nodeid;
    local_start = MPI_Wtime();
    return_value = pnga_cluster_proc_nodeid(proc);
    local_stop = MPI_Wtime();
    time_pnga_cluster_proc_nodeid += local_stop - local_start;
    return return_value;
}


Integer wnga_cluster_procid(Integer node, Integer loc_proc_id)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_cluster_procid;
    local_start = MPI_Wtime();
    return_value = pnga_cluster_procid(node, loc_proc_id);
    local_stop = MPI_Wtime();
    time_pnga_cluster_procid += local_stop - local_start;
    return return_value;
}


logical wnga_comp_patch(Integer andim, Integer *alo, Integer *ahi, Integer bndim, Integer *blo, Integer *bhi)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_comp_patch;
    local_start = MPI_Wtime();
    return_value = pnga_comp_patch(andim, alo, ahi, bndim, blo, bhi);
    local_stop = MPI_Wtime();
    time_pnga_comp_patch += local_stop - local_start;
    return return_value;
}


logical wnga_compare_distr(Integer g_a, Integer g_b)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_compare_distr;
    local_start = MPI_Wtime();
    return_value = pnga_compare_distr(g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_compare_distr += local_stop - local_start;
    return return_value;
}


void wnga_copy(Integer g_a, Integer g_b)
{
    double local_start, local_stop;
    ++count_pnga_copy;
    local_start = MPI_Wtime();
    pnga_copy(g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_copy += local_stop - local_start;
}


void wnga_copy_patch(char *trans, Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi)
{
    double local_start, local_stop;
    ++count_pnga_copy_patch;
    local_start = MPI_Wtime();
    pnga_copy_patch(trans, g_a, alo, ahi, g_b, blo, bhi);
    local_stop = MPI_Wtime();
    time_pnga_copy_patch += local_stop - local_start;
}


void wnga_copy_patch_dp(char *t_a, Integer g_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi)
{
    double local_start, local_stop;
    ++count_pnga_copy_patch_dp;
    local_start = MPI_Wtime();
    pnga_copy_patch_dp(t_a, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi);
    local_stop = MPI_Wtime();
    time_pnga_copy_patch_dp += local_stop - local_start;
}


logical wnga_create(Integer type, Integer ndim, Integer *dims, char *name, Integer *chunk, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create;
    local_start = MPI_Wtime();
    return_value = pnga_create(type, ndim, dims, name, chunk, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create += local_stop - local_start;
    return return_value;
}


logical wnga_create_bin_range(Integer g_bin, Integer g_cnt, Integer g_off, Integer *g_range)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_bin_range;
    local_start = MPI_Wtime();
    return_value = pnga_create_bin_range(g_bin, g_cnt, g_off, g_range);
    local_stop = MPI_Wtime();
    time_pnga_create_bin_range += local_stop - local_start;
    return return_value;
}


logical wnga_create_config(Integer type, Integer ndim, Integer *dims, char *name, Integer *chunk, Integer p_handle, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_config;
    local_start = MPI_Wtime();
    return_value = pnga_create_config(type, ndim, dims, name, chunk, p_handle, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create_config += local_stop - local_start;
    return return_value;
}


logical wnga_create_ghosts(Integer type, Integer ndim, Integer *dims, Integer *width, char *name, Integer *chunk, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_create_ghosts(type, ndim, dims, width, name, chunk, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_create_ghosts_config(Integer type, Integer ndim, Integer *dims, Integer *width, char *name, Integer *chunk, Integer p_handle, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_ghosts_config;
    local_start = MPI_Wtime();
    return_value = pnga_create_ghosts_config(type, ndim, dims, width, name, chunk, p_handle, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create_ghosts_config += local_stop - local_start;
    return return_value;
}


logical wnga_create_ghosts_irreg(Integer type, Integer ndim, Integer *dims, Integer *width, char *name, Integer *map, Integer *block, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_ghosts_irreg;
    local_start = MPI_Wtime();
    return_value = pnga_create_ghosts_irreg(type, ndim, dims, width, name, map, block, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create_ghosts_irreg += local_stop - local_start;
    return return_value;
}


logical wnga_create_ghosts_irreg_config(Integer type, Integer ndim, Integer *dims, Integer *width, char *name, Integer *map, Integer *block, Integer p_handle, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_ghosts_irreg_config;
    local_start = MPI_Wtime();
    return_value = pnga_create_ghosts_irreg_config(type, ndim, dims, width, name, map, block, p_handle, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create_ghosts_irreg_config += local_stop - local_start;
    return return_value;
}


Integer wnga_create_handle()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_create_handle;
    local_start = MPI_Wtime();
    return_value = pnga_create_handle();
    local_stop = MPI_Wtime();
    time_pnga_create_handle += local_stop - local_start;
    return return_value;
}


logical wnga_create_irreg(Integer type, Integer ndim, Integer *dims, char *name, Integer *map, Integer *block, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_irreg;
    local_start = MPI_Wtime();
    return_value = pnga_create_irreg(type, ndim, dims, name, map, block, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create_irreg += local_stop - local_start;
    return return_value;
}


logical wnga_create_irreg_config(Integer type, Integer ndim, Integer *dims, char *name, Integer *map, Integer *block, Integer p_handle, Integer *g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_irreg_config;
    local_start = MPI_Wtime();
    return_value = pnga_create_irreg_config(type, ndim, dims, name, map, block, p_handle, g_a);
    local_stop = MPI_Wtime();
    time_pnga_create_irreg_config += local_stop - local_start;
    return return_value;
}


logical wnga_create_mutexes(Integer num)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_create_mutexes;
    local_start = MPI_Wtime();
    return_value = pnga_create_mutexes(num);
    local_stop = MPI_Wtime();
    time_pnga_create_mutexes += local_stop - local_start;
    return return_value;
}


DoublePrecision wnga_ddot_patch_dp(Integer g_a, char *t_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, char *t_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi)
{
    DoublePrecision return_value;
    double local_start, local_stop;
    ++count_pnga_ddot_patch_dp;
    local_start = MPI_Wtime();
    return_value = pnga_ddot_patch_dp(g_a, t_a, ailo, aihi, ajlo, ajhi, g_b, t_b, bilo, bihi, bjlo, bjhi);
    local_stop = MPI_Wtime();
    time_pnga_ddot_patch_dp += local_stop - local_start;
    return return_value;
}


int wnga_deregister_type(int type)
{
    int return_value;
    double local_start, local_stop;
    ++count_pnga_deregister_type;
    local_start = MPI_Wtime();
    return_value = pnga_deregister_type(type);
    local_stop = MPI_Wtime();
    time_pnga_deregister_type += local_stop - local_start;
    return return_value;
}


logical wnga_destroy(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_destroy;
    local_start = MPI_Wtime();
    return_value = pnga_destroy(g_a);
    local_stop = MPI_Wtime();
    time_pnga_destroy += local_stop - local_start;
    return return_value;
}


logical wnga_destroy_mutexes()
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_destroy_mutexes;
    local_start = MPI_Wtime();
    return_value = pnga_destroy_mutexes();
    local_stop = MPI_Wtime();
    time_pnga_destroy_mutexes += local_stop - local_start;
    return return_value;
}


void wnga_diag(Integer g_a, Integer g_s, Integer g_v, DoublePrecision *eval)
{
    double local_start, local_stop;
    ++count_pnga_diag;
    local_start = MPI_Wtime();
    pnga_diag(g_a, g_s, g_v, eval);
    local_stop = MPI_Wtime();
    time_pnga_diag += local_stop - local_start;
}


void wnga_diag_reuse(Integer reuse, Integer g_a, Integer g_s, Integer g_v, DoublePrecision *eval)
{
    double local_start, local_stop;
    ++count_pnga_diag_reuse;
    local_start = MPI_Wtime();
    pnga_diag_reuse(reuse, g_a, g_s, g_v, eval);
    local_stop = MPI_Wtime();
    time_pnga_diag_reuse += local_stop - local_start;
}


void wnga_diag_seq(Integer g_a, Integer g_s, Integer g_v, DoublePrecision *eval)
{
    double local_start, local_stop;
    ++count_pnga_diag_seq;
    local_start = MPI_Wtime();
    pnga_diag_seq(g_a, g_s, g_v, eval);
    local_stop = MPI_Wtime();
    time_pnga_diag_seq += local_stop - local_start;
}


void wnga_diag_std(Integer g_a, Integer g_v, DoublePrecision *eval)
{
    double local_start, local_stop;
    ++count_pnga_diag_std;
    local_start = MPI_Wtime();
    pnga_diag_std(g_a, g_v, eval);
    local_stop = MPI_Wtime();
    time_pnga_diag_std += local_stop - local_start;
}


void wnga_diag_std_seq(Integer g_a, Integer g_v, DoublePrecision *eval)
{
    double local_start, local_stop;
    ++count_pnga_diag_std_seq;
    local_start = MPI_Wtime();
    pnga_diag_std_seq(g_a, g_v, eval);
    local_stop = MPI_Wtime();
    time_pnga_diag_std_seq += local_stop - local_start;
}


void wnga_distribution(Integer g_a, Integer proc, Integer *lo, Integer *hi)
{
    double local_start, local_stop;
    ++count_pnga_distribution;
    local_start = MPI_Wtime();
    pnga_distribution(g_a, proc, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_distribution += local_stop - local_start;
}


void wnga_dot(int type, Integer g_a, Integer g_b, void *value)
{
    double local_start, local_stop;
    ++count_pnga_dot;
    local_start = MPI_Wtime();
    pnga_dot(type, g_a, g_b, value);
    local_stop = MPI_Wtime();
    time_pnga_dot += local_stop - local_start;
}


void wnga_dot_patch(Integer g_a, char *t_a, Integer *alo, Integer *ahi, Integer g_b, char *t_b, Integer *blo, Integer *bhi, void *retval)
{
    double local_start, local_stop;
    ++count_pnga_dot_patch;
    local_start = MPI_Wtime();
    pnga_dot_patch(g_a, t_a, alo, ahi, g_b, t_b, blo, bhi, retval);
    local_stop = MPI_Wtime();
    time_pnga_dot_patch += local_stop - local_start;
}


logical wnga_duplicate(Integer g_a, Integer *g_b, char *array_name)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_duplicate;
    local_start = MPI_Wtime();
    return_value = pnga_duplicate(g_a, g_b, array_name);
    local_stop = MPI_Wtime();
    time_pnga_duplicate += local_stop - local_start;
    return return_value;
}


void wnga_elem_divide(Integer g_a, Integer g_b, Integer g_c)
{
    double local_start, local_stop;
    ++count_pnga_elem_divide;
    local_start = MPI_Wtime();
    pnga_elem_divide(g_a, g_b, g_c);
    local_stop = MPI_Wtime();
    time_pnga_elem_divide += local_stop - local_start;
}


void wnga_elem_divide_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_elem_divide_patch;
    local_start = MPI_Wtime();
    pnga_elem_divide_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_elem_divide_patch += local_stop - local_start;
}


void wnga_elem_maximum(Integer g_a, Integer g_b, Integer g_c)
{
    double local_start, local_stop;
    ++count_pnga_elem_maximum;
    local_start = MPI_Wtime();
    pnga_elem_maximum(g_a, g_b, g_c);
    local_stop = MPI_Wtime();
    time_pnga_elem_maximum += local_stop - local_start;
}


void wnga_elem_maximum_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_elem_maximum_patch;
    local_start = MPI_Wtime();
    pnga_elem_maximum_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_elem_maximum_patch += local_stop - local_start;
}


void wnga_elem_minimum(Integer g_a, Integer g_b, Integer g_c)
{
    double local_start, local_stop;
    ++count_pnga_elem_minimum;
    local_start = MPI_Wtime();
    pnga_elem_minimum(g_a, g_b, g_c);
    local_stop = MPI_Wtime();
    time_pnga_elem_minimum += local_stop - local_start;
}


void wnga_elem_minimum_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_elem_minimum_patch;
    local_start = MPI_Wtime();
    pnga_elem_minimum_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_elem_minimum_patch += local_stop - local_start;
}


void wnga_elem_multiply(Integer g_a, Integer g_b, Integer g_c)
{
    double local_start, local_stop;
    ++count_pnga_elem_multiply;
    local_start = MPI_Wtime();
    pnga_elem_multiply(g_a, g_b, g_c);
    local_stop = MPI_Wtime();
    time_pnga_elem_multiply += local_stop - local_start;
}


void wnga_elem_multiply_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_elem_multiply_patch;
    local_start = MPI_Wtime();
    pnga_elem_multiply_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_elem_multiply_patch += local_stop - local_start;
}


void wnga_elem_step_divide_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_elem_step_divide_patch;
    local_start = MPI_Wtime();
    pnga_elem_step_divide_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_elem_step_divide_patch += local_stop - local_start;
}


void wnga_elem_stepb_divide_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_elem_stepb_divide_patch;
    local_start = MPI_Wtime();
    pnga_elem_stepb_divide_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_elem_stepb_divide_patch += local_stop - local_start;
}


void wnga_error(char *string, Integer icode)
{
    double local_start, local_stop;
    ++count_pnga_error;
    local_start = MPI_Wtime();
    pnga_error(string, icode);
    local_stop = MPI_Wtime();
    time_pnga_error += local_stop - local_start;
}


void wnga_fence()
{
    double local_start, local_stop;
    ++count_pnga_fence;
    local_start = MPI_Wtime();
    pnga_fence();
    local_stop = MPI_Wtime();
    time_pnga_fence += local_stop - local_start;
}


void wnga_fill(Integer g_a, void *val)
{
    double local_start, local_stop;
    ++count_pnga_fill;
    local_start = MPI_Wtime();
    pnga_fill(g_a, val);
    local_stop = MPI_Wtime();
    time_pnga_fill += local_stop - local_start;
}


void wnga_fill_patch(Integer g_a, Integer *lo, Integer *hi, void *val)
{
    double local_start, local_stop;
    ++count_pnga_fill_patch;
    local_start = MPI_Wtime();
    pnga_fill_patch(g_a, lo, hi, val);
    local_stop = MPI_Wtime();
    time_pnga_fill_patch += local_stop - local_start;
}


void wnga_gather(Integer g_a, void *v, Integer subscript[], Integer nv)
{
    double local_start, local_stop;
    ++count_pnga_gather;
    local_start = MPI_Wtime();
    pnga_gather(g_a, v, subscript, nv);
    local_stop = MPI_Wtime();
    time_pnga_gather += local_stop - local_start;
}


void wnga_gather2d(Integer g_a, void *v, Integer *i, Integer *j, Integer nv)
{
    double local_start, local_stop;
    ++count_pnga_gather2d;
    local_start = MPI_Wtime();
    pnga_gather2d(g_a, v, i, j, nv);
    local_stop = MPI_Wtime();
    time_pnga_gather2d += local_stop - local_start;
}


void wnga_get(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_get;
    local_start = MPI_Wtime();
    pnga_get(g_a, lo, hi, buf, ld);
    local_stop = MPI_Wtime();
    time_pnga_get += local_stop - local_start;
}


void wnga_get_block_info(Integer g_a, Integer *num_blocks, Integer *block_dims)
{
    double local_start, local_stop;
    ++count_pnga_get_block_info;
    local_start = MPI_Wtime();
    pnga_get_block_info(g_a, num_blocks, block_dims);
    local_stop = MPI_Wtime();
    time_pnga_get_block_info += local_stop - local_start;
}


logical wnga_get_debug()
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_get_debug;
    local_start = MPI_Wtime();
    return_value = pnga_get_debug();
    local_stop = MPI_Wtime();
    time_pnga_get_debug += local_stop - local_start;
    return return_value;
}


void wnga_get_diag(Integer g_a, Integer g_v)
{
    double local_start, local_stop;
    ++count_pnga_get_diag;
    local_start = MPI_Wtime();
    pnga_get_diag(g_a, g_v);
    local_stop = MPI_Wtime();
    time_pnga_get_diag += local_stop - local_start;
}


Integer wnga_get_dimension(Integer g_a)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_get_dimension;
    local_start = MPI_Wtime();
    return_value = pnga_get_dimension(g_a);
    local_stop = MPI_Wtime();
    time_pnga_get_dimension += local_stop - local_start;
    return return_value;
}


void wnga_get_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize, void *buf, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_get_field;
    local_start = MPI_Wtime();
    pnga_get_field(g_a, lo, hi, foff, fsize, buf, ld);
    local_stop = MPI_Wtime();
    time_pnga_get_field += local_stop - local_start;
}


void wnga_get_ghost_block(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_get_ghost_block;
    local_start = MPI_Wtime();
    pnga_get_ghost_block(g_a, lo, hi, buf, ld);
    local_stop = MPI_Wtime();
    time_pnga_get_ghost_block += local_stop - local_start;
}


Integer wnga_get_pgroup(Integer g_a)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_get_pgroup;
    local_start = MPI_Wtime();
    return_value = pnga_get_pgroup(g_a);
    local_stop = MPI_Wtime();
    time_pnga_get_pgroup += local_stop - local_start;
    return return_value;
}


Integer wnga_get_pgroup_size(Integer grp_id)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_get_pgroup_size;
    local_start = MPI_Wtime();
    return_value = pnga_get_pgroup_size(grp_id);
    local_stop = MPI_Wtime();
    time_pnga_get_pgroup_size += local_stop - local_start;
    return return_value;
}


void wnga_get_proc_grid(Integer g_a, Integer *dims)
{
    double local_start, local_stop;
    ++count_pnga_get_proc_grid;
    local_start = MPI_Wtime();
    pnga_get_proc_grid(g_a, dims);
    local_stop = MPI_Wtime();
    time_pnga_get_proc_grid += local_stop - local_start;
}


void wnga_get_proc_index(Integer g_a, Integer iproc, Integer *index)
{
    double local_start, local_stop;
    ++count_pnga_get_proc_index;
    local_start = MPI_Wtime();
    pnga_get_proc_index(g_a, iproc, index);
    local_stop = MPI_Wtime();
    time_pnga_get_proc_index += local_stop - local_start;
}


void wnga_ghost_barrier()
{
    double local_start, local_stop;
    ++count_pnga_ghost_barrier;
    local_start = MPI_Wtime();
    pnga_ghost_barrier();
    local_stop = MPI_Wtime();
    time_pnga_ghost_barrier += local_stop - local_start;
}


void wnga_gop(Integer type, void *x, Integer n, char *op)
{
    double local_start, local_stop;
    ++count_pnga_gop;
    local_start = MPI_Wtime();
    pnga_gop(type, x, n, op);
    local_stop = MPI_Wtime();
    time_pnga_gop += local_stop - local_start;
}


logical wnga_has_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_has_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_has_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_has_ghosts += local_stop - local_start;
    return return_value;
}


void wnga_init_fence()
{
    double local_start, local_stop;
    ++count_pnga_init_fence;
    local_start = MPI_Wtime();
    pnga_init_fence();
    local_stop = MPI_Wtime();
    time_pnga_init_fence += local_stop - local_start;
}


void wnga_inquire(Integer g_a, Integer *type, Integer *ndim, Integer *dims)
{
    double local_start, local_stop;
    ++count_pnga_inquire;
    local_start = MPI_Wtime();
    pnga_inquire(g_a, type, ndim, dims);
    local_stop = MPI_Wtime();
    time_pnga_inquire += local_stop - local_start;
}


Integer wnga_inquire_memory()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_inquire_memory;
    local_start = MPI_Wtime();
    return_value = pnga_inquire_memory();
    local_stop = MPI_Wtime();
    time_pnga_inquire_memory += local_stop - local_start;
    return return_value;
}


void wnga_inquire_name(Integer g_a, char **array_name)
{
    double local_start, local_stop;
    ++count_pnga_inquire_name;
    local_start = MPI_Wtime();
    pnga_inquire_name(g_a, array_name);
    local_stop = MPI_Wtime();
    time_pnga_inquire_name += local_stop - local_start;
}


void wnga_inquire_type(Integer g_a, Integer *type)
{
    double local_start, local_stop;
    ++count_pnga_inquire_type;
    local_start = MPI_Wtime();
    pnga_inquire_type(g_a, type);
    local_stop = MPI_Wtime();
    time_pnga_inquire_type += local_stop - local_start;
}


logical wnga_is_mirrored(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_is_mirrored;
    local_start = MPI_Wtime();
    return_value = pnga_is_mirrored(g_a);
    local_stop = MPI_Wtime();
    time_pnga_is_mirrored += local_stop - local_start;
    return return_value;
}


void wnga_list_nodeid(Integer *list, Integer nprocs)
{
    double local_start, local_stop;
    ++count_pnga_list_nodeid;
    local_start = MPI_Wtime();
    pnga_list_nodeid(list, nprocs);
    local_stop = MPI_Wtime();
    time_pnga_list_nodeid += local_stop - local_start;
}


Integer wnga_llt_solve(Integer g_a, Integer g_b)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_llt_solve;
    local_start = MPI_Wtime();
    return_value = pnga_llt_solve(g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_llt_solve += local_stop - local_start;
    return return_value;
}


logical wnga_locate(Integer g_a, Integer *subscript, Integer *owner)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_locate;
    local_start = MPI_Wtime();
    return_value = pnga_locate(g_a, subscript, owner);
    local_stop = MPI_Wtime();
    time_pnga_locate += local_stop - local_start;
    return return_value;
}


logical wnga_locate_nnodes(Integer g_a, Integer *lo, Integer *hi, Integer *np)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_locate_nnodes;
    local_start = MPI_Wtime();
    return_value = pnga_locate_nnodes(g_a, lo, hi, np);
    local_stop = MPI_Wtime();
    time_pnga_locate_nnodes += local_stop - local_start;
    return return_value;
}


Integer wnga_locate_num_blocks(Integer g_a, Integer *lo, Integer *hi)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_locate_num_blocks;
    local_start = MPI_Wtime();
    return_value = pnga_locate_num_blocks(g_a, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_locate_num_blocks += local_stop - local_start;
    return return_value;
}


logical wnga_locate_region(Integer g_a, Integer *lo, Integer *hi, Integer *map, Integer *proclist, Integer *np)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_locate_region;
    local_start = MPI_Wtime();
    return_value = pnga_locate_region(g_a, lo, hi, map, proclist, np);
    local_stop = MPI_Wtime();
    time_pnga_locate_region += local_stop - local_start;
    return return_value;
}


void wnga_lock(Integer mutex)
{
    double local_start, local_stop;
    ++count_pnga_lock;
    local_start = MPI_Wtime();
    pnga_lock(mutex);
    local_stop = MPI_Wtime();
    time_pnga_lock += local_stop - local_start;
}


void wnga_lu_solve(char *tran, Integer g_a, Integer g_b)
{
    double local_start, local_stop;
    ++count_pnga_lu_solve;
    local_start = MPI_Wtime();
    pnga_lu_solve(tran, g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_lu_solve += local_stop - local_start;
}


void wnga_lu_solve_alt(Integer tran, Integer g_a, Integer g_b)
{
    double local_start, local_stop;
    ++count_pnga_lu_solve_alt;
    local_start = MPI_Wtime();
    pnga_lu_solve_alt(tran, g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_lu_solve_alt += local_stop - local_start;
}


void wnga_lu_solve_seq(char *trans, Integer g_a, Integer g_b)
{
    double local_start, local_stop;
    ++count_pnga_lu_solve_seq;
    local_start = MPI_Wtime();
    pnga_lu_solve_seq(trans, g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_lu_solve_seq += local_stop - local_start;
}


void wnga_mask_sync(Integer begin, Integer end)
{
    double local_start, local_stop;
    ++count_pnga_mask_sync;
    local_start = MPI_Wtime();
    pnga_mask_sync(begin, end);
    local_stop = MPI_Wtime();
    time_pnga_mask_sync += local_stop - local_start;
}


void wnga_matmul(char *transa, char *transb, void *alpha, void *beta, Integer g_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi, Integer g_c, Integer cilo, Integer cihi, Integer cjlo, Integer cjhi)
{
    double local_start, local_stop;
    ++count_pnga_matmul;
    local_start = MPI_Wtime();
    pnga_matmul(transa, transb, alpha, beta, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi, g_c, cilo, cihi, cjlo, cjhi);
    local_stop = MPI_Wtime();
    time_pnga_matmul += local_stop - local_start;
}


void wnga_matmul_mirrored(char *transa, char *transb, void *alpha, void *beta, Integer g_a, Integer ailo, Integer aihi, Integer ajlo, Integer ajhi, Integer g_b, Integer bilo, Integer bihi, Integer bjlo, Integer bjhi, Integer g_c, Integer cilo, Integer cihi, Integer cjlo, Integer cjhi)
{
    double local_start, local_stop;
    ++count_pnga_matmul_mirrored;
    local_start = MPI_Wtime();
    pnga_matmul_mirrored(transa, transb, alpha, beta, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi, g_c, cilo, cihi, cjlo, cjhi);
    local_stop = MPI_Wtime();
    time_pnga_matmul_mirrored += local_stop - local_start;
}


void wnga_matmul_patch(char *transa, char *transb, void *alpha, void *beta, Integer g_a, Integer alo[], Integer ahi[], Integer g_b, Integer blo[], Integer bhi[], Integer g_c, Integer clo[], Integer chi[])
{
    double local_start, local_stop;
    ++count_pnga_matmul_patch;
    local_start = MPI_Wtime();
    pnga_matmul_patch(transa, transb, alpha, beta, g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_matmul_patch += local_stop - local_start;
}


void wnga_median(Integer g_a, Integer g_b, Integer g_c, Integer g_m)
{
    double local_start, local_stop;
    ++count_pnga_median;
    local_start = MPI_Wtime();
    pnga_median(g_a, g_b, g_c, g_m);
    local_stop = MPI_Wtime();
    time_pnga_median += local_stop - local_start;
}


void wnga_median_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi, Integer g_m, Integer *mlo, Integer *mhi)
{
    double local_start, local_stop;
    ++count_pnga_median_patch;
    local_start = MPI_Wtime();
    pnga_median_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi, g_m, mlo, mhi);
    local_stop = MPI_Wtime();
    time_pnga_median_patch += local_stop - local_start;
}


Integer wnga_memory_avail()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_memory_avail;
    local_start = MPI_Wtime();
    return_value = pnga_memory_avail();
    local_stop = MPI_Wtime();
    time_pnga_memory_avail += local_stop - local_start;
    return return_value;
}


Integer wnga_memory_avail_type(Integer datatype)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_memory_avail_type;
    local_start = MPI_Wtime();
    return_value = pnga_memory_avail_type(datatype);
    local_stop = MPI_Wtime();
    time_pnga_memory_avail_type += local_stop - local_start;
    return return_value;
}


logical wnga_memory_limited()
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_memory_limited;
    local_start = MPI_Wtime();
    return_value = pnga_memory_limited();
    local_stop = MPI_Wtime();
    time_pnga_memory_limited += local_stop - local_start;
    return return_value;
}


void wnga_merge_distr_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi)
{
    double local_start, local_stop;
    ++count_pnga_merge_distr_patch;
    local_start = MPI_Wtime();
    pnga_merge_distr_patch(g_a, alo, ahi, g_b, blo, bhi);
    local_stop = MPI_Wtime();
    time_pnga_merge_distr_patch += local_stop - local_start;
}


void wnga_merge_mirrored(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_merge_mirrored;
    local_start = MPI_Wtime();
    pnga_merge_mirrored(g_a);
    local_stop = MPI_Wtime();
    time_pnga_merge_mirrored += local_stop - local_start;
}


void wnga_msg_brdcst(Integer type, void *buffer, Integer len, Integer root)
{
    double local_start, local_stop;
    ++count_pnga_msg_brdcst;
    local_start = MPI_Wtime();
    pnga_msg_brdcst(type, buffer, len, root);
    local_stop = MPI_Wtime();
    time_pnga_msg_brdcst += local_stop - local_start;
}


void wnga_msg_pgroup_sync(Integer grp_id)
{
    double local_start, local_stop;
    ++count_pnga_msg_pgroup_sync;
    local_start = MPI_Wtime();
    pnga_msg_pgroup_sync(grp_id);
    local_stop = MPI_Wtime();
    time_pnga_msg_pgroup_sync += local_stop - local_start;
}


void wnga_msg_sync()
{
    double local_start, local_stop;
    ++count_pnga_msg_sync;
    local_start = MPI_Wtime();
    pnga_msg_sync();
    local_stop = MPI_Wtime();
    time_pnga_msg_sync += local_stop - local_start;
}


void wnga_nbacc(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, void *alpha, Integer *nbhndl)
{
    double local_start, local_stop;
    ++count_pnga_nbacc;
    local_start = MPI_Wtime();
    pnga_nbacc(g_a, lo, hi, buf, ld, alpha, nbhndl);
    local_stop = MPI_Wtime();
    time_pnga_nbacc += local_stop - local_start;
}


void wnga_nbget(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, Integer *nbhandle)
{
    double local_start, local_stop;
    ++count_pnga_nbget;
    local_start = MPI_Wtime();
    pnga_nbget(g_a, lo, hi, buf, ld, nbhandle);
    local_stop = MPI_Wtime();
    time_pnga_nbget += local_stop - local_start;
}


void wnga_nbget_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize, void *buf, Integer *ld, Integer *nbhandle)
{
    double local_start, local_stop;
    ++count_pnga_nbget_field;
    local_start = MPI_Wtime();
    pnga_nbget_field(g_a, lo, hi, foff, fsize, buf, ld, nbhandle);
    local_stop = MPI_Wtime();
    time_pnga_nbget_field += local_stop - local_start;
}


void wnga_nbget_ghost_dir(Integer g_a, Integer *mask, Integer *nbhandle)
{
    double local_start, local_stop;
    ++count_pnga_nbget_ghost_dir;
    local_start = MPI_Wtime();
    pnga_nbget_ghost_dir(g_a, mask, nbhandle);
    local_stop = MPI_Wtime();
    time_pnga_nbget_ghost_dir += local_stop - local_start;
}


void wnga_nblock(Integer g_a, Integer *nblock)
{
    double local_start, local_stop;
    ++count_pnga_nblock;
    local_start = MPI_Wtime();
    pnga_nblock(g_a, nblock);
    local_stop = MPI_Wtime();
    time_pnga_nblock += local_stop - local_start;
}


void wnga_nbput(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, Integer *nbhandle)
{
    double local_start, local_stop;
    ++count_pnga_nbput;
    local_start = MPI_Wtime();
    pnga_nbput(g_a, lo, hi, buf, ld, nbhandle);
    local_stop = MPI_Wtime();
    time_pnga_nbput += local_stop - local_start;
}


void wnga_nbput_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize, void *buf, Integer *ld, Integer *nbhandle)
{
    double local_start, local_stop;
    ++count_pnga_nbput_field;
    local_start = MPI_Wtime();
    pnga_nbput_field(g_a, lo, hi, foff, fsize, buf, ld, nbhandle);
    local_stop = MPI_Wtime();
    time_pnga_nbput_field += local_stop - local_start;
}


Integer wnga_nbtest(Integer *nbhandle)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_nbtest;
    local_start = MPI_Wtime();
    return_value = pnga_nbtest(nbhandle);
    local_stop = MPI_Wtime();
    time_pnga_nbtest += local_stop - local_start;
    return return_value;
}


void wnga_nbwait(Integer *nbhandle)
{
    double local_start, local_stop;
    ++count_pnga_nbwait;
    local_start = MPI_Wtime();
    pnga_nbwait(nbhandle);
    local_stop = MPI_Wtime();
    time_pnga_nbwait += local_stop - local_start;
}


Integer wnga_ndim(Integer g_a)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_ndim;
    local_start = MPI_Wtime();
    return_value = pnga_ndim(g_a);
    local_stop = MPI_Wtime();
    time_pnga_ndim += local_stop - local_start;
    return return_value;
}


Integer wnga_nnodes()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_nnodes;
    local_start = MPI_Wtime();
    return_value = pnga_nnodes();
    local_stop = MPI_Wtime();
    time_pnga_nnodes += local_stop - local_start;
    return return_value;
}


Integer wnga_nodeid()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_nodeid;
    local_start = MPI_Wtime();
    return_value = pnga_nodeid();
    local_stop = MPI_Wtime();
    time_pnga_nodeid += local_stop - local_start;
    return return_value;
}


void wnga_norm1(Integer g_a, double *nm)
{
    double local_start, local_stop;
    ++count_pnga_norm1;
    local_start = MPI_Wtime();
    pnga_norm1(g_a, nm);
    local_stop = MPI_Wtime();
    time_pnga_norm1 += local_stop - local_start;
}


void wnga_norm_infinity(Integer g_a, double *nm)
{
    double local_start, local_stop;
    ++count_pnga_norm_infinity;
    local_start = MPI_Wtime();
    pnga_norm_infinity(g_a, nm);
    local_stop = MPI_Wtime();
    time_pnga_norm_infinity += local_stop - local_start;
}


void wnga_pack(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi, Integer *icount)
{
    double local_start, local_stop;
    ++count_pnga_pack;
    local_start = MPI_Wtime();
    pnga_pack(g_a, g_b, g_sbit, lo, hi, icount);
    local_stop = MPI_Wtime();
    time_pnga_pack += local_stop - local_start;
}


void wnga_patch_enum(Integer g_a, Integer lo, Integer hi, void *start, void *stride)
{
    double local_start, local_stop;
    ++count_pnga_patch_enum;
    local_start = MPI_Wtime();
    pnga_patch_enum(g_a, lo, hi, start, stride);
    local_stop = MPI_Wtime();
    time_pnga_patch_enum += local_stop - local_start;
}


logical wnga_patch_intersect(Integer *lo, Integer *hi, Integer *lop, Integer *hip, Integer ndim)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_patch_intersect;
    local_start = MPI_Wtime();
    return_value = pnga_patch_intersect(lo, hi, lop, hip, ndim);
    local_stop = MPI_Wtime();
    time_pnga_patch_intersect += local_stop - local_start;
    return return_value;
}


void wnga_periodic(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, void *alpha, Integer op_code)
{
    double local_start, local_stop;
    ++count_pnga_periodic;
    local_start = MPI_Wtime();
    pnga_periodic(g_a, lo, hi, buf, ld, alpha, op_code);
    local_stop = MPI_Wtime();
    time_pnga_periodic += local_stop - local_start;
}


Integer wnga_pgroup_absolute_id(Integer grp, Integer pid)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_absolute_id;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_absolute_id(grp, pid);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_absolute_id += local_stop - local_start;
    return return_value;
}


void wnga_pgroup_brdcst(Integer grp_id, Integer type, void *buf, Integer len, Integer originator)
{
    double local_start, local_stop;
    ++count_pnga_pgroup_brdcst;
    local_start = MPI_Wtime();
    pnga_pgroup_brdcst(grp_id, type, buf, len, originator);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_brdcst += local_stop - local_start;
}


Integer wnga_pgroup_create(Integer *list, Integer count)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_create;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_create(list, count);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_create += local_stop - local_start;
    return return_value;
}


logical wnga_pgroup_destroy(Integer grp)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_destroy;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_destroy(grp);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_destroy += local_stop - local_start;
    return return_value;
}


Integer wnga_pgroup_get_default()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_get_default;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_get_default();
    local_stop = MPI_Wtime();
    time_pnga_pgroup_get_default += local_stop - local_start;
    return return_value;
}


Integer wnga_pgroup_get_mirror()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_get_mirror;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_get_mirror();
    local_stop = MPI_Wtime();
    time_pnga_pgroup_get_mirror += local_stop - local_start;
    return return_value;
}


Integer wnga_pgroup_get_world()
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_get_world;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_get_world();
    local_stop = MPI_Wtime();
    time_pnga_pgroup_get_world += local_stop - local_start;
    return return_value;
}


void wnga_pgroup_gop(Integer p_grp, Integer type, void *x, Integer n, char *op)
{
    double local_start, local_stop;
    ++count_pnga_pgroup_gop;
    local_start = MPI_Wtime();
    pnga_pgroup_gop(p_grp, type, x, n, op);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_gop += local_stop - local_start;
}


Integer wnga_pgroup_nnodes(Integer grp)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_nnodes;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_nnodes(grp);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_nnodes += local_stop - local_start;
    return return_value;
}


Integer wnga_pgroup_nodeid(Integer grp)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_nodeid;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_nodeid(grp);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_nodeid += local_stop - local_start;
    return return_value;
}


void wnga_pgroup_set_default(Integer grp)
{
    double local_start, local_stop;
    ++count_pnga_pgroup_set_default;
    local_start = MPI_Wtime();
    pnga_pgroup_set_default(grp);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_set_default += local_stop - local_start;
}


Integer wnga_pgroup_split(Integer grp, Integer grp_num)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_split;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_split(grp, grp_num);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_split += local_stop - local_start;
    return return_value;
}


Integer wnga_pgroup_split_irreg(Integer grp, Integer mycolor)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_pgroup_split_irreg;
    local_start = MPI_Wtime();
    return_value = pnga_pgroup_split_irreg(grp, mycolor);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_split_irreg += local_stop - local_start;
    return return_value;
}


void wnga_pgroup_sync(Integer grp_id)
{
    double local_start, local_stop;
    ++count_pnga_pgroup_sync;
    local_start = MPI_Wtime();
    pnga_pgroup_sync(grp_id);
    local_stop = MPI_Wtime();
    time_pnga_pgroup_sync += local_stop - local_start;
}


void wnga_print(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_print;
    local_start = MPI_Wtime();
    pnga_print(g_a);
    local_stop = MPI_Wtime();
    time_pnga_print += local_stop - local_start;
}


void wnga_print_distribution(int fstyle, Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_print_distribution;
    local_start = MPI_Wtime();
    pnga_print_distribution(fstyle, g_a);
    local_stop = MPI_Wtime();
    time_pnga_print_distribution += local_stop - local_start;
}


void wnga_print_file(FILE *file, Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_print_file;
    local_start = MPI_Wtime();
    pnga_print_file(file, g_a);
    local_stop = MPI_Wtime();
    time_pnga_print_file += local_stop - local_start;
}


void wnga_print_patch(Integer g_a, Integer *lo, Integer *hi, Integer pretty)
{
    double local_start, local_stop;
    ++count_pnga_print_patch;
    local_start = MPI_Wtime();
    pnga_print_patch(g_a, lo, hi, pretty);
    local_stop = MPI_Wtime();
    time_pnga_print_patch += local_stop - local_start;
}


void wnga_print_patch2d(Integer g_a, Integer ilo, Integer ihi, Integer jlo, Integer jhi, Integer pretty)
{
    double local_start, local_stop;
    ++count_pnga_print_patch2d;
    local_start = MPI_Wtime();
    pnga_print_patch2d(g_a, ilo, ihi, jlo, jhi, pretty);
    local_stop = MPI_Wtime();
    time_pnga_print_patch2d += local_stop - local_start;
}


void wnga_print_patch_file(FILE *file, Integer g_a, Integer *lo, Integer *hi, Integer pretty)
{
    double local_start, local_stop;
    ++count_pnga_print_patch_file;
    local_start = MPI_Wtime();
    pnga_print_patch_file(file, g_a, lo, hi, pretty);
    local_stop = MPI_Wtime();
    time_pnga_print_patch_file += local_stop - local_start;
}


void wnga_print_patch_file2d(FILE *file, Integer g_a, Integer ilo, Integer ihi, Integer jlo, Integer jhi, Integer pretty)
{
    double local_start, local_stop;
    ++count_pnga_print_patch_file2d;
    local_start = MPI_Wtime();
    pnga_print_patch_file2d(file, g_a, ilo, ihi, jlo, jhi, pretty);
    local_stop = MPI_Wtime();
    time_pnga_print_patch_file2d += local_stop - local_start;
}


void wnga_print_stats()
{
    double local_start, local_stop;
    ++count_pnga_print_stats;
    local_start = MPI_Wtime();
    pnga_print_stats();
    local_stop = MPI_Wtime();
    time_pnga_print_stats += local_stop - local_start;
}


void wnga_proc_topology(Integer g_a, Integer proc, Integer *subscript)
{
    double local_start, local_stop;
    ++count_pnga_proc_topology;
    local_start = MPI_Wtime();
    pnga_proc_topology(g_a, proc, subscript);
    local_stop = MPI_Wtime();
    time_pnga_proc_topology += local_stop - local_start;
}


void wnga_put(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_put;
    local_start = MPI_Wtime();
    pnga_put(g_a, lo, hi, buf, ld);
    local_stop = MPI_Wtime();
    time_pnga_put += local_stop - local_start;
}


void wnga_put_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize, void *buf, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_put_field;
    local_start = MPI_Wtime();
    pnga_put_field(g_a, lo, hi, foff, fsize, buf, ld);
    local_stop = MPI_Wtime();
    time_pnga_put_field += local_stop - local_start;
}


void wnga_randomize(Integer g_a, void *val)
{
    double local_start, local_stop;
    ++count_pnga_randomize;
    local_start = MPI_Wtime();
    pnga_randomize(g_a, val);
    local_stop = MPI_Wtime();
    time_pnga_randomize += local_stop - local_start;
}


Integer wnga_read_inc(Integer g_a, Integer *subscript, Integer inc)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_read_inc;
    local_start = MPI_Wtime();
    return_value = pnga_read_inc(g_a, subscript, inc);
    local_stop = MPI_Wtime();
    time_pnga_read_inc += local_stop - local_start;
    return return_value;
}


void wnga_recip(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_recip;
    local_start = MPI_Wtime();
    pnga_recip(g_a);
    local_stop = MPI_Wtime();
    time_pnga_recip += local_stop - local_start;
}


void wnga_recip_patch(Integer g_a, Integer *lo, Integer *hi)
{
    double local_start, local_stop;
    ++count_pnga_recip_patch;
    local_start = MPI_Wtime();
    pnga_recip_patch(g_a, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_recip_patch += local_stop - local_start;
}


int wnga_register_type(size_t size)
{
    int return_value;
    double local_start, local_stop;
    ++count_pnga_register_type;
    local_start = MPI_Wtime();
    return_value = pnga_register_type(size);
    local_stop = MPI_Wtime();
    time_pnga_register_type += local_stop - local_start;
    return return_value;
}


void wnga_release(Integer g_a, Integer *lo, Integer *hi)
{
    double local_start, local_stop;
    ++count_pnga_release;
    local_start = MPI_Wtime();
    pnga_release(g_a, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_release += local_stop - local_start;
}


void wnga_release_block(Integer g_a, Integer iblock)
{
    double local_start, local_stop;
    ++count_pnga_release_block;
    local_start = MPI_Wtime();
    pnga_release_block(g_a, iblock);
    local_stop = MPI_Wtime();
    time_pnga_release_block += local_stop - local_start;
}


void wnga_release_block_grid(Integer g_a, Integer *index)
{
    double local_start, local_stop;
    ++count_pnga_release_block_grid;
    local_start = MPI_Wtime();
    pnga_release_block_grid(g_a, index);
    local_stop = MPI_Wtime();
    time_pnga_release_block_grid += local_stop - local_start;
}


void wnga_release_block_segment(Integer g_a, Integer iproc)
{
    double local_start, local_stop;
    ++count_pnga_release_block_segment;
    local_start = MPI_Wtime();
    pnga_release_block_segment(g_a, iproc);
    local_stop = MPI_Wtime();
    time_pnga_release_block_segment += local_stop - local_start;
}


void wnga_release_ghost_element(Integer g_a, Integer subscript[])
{
    double local_start, local_stop;
    ++count_pnga_release_ghost_element;
    local_start = MPI_Wtime();
    pnga_release_ghost_element(g_a, subscript);
    local_stop = MPI_Wtime();
    time_pnga_release_ghost_element += local_stop - local_start;
}


void wnga_release_ghosts(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_release_ghosts;
    local_start = MPI_Wtime();
    pnga_release_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_release_ghosts += local_stop - local_start;
}


void wnga_release_update(Integer g_a, Integer *lo, Integer *hi)
{
    double local_start, local_stop;
    ++count_pnga_release_update;
    local_start = MPI_Wtime();
    pnga_release_update(g_a, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_release_update += local_stop - local_start;
}


void wnga_release_update_block(Integer g_a, Integer iblock)
{
    double local_start, local_stop;
    ++count_pnga_release_update_block;
    local_start = MPI_Wtime();
    pnga_release_update_block(g_a, iblock);
    local_stop = MPI_Wtime();
    time_pnga_release_update_block += local_stop - local_start;
}


void wnga_release_update_block_grid(Integer g_a, Integer *index)
{
    double local_start, local_stop;
    ++count_pnga_release_update_block_grid;
    local_start = MPI_Wtime();
    pnga_release_update_block_grid(g_a, index);
    local_stop = MPI_Wtime();
    time_pnga_release_update_block_grid += local_stop - local_start;
}


void wnga_release_update_block_segment(Integer g_a, Integer iproc)
{
    double local_start, local_stop;
    ++count_pnga_release_update_block_segment;
    local_start = MPI_Wtime();
    pnga_release_update_block_segment(g_a, iproc);
    local_stop = MPI_Wtime();
    time_pnga_release_update_block_segment += local_stop - local_start;
}


void wnga_release_update_ghost_element(Integer g_a, Integer subscript[])
{
    double local_start, local_stop;
    ++count_pnga_release_update_ghost_element;
    local_start = MPI_Wtime();
    pnga_release_update_ghost_element(g_a, subscript);
    local_stop = MPI_Wtime();
    time_pnga_release_update_ghost_element += local_stop - local_start;
}


void wnga_release_update_ghosts(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_release_update_ghosts;
    local_start = MPI_Wtime();
    pnga_release_update_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_release_update_ghosts += local_stop - local_start;
}


void wnga_scale(Integer g_a, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_scale;
    local_start = MPI_Wtime();
    pnga_scale(g_a, alpha);
    local_stop = MPI_Wtime();
    time_pnga_scale += local_stop - local_start;
}


void wnga_scale_cols(Integer g_a, Integer g_v)
{
    double local_start, local_stop;
    ++count_pnga_scale_cols;
    local_start = MPI_Wtime();
    pnga_scale_cols(g_a, g_v);
    local_stop = MPI_Wtime();
    time_pnga_scale_cols += local_stop - local_start;
}


void wnga_scale_patch(Integer g_a, Integer *lo, Integer *hi, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_scale_patch;
    local_start = MPI_Wtime();
    pnga_scale_patch(g_a, lo, hi, alpha);
    local_stop = MPI_Wtime();
    time_pnga_scale_patch += local_stop - local_start;
}


void wnga_scale_rows(Integer g_a, Integer g_v)
{
    double local_start, local_stop;
    ++count_pnga_scale_rows;
    local_start = MPI_Wtime();
    pnga_scale_rows(g_a, g_v);
    local_stop = MPI_Wtime();
    time_pnga_scale_rows += local_stop - local_start;
}


void wnga_scan_add(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi, Integer excl)
{
    double local_start, local_stop;
    ++count_pnga_scan_add;
    local_start = MPI_Wtime();
    pnga_scan_add(g_a, g_b, g_sbit, lo, hi, excl);
    local_stop = MPI_Wtime();
    time_pnga_scan_add += local_stop - local_start;
}


void wnga_scan_copy(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi)
{
    double local_start, local_stop;
    ++count_pnga_scan_copy;
    local_start = MPI_Wtime();
    pnga_scan_copy(g_a, g_b, g_sbit, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_scan_copy += local_stop - local_start;
}


void wnga_scatter(Integer g_a, void *v, Integer *subscript, Integer nv)
{
    double local_start, local_stop;
    ++count_pnga_scatter;
    local_start = MPI_Wtime();
    pnga_scatter(g_a, v, subscript, nv);
    local_stop = MPI_Wtime();
    time_pnga_scatter += local_stop - local_start;
}


void wnga_scatter2d(Integer g_a, void *v, Integer *i, Integer *j, Integer nv)
{
    double local_start, local_stop;
    ++count_pnga_scatter2d;
    local_start = MPI_Wtime();
    pnga_scatter2d(g_a, v, i, j, nv);
    local_stop = MPI_Wtime();
    time_pnga_scatter2d += local_stop - local_start;
}


void wnga_scatter_acc(Integer g_a, void *v, Integer subscript[], Integer nv, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_scatter_acc;
    local_start = MPI_Wtime();
    pnga_scatter_acc(g_a, v, subscript, nv, alpha);
    local_stop = MPI_Wtime();
    time_pnga_scatter_acc += local_stop - local_start;
}


void wnga_scatter_acc2d(Integer g_a, void *v, Integer *i, Integer *j, Integer nv, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_scatter_acc2d;
    local_start = MPI_Wtime();
    pnga_scatter_acc2d(g_a, v, i, j, nv, alpha);
    local_stop = MPI_Wtime();
    time_pnga_scatter_acc2d += local_stop - local_start;
}


void wnga_select_elem(Integer g_a, char *op, void *val, Integer *subscript)
{
    double local_start, local_stop;
    ++count_pnga_select_elem;
    local_start = MPI_Wtime();
    pnga_select_elem(g_a, op, val, subscript);
    local_stop = MPI_Wtime();
    time_pnga_select_elem += local_stop - local_start;
}


void wnga_set_array_name(Integer g_a, char *array_name)
{
    double local_start, local_stop;
    ++count_pnga_set_array_name;
    local_start = MPI_Wtime();
    pnga_set_array_name(g_a, array_name);
    local_stop = MPI_Wtime();
    time_pnga_set_array_name += local_stop - local_start;
}


void wnga_set_block_cyclic(Integer g_a, Integer *dims)
{
    double local_start, local_stop;
    ++count_pnga_set_block_cyclic;
    local_start = MPI_Wtime();
    pnga_set_block_cyclic(g_a, dims);
    local_stop = MPI_Wtime();
    time_pnga_set_block_cyclic += local_stop - local_start;
}


void wnga_set_block_cyclic_proc_grid(Integer g_a, Integer *dims, Integer *proc_grid)
{
    double local_start, local_stop;
    ++count_pnga_set_block_cyclic_proc_grid;
    local_start = MPI_Wtime();
    pnga_set_block_cyclic_proc_grid(g_a, dims, proc_grid);
    local_stop = MPI_Wtime();
    time_pnga_set_block_cyclic_proc_grid += local_stop - local_start;
}


void wnga_set_chunk(Integer g_a, Integer *chunk)
{
    double local_start, local_stop;
    ++count_pnga_set_chunk;
    local_start = MPI_Wtime();
    pnga_set_chunk(g_a, chunk);
    local_stop = MPI_Wtime();
    time_pnga_set_chunk += local_stop - local_start;
}


void wnga_set_data(Integer g_a, Integer ndim, Integer *dims, Integer type)
{
    double local_start, local_stop;
    ++count_pnga_set_data;
    local_start = MPI_Wtime();
    pnga_set_data(g_a, ndim, dims, type);
    local_stop = MPI_Wtime();
    time_pnga_set_data += local_stop - local_start;
}


void wnga_set_debug(logical flag)
{
    double local_start, local_stop;
    ++count_pnga_set_debug;
    local_start = MPI_Wtime();
    pnga_set_debug(flag);
    local_stop = MPI_Wtime();
    time_pnga_set_debug += local_stop - local_start;
}


void wnga_set_diagonal(Integer g_a, Integer g_v)
{
    double local_start, local_stop;
    ++count_pnga_set_diagonal;
    local_start = MPI_Wtime();
    pnga_set_diagonal(g_a, g_v);
    local_stop = MPI_Wtime();
    time_pnga_set_diagonal += local_stop - local_start;
}


void wnga_set_ghost_corner_flag(Integer g_a, logical flag)
{
    double local_start, local_stop;
    ++count_pnga_set_ghost_corner_flag;
    local_start = MPI_Wtime();
    pnga_set_ghost_corner_flag(g_a, flag);
    local_stop = MPI_Wtime();
    time_pnga_set_ghost_corner_flag += local_stop - local_start;
}


logical wnga_set_ghost_info(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_set_ghost_info;
    local_start = MPI_Wtime();
    return_value = pnga_set_ghost_info(g_a);
    local_stop = MPI_Wtime();
    time_pnga_set_ghost_info += local_stop - local_start;
    return return_value;
}


void wnga_set_ghosts(Integer g_a, Integer *width)
{
    double local_start, local_stop;
    ++count_pnga_set_ghosts;
    local_start = MPI_Wtime();
    pnga_set_ghosts(g_a, width);
    local_stop = MPI_Wtime();
    time_pnga_set_ghosts += local_stop - local_start;
}


void wnga_set_irreg_distr(Integer g_a, Integer *map, Integer *block)
{
    double local_start, local_stop;
    ++count_pnga_set_irreg_distr;
    local_start = MPI_Wtime();
    pnga_set_irreg_distr(g_a, map, block);
    local_stop = MPI_Wtime();
    time_pnga_set_irreg_distr += local_stop - local_start;
}


void wnga_set_irreg_flag(Integer g_a, logical flag)
{
    double local_start, local_stop;
    ++count_pnga_set_irreg_flag;
    local_start = MPI_Wtime();
    pnga_set_irreg_flag(g_a, flag);
    local_stop = MPI_Wtime();
    time_pnga_set_irreg_flag += local_stop - local_start;
}


void wnga_set_memory_limit(Integer mem_limit)
{
    double local_start, local_stop;
    ++count_pnga_set_memory_limit;
    local_start = MPI_Wtime();
    pnga_set_memory_limit(mem_limit);
    local_stop = MPI_Wtime();
    time_pnga_set_memory_limit += local_stop - local_start;
}


void wnga_set_pgroup(Integer g_a, Integer p_handle)
{
    double local_start, local_stop;
    ++count_pnga_set_pgroup;
    local_start = MPI_Wtime();
    pnga_set_pgroup(g_a, p_handle);
    local_stop = MPI_Wtime();
    time_pnga_set_pgroup += local_stop - local_start;
}


void wnga_set_restricted(Integer g_a, Integer *list, Integer size)
{
    double local_start, local_stop;
    ++count_pnga_set_restricted;
    local_start = MPI_Wtime();
    pnga_set_restricted(g_a, list, size);
    local_stop = MPI_Wtime();
    time_pnga_set_restricted += local_stop - local_start;
}


void wnga_set_restricted_range(Integer g_a, Integer lo_proc, Integer hi_proc)
{
    double local_start, local_stop;
    ++count_pnga_set_restricted_range;
    local_start = MPI_Wtime();
    pnga_set_restricted_range(g_a, lo_proc, hi_proc);
    local_stop = MPI_Wtime();
    time_pnga_set_restricted_range += local_stop - local_start;
}


logical wnga_set_update4_info(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_set_update4_info;
    local_start = MPI_Wtime();
    return_value = pnga_set_update4_info(g_a);
    local_stop = MPI_Wtime();
    time_pnga_set_update4_info += local_stop - local_start;
    return return_value;
}


logical wnga_set_update5_info(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_set_update5_info;
    local_start = MPI_Wtime();
    return_value = pnga_set_update5_info(g_a);
    local_stop = MPI_Wtime();
    time_pnga_set_update5_info += local_stop - local_start;
    return return_value;
}


void wnga_shift_diagonal(Integer g_a, void *c)
{
    double local_start, local_stop;
    ++count_pnga_shift_diagonal;
    local_start = MPI_Wtime();
    pnga_shift_diagonal(g_a, c);
    local_stop = MPI_Wtime();
    time_pnga_shift_diagonal += local_stop - local_start;
}


Integer wnga_solve(Integer g_a, Integer g_b)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_solve;
    local_start = MPI_Wtime();
    return_value = pnga_solve(g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_solve += local_stop - local_start;
    return return_value;
}


Integer wnga_spd_invert(Integer g_a)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_spd_invert;
    local_start = MPI_Wtime();
    return_value = pnga_spd_invert(g_a);
    local_stop = MPI_Wtime();
    time_pnga_spd_invert += local_stop - local_start;
    return return_value;
}


void wnga_step_bound_info(Integer g_xx, Integer g_vv, Integer g_xxll, Integer g_xxuu, void *boundmin, void *wolfemin, void *boundmax)
{
    double local_start, local_stop;
    ++count_pnga_step_bound_info;
    local_start = MPI_Wtime();
    pnga_step_bound_info(g_xx, g_vv, g_xxll, g_xxuu, boundmin, wolfemin, boundmax);
    local_stop = MPI_Wtime();
    time_pnga_step_bound_info += local_stop - local_start;
}


void wnga_step_bound_info_patch(Integer g_xx, Integer *xxlo, Integer *xxhi, Integer g_vv, Integer *vvlo, Integer *vvhi, Integer g_xxll, Integer *xxlllo, Integer *xxllhi, Integer g_xxuu, Integer *xxuulo, Integer *xxuuhi, void *boundmin, void *wolfemin, void *boundmax)
{
    double local_start, local_stop;
    ++count_pnga_step_bound_info_patch;
    local_start = MPI_Wtime();
    pnga_step_bound_info_patch(g_xx, xxlo, xxhi, g_vv, vvlo, vvhi, g_xxll, xxlllo, xxllhi, g_xxuu, xxuulo, xxuuhi, boundmin, wolfemin, boundmax);
    local_stop = MPI_Wtime();
    time_pnga_step_bound_info_patch += local_stop - local_start;
}


void wnga_step_mask_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, Integer g_c, Integer *clo, Integer *chi)
{
    double local_start, local_stop;
    ++count_pnga_step_mask_patch;
    local_start = MPI_Wtime();
    pnga_step_mask_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi);
    local_stop = MPI_Wtime();
    time_pnga_step_mask_patch += local_stop - local_start;
}


void wnga_step_max(Integer g_a, Integer g_b, void *retval)
{
    double local_start, local_stop;
    ++count_pnga_step_max;
    local_start = MPI_Wtime();
    pnga_step_max(g_a, g_b, retval);
    local_stop = MPI_Wtime();
    time_pnga_step_max += local_stop - local_start;
}


void wnga_step_max_patch(Integer g_a, Integer *alo, Integer *ahi, Integer g_b, Integer *blo, Integer *bhi, void *result)
{
    double local_start, local_stop;
    ++count_pnga_step_max_patch;
    local_start = MPI_Wtime();
    pnga_step_max_patch(g_a, alo, ahi, g_b, blo, bhi, result);
    local_stop = MPI_Wtime();
    time_pnga_step_max_patch += local_stop - local_start;
}


void wnga_strided_acc(Integer g_a, Integer *lo, Integer *hi, Integer *skip, void *buf, Integer *ld, void *alpha)
{
    double local_start, local_stop;
    ++count_pnga_strided_acc;
    local_start = MPI_Wtime();
    pnga_strided_acc(g_a, lo, hi, skip, buf, ld, alpha);
    local_stop = MPI_Wtime();
    time_pnga_strided_acc += local_stop - local_start;
}


void wnga_strided_get(Integer g_a, Integer *lo, Integer *hi, Integer *skip, void *buf, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_strided_get;
    local_start = MPI_Wtime();
    pnga_strided_get(g_a, lo, hi, skip, buf, ld);
    local_stop = MPI_Wtime();
    time_pnga_strided_get += local_stop - local_start;
}


void wnga_strided_put(Integer g_a, Integer *lo, Integer *hi, Integer *skip, void *buf, Integer *ld)
{
    double local_start, local_stop;
    ++count_pnga_strided_put;
    local_start = MPI_Wtime();
    pnga_strided_put(g_a, lo, hi, skip, buf, ld);
    local_stop = MPI_Wtime();
    time_pnga_strided_put += local_stop - local_start;
}


void wnga_summarize(Integer verbose)
{
    double local_start, local_stop;
    ++count_pnga_summarize;
    local_start = MPI_Wtime();
    pnga_summarize(verbose);
    local_stop = MPI_Wtime();
    time_pnga_summarize += local_stop - local_start;
}


void wnga_symmetrize(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_symmetrize;
    local_start = MPI_Wtime();
    pnga_symmetrize(g_a);
    local_stop = MPI_Wtime();
    time_pnga_symmetrize += local_stop - local_start;
}


void wnga_sync()
{
    double local_start, local_stop;
    ++count_pnga_sync;
    local_start = MPI_Wtime();
    pnga_sync();
    local_stop = MPI_Wtime();
    time_pnga_sync += local_stop - local_start;
}


double wnga_timer()
{
    double return_value;
    double local_start, local_stop;
    ++count_pnga_timer;
    local_start = MPI_Wtime();
    return_value = pnga_timer();
    local_stop = MPI_Wtime();
    time_pnga_timer += local_stop - local_start;
    return return_value;
}


Integer wnga_total_blocks(Integer g_a)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_total_blocks;
    local_start = MPI_Wtime();
    return_value = pnga_total_blocks(g_a);
    local_stop = MPI_Wtime();
    time_pnga_total_blocks += local_stop - local_start;
    return return_value;
}


void wnga_transpose(Integer g_a, Integer g_b)
{
    double local_start, local_stop;
    ++count_pnga_transpose;
    local_start = MPI_Wtime();
    pnga_transpose(g_a, g_b);
    local_stop = MPI_Wtime();
    time_pnga_transpose += local_stop - local_start;
}


Integer wnga_type_c2f(Integer type)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_type_c2f;
    local_start = MPI_Wtime();
    return_value = pnga_type_c2f(type);
    local_stop = MPI_Wtime();
    time_pnga_type_c2f += local_stop - local_start;
    return return_value;
}


Integer wnga_type_f2c(Integer type)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_type_f2c;
    local_start = MPI_Wtime();
    return_value = pnga_type_f2c(type);
    local_stop = MPI_Wtime();
    time_pnga_type_f2c += local_stop - local_start;
    return return_value;
}


void wnga_unlock(Integer mutex)
{
    double local_start, local_stop;
    ++count_pnga_unlock;
    local_start = MPI_Wtime();
    pnga_unlock(mutex);
    local_stop = MPI_Wtime();
    time_pnga_unlock += local_stop - local_start;
}


void wnga_unpack(Integer g_a, Integer g_b, Integer g_sbit, Integer lo, Integer hi, Integer *icount)
{
    double local_start, local_stop;
    ++count_pnga_unpack;
    local_start = MPI_Wtime();
    pnga_unpack(g_a, g_b, g_sbit, lo, hi, icount);
    local_stop = MPI_Wtime();
    time_pnga_unpack += local_stop - local_start;
}


void wnga_update1_ghosts(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_update1_ghosts;
    local_start = MPI_Wtime();
    pnga_update1_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update1_ghosts += local_stop - local_start;
}


logical wnga_update2_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update2_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update2_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update2_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update3_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update3_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update3_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update3_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update44_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update44_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update44_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update44_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update4_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update4_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update4_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update4_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update55_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update55_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update55_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update55_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update5_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update5_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update5_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update5_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update6_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update6_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update6_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update6_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update7_ghosts(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update7_ghosts;
    local_start = MPI_Wtime();
    return_value = pnga_update7_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update7_ghosts += local_stop - local_start;
    return return_value;
}


logical wnga_update_ghost_dir(Integer g_a, Integer pdim, Integer pdir, logical pflag)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_update_ghost_dir;
    local_start = MPI_Wtime();
    return_value = pnga_update_ghost_dir(g_a, pdim, pdir, pflag);
    local_stop = MPI_Wtime();
    time_pnga_update_ghost_dir += local_stop - local_start;
    return return_value;
}


void wnga_update_ghosts(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_update_ghosts;
    local_start = MPI_Wtime();
    pnga_update_ghosts(g_a);
    local_stop = MPI_Wtime();
    time_pnga_update_ghosts += local_stop - local_start;
}


logical wnga_uses_ma()
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_uses_ma;
    local_start = MPI_Wtime();
    return_value = pnga_uses_ma();
    local_stop = MPI_Wtime();
    time_pnga_uses_ma += local_stop - local_start;
    return return_value;
}


logical wnga_uses_proc_grid(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_uses_proc_grid;
    local_start = MPI_Wtime();
    return_value = pnga_uses_proc_grid(g_a);
    local_stop = MPI_Wtime();
    time_pnga_uses_proc_grid += local_stop - local_start;
    return return_value;
}


logical wnga_valid_handle(Integer g_a)
{
    logical return_value;
    double local_start, local_stop;
    ++count_pnga_valid_handle;
    local_start = MPI_Wtime();
    return_value = pnga_valid_handle(g_a);
    local_stop = MPI_Wtime();
    time_pnga_valid_handle += local_stop - local_start;
    return return_value;
}


Integer wnga_verify_handle(Integer g_a)
{
    Integer return_value;
    double local_start, local_stop;
    ++count_pnga_verify_handle;
    local_start = MPI_Wtime();
    return_value = pnga_verify_handle(g_a);
    local_stop = MPI_Wtime();
    time_pnga_verify_handle += local_stop - local_start;
    return return_value;
}


DoublePrecision wnga_wtime()
{
    DoublePrecision return_value;
    double local_start, local_stop;
    ++count_pnga_wtime;
    local_start = MPI_Wtime();
    return_value = pnga_wtime();
    local_stop = MPI_Wtime();
    time_pnga_wtime += local_stop - local_start;
    return return_value;
}


void wnga_zero(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_zero;
    local_start = MPI_Wtime();
    pnga_zero(g_a);
    local_stop = MPI_Wtime();
    time_pnga_zero += local_stop - local_start;
}


void wnga_zero_diagonal(Integer g_a)
{
    double local_start, local_stop;
    ++count_pnga_zero_diagonal;
    local_start = MPI_Wtime();
    pnga_zero_diagonal(g_a);
    local_stop = MPI_Wtime();
    time_pnga_zero_diagonal += local_stop - local_start;
}


void wnga_zero_patch(Integer g_a, Integer *lo, Integer *hi)
{
    double local_start, local_stop;
    ++count_pnga_zero_patch;
    local_start = MPI_Wtime();
    pnga_zero_patch(g_a, lo, hi);
    local_stop = MPI_Wtime();
    time_pnga_zero_patch += local_stop - local_start;
}

void wnga_initialize()
{
    ++count_pnga_initialize;
    pnga_initialize();
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
}

void wnga_initialize_ltd(Integer limit)
{
    ++count_pnga_initialize_ltd;
    pnga_initialize_ltd(limit);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
}

void wnga_terminate()
{
    ++count_pnga_terminate;
    pnga_terminate();
    /* don't dump info if terminate more than once */
    if (1 == count_pnga_terminate) {

        double recvbuf = 0.0;

        MPI_Reduce(&time_pnga_abs_value, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_abs_value,%ld,%lf\n", count_pnga_abs_value, recvbuf);
        }

        MPI_Reduce(&time_pnga_abs_value_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_abs_value_patch,%ld,%lf\n", count_pnga_abs_value_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_acc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_acc,%ld,%lf\n", count_pnga_acc, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_block_grid_idx, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_block_grid_idx,%ld,%lf\n", count_pnga_access_block_grid_idx, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_block_grid_ptr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_block_grid_ptr,%ld,%lf\n", count_pnga_access_block_grid_ptr, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_block_idx, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_block_idx,%ld,%lf\n", count_pnga_access_block_idx, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_block_ptr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_block_ptr,%ld,%lf\n", count_pnga_access_block_ptr, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_block_segment_idx, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_block_segment_idx,%ld,%lf\n", count_pnga_access_block_segment_idx, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_block_segment_ptr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_block_segment_ptr,%ld,%lf\n", count_pnga_access_block_segment_ptr, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_ghost_element, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_ghost_element,%ld,%lf\n", count_pnga_access_ghost_element, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_ghost_element_ptr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_ghost_element_ptr,%ld,%lf\n", count_pnga_access_ghost_element_ptr, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_ghost_ptr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_ghost_ptr,%ld,%lf\n", count_pnga_access_ghost_ptr, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_ghosts,%ld,%lf\n", count_pnga_access_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_idx, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_idx,%ld,%lf\n", count_pnga_access_idx, recvbuf);
        }

        MPI_Reduce(&time_pnga_access_ptr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_access_ptr,%ld,%lf\n", count_pnga_access_ptr, recvbuf);
        }

        MPI_Reduce(&time_pnga_add, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_add,%ld,%lf\n", count_pnga_add, recvbuf);
        }

        MPI_Reduce(&time_pnga_add_constant, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_add_constant,%ld,%lf\n", count_pnga_add_constant, recvbuf);
        }

        MPI_Reduce(&time_pnga_add_constant_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_add_constant_patch,%ld,%lf\n", count_pnga_add_constant_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_add_diagonal, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_add_diagonal,%ld,%lf\n", count_pnga_add_diagonal, recvbuf);
        }

        MPI_Reduce(&time_pnga_add_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_add_patch,%ld,%lf\n", count_pnga_add_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_allocate, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_allocate,%ld,%lf\n", count_pnga_allocate, recvbuf);
        }

        MPI_Reduce(&time_pnga_bin_index, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_bin_index,%ld,%lf\n", count_pnga_bin_index, recvbuf);
        }

        MPI_Reduce(&time_pnga_bin_sorter, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_bin_sorter,%ld,%lf\n", count_pnga_bin_sorter, recvbuf);
        }

        MPI_Reduce(&time_pnga_brdcst, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_brdcst,%ld,%lf\n", count_pnga_brdcst, recvbuf);
        }

        MPI_Reduce(&time_pnga_check_handle, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_check_handle,%ld,%lf\n", count_pnga_check_handle, recvbuf);
        }

        MPI_Reduce(&time_pnga_cluster_nnodes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_cluster_nnodes,%ld,%lf\n", count_pnga_cluster_nnodes, recvbuf);
        }

        MPI_Reduce(&time_pnga_cluster_nodeid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_cluster_nodeid,%ld,%lf\n", count_pnga_cluster_nodeid, recvbuf);
        }

        MPI_Reduce(&time_pnga_cluster_nprocs, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_cluster_nprocs,%ld,%lf\n", count_pnga_cluster_nprocs, recvbuf);
        }

        MPI_Reduce(&time_pnga_cluster_proc_nodeid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_cluster_proc_nodeid,%ld,%lf\n", count_pnga_cluster_proc_nodeid, recvbuf);
        }

        MPI_Reduce(&time_pnga_cluster_procid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_cluster_procid,%ld,%lf\n", count_pnga_cluster_procid, recvbuf);
        }

        MPI_Reduce(&time_pnga_comp_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_comp_patch,%ld,%lf\n", count_pnga_comp_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_compare_distr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_compare_distr,%ld,%lf\n", count_pnga_compare_distr, recvbuf);
        }

        MPI_Reduce(&time_pnga_copy, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_copy,%ld,%lf\n", count_pnga_copy, recvbuf);
        }

        MPI_Reduce(&time_pnga_copy_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_copy_patch,%ld,%lf\n", count_pnga_copy_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_copy_patch_dp, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_copy_patch_dp,%ld,%lf\n", count_pnga_copy_patch_dp, recvbuf);
        }

        MPI_Reduce(&time_pnga_create, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create,%ld,%lf\n", count_pnga_create, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_bin_range, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_bin_range,%ld,%lf\n", count_pnga_create_bin_range, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_config, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_config,%ld,%lf\n", count_pnga_create_config, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_ghosts,%ld,%lf\n", count_pnga_create_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_ghosts_config, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_ghosts_config,%ld,%lf\n", count_pnga_create_ghosts_config, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_ghosts_irreg, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_ghosts_irreg,%ld,%lf\n", count_pnga_create_ghosts_irreg, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_ghosts_irreg_config, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_ghosts_irreg_config,%ld,%lf\n", count_pnga_create_ghosts_irreg_config, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_handle, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_handle,%ld,%lf\n", count_pnga_create_handle, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_irreg, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_irreg,%ld,%lf\n", count_pnga_create_irreg, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_irreg_config, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_irreg_config,%ld,%lf\n", count_pnga_create_irreg_config, recvbuf);
        }

        MPI_Reduce(&time_pnga_create_mutexes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_create_mutexes,%ld,%lf\n", count_pnga_create_mutexes, recvbuf);
        }

        MPI_Reduce(&time_pnga_ddot_patch_dp, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_ddot_patch_dp,%ld,%lf\n", count_pnga_ddot_patch_dp, recvbuf);
        }

        MPI_Reduce(&time_pnga_deregister_type, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_deregister_type,%ld,%lf\n", count_pnga_deregister_type, recvbuf);
        }

        MPI_Reduce(&time_pnga_destroy, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_destroy,%ld,%lf\n", count_pnga_destroy, recvbuf);
        }

        MPI_Reduce(&time_pnga_destroy_mutexes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_destroy_mutexes,%ld,%lf\n", count_pnga_destroy_mutexes, recvbuf);
        }

        MPI_Reduce(&time_pnga_diag, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_diag,%ld,%lf\n", count_pnga_diag, recvbuf);
        }

        MPI_Reduce(&time_pnga_diag_reuse, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_diag_reuse,%ld,%lf\n", count_pnga_diag_reuse, recvbuf);
        }

        MPI_Reduce(&time_pnga_diag_seq, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_diag_seq,%ld,%lf\n", count_pnga_diag_seq, recvbuf);
        }

        MPI_Reduce(&time_pnga_diag_std, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_diag_std,%ld,%lf\n", count_pnga_diag_std, recvbuf);
        }

        MPI_Reduce(&time_pnga_diag_std_seq, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_diag_std_seq,%ld,%lf\n", count_pnga_diag_std_seq, recvbuf);
        }

        MPI_Reduce(&time_pnga_distribution, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_distribution,%ld,%lf\n", count_pnga_distribution, recvbuf);
        }

        MPI_Reduce(&time_pnga_dot, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_dot,%ld,%lf\n", count_pnga_dot, recvbuf);
        }

        MPI_Reduce(&time_pnga_dot_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_dot_patch,%ld,%lf\n", count_pnga_dot_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_duplicate, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_duplicate,%ld,%lf\n", count_pnga_duplicate, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_divide, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_divide,%ld,%lf\n", count_pnga_elem_divide, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_divide_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_divide_patch,%ld,%lf\n", count_pnga_elem_divide_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_maximum, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_maximum,%ld,%lf\n", count_pnga_elem_maximum, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_maximum_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_maximum_patch,%ld,%lf\n", count_pnga_elem_maximum_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_minimum, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_minimum,%ld,%lf\n", count_pnga_elem_minimum, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_minimum_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_minimum_patch,%ld,%lf\n", count_pnga_elem_minimum_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_multiply, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_multiply,%ld,%lf\n", count_pnga_elem_multiply, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_multiply_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_multiply_patch,%ld,%lf\n", count_pnga_elem_multiply_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_step_divide_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_step_divide_patch,%ld,%lf\n", count_pnga_elem_step_divide_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_elem_stepb_divide_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_elem_stepb_divide_patch,%ld,%lf\n", count_pnga_elem_stepb_divide_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_error, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_error,%ld,%lf\n", count_pnga_error, recvbuf);
        }

        MPI_Reduce(&time_pnga_fence, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_fence,%ld,%lf\n", count_pnga_fence, recvbuf);
        }

        MPI_Reduce(&time_pnga_fill, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_fill,%ld,%lf\n", count_pnga_fill, recvbuf);
        }

        MPI_Reduce(&time_pnga_fill_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_fill_patch,%ld,%lf\n", count_pnga_fill_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_gather, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_gather,%ld,%lf\n", count_pnga_gather, recvbuf);
        }

        MPI_Reduce(&time_pnga_gather2d, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_gather2d,%ld,%lf\n", count_pnga_gather2d, recvbuf);
        }

        MPI_Reduce(&time_pnga_get, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get,%ld,%lf\n", count_pnga_get, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_block_info, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_block_info,%ld,%lf\n", count_pnga_get_block_info, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_debug, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_debug,%ld,%lf\n", count_pnga_get_debug, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_diag, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_diag,%ld,%lf\n", count_pnga_get_diag, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_dimension, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_dimension,%ld,%lf\n", count_pnga_get_dimension, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_field, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_field,%ld,%lf\n", count_pnga_get_field, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_ghost_block, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_ghost_block,%ld,%lf\n", count_pnga_get_ghost_block, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_pgroup, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_pgroup,%ld,%lf\n", count_pnga_get_pgroup, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_pgroup_size, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_pgroup_size,%ld,%lf\n", count_pnga_get_pgroup_size, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_proc_grid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_proc_grid,%ld,%lf\n", count_pnga_get_proc_grid, recvbuf);
        }

        MPI_Reduce(&time_pnga_get_proc_index, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_get_proc_index,%ld,%lf\n", count_pnga_get_proc_index, recvbuf);
        }

        MPI_Reduce(&time_pnga_ghost_barrier, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_ghost_barrier,%ld,%lf\n", count_pnga_ghost_barrier, recvbuf);
        }

        MPI_Reduce(&time_pnga_gop, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_gop,%ld,%lf\n", count_pnga_gop, recvbuf);
        }

        MPI_Reduce(&time_pnga_has_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_has_ghosts,%ld,%lf\n", count_pnga_has_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_init_fence, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_init_fence,%ld,%lf\n", count_pnga_init_fence, recvbuf);
        }

        MPI_Reduce(&time_pnga_initialize, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_initialize,%ld,%lf\n", count_pnga_initialize, recvbuf);
        }

        MPI_Reduce(&time_pnga_initialize_ltd, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_initialize_ltd,%ld,%lf\n", count_pnga_initialize_ltd, recvbuf);
        }

        MPI_Reduce(&time_pnga_inquire, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_inquire,%ld,%lf\n", count_pnga_inquire, recvbuf);
        }

        MPI_Reduce(&time_pnga_inquire_memory, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_inquire_memory,%ld,%lf\n", count_pnga_inquire_memory, recvbuf);
        }

        MPI_Reduce(&time_pnga_inquire_name, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_inquire_name,%ld,%lf\n", count_pnga_inquire_name, recvbuf);
        }

        MPI_Reduce(&time_pnga_inquire_type, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_inquire_type,%ld,%lf\n", count_pnga_inquire_type, recvbuf);
        }

        MPI_Reduce(&time_pnga_is_mirrored, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_is_mirrored,%ld,%lf\n", count_pnga_is_mirrored, recvbuf);
        }

        MPI_Reduce(&time_pnga_list_nodeid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_list_nodeid,%ld,%lf\n", count_pnga_list_nodeid, recvbuf);
        }

        MPI_Reduce(&time_pnga_llt_solve, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_llt_solve,%ld,%lf\n", count_pnga_llt_solve, recvbuf);
        }

        MPI_Reduce(&time_pnga_locate, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_locate,%ld,%lf\n", count_pnga_locate, recvbuf);
        }

        MPI_Reduce(&time_pnga_locate_nnodes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_locate_nnodes,%ld,%lf\n", count_pnga_locate_nnodes, recvbuf);
        }

        MPI_Reduce(&time_pnga_locate_num_blocks, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_locate_num_blocks,%ld,%lf\n", count_pnga_locate_num_blocks, recvbuf);
        }

        MPI_Reduce(&time_pnga_locate_region, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_locate_region,%ld,%lf\n", count_pnga_locate_region, recvbuf);
        }

        MPI_Reduce(&time_pnga_lock, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_lock,%ld,%lf\n", count_pnga_lock, recvbuf);
        }

        MPI_Reduce(&time_pnga_lu_solve, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_lu_solve,%ld,%lf\n", count_pnga_lu_solve, recvbuf);
        }

        MPI_Reduce(&time_pnga_lu_solve_alt, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_lu_solve_alt,%ld,%lf\n", count_pnga_lu_solve_alt, recvbuf);
        }

        MPI_Reduce(&time_pnga_lu_solve_seq, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_lu_solve_seq,%ld,%lf\n", count_pnga_lu_solve_seq, recvbuf);
        }

        MPI_Reduce(&time_pnga_mask_sync, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_mask_sync,%ld,%lf\n", count_pnga_mask_sync, recvbuf);
        }

        MPI_Reduce(&time_pnga_matmul, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_matmul,%ld,%lf\n", count_pnga_matmul, recvbuf);
        }

        MPI_Reduce(&time_pnga_matmul_mirrored, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_matmul_mirrored,%ld,%lf\n", count_pnga_matmul_mirrored, recvbuf);
        }

        MPI_Reduce(&time_pnga_matmul_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_matmul_patch,%ld,%lf\n", count_pnga_matmul_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_median, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_median,%ld,%lf\n", count_pnga_median, recvbuf);
        }

        MPI_Reduce(&time_pnga_median_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_median_patch,%ld,%lf\n", count_pnga_median_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_memory_avail, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_memory_avail,%ld,%lf\n", count_pnga_memory_avail, recvbuf);
        }

        MPI_Reduce(&time_pnga_memory_avail_type, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_memory_avail_type,%ld,%lf\n", count_pnga_memory_avail_type, recvbuf);
        }

        MPI_Reduce(&time_pnga_memory_limited, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_memory_limited,%ld,%lf\n", count_pnga_memory_limited, recvbuf);
        }

        MPI_Reduce(&time_pnga_merge_distr_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_merge_distr_patch,%ld,%lf\n", count_pnga_merge_distr_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_merge_mirrored, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_merge_mirrored,%ld,%lf\n", count_pnga_merge_mirrored, recvbuf);
        }

        MPI_Reduce(&time_pnga_msg_brdcst, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_msg_brdcst,%ld,%lf\n", count_pnga_msg_brdcst, recvbuf);
        }

        MPI_Reduce(&time_pnga_msg_pgroup_sync, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_msg_pgroup_sync,%ld,%lf\n", count_pnga_msg_pgroup_sync, recvbuf);
        }

        MPI_Reduce(&time_pnga_msg_sync, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_msg_sync,%ld,%lf\n", count_pnga_msg_sync, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbacc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbacc,%ld,%lf\n", count_pnga_nbacc, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbget, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbget,%ld,%lf\n", count_pnga_nbget, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbget_field, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbget_field,%ld,%lf\n", count_pnga_nbget_field, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbget_ghost_dir, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbget_ghost_dir,%ld,%lf\n", count_pnga_nbget_ghost_dir, recvbuf);
        }

        MPI_Reduce(&time_pnga_nblock, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nblock,%ld,%lf\n", count_pnga_nblock, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbput, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbput,%ld,%lf\n", count_pnga_nbput, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbput_field, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbput_field,%ld,%lf\n", count_pnga_nbput_field, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbtest, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbtest,%ld,%lf\n", count_pnga_nbtest, recvbuf);
        }

        MPI_Reduce(&time_pnga_nbwait, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nbwait,%ld,%lf\n", count_pnga_nbwait, recvbuf);
        }

        MPI_Reduce(&time_pnga_ndim, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_ndim,%ld,%lf\n", count_pnga_ndim, recvbuf);
        }

        MPI_Reduce(&time_pnga_nnodes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nnodes,%ld,%lf\n", count_pnga_nnodes, recvbuf);
        }

        MPI_Reduce(&time_pnga_nodeid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_nodeid,%ld,%lf\n", count_pnga_nodeid, recvbuf);
        }

        MPI_Reduce(&time_pnga_norm1, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_norm1,%ld,%lf\n", count_pnga_norm1, recvbuf);
        }

        MPI_Reduce(&time_pnga_norm_infinity, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_norm_infinity,%ld,%lf\n", count_pnga_norm_infinity, recvbuf);
        }

        MPI_Reduce(&time_pnga_pack, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pack,%ld,%lf\n", count_pnga_pack, recvbuf);
        }

        MPI_Reduce(&time_pnga_patch_enum, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_patch_enum,%ld,%lf\n", count_pnga_patch_enum, recvbuf);
        }

        MPI_Reduce(&time_pnga_patch_intersect, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_patch_intersect,%ld,%lf\n", count_pnga_patch_intersect, recvbuf);
        }

        MPI_Reduce(&time_pnga_periodic, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_periodic,%ld,%lf\n", count_pnga_periodic, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_absolute_id, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_absolute_id,%ld,%lf\n", count_pnga_pgroup_absolute_id, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_brdcst, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_brdcst,%ld,%lf\n", count_pnga_pgroup_brdcst, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_create, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_create,%ld,%lf\n", count_pnga_pgroup_create, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_destroy, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_destroy,%ld,%lf\n", count_pnga_pgroup_destroy, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_get_default, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_get_default,%ld,%lf\n", count_pnga_pgroup_get_default, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_get_mirror, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_get_mirror,%ld,%lf\n", count_pnga_pgroup_get_mirror, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_get_world, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_get_world,%ld,%lf\n", count_pnga_pgroup_get_world, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_gop, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_gop,%ld,%lf\n", count_pnga_pgroup_gop, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_nnodes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_nnodes,%ld,%lf\n", count_pnga_pgroup_nnodes, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_nodeid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_nodeid,%ld,%lf\n", count_pnga_pgroup_nodeid, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_set_default, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_set_default,%ld,%lf\n", count_pnga_pgroup_set_default, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_split, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_split,%ld,%lf\n", count_pnga_pgroup_split, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_split_irreg, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_split_irreg,%ld,%lf\n", count_pnga_pgroup_split_irreg, recvbuf);
        }

        MPI_Reduce(&time_pnga_pgroup_sync, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_pgroup_sync,%ld,%lf\n", count_pnga_pgroup_sync, recvbuf);
        }

        MPI_Reduce(&time_pnga_print, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print,%ld,%lf\n", count_pnga_print, recvbuf);
        }

        MPI_Reduce(&time_pnga_print_distribution, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print_distribution,%ld,%lf\n", count_pnga_print_distribution, recvbuf);
        }

        MPI_Reduce(&time_pnga_print_file, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print_file,%ld,%lf\n", count_pnga_print_file, recvbuf);
        }

        MPI_Reduce(&time_pnga_print_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print_patch,%ld,%lf\n", count_pnga_print_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_print_patch2d, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print_patch2d,%ld,%lf\n", count_pnga_print_patch2d, recvbuf);
        }

        MPI_Reduce(&time_pnga_print_patch_file, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print_patch_file,%ld,%lf\n", count_pnga_print_patch_file, recvbuf);
        }

        MPI_Reduce(&time_pnga_print_patch_file2d, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print_patch_file2d,%ld,%lf\n", count_pnga_print_patch_file2d, recvbuf);
        }

        MPI_Reduce(&time_pnga_print_stats, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_print_stats,%ld,%lf\n", count_pnga_print_stats, recvbuf);
        }

        MPI_Reduce(&time_pnga_proc_topology, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_proc_topology,%ld,%lf\n", count_pnga_proc_topology, recvbuf);
        }

        MPI_Reduce(&time_pnga_put, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_put,%ld,%lf\n", count_pnga_put, recvbuf);
        }

        MPI_Reduce(&time_pnga_put_field, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_put_field,%ld,%lf\n", count_pnga_put_field, recvbuf);
        }

        MPI_Reduce(&time_pnga_randomize, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_randomize,%ld,%lf\n", count_pnga_randomize, recvbuf);
        }

        MPI_Reduce(&time_pnga_read_inc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_read_inc,%ld,%lf\n", count_pnga_read_inc, recvbuf);
        }

        MPI_Reduce(&time_pnga_recip, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_recip,%ld,%lf\n", count_pnga_recip, recvbuf);
        }

        MPI_Reduce(&time_pnga_recip_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_recip_patch,%ld,%lf\n", count_pnga_recip_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_register_type, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_register_type,%ld,%lf\n", count_pnga_register_type, recvbuf);
        }

        MPI_Reduce(&time_pnga_release, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release,%ld,%lf\n", count_pnga_release, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_block, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_block,%ld,%lf\n", count_pnga_release_block, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_block_grid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_block_grid,%ld,%lf\n", count_pnga_release_block_grid, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_block_segment, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_block_segment,%ld,%lf\n", count_pnga_release_block_segment, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_ghost_element, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_ghost_element,%ld,%lf\n", count_pnga_release_ghost_element, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_ghosts,%ld,%lf\n", count_pnga_release_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_update, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_update,%ld,%lf\n", count_pnga_release_update, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_update_block, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_update_block,%ld,%lf\n", count_pnga_release_update_block, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_update_block_grid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_update_block_grid,%ld,%lf\n", count_pnga_release_update_block_grid, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_update_block_segment, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_update_block_segment,%ld,%lf\n", count_pnga_release_update_block_segment, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_update_ghost_element, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_update_ghost_element,%ld,%lf\n", count_pnga_release_update_ghost_element, recvbuf);
        }

        MPI_Reduce(&time_pnga_release_update_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_release_update_ghosts,%ld,%lf\n", count_pnga_release_update_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_scale, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scale,%ld,%lf\n", count_pnga_scale, recvbuf);
        }

        MPI_Reduce(&time_pnga_scale_cols, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scale_cols,%ld,%lf\n", count_pnga_scale_cols, recvbuf);
        }

        MPI_Reduce(&time_pnga_scale_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scale_patch,%ld,%lf\n", count_pnga_scale_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_scale_rows, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scale_rows,%ld,%lf\n", count_pnga_scale_rows, recvbuf);
        }

        MPI_Reduce(&time_pnga_scan_add, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scan_add,%ld,%lf\n", count_pnga_scan_add, recvbuf);
        }

        MPI_Reduce(&time_pnga_scan_copy, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scan_copy,%ld,%lf\n", count_pnga_scan_copy, recvbuf);
        }

        MPI_Reduce(&time_pnga_scatter, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scatter,%ld,%lf\n", count_pnga_scatter, recvbuf);
        }

        MPI_Reduce(&time_pnga_scatter2d, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scatter2d,%ld,%lf\n", count_pnga_scatter2d, recvbuf);
        }

        MPI_Reduce(&time_pnga_scatter_acc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scatter_acc,%ld,%lf\n", count_pnga_scatter_acc, recvbuf);
        }

        MPI_Reduce(&time_pnga_scatter_acc2d, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_scatter_acc2d,%ld,%lf\n", count_pnga_scatter_acc2d, recvbuf);
        }

        MPI_Reduce(&time_pnga_select_elem, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_select_elem,%ld,%lf\n", count_pnga_select_elem, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_array_name, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_array_name,%ld,%lf\n", count_pnga_set_array_name, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_block_cyclic, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_block_cyclic,%ld,%lf\n", count_pnga_set_block_cyclic, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_block_cyclic_proc_grid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_block_cyclic_proc_grid,%ld,%lf\n", count_pnga_set_block_cyclic_proc_grid, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_chunk, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_chunk,%ld,%lf\n", count_pnga_set_chunk, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_data, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_data,%ld,%lf\n", count_pnga_set_data, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_debug, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_debug,%ld,%lf\n", count_pnga_set_debug, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_diagonal, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_diagonal,%ld,%lf\n", count_pnga_set_diagonal, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_ghost_corner_flag, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_ghost_corner_flag,%ld,%lf\n", count_pnga_set_ghost_corner_flag, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_ghost_info, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_ghost_info,%ld,%lf\n", count_pnga_set_ghost_info, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_ghosts,%ld,%lf\n", count_pnga_set_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_irreg_distr, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_irreg_distr,%ld,%lf\n", count_pnga_set_irreg_distr, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_irreg_flag, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_irreg_flag,%ld,%lf\n", count_pnga_set_irreg_flag, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_memory_limit, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_memory_limit,%ld,%lf\n", count_pnga_set_memory_limit, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_pgroup, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_pgroup,%ld,%lf\n", count_pnga_set_pgroup, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_restricted, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_restricted,%ld,%lf\n", count_pnga_set_restricted, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_restricted_range, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_restricted_range,%ld,%lf\n", count_pnga_set_restricted_range, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_update4_info, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_update4_info,%ld,%lf\n", count_pnga_set_update4_info, recvbuf);
        }

        MPI_Reduce(&time_pnga_set_update5_info, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_set_update5_info,%ld,%lf\n", count_pnga_set_update5_info, recvbuf);
        }

        MPI_Reduce(&time_pnga_shift_diagonal, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_shift_diagonal,%ld,%lf\n", count_pnga_shift_diagonal, recvbuf);
        }

        MPI_Reduce(&time_pnga_solve, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_solve,%ld,%lf\n", count_pnga_solve, recvbuf);
        }

        MPI_Reduce(&time_pnga_spd_invert, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_spd_invert,%ld,%lf\n", count_pnga_spd_invert, recvbuf);
        }

        MPI_Reduce(&time_pnga_step_bound_info, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_step_bound_info,%ld,%lf\n", count_pnga_step_bound_info, recvbuf);
        }

        MPI_Reduce(&time_pnga_step_bound_info_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_step_bound_info_patch,%ld,%lf\n", count_pnga_step_bound_info_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_step_mask_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_step_mask_patch,%ld,%lf\n", count_pnga_step_mask_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_step_max, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_step_max,%ld,%lf\n", count_pnga_step_max, recvbuf);
        }

        MPI_Reduce(&time_pnga_step_max_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_step_max_patch,%ld,%lf\n", count_pnga_step_max_patch, recvbuf);
        }

        MPI_Reduce(&time_pnga_strided_acc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_strided_acc,%ld,%lf\n", count_pnga_strided_acc, recvbuf);
        }

        MPI_Reduce(&time_pnga_strided_get, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_strided_get,%ld,%lf\n", count_pnga_strided_get, recvbuf);
        }

        MPI_Reduce(&time_pnga_strided_put, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_strided_put,%ld,%lf\n", count_pnga_strided_put, recvbuf);
        }

        MPI_Reduce(&time_pnga_summarize, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_summarize,%ld,%lf\n", count_pnga_summarize, recvbuf);
        }

        MPI_Reduce(&time_pnga_symmetrize, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_symmetrize,%ld,%lf\n", count_pnga_symmetrize, recvbuf);
        }

        MPI_Reduce(&time_pnga_sync, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_sync,%ld,%lf\n", count_pnga_sync, recvbuf);
        }

        MPI_Reduce(&time_pnga_terminate, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_terminate,%ld,%lf\n", count_pnga_terminate, recvbuf);
        }

        MPI_Reduce(&time_pnga_timer, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_timer,%ld,%lf\n", count_pnga_timer, recvbuf);
        }

        MPI_Reduce(&time_pnga_total_blocks, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_total_blocks,%ld,%lf\n", count_pnga_total_blocks, recvbuf);
        }

        MPI_Reduce(&time_pnga_transpose, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_transpose,%ld,%lf\n", count_pnga_transpose, recvbuf);
        }

        MPI_Reduce(&time_pnga_type_c2f, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_type_c2f,%ld,%lf\n", count_pnga_type_c2f, recvbuf);
        }

        MPI_Reduce(&time_pnga_type_f2c, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_type_f2c,%ld,%lf\n", count_pnga_type_f2c, recvbuf);
        }

        MPI_Reduce(&time_pnga_unlock, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_unlock,%ld,%lf\n", count_pnga_unlock, recvbuf);
        }

        MPI_Reduce(&time_pnga_unpack, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_unpack,%ld,%lf\n", count_pnga_unpack, recvbuf);
        }

        MPI_Reduce(&time_pnga_update1_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update1_ghosts,%ld,%lf\n", count_pnga_update1_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update2_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update2_ghosts,%ld,%lf\n", count_pnga_update2_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update3_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update3_ghosts,%ld,%lf\n", count_pnga_update3_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update44_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update44_ghosts,%ld,%lf\n", count_pnga_update44_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update4_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update4_ghosts,%ld,%lf\n", count_pnga_update4_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update55_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update55_ghosts,%ld,%lf\n", count_pnga_update55_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update5_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update5_ghosts,%ld,%lf\n", count_pnga_update5_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update6_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update6_ghosts,%ld,%lf\n", count_pnga_update6_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update7_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update7_ghosts,%ld,%lf\n", count_pnga_update7_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_update_ghost_dir, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update_ghost_dir,%ld,%lf\n", count_pnga_update_ghost_dir, recvbuf);
        }

        MPI_Reduce(&time_pnga_update_ghosts, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_update_ghosts,%ld,%lf\n", count_pnga_update_ghosts, recvbuf);
        }

        MPI_Reduce(&time_pnga_uses_ma, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_uses_ma,%ld,%lf\n", count_pnga_uses_ma, recvbuf);
        }

        MPI_Reduce(&time_pnga_uses_proc_grid, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_uses_proc_grid,%ld,%lf\n", count_pnga_uses_proc_grid, recvbuf);
        }

        MPI_Reduce(&time_pnga_valid_handle, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_valid_handle,%ld,%lf\n", count_pnga_valid_handle, recvbuf);
        }

        MPI_Reduce(&time_pnga_verify_handle, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_verify_handle,%ld,%lf\n", count_pnga_verify_handle, recvbuf);
        }

        MPI_Reduce(&time_pnga_wtime, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_wtime,%ld,%lf\n", count_pnga_wtime, recvbuf);
        }

        MPI_Reduce(&time_pnga_zero, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_zero,%ld,%lf\n", count_pnga_zero, recvbuf);
        }

        MPI_Reduce(&time_pnga_zero_diagonal, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_zero_diagonal,%ld,%lf\n", count_pnga_zero_diagonal, recvbuf);
        }

        MPI_Reduce(&time_pnga_zero_patch, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("pnga_zero_patch,%ld,%lf\n", count_pnga_zero_patch, recvbuf);
        }

    }
}

