/* FILE GENERATED AUTOMATICALLY, DO NOT EDIT */

#include <interfaces/highs_c_api.h>
#include "HiGHSMILPSolver.h"

using namespace SMSpp_di_unipi_it;

const std::array< std::string, HiGHS_NUM_INT_PARS > HiGHSMILPSolver::SMSpp_to_HiGHS_int_pars{
 "random_seed",
 "threads",
 "highs_debug_level",
 "highs_analysis_level",
 "simplex_strategy",
 "simplex_scale_strategy",
 "simplex_crash_strategy",
 "simplex_dual_edge_weight_strategy",
 "simplex_primal_edge_weight_strategy",
 "simplex_iteration_limit",
 "simplex_update_limit",
 "simplex_min_concurrency",
 "simplex_max_concurrency",
 "output_flag",
 "log_to_console",
 "write_solution_to_file",
 "write_solution_style",
 "glpsol_cost_row_location",
 "icrash",
 "icrash_dualize",
 "icrash_iterations",
 "icrash_approx_iter",
 "icrash_exact",
 "icrash_breakpoints",
 "write_model_to_file",
 "mip_detect_symmetry",
 "mip_max_nodes",
 "mip_max_stall_nodes",
 "mip_improving_solution_save",
 "mip_improving_solution_report_sparse",
 "mip_max_leaves",
 "mip_max_improving_sols",
 "mip_lp_age_limit",
 "mip_pool_age_limit",
 "mip_pool_soft_limit",
 "mip_pscost_minreliable",
 "mip_min_cliquetable_entries_for_parallelism",
 "mip_report_level",
 "ipm_iteration_limit",
 "log_dev_level",
 "solve_relaxation",
 "allow_unbounded_or_infeasible",
 "use_implied_bounds_from_presolve",
 "lp_presolve_requires_basis_postsolve",
 "mps_parser_type_free",
 "keep_n_rows",
 "cost_scale_factor",
 "allowed_matrix_scale_factor",
 "allowed_cost_scale_factor",
 "ipx_dualize_strategy",
 "simplex_dualize_strategy",
 "simplex_permute_strategy",
 "max_dual_simplex_cleanup_level",
 "max_dual_simplex_phase1_cleanup_level",
 "simplex_price_strategy",
 "simplex_unscaled_solution_strategy",
 "simplex_initial_condition_check",
 "no_unnecessary_rebuild_refactor",
 "presolve_reduction_limit",
 "presolve_rule_off",
 "presolve_rule_logging",
 "presolve_substitution_maxfillin",
 "use_original_HFactor_logic",
 "less_infeasible_DSE_check",
 "less_infeasible_DSE_choose_row",
};

const std::array< std::string, HiGHS_NUM_DBL_PARS > HiGHSMILPSolver::SMSpp_to_HiGHS_dbl_pars{
 "time_limit",
 "infinite_cost",
 "infinite_bound",
 "small_matrix_value",
 "large_matrix_value",
 "primal_feasibility_tolerance",
 "dual_feasibility_tolerance",
 "ipm_optimality_tolerance",
 "objective_bound",
 "objective_target",
 "icrash_starting_weight",
 "mip_feasibility_tolerance",
 "mip_heuristic_effort",
 "mip_rel_gap",
 "mip_abs_gap",
 "mip_min_logging_interval",
 "simplex_initial_condition_tolerance",
 "rebuild_refactor_solution_error_tolerance",
 "dual_steepest_edge_weight_error_tolerance",
 "dual_steepest_edge_weight_log_error_threshold",
 "dual_simplex_cost_perturbation_multiplier",
 "primal_simplex_bound_perturbation_multiplier",
 "dual_simplex_pivot_growth_tolerance",
 "presolve_pivot_threshold",
 "factor_pivot_threshold",
 "factor_pivot_tolerance",
 "start_crossover_tolerance",
};

const std::array< std::string, HiGHS_NUM_STR_PARS > HiGHSMILPSolver::SMSpp_to_HiGHS_str_pars{
 "presolve",
 "solver",
 "parallel",
 "run_crossover",
 "ranging",
 "solution_file",
 "log_file",
 "icrash_strategy",
 "write_model_file",
 "mip_improving_solution_file",
};

const std::array< std::pair< std::string, int >, HiGHS_NUM_INT_PARS >
 HiGHSMILPSolver::HiGHS_to_SMSpp_int_pars{
 {
  { "random_seed", intFirstHiGHSPar + 0 },
  { "threads", intFirstHiGHSPar + 1 },
  { "highs_debug_level", intFirstHiGHSPar + 2 },
  { "highs_analysis_level", intFirstHiGHSPar + 3 },
  { "simplex_strategy", intFirstHiGHSPar + 4 },
  { "simplex_scale_strategy", intFirstHiGHSPar + 5 },
  { "simplex_crash_strategy", intFirstHiGHSPar + 6 },
  { "simplex_dual_edge_weight_strategy", intFirstHiGHSPar + 7 },
  { "simplex_primal_edge_weight_strategy", intFirstHiGHSPar + 8 },
  { "simplex_iteration_limit", intFirstHiGHSPar + 9 },
  { "simplex_update_limit", intFirstHiGHSPar + 10 },
  { "simplex_min_concurrency", intFirstHiGHSPar + 11 },
  { "simplex_max_concurrency", intFirstHiGHSPar + 12 },
  { "output_flag", intFirstHiGHSPar + 13 },
  { "log_to_console", intFirstHiGHSPar + 14 },
  { "write_solution_to_file", intFirstHiGHSPar + 15 },
  { "write_solution_style", intFirstHiGHSPar + 16 },
  { "glpsol_cost_row_location", intFirstHiGHSPar + 17 },
  { "icrash", intFirstHiGHSPar + 18 },
  { "icrash_dualize", intFirstHiGHSPar + 19 },
  { "icrash_iterations", intFirstHiGHSPar + 20 },
  { "icrash_approx_iter", intFirstHiGHSPar + 21 },
  { "icrash_exact", intFirstHiGHSPar + 22 },
  { "icrash_breakpoints", intFirstHiGHSPar + 23 },
  { "write_model_to_file", intFirstHiGHSPar + 24 },
  { "mip_detect_symmetry", intFirstHiGHSPar + 25 },
  { "mip_max_nodes", intFirstHiGHSPar + 26 },
  { "mip_max_stall_nodes", intFirstHiGHSPar + 27 },
  { "mip_improving_solution_save", intFirstHiGHSPar + 28 },
  { "mip_improving_solution_report_sparse", intFirstHiGHSPar + 29 },
  { "mip_max_leaves", intFirstHiGHSPar + 30 },
  { "mip_max_improving_sols", intFirstHiGHSPar + 31 },
  { "mip_lp_age_limit", intFirstHiGHSPar + 32 },
  { "mip_pool_age_limit", intFirstHiGHSPar + 33 },
  { "mip_pool_soft_limit", intFirstHiGHSPar + 34 },
  { "mip_pscost_minreliable", intFirstHiGHSPar + 35 },
  { "mip_min_cliquetable_entries_for_parallelism", intFirstHiGHSPar + 36 },
  { "mip_report_level", intFirstHiGHSPar + 37 },
  { "ipm_iteration_limit", intFirstHiGHSPar + 38 },
  { "log_dev_level", intFirstHiGHSPar + 39 },
  { "solve_relaxation", intFirstHiGHSPar + 40 },
  { "allow_unbounded_or_infeasible", intFirstHiGHSPar + 41 },
  { "use_implied_bounds_from_presolve", intFirstHiGHSPar + 42 },
  { "lp_presolve_requires_basis_postsolve", intFirstHiGHSPar + 43 },
  { "mps_parser_type_free", intFirstHiGHSPar + 44 },
  { "keep_n_rows", intFirstHiGHSPar + 45 },
  { "cost_scale_factor", intFirstHiGHSPar + 46 },
  { "allowed_matrix_scale_factor", intFirstHiGHSPar + 47 },
  { "allowed_cost_scale_factor", intFirstHiGHSPar + 48 },
  { "ipx_dualize_strategy", intFirstHiGHSPar + 49 },
  { "simplex_dualize_strategy", intFirstHiGHSPar + 50 },
  { "simplex_permute_strategy", intFirstHiGHSPar + 51 },
  { "max_dual_simplex_cleanup_level", intFirstHiGHSPar + 52 },
  { "max_dual_simplex_phase1_cleanup_level", intFirstHiGHSPar + 53 },
  { "simplex_price_strategy", intFirstHiGHSPar + 54 },
  { "simplex_unscaled_solution_strategy", intFirstHiGHSPar + 55 },
  { "simplex_initial_condition_check", intFirstHiGHSPar + 56 },
  { "no_unnecessary_rebuild_refactor", intFirstHiGHSPar + 57 },
  { "presolve_reduction_limit", intFirstHiGHSPar + 58 },
  { "presolve_rule_off", intFirstHiGHSPar + 59 },
  { "presolve_rule_logging", intFirstHiGHSPar + 60 },
  { "presolve_substitution_maxfillin", intFirstHiGHSPar + 61 },
  { "use_original_HFactor_logic", intFirstHiGHSPar + 62 },
  { "less_infeasible_DSE_check", intFirstHiGHSPar + 63 },
  { "less_infeasible_DSE_choose_row", intFirstHiGHSPar + 64 },
 }
};

const std::array< std::pair< std::string, int >, HiGHS_NUM_DBL_PARS >
 HiGHSMILPSolver::HiGHS_to_SMSpp_dbl_pars{
 {
  { "time_limit", dblFirstHiGHSPar + 0 },
  { "infinite_cost", dblFirstHiGHSPar + 1 },
  { "infinite_bound", dblFirstHiGHSPar + 2 },
  { "small_matrix_value", dblFirstHiGHSPar + 3 },
  { "large_matrix_value", dblFirstHiGHSPar + 4 },
  { "primal_feasibility_tolerance", dblFirstHiGHSPar + 5 },
  { "dual_feasibility_tolerance", dblFirstHiGHSPar + 6 },
  { "ipm_optimality_tolerance", dblFirstHiGHSPar + 7 },
  { "objective_bound", dblFirstHiGHSPar + 8 },
  { "objective_target", dblFirstHiGHSPar + 9 },
  { "icrash_starting_weight", dblFirstHiGHSPar + 10 },
  { "mip_feasibility_tolerance", dblFirstHiGHSPar + 11 },
  { "mip_heuristic_effort", dblFirstHiGHSPar + 12 },
  { "mip_rel_gap", dblFirstHiGHSPar + 13 },
  { "mip_abs_gap", dblFirstHiGHSPar + 14 },
  { "mip_min_logging_interval", dblFirstHiGHSPar + 15 },
  { "simplex_initial_condition_tolerance", dblFirstHiGHSPar + 16 },
  { "rebuild_refactor_solution_error_tolerance", dblFirstHiGHSPar + 17 },
  { "dual_steepest_edge_weight_error_tolerance", dblFirstHiGHSPar + 18 },
  { "dual_steepest_edge_weight_log_error_threshold", dblFirstHiGHSPar + 19 },
  { "dual_simplex_cost_perturbation_multiplier", dblFirstHiGHSPar + 20 },
  { "primal_simplex_bound_perturbation_multiplier", dblFirstHiGHSPar + 21 },
  { "dual_simplex_pivot_growth_tolerance", dblFirstHiGHSPar + 22 },
  { "presolve_pivot_threshold", dblFirstHiGHSPar + 23 },
  { "factor_pivot_threshold", dblFirstHiGHSPar + 24 },
  { "factor_pivot_tolerance", dblFirstHiGHSPar + 25 },
  { "start_crossover_tolerance", dblFirstHiGHSPar + 26 },
 }
};

const std::array< std::pair< std::string, int >, HiGHS_NUM_STR_PARS >
 HiGHSMILPSolver::HiGHS_to_SMSpp_str_pars{
 {
  { "presolve", strFirstHiGHSPar + 0 },
  { "solver", strFirstHiGHSPar + 1 },
  { "parallel", strFirstHiGHSPar + 2 },
  { "run_crossover", strFirstHiGHSPar + 3 },
  { "ranging", strFirstHiGHSPar + 4 },
  { "solution_file", strFirstHiGHSPar + 5 },
  { "log_file", strFirstHiGHSPar + 6 },
  { "icrash_strategy", strFirstHiGHSPar + 7 },
  { "write_model_file", strFirstHiGHSPar + 8 },
  { "mip_improving_solution_file", strFirstHiGHSPar + 9 },
 }
};