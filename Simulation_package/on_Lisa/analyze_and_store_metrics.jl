#sourcing

#args: path to load from, file to load, where to store

include(ARGS[1]*"single_run_analysis_functions.jl");  
include(ARGS[1]*"all_functions_5d.jl");  
include(ARGS[1]*"automated_comparison_functions.jl");  

using JLD2
population_states = load_object(ARGS[2])


population_states = reduce_population_states_to_connected_to_pst_to_3(population_states, max_n, s5tn_v1_UInt,s5tn_v2_UInt)
PST_DF, DF_cond = DF_cond_from_population_states_3(population_states,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
save_top_x_by_k_3_version(ARGS[3]*ARGS[4]*"_top_20_pure.txt", DF_cond,20,1)
save_top_x_by_k_3_version(ARGS[3]*ARGS[4]*"_top_20_mixed.txt", DF_cond,20,2)
