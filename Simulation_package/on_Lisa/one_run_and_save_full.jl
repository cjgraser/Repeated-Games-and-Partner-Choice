#no_leaves = 1; 
#trimming = -1; --> these are parsed as arguments





a5 =ARGS[5]
bi = parse(Int64, ARGS[1]);
betai = parse(Int64, ARGS[2]);
no_leaves= parse(Int64, ARGS[3]);
trimming = parse(Int64, ARGS[4]); #args 5 and 6 are the source direcotry for iterators_and_parameters and the destination directory to store to
include(ARGS[5]*"iterators_and_parameters.jl");  #this takes shortly compared to the whole, so doesn't matter to import multiple times


match_types, population_states, errors, leaves = execute_in_function(bi,betai,no_leaves, trimming)

if trimming ==1
    trimv = "_trimming"
elseif trimming  ==0
    trimv = "_reattach"
elseif trimming ==-1
    trimv = "_original_trim"
end

if no_leaves ==1
    no_leavev = "_no_leaves"
else
    no_leavev = "_leaves"
end


name_id = "STORED_b_"*string(bv[bi])*"_beta_"*string(betav[betai])*no_leavev*trimv;

using JLD2
save_object(ARGS[6]*name_id*"_match_types.jld2",match_types)
save_object(ARGS[6]*name_id*"_population_states.jld2",population_states)
save_object(ARGS[6]*name_id*"_leaves.jld2",leaves)
#save_object(ARGS[6]*name_id*"_betters.jld2",betters)



