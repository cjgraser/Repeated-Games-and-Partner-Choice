include("functions_to_fill_markov_matrix.jl")




b = 2; c = 1; PM = zeros(2,2);PM[1,:] = [b-c -c];PM[2,:] = [b 0];PM = PM+2*c*ones(2,2);
delta = 0.8;tolerance = 0.01;
#
EB11,EB12,P11,P12  =generate_P_and_EB_matrices(Strategy_number_v ,PM,delta,tolerance,max_n_markov)


N=100;

Big_Markov = leave_markov_fill(Strategy_number_v,N, EB11,EB12,P11,P12)


Big_Markov =weights_Markov_by_uniform_Muataion_kernel!(Big_Markov,Strategy_number_v);



using DiscreteMarkovChains

chain = DiscreteMarkovChain(Big_Markov)

sum(stationary_distribution(chain))
s = stationary_distribution(chain)

using Plots


clean_to_permitted_without_leaves = function(Strategy_number_v,max_n_markov)
    cond_v = zeros(Int8,length(Strategy_number_v))
    for i = 1:length(Strategy_number_v)
        stra = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
        stra = reduce_to_connected(stra)  #including leave state disconnections
        if all(stra[1,stra[2,:].!=0].!=0)
            cond_v[i] =1
        end
    end
    Strategy_number_v_new =Strategy_number_v[cond_v.==1]
    return Strategy_number_v_new, cond_v
end

strategy_number_v_no_leaves,cond_v = clean_to_permitted_without_leaves(Strategy_number_v,max_n_markov)


Big_Markov_no_leaves = Big_Markov[(cond_v.==1),(cond_v.==1)];

Big_Markov_no_leaves_corrected = deepcopy(Big_Markov_no_leaves);
for i in 1:size(Big_Markov_no_leaves,1)
    Big_Markov_no_leaves_corrected[i,:] =     Big_Markov_no_leaves_corrected[i,:].*(length(Strategy_number_v)/length(strategy_number_v_no_leaves))
    Big_Markov_no_leaves_corrected[i,i] = 1-sum(Big_Markov_no_leaves_corrected[i,:])+Big_Markov_no_leaves_corrected[i,i]
end

index_mapping_no_leaves_to_leaves = collect(1:length(Strategy_number_v))[cond_v.==1]; #...no_leaves[i] is the ...leaves[index_mapping...[i]] index

