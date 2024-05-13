strategy_to_p_adic = function(strategy,max_n_markov)
    strategy[2:3,:] = strategy[2:3,:].-ones(Int64,2,max_n_markov)
    strategy = vec(strategy);
    number = 0;

    for    i = 1:length(strategy)
        number = number+ strategy[i]*max_n_markov^(i-1)
    end
    return number
end

max_n_markov = 4; # how many states at max

p_adic_to_max_n_state_strategy = function(number,max_n_markov)
    strategy = zeros(Int64,3,max_n_markov) + reshape(num_to_p_adic_v(number,3*max_n_markov,max_n_markov),3,max_n_markov)
    strategy[2:3,:] = strategy[2:3,:].+ones(Int64,2,max_n_markov)
    return strategy
end


reduce_to_connected = function(strategy)  
    indices_stay, indices_leave, reachable_stay, reachable_leave, n  = identify_disconnections(strategy);

    if (length(reachable_stay)+length(reachable_leave))<n
        d_indices = setdiff(collect(1:n), union(reachable_stay,reachable_leave))
        #now the idm'th column is deleted, and all numbers in strategy[2:end,:]>idm are reduced by one
        for idm = sort(d_indices, rev=true)
            if idm==size(strategy,2)
                strategy[:,idm].=0;
            else
                strategy[:,idm:end-1] = strategy[:,idm+1:end];
                strategy[:,size(strategy,2)].=0;
            end
            strategy[2:end,:]  =strategy[2:end,:].- (strategy[2:end,:].>idm);
        end
    end
    return strategy
end


num_to_p_adic_v = function(n,k,p)
    bin = zeros(Int64,k)
    for i = k:-1:1
        bin[i]=div(n,p^(i-1))
        n = mod(n,p^(i-1))
    end
    return bin
end

Strategy_number_v =collect(0:((max_n_markov)^(3*max_n_markov)-1));



clean_to_permitted = function(Strategy_number_v,max_n)
    cond_v = zeros(Int8,length(Strategy_number_v))
    for i = 1:length(Strategy_number_v)
        stra = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n)
        if stra[1,1]!=0 && all(stra[1,:].<3)
            cond_v[i] =1
        end
    end
    Strategy_number_v_new =Strategy_number_v[cond_v.==1]
    return Strategy_number_v_new
end



clean_to_unique = function(Strategy_number_v,max_n_markov)
    strategy_matrix = zeros(UInt64,length(Strategy_number_v),3);
    for i = 1:length(Strategy_number_v)
        println(i)
        stra = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
        stra = reduce_to_connected(stra);
        stra[2:end,stra[1,:].==0].=1;
        stra = hcat(stra,zeros(Int64,3,10-max_n_markov));
        stra = one_strategy_same_state_cleaning(stra,max_n);
        strategy_matrix[i,:] =  strat_3_to_number(stra,max_n,s5tn_v1,s5tn_v2)
    end

    SM = DataFrame(strategy_matrix,:auto);
    Strategy_number_v_new = Strategy_number_v[nonunique(SM).==0]
    return Strategy_number_v_new
end

using DataFrames


## take the same machinery that reduces a strategy to a unique representation 
## reduce everything to this unique representation in a dataframe 
## then simply take unique rows of that dataframe

# LOAD IN THESE TWO SCRIPTS FROM APPROPRIATE PATH
#include(".../single_run_analysis_functions.jl");
#include(".../all_functions_5d.jl");

number_to_strat_3_with_same_state_cleaning = function(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
    strategy_r = zeros(UInt128, 3, max_n);
    for i = max_n:-1:1
        strategy_r[1,i] = div(s[1],s5tn_v1_UInt[i]);
        s[1] = mod(s[1],s5tn_v1_UInt[i]);
        for j = 2:3
        strategy_r[j,i] = div(s[j],s5tn_v2_UInt[i]);
        s[j] = mod(s[j],s5tn_v2_UInt[i]);
        end
    end
    strategy_r = convert(Array{Int64},strategy_r)
    strategy_r = condense_strat_3(strategy_r, max_n)
    strategy_r = sort_row_ascend_3_version(strategy_r,max_n)
    
    strategy_r = same_state_claning(strategy_r)

    strategy_r = remove_equivalent_cycles(strategy_r) 


    return strategy_r
end

same_state_claning= function(strategy)
    d_indices = [];#
    k_indices = []; #everything pointing to the ith element of d_indices must point to k_indices[i]
    decimal_representation = [1 10 100]*strategy;
    i=2
    while i<= length(decimal_representation) && decimal_representation[i]!=0
        if any(decimal_representation[1:(i-1)].==decimal_representation[i])
        append!(d_indices,i)
        append!(k_indices,findfirst((decimal_representation[1:(i-1)].==decimal_representation[i])))
        end 
        i = i+1;
    end


    for i = 1:length(d_indices)
        strategy[2:end,:]= strategy[2:end,:].-(d_indices[i]-k_indices[i])*(strategy[2:end,:].==d_indices[i]);
    end
            for idm = sort(d_indices, rev=true)
                if idm==size(strategy,2)
                    strategy[:,idm].=0;
                else
                    strategy[:,idm:end-1] = strategy[:,idm+1:end];
                    strategy[:,size(strategy,2)].=0;
                end
                strategy[2:end,:]  =strategy[2:end,:].- (strategy[2:end,:].>idm);
            end

    strategy_new = hcat(strategy,zeros(Int64,3,10-size(strategy,2)))
    return strategy_new
end






one_strategy_same_state_cleaning = function(s,max_n)

    strategy_r = deepcopy(s)
    strategy_r = condense_strat_3(strategy_r, max_n)
    
    
    
    strategy_r = same_state_claning(strategy_r)
    strategy_r = remove_equivalent_cycles(strategy_r) 
    strategy_r = sort_row_ascend_3_version(strategy_r,max_n)

    return strategy_r
end





remove_equivalent_cycles = function(strategy)
    column_sets = find_column_sets(strategy)

    while length(column_sets)>0
        d_indices = column_sets[end][2:end]
        k_index = column_sets[end][1]    
        for i = 1:length(d_indices)
            strategy[2:end,:]= strategy[2:end,:].-(d_indices[i]-k_index)*(strategy[2:end,:].==d_indices[i]);
        end
                for idm = sort(d_indices, rev=true)
                    if idm==size(strategy,2)
                        strategy[:,idm].=0;
                    else
                        strategy[:,idm:end-1] = strategy[:,idm+1:end];
                        strategy[:,size(strategy,2)].=0;
                    end
                    strategy[2:end,:]  =strategy[2:end,:].- (strategy[2:end,:].>idm);
                end
    
        column_sets = find_column_sets(strategy)
    end

    strategy_new = hcat(strategy,zeros(Int64,3,10-size(strategy,2)))

    return strategy_new 
end


function all_in_vector(matrix, vector)
    all(elem in vector for elem in Iterators.flatten(matrix))
end

using Combinatorics

function find_column_sets(matrix)
    n, m = size(matrix)
    column_sets = []
    for size in 2:m # iterate over all possible set sizes
        for col_set in combinations(1:m, size) # generate all possible column sets of the given size
            #println(col_set)
            if allsame(matrix[1, col_set]) &&all_in_vector(matrix[2:end, col_set],col_set)
                push!(column_sets, col_set) # add the current set to the list
            end
        end
    end
    return column_sets
end

