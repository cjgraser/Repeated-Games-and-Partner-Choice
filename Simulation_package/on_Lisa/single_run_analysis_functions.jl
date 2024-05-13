using DataFrames
using Statistics
repeats = 3
max_n = 10
base_s5tn_v1 = 3;
s5tn_v1 = base_s5tn_v1.^collect(0:(max_n-1)); 
s5tn_v2 = (max_n+1).^collect(0:(max_n-1)); 
s5tn_v1_UInt = convert(Vector{UInt128}, s5tn_v1)
s5tn_v2_UInt = convert(Vector{UInt128}, s5tn_v2)

number_to_strat_3 = function(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
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

    return strategy_r
end

number_to_strat_3_no_cleaning = function(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
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
    return strategy_r
end

number_to_strat_5 = function(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
    strategy_r = zeros(UInt128, 5, max_n);
    for i = max_n:-1:1
        strategy_r[1,i] = div(s[1],s5tn_v1_UInt[i]);
        s[1] = mod(s[1],s5tn_v1_UInt[i]);
        for j = 2:5
        strategy_r[j,i] = div(s[j],s5tn_v2_UInt[i]);
        s[j] = mod(s[j],s5tn_v2_UInt[i]);
        end
    end
    strategy_r = convert(Array{Int64},strategy_r)
    strategy_r = condense_strat(strategy_r, max_n)
    strategy_r = sort_row_ascend(strategy_r,max_n)
    return strategy_r
end



number_to_strat_5_no_cleaning = function(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
    strategy_r = zeros(UInt128, 5, max_n);
    for i = max_n:-1:1
        strategy_r[1,i] = div(s[1],s5tn_v1_UInt[i]);
        s[1] = mod(s[1],s5tn_v1_UInt[i]);
        for j = 2:5
        strategy_r[j,i] = div(s[j],s5tn_v2_UInt[i]);
        s[j] = mod(s[j],s5tn_v2_UInt[i]);
        end
    end
    strategy_r = convert(Array{Int64},strategy_r)
    return strategy_r
end

allsame(x) = all(y -> y == first(x), x)


condense_strat_3 = function(strategy, max_n)
    #find n 
    if strategy[2,end]!=0 
        n = max_n;
    else
        n =-1+ findfirst( x -> strategy[2,x] ==0, 1:max_n);
    end
    #
    if allsame(strategy[1,1:n])
        strategy[:,1] = [strategy[1,1] , 1 , 1]
        strategy[:,2:end].=0;
    end
    return(strategy)
end


condense_strat = function(strategy, max_n)
    #find n 
    if strategy[2,end]!=0 
        n = max_n;
    else
        n =-1+ findfirst( x -> strategy[2,x] ==0, 1:max_n);
    end
    #
    if allsame(strategy[1,1:n])
        strategy[:,1] = [strategy[1,1] , 1 , 1, 1, 1]
        strategy[:,2:end].=0;
    end
    return(strategy)
end

condense_data_frame = function(DF,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)

    for i = 1:nrow(DF)
        DF = inner_loop_simplifying(DF,i,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    end
 
    for i = 1:nrow(DF)
        for j = repeats:-1:2
            for k = (j-1):-1:1 
                if collect(DF[i,(5*(j-1)+1):5*j])==collect(DF[i,(5*(k-1)+1):5*k])
                    DF[i,(5*(k-1)+1):5*(repeats-1)] = collect(DF[i,(5*(k)+1):5*(repeats)]);
                    DF[i,(5*(repeats-1)+1):5*repeats].=0 
                end
            end 
        end
        println(i/nrow(DF))
    end
    return DF
end


inner_loop_simplifying = function(DF,i,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    for j = 1:repeats
        s = collect(DF[i,(5*(j-1)+1):5*j]);
        if sum(s.!=0)>0
        s = convert(Array{UInt128}, s);
        stra = number_to_strat_5(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt);  
        stra = strat_5_to_number(stra,max_n,s5tn_v1,s5tn_v2);
        DF[i,(5*(j-1)+1):5*j] = stra;
        end
    end
    return DF
end




condense_data_frame_3 = function(DF,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    for i = 1:nrow(DF)
        DF = inner_loop_simplifying_3(DF,i,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    end
    for i = 1:nrow(DF)
        for j = repeats:-1:2
            for k = (j-1):-1:1 
                if collect(DF[i,(3*(j-1)+1):3*j])==collect(DF[i,(3*(k-1)+1):3*k])
                    DF[i,(3*(k-1)+1):3*(repeats-1)] = collect(DF[i,(3*(k)+1):3*(repeats)]);
                    DF[i,(3*(repeats-1)+1):3*repeats].=0 
                end
            end 
        end
        println(i/nrow(DF))
    end
    return DF
end


inner_loop_simplifying_3 = function(DF,i,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    for j = 1:repeats
        s = collect(DF[i,(3*(j-1)+1):3*j]);
        if sum(s.!=0)>0
        s = convert(Array{UInt128}, s);
        stra = number_to_strat_3(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt); 
        stra = strat_3_to_number(stra,max_n,s5tn_v1,s5tn_v2);
        DF[i,(3*(j-1)+1):3*j] = stra;
        end
    end
    return DF
end




sort_row_ascend_3_version =function(strategy,max_n)
    if strategy[2,end]!=0 
        n = max_n;
    else
        n =-1+ findfirst( x -> strategy[2,x] ==0, 1:max_n);
    end
    stop = 0;
    iterations = 0; 
    while stop ==0 && iterations< ((n-1)*(n-2))/2
    stop = 1
        for i=2:(n-1)
            iterations = iterations+1;
            if compare_two_3_version(strategy,i,i+1)
                stop = 0
                strategy = switch_two_states(strategy,i,i+1) 
            end 
        end
    end
    return strategy
    end
    



sort_row_ascend =function(strategy,max_n)
    if strategy[2,end]!=0 
        n = max_n;
    else
        n =-1+ findfirst( x -> strategy[2,x] ==0, 1:max_n);
    end
    stop = 0;
    iterations = 0; 
    while stop ==0 && iterations< ((n-1)*(n-2))/2
    stop = 1
        for i=2:(n-1)
            iterations = iterations+1;
            if compare_two(strategy,i,i+1)
                stop = 0
                strategy = switch_two_states(strategy,i,i+1) 
            end 
        end
    end
    return strategy
    end
    
    
compare_two = function(strategy, i,j)
        #switch = false #j is the one that is later than i. So switches happen if j has the lower number than i
        if strategy[1,i] > strategy[1,j]
            switch = true
        elseif strategy[1,i] == strategy[1,j] && strategy[2,i] > strategy[2,j]
            switch = true
        elseif strategy[1,i] == strategy[1,j] && strategy[2,i] == strategy[2,j] && strategy[3,i] > strategy[3,j]
            switch = true
        elseif strategy[1,i] == strategy[1,j] && strategy[2,i] == strategy[2,j] && strategy[3,i] == strategy[3,j] && strategy[4,i] > strategy[4,j]
            switch = true
        elseif strategy[1,i] == strategy[1,j] && strategy[2,i] == strategy[2,j] && strategy[3,i] == strategy[3,j] && strategy[4,i] == strategy[4,j] && strategy[5,i] > strategy[5,j]
            switch = true
        else 
            switch = false
        end
        return switch
    end



compare_two_3_version = function(strategy, i,j)
        #switch = false #j is the one that is later than i. So switches happen if j has the lower number than i
        if strategy[1,i] > strategy[1,j]
            switch = true
        elseif strategy[1,i] == strategy[1,j] && strategy[2,i] > strategy[2,j]
            switch = true
        elseif strategy[1,i] == strategy[1,j] && strategy[2,i] == strategy[2,j] && strategy[3,i] > strategy[3,j]
            switch = true
        else 
            switch = false
        end
        return switch
end
    
switch_two_states = function(strategy,i,j)
    strategy[:,[i,j]]=strategy[:,[j,i]]
    strategy[2:end,:]  =strategy[2:end,:].+(strategy[2:end,:].==j).*(i-j).+(strategy[2:end,:].==i).*(j-i);
    return strategy
end
    
    
    
    
unreached_neutral_check = function(strategy1,strategy2)
state1 = 1
state2 = 1
exit_state = 0 
exit_other = 0
counts = 0
while exit_state ==0 &&  exit_other == 0
    counts = counts+1
    if strategy1[1,state1] != strategy2[1,state2]  # they must always behave the same
        exit_state = 1
        if strategy1[1,state1] ==0
            exit_other = 1 #just exit here because the game ends
        end        
    else 
        state1_new = strategy1[strategy2[1,state2]+1,state1]
        state2 = strategy2[strategy1[1,state1]+1,state2]
        state1 = state1_new
    end

    if counts>(size(strategy,2))*(size(strategy2,2)-1) #this is a rough proxy
        exit_other=2
    end


end

return exit_state, exit_other
end


invading_check = function(strategy1, strategy2, ipayoff2b2)
    #the answer that this one gives is: what does delta have to be s.t. 2 can invade
    # sum payoffs until convergence, or leave 
    # USE THIS FUNCITON, if strategy1 and strategy 2 are different, as determined by 'unreached_neutral_check', (and not on a leave state)
    #...meaning they behave differently on some state, before a leave occurs (or no leaves occur ever)
    payoff_vector_1 = []
    payoff_vector_2 = []
    state1 = 1
    state2 = 1

    leave_happened = 0
    count = 0
    while leave_happened ==0 && count<30

        if     strategy1[1,state1]==0 || strategy2[1,state2]==0
            leave_happened = 1
        else
        payoff_vector_1 = append!(payoff_vector_1, ipayoff2b2[strategy1[1,state1], strategy2[1,state2]])
        payoff_vector_2 = append!(payoff_vector_2, ipayoff2b2[strategy2[1,state2], strategy1[1,state1]])
        
        state1_new = strategy1[strategy2[1,state2]+1,state1]
        state2 = strategy2[strategy1[1,state1]+1,state2]
        state1 = state1_new
        end
        
        count = count+1 
        #one good criterion for convergence is if the states stay the same. sufficient but not necessary
        # i'd say, let it run for max 30 periods.... that's a good enough proxy
    end 

    #Then, i have a payoff vector: 
    payoff_diffs_delta = zeros(1,length(payoff_vector_1))
    i=1
    for delta = 0.7:0.05:1
        payoff_diffs_delta[1,i] = sum([delta^x for x in 0:(length(payoff_vector_1)-1) ].* (payoff_vector_1 - payoff_vector_2))
        i=i+1
    end     

    return payoff_diffs_delta 
end

neutral_invading_check_against_itself = function(strategy1)
    payoff_vector_1 = []
    state1 = 1
    leave_happened = 0
    count = 0
    while leave_happened ==0 && count<30

        if     strategy1[1,state1]==0 
            leave_happened = 1
        else
        payoff_vector_1 = append!(payoff_vector_1, ipayoff2b2[strategy1[1,state1], strategy1[state1]])
        state1 = strategy1[strategy1[1,state1]+1,state1]
        end
        
        count = count+1 
    end 
    payoffs_delta = zeros(1,length(payoff_vector_1))
    i=1
    for delta = 0.7:0.05:1
        payoffs_delta[1,i] = sum([delta^x for x in 0:(length(payoff_vector_1)-1) ].* (payoff_vector_1))
        i=i+1
    end     
    return payoffs_delta

end




##########

prepare_dataframe = function(population_states,number_individuals_threshold,number_of_strategies_maximum,three_or_five,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    pst_matrix = zeros(length(population_states),three_or_five*number_of_strategies_maximum)
    for i = 1:length(population_states)
        pst = population_states[i]
        pst = sort!(pst, [:nrow], rev=true)
        j=1;
        while j<=size(pst,1)&& pst[j,6]>number_individuals_threshold && j<=number_of_strategies_maximum 
            pst_matrix[i,(three_or_five*(j-1)+1):three_or_five*j].= collect(pst[j,1:three_or_five]);
            j= j+1
        end
    end
    PST_DF = DataFrame(pst_matrix, :auto)
    DF_cond = condense_data_frame(PST_DF,max_n,number_of_strategies_maximum,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    return PST_DF,DF_cond
end



### 
pic_and_strats_5 = function(ipath,PST_DF,DF_cond)
    open(ipath, "w") do file
    iprev = 0;
    for i = 1:ntime-1
        if collect(DF_cond[i,1:end-1]) != collect(DF_cond[i+1,1:end-1])
            j=1;
            s1 = collect(PST_DF[i,(5*(j-1)+1):5*j]);
            s1 = convert(Array{UInt128}, s1)
            s2 = collect(PST_DF[i+1,(5*(j-1)+1):5*j]);
            s2 = convert(Array{UInt128}, s2)
#            println(i)
                write(file, "\n")
                write(file, "$iprev to $i:  \n")
            for j = 1:3
                s = collect(PST_DF[i,(5*(j-1)+1):5*j]);
                s = convert(Array{UInt128}, s)
                stra = number_to_strat_5(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
                if stra[1,1]!=0 
                    stra1 = stra[1,:]
                    stra2 = stra[2,:]
                    #stra1 = action_to_symbol(stra1)
                    stra3 = stra[3,:]   
                    stra4 = stra[4,:]
                    stra5 = stra[5,:]
                        write(file, "$stra1 Strategy\n")
                        write(file, "$stra2 CC\n")
                        write(file, "$stra3 CD\n")
                        write(file, "$stra4 DC\n")
                        write(file, "$stra5 DD\n")
                        write(file, "\n")
                end
            end
            iprev = i;
        end
    end
    end
#    return 
end


######



action_to_symbol = function(v1)
    v2 = Array{String}(undef,size(v1,1)) 
    for i = 1:size(v1,1)
        if v1[i]==0
            v2[i] = "X"
        elseif v1[i] ==1
            v2[i] = "C"
        elseif v1[i]==2
            v2[i] = "D"
        end
    end
    return v2
end



strat_3_to_number = function(strat_,max_n,s5tn_v1,s5tn_v2)
    strat = copy(strat_)
    s = zeros(UInt64, 3);
    s[1]=strat[1,:]'*s5tn_v1
    s[2]=strat[2,:]'*s5tn_v2
    s[3]=strat[3,:]'*s5tn_v2
    return s
end


save_top_x_by_k = function(storepath, df2,x,k) #order looking at the k most prevalent strategies (k ={1,2,3}), and store the x most 
    if k ==1
        df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5]), nrow),[:nrow],rev=true)
    elseif k==2
        df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10]), nrow),[:nrow],rev=true)
    else #any number not 1 or 2 will be interpreted as 3
        df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14, :x15]), nrow),[:nrow],rev=true)
    end

    open(storepath, "w") do file
    for i = 1:x
    write(file,"\nrank $i ")
    dd = df2[i,end]
    write(file, "\ncount: $dd \n")
    for j = 1:k
        s = collect(df2[i,(5*(j-1)+1):5*j]);
        s = convert(Array{UInt128}, s)
        stra = number_to_strat_5(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
        if stra[1,1]!=0 
        stra1 = stra[1,:]
        stra2 = stra[2,:]
        stra3 = stra[3,:]   
        stra4 = stra[4,:]
        stra5 = stra[5,:]
            write(file, "$stra1 Strategy\n")
            write(file, "$stra2 CC\n")
            write(file, "$stra3 CD\n")
            write(file, "$stra4 DC\n")
            write(file, "$stra5 DD\n")
            write(file, "\n")
        end
    end
    end
    end

end