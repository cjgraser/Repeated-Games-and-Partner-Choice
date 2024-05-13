using NLsolve
include("../../TWO STATES/Full_functions_relatedness.jl")
include("general_markov.jl")

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


clean_to_3_cd_state = function(Strategy_number_v,max_n_markov)
    cond_v = zeros(Int8,length(Strategy_number_v))
    for i = 1:length(Strategy_number_v)
        stra = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
    
        if sum(stra[1,:].!=0)<=3
            cond_v[i] =1
        end
    end
    Strategy_number_v_new =Strategy_number_v[cond_v.==1]
    return Strategy_number_v_new, cond_v
end


clean_to_3_cd_state_permitted_without_leaves = function(Strategy_number_v,max_n_markov)
    cond_v = zeros(Int8,length(Strategy_number_v))
    for i = 1:length(Strategy_number_v)
        stra = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
        stra = reduce_to_connected(stra) 
        if all(stra[1,stra[2,:].!=0].!=0) && sum(stra[1,:].!=0)<=3
            cond_v[i] =1
        end
    end
    Strategy_number_v_new =Strategy_number_v[cond_v.==1]
    return Strategy_number_v_new, cond_v
end


weights_Markov_by_uniform_Muataion_kernel! =function(Big_Markov,Strategy_number_v)
    Weight_Matrix = zeros(length(Strategy_number_v),length(Strategy_number_v))
    
    for i = 1:length(Strategy_number_v)
        println("i = $i")
        for j = 1:length(Strategy_number_v)
            if i!=j
            Weight_Matrix[i,j] = 1/length(Strategy_number_v)
            else
                Weight_Matrix[i,j] =0
            end
        end    
    end

    Big_Markov = Big_Markov.*Weight_Matrix

    for i = 1:length(Strategy_number_v)
        Big_Markov[i,i] = 1 - sum(Big_Markov[i,:])
    end

    return Big_Markov
end


generate_P_and_EB_matrices = function(Strategy_number_v ,PM,delta,tolerance,max_n_markov)
    EB11 = zeros(1,length(Strategy_number_v));
    P11 = zeros(1,length(Strategy_number_v));
    EB12 = zeros(length(Strategy_number_v),length(Strategy_number_v));
    P12 = zeros(length(Strategy_number_v),length(Strategy_number_v));

    println("building single")
    Threads.@threads    for i = 1:length(Strategy_number_v)
        println(i)
        strategy1 = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
        EB11[i],P11[i],x =get_ebs_and_ps(strategy1,strategy1,delta,PM,tolerance) #--> do this only for one pair!! and create big lookup table
        if mod(i,1000)==0
            println(i/length(Strategy_number_v))
        end
    end


    for i = 1:(length(Strategy_number_v)-1)
        strategy1 = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
        Threads.@threads for j = i+1:length(Strategy_number_v)
            strategy2 = p_adic_to_max_n_state_strategy(Strategy_number_v[j],max_n_markov)
            eb12,p12,p21 =get_ebs_and_ps(strategy1,strategy2,delta,PM,tolerance) #--> do this only for one pair!! and create big lookup table
            EB12[i,j] = eb12;
            EB12[j,i] = eb12;
            P12[i,j] = p12;
            P12[j,i] = p21;
        end
        if mod(i,10)==0
            println(i/length(Strategy_number_v))
        end
    end
    return EB11,EB12,P11,P12 
end


inv_p = function(delta,eb11,eb12,eb22,p11,p12,p21,p22,N)

    p1_v, p2_v = long_pic(delta, eb11, eb12, eb22,p11,p12,p21,p22,N);
    
    payoffs_1 = [p1_v[i+1] for i = 1:N-1];
    payoffs_2 = [p2_v[i+1] for i = 1:N-1];
    payoff_ratios = payoffs_2./ payoffs_1    # = T^-/T^+
    
    den =  1 ; #denominator 

    for j = 1:(N-1)
        den = den+prod(payoff_ratios[1:j])
    end

    inv_prob_ij = 1/den

    ### and now flip it, to only calculate it once!
    payoff_ratios = payoff_ratios.^(-1)
    den =  1 ; 

    for j = 1:(N-1)
        den = den+prod(payoff_ratios[1:j])
    end

    inv_prob_ji = 1/den


    return inv_prob_ij, inv_prob_ji
end



long_pic = function(delta, eb_11, eb_12, eb_22,p11,p12,p21,p22,N)
    x1_v = collect(0:(1/N):1)
    p1_v = zeros(length(x1_v))
    p2_v = zeros(length(x1_v))
    
    x11_pre = 1.0;
    for i = 2:(length(x1_v)-1)
        b11 = (1- ((1-eb_11)*delta))
        b12 = (1- ((1-eb_12)*delta))
        b22 = (1- ((1-eb_22)*delta))
        x1 = x1_v[i]
        function f!(F, x11)
            F[1] = (((x11[1]*b11 + b12* (x1 - x11[1]) )^2)/(x11[1]*b11 + 2*b12*(x1 - x11[1]) + b22*(1+x11[1] - 2*x1)))-(x11[1]*b11)
        end
        x11_sol_1  = nlsolve(f!, [x11_pre])
        x11_sol_2  = nlsolve(f!, [0.0])

        x11_1 = x11_sol_1.zero[1]
        x11_2 = x11_sol_2.zero[1]

        x12_1 = 2*(x1-x11_1)
        x12_2 = 2*(x1-x11_2)

        if x11_1>=0 && x12_1>=0
            x11 = x11_1
            x12 = x12_1
        elseif x11_2>=0 && x12_2>=0
            x11 = x11_2
            x12 = x12_2
        else
            println("NO SOLUTION")
        end
        x11_pre = x11;
        p1_v[i] = (1/x1)*(p11*(x11)+(1/2)*p12*(x12))
        p2_v[i] = (1/(1-x1))*(p22*(1-x11-x12)+(1/2)*p21*(x12))
    end    
    return p1_v, p2_v
end

long_pic_with_stored_frequencies = function(delta, eb_11, eb_12, eb_22,p11,p12,p21,p22,N)
    x1_v = collect(0:(1/N):1)
    p1_v = zeros(length(x1_v))
    p2_v = zeros(length(x1_v))
    FREQS = zeros(3,length(x1_v)) 
    x11_pre = 1.0;
    for i = 2:(length(x1_v)-1)
        b11 = (1- ((1-eb_11)*delta))
        b12 = (1- ((1-eb_12)*delta))
        b22 = (1- ((1-eb_22)*delta))
        x1 = x1_v[i]
        function f!(F, x11)
            F[1] = (((x11[1]*b11 + b12* (x1 - x11[1]) )^2)/(x11[1]*b11 + 2*b12*(x1 - x11[1]) + b22*(1+x11[1] - 2*x1)))-(x11[1]*b11)
        end
        x11_sol_1  = nlsolve(f!, [x11_pre])
        x11_sol_2  = nlsolve(f!, [0.0])

        x11_1 = x11_sol_1.zero[1]
        x11_2 = x11_sol_2.zero[1]

        x12_1 = 2*(x1-x11_1)
        x12_2 = 2*(x1-x11_2)

        if x11_1>=0 && x12_1>=0
            x11 = x11_1
            x12 = x12_1
        elseif x11_2>=0 && x12_2>=0
            x11 = x11_2
            x12 = x12_2
        else
            println("NO SOLUTION")
        end
        x11_pre = x11;
        p1_v[i] = (1/x1)*(p11*(x11)+(1/2)*p12*(x12))
        p2_v[i] = (1/(1-x1))*(p22*(1-x11-x12)+(1/2)*p21*(x12))
        FREQS[1,i] = x11
        FREQS[2,i] = x12
        FREQS[3,i] = 1-x12-x11

    end    
    return p1_v, p2_v, FREQS
end

get_ebs_and_ps= function(strategy1,strategy2,delta,PM,tolerance) #--> do this only for one pair!! and create big lookup table

    # take existing machinery for strategy-pairs that don't break up endogenously
    # so the p1_p2_eb_analysis_pair function only calculates p if eb>0
    println(cut_strategy)
    cut_strategy1,n1  = cut_strategy(strategy1)
    cut_strategy2,n2  = cut_strategy(strategy2)

    p12,p21,eb12 = p1_p2_eb_analysis_pair(cut_strategy1, cut_strategy2,PM,delta)

    if eb12==0
        p12,p21=no_errors_AVERAGE_payoff_strategy1_strategy2(cut_strategy1, cut_strategy2, 0, delta, tolerance, PM)
    end        
    return eb12,p12,p21
end 


no_errors_AVERAGE_payoff_strategy1_strategy2 = function(strategy1, strategy2, error_rate, delta, tolerance, PM)
    cut_strategy1 = strategy1;
    n1 = size(cut_strategy1,2);
    cut_strategy2 = strategy2;
    n2 = size(cut_strategy2,2);

    ## determine how many powers will be calculated
     hp = maximum(PM)
    n_rounds = find_number_periods(hp,delta, tolerance)

    #### set up the transition matrix
        ## start by determining the actions in each state 
        # indices of markov matrix states are 1,2,3,... --> (1,1),(2,1),(3,1),...,(n1,1),(1,2),... 
    actions = actions_vector(n1,n2,cut_strategy1,cut_strategy2)
    actions[1,actions[1,:].==0].=1
    actions[2,actions[2,:].==0].=1
#states wich 0 are never reached (otherwise this functino wouldn't be reached) so i can assign whatever action to these. Zero becomes problematic for the downstream function, so just call it 1
    # --> 
        ## set up a markov matrix from it
    MARKOV = no_errors_full_markov_from_actions_and_strategies(cut_strategy1,cut_strategy2,actions,n1,n2,error_rate)



    ##### Payoffs per state
    payoffs = payoffs_from_actions(actions,PM)
    ###### now calculate distributions and resulting payoffs
    distr = hcat(1,zeros(1,n1*n2-1))
    p1 = 0;
    p2 = 0;
    den = 0;
    for t = 1:n_rounds
        p1 = p1 + delta^(t-1)*(distr*payoffs[1,:])[1]
        p2 = p2 + delta^(t-1)*(distr*payoffs[2,:])[1]
        distr = distr*MARKOV
        den = den + delta^(t-1);
    end    
    p1 = p1/den;
    p2 = p2/den;
    return p1,p2
end



p1_p2_eb_analysis_pair = function(strategy1, strategy2, PM, delta)
    eb = 0;
    p1 = 0;
    p2 = 0;
    ends_in_leave = 0;
    rounds_before_leave = 0;
    state1 = 1;
    state2 =1;
    statecombs = ones(2,1); 
    if any(strategy1[1,:].==0) || any(strategy2[1,:].==0)
        while  ends_in_leave ==0 && (rounds_before_leave == 0 || any((statecombs[1,1:end-1].==state1).*(statecombs[2,1:end-1].==state2))==false) # check if this state combination has never been reached... otherwise we have seen all states
            state1_new = strategy1[1+strategy2[1,state2],state1];
            state2_new = strategy2[1+strategy1[1,state1],state2];
            state1 = state1_new;
            state2 = state2_new;
            rounds_before_leave = rounds_before_leave+1;
            statecombs = hcat(statecombs,[state1, state2]);
            if strategy1[1,state1]==0 || strategy2[1,state2]==0
                ends_in_leave=1;
            end
        end
    end



    actions_until_leave = hcat(strategy1[1,convert(Array{Int64},statecombs[1,1:end-1])],strategy2[1,convert(Array{Int64},statecombs[2,1:end-1])])'; #resident 1,:  entrant 2,:

    if ends_in_leave==1
        den = 1;
        p1 = PM[actions_until_leave[1,1],actions_until_leave[2,1]];
        p2 = PM[actions_until_leave[2,1],actions_until_leave[1,1]];
        for i = 2:rounds_before_leave
            den = den+delta^(i-1);
            p1 = p1 + delta^(i-1)*PM[actions_until_leave[1,i],actions_until_leave[2,i]];
            p2 = p2 + delta^(i-1)*PM[actions_until_leave[2,i],actions_until_leave[1,i]];
        end
        eb = 1/den
        p1 = p1/den
        p2 = p2/den
    end


    return p1,p2,eb
end 



no_errors_deterministic_markov_from_actions_and_strategies = function(strategy1,strategy2,actions,n1,n2)
    M =zeros(size(actions,2),size(actions,2))
    i = 1
    for s2 = 1:n2
        for s1 = 1:n1
            s1_new = strategy1[1+actions[2,i],s1]#state that 1 moves to
            s2_new = strategy2[1+actions[1,i],s2]#state that 2 moves to 
            M[i, n1*(s2_new-1)+s1_new] = 1
            i = i+1
        end
    end

    return M
end


no_errors_full_markov_from_actions_and_strategies = function(strategy1,strategy2,actions,n1,n2,error_rate)
    M00 = no_errors_deterministic_markov_from_actions_and_strategies(strategy1,strategy2,actions,n1,n2)
    return M00
end



leave_markov_fill = function(Strategy_number_v,N,delta, EB11,EB12,P11,P12)
    Big_Markov = zeros(length(Strategy_number_v),length(Strategy_number_v))

    for i = 1:(length(Strategy_number_v)-1)
        println("i = $i")
        eb11 = EB11[i] 
        p11 = P11[i]
        Threads.@threads for j = (i+1):length(Strategy_number_v)
            if j!=i
                eb12 = EB12[i,j]
                eb22 = EB11[j]
                p12 = P12[i,j] 
                p21 = P12[j,i] 
                p22 = P11[j]   
                Big_Markov[j,i] , Big_Markov[i,j] = inv_p(delta,eb11,eb12,eb22,p11,p12,p21,p22,N)  #this is still a bit excessive on passing objects
            end
        end    
    end

    for i = 1:length(Strategy_number_v)
        Big_Markov[i,i] = 1 - sum(Big_Markov[i,:])
    end

    return Big_Markov
end
