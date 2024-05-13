max_n = 10
base_s5tn_v1 = 3;
s5tn_v1 = base_s5tn_v1.^collect(0:(max_n-1)); #for the first line of a strategy
s5tn_v2 = (max_n+1).^collect(0:(max_n-1)); # for all other lines of a strategy 
s5tn_v1_UInt = convert(Vector{UInt128}, s5tn_v1);#for the first line of a strategy
s5tn_v2_UInt = convert(Vector{UInt128}, s5tn_v2); # for all other lines of a strategy 




get_n = function(strategy) # returns the number of states of a strategy 
    if strategy[2,end]!=0 
        n = size(strategy,2);
    else
        n =-1+ findfirst( x -> strategy[2,x] ==0, 1:size(strategy,2));
    end
    return n
end

cut_strategy = function(strategy) #reduces the strategies from a 3 (5) by max_n matrix to a 3 (5) by n matrix, (n = number of states) 
    n = get_n(strategy)
    red_strategy = strategy[:,1:n]
    return red_strategy,n
end

convert_strat_5_to_3 = function(strategy5)
    strategy3 = zeros(Int64,3,size(strategy5,2))
    strategy3[1,:] = strategy5[1,:]
    for i = 1:size(strategy5,2)
        if strategy5[1,i]==1
            strategy3[2,i] = strategy5[2,i]
            strategy3[3,i] = strategy5[3,i]
        elseif strategy5[1,i]==2 
            strategy3[2,i] = strategy5[4,i] 
            strategy3[3,i] = strategy5[5,i]
        else  ## in case it's a leave --> conversion to 3strat isn't clear --> adapt for specific use-case 
        end
    end
    return strategy3
end



number_to_strat_3_basic = function(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt) # converts from triple of numbers to matrix
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





number_to_strat_5_basic = function(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt) # converts from quintuple of numbers to matrix
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



strat_3_to_number = function(strat_,max_n,s5tn_v1,s5tn_v2) # converts from matrix to triple of numbers 
    strat = copy(strat_)
    s = zeros(UInt64, 3);
    s[1]=strat[1,:]'*s5tn_v1
    s[2]=strat[2,:]'*s5tn_v2
    s[3]=strat[3,:]'*s5tn_v2
    return s
end



strat_5_to_number = function(strat_,max_n,s5tn_v1,s5tn_v2) # converts from matrix to quintuple of numbers 
    strat = copy(strat_)
    s = zeros(UInt64, 5);
    s[1]=strat[1,:]'*s5tn_v1
    s[2]=strat[2,:]'*s5tn_v2
    s[3]=strat[3,:]'*s5tn_v2
    s[4]=strat[4,:]'*s5tn_v2
    s[5]=strat[5,:]'*s5tn_v2
    return s
end