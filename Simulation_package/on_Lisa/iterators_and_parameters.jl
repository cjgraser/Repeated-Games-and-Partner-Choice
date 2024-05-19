imax_n = 10;istore_interval = 1000;inperiods = 100000000;iN = 100; imutation_rate = 0.05;
inperiods = 10000
a = 0; c = 1; gamma = 0.9502; eps = 0; w =0.9;
ntime = convert(Int128,inperiods/istore_interval);
#picture of: beta, and of b
bv = collect(1:(1/2):6);
betav = collect(0.99:-0.01:0.09);

include(a5*"all_functions_5d.jl"); #this takes shortly compared to the whole, so doesn't matter to import multiple times

execute_in_function = function(bi, betai, no_leaves, trimming)
    b = bv[bi];
    beta = betav[betai]
    ipayoff2b2 =zeros(2, 2); ipayoff2b2[1,:] = [b-c -c]; ipayoff2b2[2,:] = [b 0]; ipayoff2b2 = ipayoff2b2.+(c+a);ipayoff2b2 =ipayoff2b2./maximum(ipayoff2b2);
    match_types, population_states, errors, leaves = simulation(ipayoff2b2,imax_n,istore_interval,inperiods,iN,imutation_rate,eps,gamma,beta,w,no_leaves,trimming);
    return match_types, population_states, errors, leaves
end









