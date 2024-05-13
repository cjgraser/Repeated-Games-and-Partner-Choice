#import Pkg;
#Pkg.add("Distributions")
#Pkg.add("Random")
#Pkg.add("StatsBase")
#Pkg.add("DataFrames")
#Pkg.add("FreqTables")
using Distributions
using Random
using StatsBase
using DataFrames
using FreqTables


## in the current state, what are the intended actions, i.e. the actions if no error happens --> error happens on C-D, not on leave  
intended_actions(PS,state_v,index_shifts) = PS[state_v+index_shifts]; ## PS = Population[:,1,:]' #THE TRANSPOSE IS IMPORTANT


## gather all actions for the entire population, with the possibility of errors
pop_actions = function(N,epsilon,Pop,state_v,index_shifts)
    i_actions = intended_actions(Pop[:,1,:]', state_v,index_shifts); #which is a vector of 0,1,2
    play_errors = rand(N).<epsilon; 
    # errors can only act on the fields that are 1 or 2
    actions = i_actions  + ((i_actions.==1)-(i_actions.==2)).*play_errors;
    return actions
end



## gather all payoffs for the population 
pop_payoffs  = function(Pairs,actions,N,payoff2b2,w)
    payoffs = zeros(1,N);
    for i in 1:size(Pairs,2)
        payoffs[Pairs[1,i]] = payoff2b2[actions[Pairs[1,i]], actions[Pairs[2,i]]];
        payoffs[Pairs[2,i]] = payoff2b2[actions[Pairs[2,i]], actions[Pairs[1,i]]];
    end
    payoffs = (1-w).*ones(1,N) + (w).*payoffs;
    return payoffs
end


## randomly break up some pairs
pop_breaks = function(Pairs,beta)
    r = rand(size(Pairs,2));
    Pairs = Pairs[:, (r.<(beta^2))]; # the events are i.i.d.
    return Pairs
end

## given the actions that were played, perform the state transitions and gather for whole population 
state_transitions_pop5 = function(actions,Pairs, state_v_pre, Population) #use this only on those pairs that didn't break, so actions(pairs) only contain 1's and 2's
    state_v = copy(state_v_pre)
    for i in 1:size(Pairs,2)
    state_v[Pairs[1,i]] = Population[Pairs[1,i],2*(actions[Pairs[1,i]]-1)+actions[Pairs[2,i]]+1,state_v_pre[Pairs[1,i]]];
    state_v[Pairs[2,i]] = Population[Pairs[2,i],2*(actions[Pairs[2,i]]-1)+actions[Pairs[1,i]]+1,state_v_pre[Pairs[2,i]]];
    end
    return state_v
end

## break up pairs in which the current action for at least one partner is 'leave'
break_from_leave  = function(Pairs,state_v,Population,index_shifts)
    ## part that deletes the pairs that have a leaving individual in them
    IA = intended_actions(Population[:,1,:]',state_v,index_shifts);
    t_pre = IA[Pairs];
    t_pre = vec(sum((t_pre.==0),dims=1));
    Pairs = Pairs[:,findall(t_pre.==0)];
return Pairs
end

## use random number in normalized, stacked vector to map draw from a uniform [0,1] to an event, and return a random number
## where in a given interval a randomly drawn number lies is uncorrelated with whether or not it lies in this interval. So this exact location can be used as a new random number.
## (the rationale for doing this, is that drawing random numbers takes time, and that we have few events 
## -- all with probabilities that are large compared to the precision of the random number, so re-using the random number is computationally cheaper, and inconsequential for the result)
chooseV_recycle = function(v,r) #v must already represent a cdf, with first element 0, and last 1.
    index =  findfirst( x -> v[x] > r, 1:length(v));
    r_new = (r-v[index-1])/v[index];
    index = index-1;
     return index, r_new
end

## same without returning a random number
chooseV = function(v,r) #v must already represent a cdf
   index =  findfirst( x -> v[x] > r, 1:length(v));
    return index
end


#perform a single  mutation step on a strategy 
mutation_basic_5d = function(mtd,rn,strategy,n) # 5d indicates that the strategy is 5 by max_n -- there is an equivalent version for 3 by max_n matrices

    if rn < mtd[1] #arrow mutation (as described in thesis)
        if n!=1
            rn = rn/mtd[1];
            # which arrow mutates is determined by some chooseV
            # chosen u.a.r. index of mutation (idm) (this is a linear index --> below, transform this linear index into cartesian)
            idm,rn=chooseV_recycle(collect(0:(1/(4*n)):1),rn); 
            idmrow = mod(idm,4)+2; # first part takes values 0,1,2,3 --> map to 2,3,4,5
            idmcol = convert(Int64,floor((idm-1)/4)+1);  # like this, it takes values 1,...,n //  example n=2 --> idm = 1,2,3,4 map to 1 --> 5,6,7,8 map to 2 
            push = convert(Int64,ceil((n-1)*rn)); # we are in Z modulo n. we move between 1 and n-1 steps --> that puts us to a proper different spot in Z modulo n 
            # strategy[idmrow,idmcol]+push-1 --> subtract one to move from 1-N to Z modulo n, and then add one to map back 
            strategy[idmrow,idmcol] = 1+mod(strategy[idmrow,idmcol]+push-1,n);# could be +1, +2,... modn (can only point to existing states)
        end
    elseif rn<mtd[2] #state-output mutation (as described in thesis)
        rn = (rn-mtd[1])/(mtd[2]-mtd[1]);
        idm,rn=chooseV_recycle(collect(0:(1/(n)):1),rn) #which state's action is chosen 
        if idm !=1 #only really go through with this, if the mutation doesn't change the first state to a zero
        strategy[1,idm] = mod(strategy[1,idm]+ceil(2*rn),3) #could be 0,1,2 --except the one that already is specified there // all permitted
        else
        strategy[1,idm] = mod(strategy[1,idm],2)+1; #because 0 (leave) is not permitted in the first state. First, map back to {0,1} then add one, modulo it by 2, then map back to {1,2}, so it becomes  mod(strategy[1,idm],2)+1-1, which cancels
        end
    elseif rn< mtd[3] #state-deletion  (as described in thesis)
        if n>1
        rn = (rn-mtd[2])/(mtd[3]-mtd[2]);
        #let's make it simple and make sure that state 1 is not deleted
        idm,rn=chooseV_recycle(collect(0:(1/(n-1)):1),rn);#which state is chosen to be deleted
        idm = idm+1 #because we excluded the first state
        # redirect all pointers randomly 
        point_to_v = union(1:(idm-1),(idm+1):n); #possible states to point to 
        for i in union(1:(idm-1),(idm+1):n)
            for j in 2:5
                if strategy[j,i]==idm
                        strategy[j,i] = rand(point_to_v); 
                end
            end
        end         
        #now the idm'th column is deleted, and all numbers in strategy[2:end,:]>idm are reduced by one
            if idm==size(strategy,2)
                strategy[:,idm].=0;
            else
                strategy[:,idm:end-1] = strategy[:,idm+1:end];                
            end
            strategy[2:end,:]  =strategy[2:end,:].- (strategy[2:end,:].>idm);
        end
    else #state-addition  (as described in thesis)
        rn = (rn-mtd[3])/(1- mtd[3])
        #add state: this is the most involved one
        n = n+1;
        #need to draw 4 numbers 
        index1,rn = chooseV_recycle(collect(0:1/n:1),rn);
        index2,rn = chooseV_recycle(collect(0:1/n:1),rn);
        index3,rn = chooseV_recycle(collect(0:1/n:1),rn);
        index4,rn = chooseV_recycle(collect(0:1/n:1),rn);
        strategy[2,n]=index1;
        strategy[3,n]=index2;
        strategy[4,n]=index3;
        strategy[5,n]=index4; 
        strategy[1,n] = mod(ceil(3*rn),3); #lastly we draw an action u.a.r. 
        # now we have wiring from the state. but we also need wiring to the state
        # we randomly draw one element that links there
        # we're still exploiting the same random number
        rn = mod(rn,1/3);
        idm = chooseV(0:1/(4*(n-1)):1,rn);
        idmrow = mod(idm,4)+2; # first part takes values 0,1,2,3 --> map to 2,3,4,5
        idmcol = convert(Int64,floor((idm-1)/4)+1);  # like this, it takes values 1,...,n //  example n=2 --> idm = 1,2,3,4 map to 1 --> 5,6,7,8 map to 2
        strategy[idmrow,idmcol] = n;    
    end
    return strategy
end



## perform a selection step 
Moran_selection =  function(payoffs,N, RI,Pairs) #RI the number of replacing individuals
    # returns the vector of indices that are new, and a reference to what strategy each of these new individuals play, in reference to the old population
    payoffs = payoffs/sum(payoffs);
    for i in 2:N
        payoffs[i] = payoffs[i]+ payoffs[i-1];  
    end 
    nr  = rand(RI);#new randowm draw
    repstr = zeros(Int64,RI);
    for i in 1:RI
        repstr[i] = findfirst( x -> payoffs[x] > nr[i], 1:N); ##
    end
    newinds = vec(Pairs[:, shuffle(1:div(N,2))[1:div(RI,2)]]); #the ones that died for this
    return newinds, repstr
end 

## perform selection step and gather results and implications for pairs
pop_selection = function(payoffs,Pairs,Population,state_v, N, RI)
    ## new_inds are the indivuals that are being replaced
    ## repstr are the reproducing strategies --> so in the population object the new_inds will have their strategy updated to the strategies in the repstr object
    new_inds, repstr  = Moran_selection(payoffs,N, RI,Pairs); #the new individuals
    #states are set to one
    state_v[new_inds] .= 1;
    Population[new_inds,:,:] = Population[repstr,:,:]; 
    ## part that deletes the pairs that have a new individual in them
    a = in(new_inds).(Pairs[1,:]);
    b = in(new_inds).(Pairs[2,:]);
    changed = (a.+b);
    Pairs = Pairs[:,(changed.==0)]
    return Pairs, Population, state_v 
end


## perform rematching, and set states of rematched ones to one
pop_matching = function(Pairs,otN, state_v)  
    singles = setdiff(otN,vec(Pairs));
    state_v[singles] .= 1;
    singles = shuffle(singles);
    Pairs = hcat(Pairs,reshape(singles, 2,:));
    return Pairs, state_v
end


## this function does a full simulation
pop_simulation = function(Population,Strategy_stored_Blank, index_shifts, Pairs, otN, beta, gamma, N,  epsilon,payoff2b2, nperiods, store_interval, mutation_rate, match_types, population_states,max_n,leaves,w,no_leave,trimming)
    # define some very frequently used objects (most are the same as in the data_types_back_and_forths.jl file)
    base_s5tn_v1 = 3;
    s5tn_v1 = base_s5tn_v1.^collect(0:(max_n-1)); 
    s5tn_v2 = (max_n+1).^collect(0:(max_n-1)); # because it could be empty (zero)

    #initialize the quantities that are kept track of during the simulation 
    

    state_v = ones(Int64,N); #what states the automata are in 
    Pairs_Pre = Pairs; #the matching
    actions = zeros(Int64,N); #the actions that are played in the current round 
    errors = zeros(convert(Int64,(nperiods/store_interval))); # this is simply a tool to keep track of whether max_n is chosen appropriately --> if strategies mutate to something larger than max_n the mutation is aborted and this error variable is increased by one

    # we take a snapshot of the population every store_interval stage-games. all steps in between two snapshots are in the "store_interval_runs()" function
        for t in 1:convert(Int64,(nperiods/store_interval))
            Pairs_Pre, Population, Pairs, state_v, actions, errors_i, leave_count =store_interval_runs(store_interval, Population, Pairs, state_v, otN, beta, gamma, N, max_n , epsilon,payoff2b2, index_shifts, mutation_rate,errors,w,no_leave,trimming)
            errors[t] =errors_i;
            match_types[:,t] = pair_types(Pairs_Pre,actions); #because actions are recoreded at the beginning of period, and not later 
            population_states[t] = store_groupcounts5(Population,max_n,N,Strategy_stored_Blank,base_s5tn_v1,s5tn_v1,s5tn_v2);
            leaves[t] = leave_count;
            if t%1000 ==0
                println(t/(nperiods/store_interval))
            end
        end
        errors = sum(errors)

        # the leaves obect counts how many individuals play leave at each snapshot 
    return match_types, population_states, errors, leaves#... all the metrics that are stored, over long time
end


## inner loop of the above function --> runs the one_period function store_interval times (one_period is one stage game) 
store_interval_runs = function(store_interval, Population, Pairs, state_v, otN, beta, gamma, N, max_n , epsilon,payoff2b2, index_shifts, mutation_rate,errors,w,no_leave,trimming)
    Pairs_Pre = Pairs
    actions = zeros(Int64,N)
    leave_count = 0
    errors = zeros(convert(Int64,store_interval));
        for st in 1:store_interval
            Pairs_Pre = Pairs;
            Population, Pairs, state_v, actions_i, errors_i, leaves_i =one_period(Population, Pairs, state_v, otN, beta, gamma, N, max_n , epsilon,payoff2b2, index_shifts, mutation_rate,w,no_leave,trimming);
            errors[st] = errors_i;
            actions = actions_i;
            leave_count = leaves_i;
        end
        errors = sum(errors)
        return Pairs_Pre, Population, Pairs, state_v, actions, errors, leave_count
end    
    


## all things that happen in one stage game and before the next stage game
one_period = function(Population, Pairs, state_v, otN, beta, gamma, N, max_n , epsilon,payoff2b2, index_shifts, mutation_rate,w,no_leave,trimming)
    count_errors = 0;
    #1) what actions are played
    actions = pop_actions(N,epsilon,Population,state_v,index_shifts);
    #2) what payoffs result
    payoffs = pop_payoffs(Pairs,actions,N,payoff2b2,w);
    state_v  = state_transitions_pop5(actions,Pairs,state_v,Population);
    #3) who reproduces and who gets replaced  
    RI = rand(Binomial(N,(1-gamma)/2));#draw how many individuals are replaces --> divided by two, because the pop_selection function replaces both individuals in a pair at the same time, to minimize necessary break-ups per selection event
    if RI>0
        Pairs,Population,state_v= pop_selection(payoffs,Pairs,Population,state_v, N, 2*RI);
    end
    #4) leave-events -- if no_leave==1, nothing happens, because no state-output is zero
    count_pre = size(Pairs,2);
    Pairs = break_from_leave(Pairs,state_v,Population,index_shifts);
    leaves_i = count_pre - size(Pairs,2);

    #5) exogenous break-ups
    Pairs = pop_breaks(Pairs,beta); #the state_v of the broken ones is not updated yet, but will be in the matching function

    #6) mutations #mutation_rate^2 is so close to zero that at max one mutation happens
    rn  = rand(1)[1];#draw a random number
    if rn<mutation_rate
        rn = rn/mutation_rate; 
        mut_ind = convert(Int64, ceil(rn*N)); #randomly choose an individual that mutates
        
        # n = length of this strategy
        if Population[mut_ind,2,end]!=0 
            n = max_n;
        else
            n =-1+ findfirst( x -> Population[mut_ind,2,x] ==0, 1:max_n);
        end


        # weight_increase = likelihood(adding state)/(likelihood(adding state)+likelihood(deleting state))
        # this is set to zero if n==max_n
        weight_increase = 0.45;

        if n == max_n
            count_errors = 1;
            weight_increase= 0;
        end
        
        decay = 0.7; #number of mutations follows an exponential distribution with this decay parameter
        new_strategy = exp_decay_mutation!(rn,decay,Population[mut_ind,:,:],weight_increase,max_n,no_leave);


        # the different versions of the mutation kernel 
        if trimming == 1
            Population[mut_ind,:,:]  = reduce_to_connected_graph(new_strategy);
        elseif trimming ==0
            Population[mut_ind,:,:]  = random_reconnect_graph(new_strategy);
        else
            Population[mut_ind,:,:]  = new_strategy;
        end
        # mut_ind is eliminated from Pairs, if it is in a pair
        Pairs = Pairs[:,Pairs[1,:].!= mut_ind];
        Pairs = Pairs[:,Pairs[2,:].!= mut_ind];
    end

    #5) re-matching 
    Pairs,state_v = pop_matching(Pairs,otN,state_v);

    return Population, Pairs, state_v, actions, count_errors, leaves_i
end



## exponential decay mutation 
exp_decay_mutation! = function(rn, decay, strategy,weight_increase,max_n,no_leave)
     for i = 1:draw_exponential_decay(decay)
        rn = rand(1)[1]
        if strategy[2,end]!=0 
            n = max_n;
        else
            n =-1+ findfirst( x -> strategy[2,x] ==0, 1:max_n);
        end
        # THESE ARE THE WEIGHTS OF THE DIFFERENT TYPES OF MUTATIONS

        mtd = [0.2 0.4 0.4+(0.6*(1-weight_increase)) 1];
        
        
        if n == max_n
            mtd  = [0.2 0.4 1]
        end

        if no_leave == 0
            strategy = mutation_basic_5d(mtd,rn,strategy,n);
        else
            strategy = mutation_basic_5d_no_leave(mtd,rn,strategy,n);
        end
    end
    return strategy
end


draw_exponential_decay = function(decay)
    r = rand(1)[1]
    i = 1;
    while r>1-((decay)^i)
        i = i+1;
    end
    return i
end

store_groupcounts5 = function(population,max_n,N,Strategy_stored_Blank,base_s5tn_v1,s5tn_v1,s5tn_v2)
    # groups all individuals with identical strategy to reduce required storage space 
    Strategy_stored = deepcopy(Strategy_stored_Blank);
    for i in 1:N
        Strategy_stored[i,:]=strat_5_to_number(population[i,:,:],max_n,s5tn_v1,s5tn_v2);
    end
    return combine(groupby(Strategy_stored, [:x1, :x2, :x3, :x4, :x5]), nrow)
end


pair_types =function(Pairs,actions)
    # how many pairs play CC, CD or DD
    a_i(i) = actions[i];
    t = zeros(3);
    t_pre=a_i.(Pairs);
    t_pre = sum(t_pre,dims=1);
    t[1] = sum(t_pre.==2);
    t[2] = sum(t_pre.==3);
    t[3] = sum(t_pre.==4);
    return t
end 

strat_5_to_number = function(strat_,max_n,s5tn_v1,s5tn_v2)
# as in data_types_back_and_forths
    strat = deepcopy(strat_) 
    s = zeros(UInt64, 5);
    s[1]=strat[1,:]'*s5tn_v1
    s[2]=strat[2,:]'*s5tn_v2
    s[3]=strat[3,:]'*s5tn_v2
    s[4]=strat[4,:]'*s5tn_v2
    s[5]=strat[5,:]'*s5tn_v2
    return s
end


simulation = function(payoff2b2,max_n,store_interval,nperiods,N,mutation_rate,epsilon,gamma,beta,w,no_leave,trimming)
    # This is the wrapper for the pop_simulation function 
    # all the things that have to be declared once 
    Strategy_stored_Blank = DataFrame(zeros(Int64,N,5),:auto);
    otN = collect(1:N);
    index_shifts = [ (x-1)*max_n for x in 1:N ]; 
    
    # here things that are stored 
    ntime = convert(Int128,nperiods/store_interval)
    match_types = zeros(3,ntime);
    population_states = Array{Any}(undef,ntime);
    leaves = zeros(ntime);

    #initialize
    Population = zeros(Int64,N,5,max_n); #start with allC
    Population[:,:,1].= 1; 
    Pairs = reshape(collect(1:N), (2,Int(N/2)));
    match_types, population_states, errors , leaves= pop_simulation(Population,Strategy_stored_Blank, index_shifts, Pairs, otN, beta, gamma, N, epsilon,payoff2b2, nperiods, store_interval, mutation_rate, match_types, population_states,max_n, leaves,w,no_leave,trimming);
        return match_types, population_states, errors, leaves
end






reduce_to_connected_graph = function(strategy)     # delets states in an automaton that cannot be reached 
    reachable ,n = identify_disconnections_graph(strategy);
    # now, carefully delete all columns that are disconnected, and change indices accordingly (nothing has pointed to these columns previously)
    if (length(reachable))<n
        #indices that need to be deleted 
        d_indices = setdiff(collect(1:n), reachable)
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
    return (strategy)
end



identify_disconnections_graph  = function(strategy) #identifies components of the graph
    #works on 3 and 5 strategies
        if strategy[2,end]!=0 
            n = size(strategy,2);
        else
            n =-1+ findfirst( x -> strategy[2,x] ==0, 1:size(strategy,2)); # n states are defined 
        end
        reachable = unique(union(1,unique(strategy[2:end,1]))) #those are reachable from state 1
        i=1; #the index in the reachable_stay vector 
        while (length(reachable))<n && (i<length(reachable))  #if we haven't identified all states and if and we still have reachable states to consider
            i = i+1;
            a = strategy[2:end,reachable[i]]; 
            ap = setdiff(a, reachable);
            reachable = append!(reachable,ap);
        end
        return reachable ,n
end



random_reconnect_graph  = function(strategy)
    reachable ,n = identify_disconnections_graph(strategy);
    strategy  = random_reconnect_graph_loop(strategy, reachable, n);
    return (strategy)
end

random_reconnect_graph_loop = function(strategy, reachable, n)
    while length(reachable)<n
        strategy,reachable = random_reconnection_step(strategy, reachable, n)
    end
    return strategy
end

random_reconnection_step = function(strategy, reachable, n)
    # we have an unconnected graph (otherwise the function wouldn't be called)
    d_indices = setdiff(collect(1:n), reachable)
    d_indices = shuffle(d_indices)
    k = length(reachable)
    strategy_new = deepcopy(strategy)
    arrow = shuffle(collect(1:((size(strategy,1)-1)*k)))[1] #choose a random arrow to point there
    @views strategy_new[2:end, reachable][arrow] = d_indices[1]; #the pointing has happened
    reachable_new ,n = identify_disconnections_graph(strategy_new)
    if length(setdiff(reachable,reachable_new))==0 #if there is nothing in reachable that isn't in reachable_new
        return strategy_new, reachable_new
    else
        return strategy, reachable
    end
end



mutation_basic_5d_no_leave = function(mtd,rn,strategy,n) #same as the one with leaving, just that no state output can be changed into a 0
    if rn < mtd[1] #arrow mutation
        if n!=1
            rn = rn/mtd[1];
            idm,rn=chooseV_recycle(collect(0:(1/(4*n)):1),rn); 
            idmrow = mod(idm,4)+2; 
            idmcol = convert(Int64,floor((idm-1)/4)+1);
            push = convert(Int64,ceil((n-1)*rn)); 
            strategy[idmrow,idmcol] = 1+mod(strategy[idmrow,idmcol]+push-1,n);
        end
    elseif rn<mtd[2]
        rn = (rn-mtd[1])/(mtd[2]-mtd[1]);
        idm,rn=chooseV_recycle(collect(0:(1/(n)):1),rn)  
        #####################version without leaves ##################
        strategy[1,idm] = mod(strategy[1,idm],2)+1;
    elseif rn< mtd[3]
        if n>1
        rn = (rn-mtd[2])/(mtd[3]-mtd[2]);
        idm,rn=chooseV_recycle(collect(0:(1/(n-1)):1),rn);
        idm = idm+1 

        point_to_v = union(1:(idm-1),(idm+1):n);  
        for i in union(1:(idm-1),(idm+1):n)
            for j in 2:5
                if strategy[j,i]==idm
                        strategy[j,i] = rand(point_to_v); 
                end
            end
        end    
            if idm==size(strategy,2)
                strategy[:,idm].=0;
            else
                strategy[:,idm:end-1] = strategy[:,idm+1:end];                
            end
            strategy[2:end,:]  =strategy[2:end,:].- (strategy[2:end,:].>idm);
        end
    else
        rn = (rn-mtd[3])/(1- mtd[3])
        n = n+1;
        index1,rn = chooseV_recycle(collect(0:1/n:1),rn);
        index2,rn = chooseV_recycle(collect(0:1/n:1),rn);
        index3,rn = chooseV_recycle(collect(0:1/n:1),rn);
        index4,rn = chooseV_recycle(collect(0:1/n:1),rn);
        strategy[2,n]=index1;
        strategy[3,n]=index2;
        strategy[4,n]=index3;
        strategy[5,n]=index4; 
        ################ adapted in no-leave-simulation###############
        strategy[1,n] =1+ mod(ceil(2*rn),2); 
        rn = mod(rn,1/3);
        idm = chooseV(0:1/(4*(n-1)):1,rn);
        idmrow = mod(idm,4)+2;
        idmcol = convert(Int64,floor((idm-1)/4)+1);
        strategy[idmrow,idmcol] = n;    
    end
    return strategy
end

strat_3_to_number = function(strat_,max_n,s5tn_v1,s5tn_v2)
# same as in data_types_back_and_forths.jl file
    strat = copy(strat_)
    s = zeros(UInt64, 3);
    s[1]=strat[1,:]'*s5tn_v1
    s[2]=strat[2,:]'*s5tn_v2
    s[3]=strat[3,:]'*s5tn_v2
    return s
end