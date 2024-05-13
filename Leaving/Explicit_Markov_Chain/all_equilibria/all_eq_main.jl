## until length 5 

max_n_markov = 4

Strategy_number_v =collect(0:((max_n_markov)^(3*max_n_markov)-1));

include("../general_markov.jl")


Strategy_number_v = clean_to_unique(Strategy_number_v,max_n_markov)
Strategy_number_v = clean_to_permitted(Strategy_number_v,max_n_markov)

Strategy_number_v = load_object("Strategy_number_v_unique_4.jld2")

include("../../../Eq_Check_Module/new_best_responder_functions.jl")
include("../functions_to_fill_markov_matrix.jl")

Leaving_Payoff_gap_eq_v = zeros(length(Strategy_number_v))
Leaving_eq_v= zeros(length(Strategy_number_v))


minx1_initial = 0.001
maxx1_initial = 0.999
gridsteps = 100
tolerance_p = 10^(-7)
tolerance = 10^(-6)
be_p_v = collect(0.1:0.1:0.9)
delta = 0.9
PM = [2 0; 3 1]
leaving = 1

for i = 1:length(Strategy_number_v)
    println(i)
    stra = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
    stra = reduce_to_connected(stra)
    stra,x = cut_strategy(stra)
    pop = [stra]
    Leaving_eq_v[i],Leaving_Payoff_gap_eq_v[i],BE_dump = payoff_gap_eq_assessment_3_version(pop,PM,delta,tolerance,tolerance_p,be_p_v,leaving,minx1_initial,maxx1_initial,gridsteps)
end



no_leave_indices = zeros(length(Strategy_number_v))
No_Leaving_Payoff_gap_eq_v = zeros(length(Strategy_number_v))
No_Leaving_eq_v= zeros(length(Strategy_number_v))


leaving =0 
for i = 1:length(Strategy_number_v)
    println(i)
    stra = p_adic_to_max_n_state_strategy(Strategy_number_v[i],max_n_markov)
    stra = reduce_to_connected(stra)
    stra,x = cut_strategy(stra)
    if all(stra[1,:].!=0)
    no_leave_indices[i]=1
    pop = [stra]
    No_Leaving_eq_v[i],No_Leaving_Payoff_gap_eq_v[i],BE_dump = payoff_gap_eq_assessment_3_version(pop,PM,delta,tolerance,tolerance_p,be_p_v,leaving,minx1_initial,maxx1_initial,gridsteps)
    end
end


S_eq_No_L = Strategy_number_v[No_Leaving_eq_v.==1]
stra = p_adic_to_max_n_state_strategy(Strategy_number_v[end-1],max_n_markov)


remove_equivalent_cycles(p_adic_to_max_n_state_strategy(S_eq_No_L[end],max_n_markov)   ) #takes in an already cleaned strategy 
remove_equivalent_cycles(p_adic_to_max_n_state_strategy(S_eq_No_L[end-1],max_n_markov))    #takes in an already cleaned strategy 


S_eq_No_L = Strategy_number_v[No_Leaving_eq_v.==1]
S_eq_L = Strategy_number_v[Leaving_eq_v.==1]


all_eq_to_txt = function(storepath,S_eq,max_n_markov)
    S_eq, Lettervec = sort_the_Seq(S_eq,max_n_markov)

    open(storepath, "w") do file
    write(file,"ALL EQUILIBRIA WITH UP TO 4 STATES \n")

    write(file,"\n")
    lv = unique(Lettervec)
    write(file," $lv")
    write(file,"\n")
    write(file,"Strategies \n")

    for i = 1:length(S_eq)
        write(file,"strategy $i \n")
        l = Lettervec[i]
        write(file,"$l \n")
    stra = p_adic_to_max_n_state_strategy(S_eq[i],max_n_markov)
        if stra[1,1]!=0 
        stra1 = stra[1,:]
        stra2 = stra[2,:]
        stra3 = stra[3,:]   
            write(file, "$stra1 Strategy\n")
            write(file, "$stra2 C\n")
            write(file, "$stra3 D\n")
            write(file, "\n")
        end
    end
    end
end



sort_the_Seq = function(Seq,max_n_markov)

    Lettervec =  fill("0", length(Seq))
    
    for stra_i = 1:length(Seq)
        #println(stra_i)
        stra = p_adic_to_max_n_state_strategy(Seq[stra_i],max_n_markov)
        stateseq = [1]
        state = 1
        new_state = stra[1+stra[1,state],state]

        while all(stateseq.!=new_state)
            state = new_state
            stateseq = append!(stateseq,state)
            if stra[1,state]==0
                new_state = state    
            else
                new_state = stra[1+stra[1,state],state]    
            end
        end

        actions = stra[1,stateseq]
        cyclefrom = findfirst(stateseq.==new_state)
        in_letters = []
        for i = 1:(cyclefrom-1)
            if actions[i]==1
                in_letters = append!(in_letters,"c")
            elseif actions[i]==2
                in_letters = append!(in_letters,"d")
            elseif actions[i]==0
                in_letters = append!(in_letters,"L")
            end
        end
        for i = cyclefrom:length(actions)
            if actions[i]==1
                in_letters = append!(in_letters,"C")
            elseif actions[i]==2
                in_letters = append!(in_letters,"D")
            elseif actions[i]==0
                in_letters = append!(in_letters,"L")
            end
        end
        Lettervec[stra_i]= join(in_letters)
    end

    p = sortperm(Lettervec)
    sorted_Seq= Seq[p]
    Lettervec= Lettervec[p]
    return sorted_Seq, Lettervec
end