DF_cond_from_population_states = function(population_states,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    pst = population_states[1]
    pst_matrix = zeros(length(population_states),15)

    for i = 1:length(population_states)
        pst = population_states[i]
        pst = sort!(pst, [:nrow], rev=true)
        j=1;
        while j<=size(pst,1)&& pst[j,6]>10 && j<4 
            pst_matrix[i,(5*(j-1)+1):5*j].= collect(pst[j,1:5]);
            j= j+1
        end
    end

    PST_DF = DataFrame(pst_matrix, :auto)

    DF_cond = condense_data_frame(PST_DF,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    return PST_DF, DF_cond
end



DF_cond_from_population_states_3 = function(population_states,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    pst = population_states[1]
    pst_matrix = zeros(length(population_states),9)

    for i = 1:length(population_states)
        pst = population_states[i]
        pst = sort!(pst, [:nrow], rev=true)
        j=1;
        while j<=size(pst,1)&& pst[j,4]>10 && j<4  ##
            pst_matrix[i,(3*(j-1)+1):3*j].= collect(pst[j,1:3]);
            j= j+1
        end
    end

    PST_DF = DataFrame(pst_matrix, :auto)

    DF_cond = condense_data_frame_3(PST_DF,max_n,repeats,s5tn_v1_UInt,s5tn_v2_UInt,s5tn_v1,s5tn_v2)
    return PST_DF, DF_cond
end



reduce_population_states_to_connected_to_pst_to_3 = function(population_states_5, max_n, s5tn_v1_UInt,s5tn_v2_UInt)
    ## first reduce everything to 3's1
    population_states = deepcopy(population_states_5)
    for i = 1:length(population_states)
    population_states[i] = df_conversion_5_to_3_with_nrow(population_states_5[i])
    end
    ## then reduce everything to connected strategies
    for i = 1:length(population_states)
        pst = deepcopy(population_states[i])
        pst = sort!(pst, [:nrow], rev=true)
        for j = 1:size(population_states[i],1)
            strategy_explicit = number_to_strat_3(collect(pst[j,1:3]),max_n, s5tn_v1_UInt,s5tn_v2_UInt);
            stra_red = reduce_to_connected(strategy_explicit) ## 'graph' doesn't consider the states --> so it would keep those
            if sum(strategy_explicit[2:3,1])==2 && strategy_explicit[1,2]!=0
                println("i = $i, j = $j")
            end
            pst[j,1:3] = strat_3_to_number(stra_red,max_n,s5tn_v1,s5tn_v2);
        end

        population_states[i] = combine(groupby(pst, [:x1, :x2, :x3]), :nrow .=> sum);
        if names(population_states[i])[4] !="nrow"
            rename!(population_states[i], :nrow_sum => :nrow)
        end
    end
    return population_states
end

df_conversion_5_to_3_with_nrow = function(df)
    df_3 =DataFrame(zeros(Int64,size(df,1),1+convert(Int64,3*floor(size(df,2)/5))),:auto);
    df_3[:,end] = df[:,end];

    for i = 1:size(df_3,1)
        for j = 1:convert(Int64,floor(size(df,2)/5))
            s = collect(df[i,(5*(j-1)+1):5*j]);
            s = convert(Array{UInt128}, s);
            stra = number_to_strat_5_no_cleaning(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt);
            if stra[1,1]!=0
                stra_new = convert_strat_5_to_3(stra)
                stra_final = strat_3_to_number(stra_new,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
                df_3[i,(3*(j-1)+1):3*j] = stra_final
            end
        end
    end

    rename!(df_3, :x4 => :nrow)
    return df_3
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
        elseif  strategy5[1,i]==0 && strategy5[2,i]!=0## in case it's a leave --> conversion to 3strat isn't clear --> keep both around in a way 
            strategy3[2,i] = 1#strategy5[4,i] #[2:3,]
            strategy3[3,i] = 1#strategy5[5,i]#[2:3,]

        end
    end
    return strategy3
end



reduce_population_states_to_connected_to_pst = function(population_states, max_n, s5tn_v1_UInt,s5tn_v2_UInt)
    for i = 1:length(population_states)
        pst = deepcopy(population_states[i])
        pst = sort!(pst, [:nrow], rev=true)
        for j = 1:size(population_states[i],1)
            strategy_explicit = number_to_strat_5(collect(pst[j,1:5]),max_n, s5tn_v1_UInt,s5tn_v2_UInt);
            stra_red = reduce_to_connected_graph(strategy_explicit)
            if sum(strategy_explicit[2:5,1])==4 && strategy_explicit[1,2]!=0
                println("i = $i, j = $j")
            end
            pst[j,1:5] = strat_5_to_number(stra_red,max_n,s5tn_v1,s5tn_v2);
        end

        population_states[i] = combine(groupby(pst, [:x1, :x2, :x3, :x4, :x5]), :nrow .=> sum);
        if names(population_states[i])[6] !="nrow"
            rename!(population_states[i], :nrow_sum => :nrow)
        end
    end
    return population_states
end




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

save_top_x_by_k = function(storepath, DF_cond,x,k) #order looking at the k most prevalent strategies (k ={1,2,3}), and store the x most 
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
        #display(stra)
        stra1 = stra[1,:]
        #display(stra1)
        #display(stra1)
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
    end
    end

end



save_top_x_by_k_3_version = function(storepath, DF_cond,x,k) #order looking at the k most prevalent strategies (k ={1,2,3}), and store the x most 
    if k ==1
        df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, ]), nrow),[:nrow],rev=true)
    elseif k==2
        df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5, :x6]), nrow),[:nrow],rev=true)
    else #any number not 1 or 2 will be interpreted as 3
        df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9]), nrow),[:nrow],rev=true)
    end

    open(storepath, "w") do file
    for i = 1:x
    write(file,"\nrank $i ")
    dd = df2[i,end]
    write(file, "\ncount: $dd \n")
    for j = 1:k
        s = collect(df2[i,(3*(j-1)+1):3*j]);
        s = convert(Array{UInt128}, s)
        stra = number_to_strat_3(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
        if stra[1,1]!=0 
        #display(stra)
        stra1 = stra[1,:]
        #display(stra1)
        #display(stra1)
        stra2 = stra[2,:]
        #stra1 = action_to_symbol(stra1)
        stra3 = stra[3,:]   
#        stra4 = stra[4,:]
#        stra5 = stra[5,:]
            write(file, "$stra1 Strategy\n")
            write(file, "$stra2 C\n")
            write(file, "$stra3 D\n")
 #           write(file, "$stra4 DC\n")
 #           write(file, "$stra5 DD\n")
            write(file, "\n")
        end
    end
    end
    end

end


DF_to_single_pop_and_durations_pures = function(DF_cond)
    df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5]), nrow),[:nrow],rev=true)
    population_states = Array{Any}(undef,0);
    durations = []
    for i = 1:size(df2,1)
        durations = append!(durations,df2[i,6])
        population_states = append!(population_states,df2[i,1:5])
    end

    return population_states, durations
end



find_occurances_c_i_d_i = function(i,DF_cond,max_n)
    c_i =zeros(Int64,3,max_n); 
    d_i = zeros(Int64,3,max_n); 
    d_i[1,1:(1+i)].=1
    d_i[:,(i+2)]= [0 1 1]
    c_i[2:3,1] = 
    df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5, :x6]), nrow),[:nrow],rev=true)

    return
end


find_occurances_c_0_d_0 = function(DF_cond,max_n,s5tn_v1,s5tn_v2)
    c_0 =zeros(Int64,3,max_n); 
    d_0 = zeros(Int64,3,max_n); 
    c_0[:,1] = [1 1 2] 
    c_0[:,2] = [0 1 1] 
    c_0_number =strat_3_to_number(c_0,max_n,s5tn_v1,s5tn_v2)

    df2 = sort!(combine(groupby(DF_cond, [:x1, :x2, :x3, :x4, :x5, :x6]), nrow),[:nrow],rev=true)
    c_0_indices = ((df2[:,:x1].==c_0_number[1])+ (df2[:,:x2].==c_0_number[2]) + (df2[:,:x3].==c_0_number[3])).==3
    df3 = deepcopy(df2[c_0_indices,:])
    keep_indices = zeros(Int64,nrow(df3))
    for i = 2:nrow(df3)
        #println(i)
        strategy_explicit = number_to_strat_3(collect(df3[i,4:6]),max_n, s5tn_v1_UInt,s5tn_v2_UInt)
        if strategy_explicit[1,1]==2 && strategy_explicit[1,strategy_explicit[3,1]]==0 
            keep_indices[i]=1
        end
    end

    absolute = sum(collect(df3[(keep_indices.==1),7]))


    proportion = absolute/sum(collect(df2[:,7]))

    return absolute, proportion
end
df3[1,:]

