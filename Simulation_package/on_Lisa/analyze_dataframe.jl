convert_strat_5_to_3 = function(strategy5)
    strategy3 = zeros(Int64,3,size(strategy5,2))
    strategy3[1,:] = strategy5[1,:]
    for i = 1:size(strategy5,2)
        if strategy5[1,i]==1
            strategy3[2,i] = strategy5[2,i] #[2:3,]
            strategy3[3,i] = strategy5[3,i]#[2:3,]
        elseif strategy5[1,i]==2 #also encompa
            strategy3[2,i] = strategy5[4,i] #[2:3,]
            strategy3[3,i] = strategy5[5,i]#[2:3,]
        else  ## in case it's a leave --> conversion to 3strat isn't clear --> keep both around in a way 
            #strategy3[2,i] = 10*strategy5[2,i]+strategy5[3,i] #[2:3,]
            #strategy3[3,i] = 10*strategy5[4,i]+strategy5[5,i]#[2:3,]
        end
    end
    return strategy3
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

df_conversion_5_to_3 = function(df)
    df_3 =DataFrame(zeros(Int64,size(df,1),convert(Int64,3*floor(size(df,2)/5))),:auto);
    for i = 1:size(df_3,1)
        #print(i)
        for j = 1:convert(Int64,floor(size(df,2)/5))
            s = collect(df[i,(5*(j-1)+1):5*j]);
            s = convert(Array{UInt128}, s);
            stra = number_to_strat_5(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt);
            if stra[1,1]!=0
            stra_new = convert_strat_5_to_3(stra);
            stra_final = strat_3_to_number(stra_new,0,s5tn_v1,s5tn_v2);
            df_3[i,(3*(j-1)+1):3*j] = stra_final
            end
        end
    end
    return df_3
end



detect_terminal_cycle = function(strategy3)
    #display(strategy3)
    leave = 1;
    c_tc = []; #cooperation in terminal cycle
    state = 1;
    #strategy playing against itself --> implies that both are always gonna be in the same state, and when a state is revisited, this signifies a cycle
    statelist = [1];
    statelist_pre = [];
    reached = 0; #did we reach a cycle (can be a 1-cycle obviously)
    while reached ==0
        state = strategy3[strategy3[1,state]+1,state]
        if strategy3[1,state] ==0
            reached = -1
            append!(statelist,state)
        else
            
            if  sum(statelist.==state)==0
                append!(statelist,state)
            else
                i = findfirst(statelist.==state)
                statelist_pre = statelist[1:(i-1)]
                statelist = statelist[i:end]
                reached = 1
            end
        end
    end 
    if reached ==1
        leave = 0
        c_tc = 2 - mean(strategy3[1,statelist]); #cooperation level in the terminal cycle
    end
    return  leave,c_tc, statelist, statelist_pre
 end


df_detect_terminal_cycle = function(df,nstrat)
    terminal_cycle = zeros(size(df,1),nstrat)
    terminal_cycle_value = zeros(size(df,1),nstrat)
    for i = 1:size(df,1)
        for j = 1:nstrat
            s = collect(PST_DF[i,(3*(j-1)+1):3*j]);
            s = convert(Array{UInt128}, s)
            stra = number_to_strat_3(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
            if stra[1,1]!=0 
                #println(i)
                leave,c_tc, statelist, statelist_pre = detect_terminal_cycle(stra)
                if leave !=1
                    terminal_cycle[i] = 1 
                    terminal_cycle_value[i] = c_tc
                end
            end
        end
    end

    return terminal_cycle, terminal_cycle_value
end






detect_c_i = function(strategy3)
    ci = 0; 
    i = 0;
    leave,c_tc, statelist, statelist_pre = detect_terminal_cycle(strategy3)
    if c_tc ==1 #cooperate fully in terminal cycle
        if length(statelist_pre)>0 #need a non-empty runup
            if prod(stra[1,statelist_pre].==2) ==1 # defect in all runup states
                if  prod(stra[1,stra[2,statelist_pre]].==0)==1 #leave in all cooperations in the handshake phase
                    if  prod(stra[1,stra[3,statelist]].==0)==1 #leave in all defections in the terminal cycle
                        i = length(statelist_pre)
                        ci = 1 #for true
                    end
                end
            end
        end
    end
    return  ci,i
end






###############


# go through --> elicit all properties, return all properties in a pst_matrix
df_stay_space_analysis = function(df,nstrat)
    # this looks at strategies that don't leave themselves
    stayspace = zeros(size(df,1),nstrat,14); # [:,:,1] runup length, [:,:,2] terminal cycle length
    # [:,:,3] runup C - proportion, [:,:,4] terminal cycle C - proportion,
    # [:,:,5] runup leave_proportion (if opponent behaves differently), # [:,:,6] terminal cycle leave proportion if opponent behaves differently
    # [:,:,7] runup defection proportion (if opponent behaves differently), # [:,:,8] terminal cycle defection proportion if opponent behaves differently
    # [:,:,9] is it c_i --> if yes, i, if not, 0
    # [:,:,10] is it like c_i in the terminal cycle (but might punish before by defecting)--> if yes, i, if not, 0
    # [:,:,11] is it like c_i in the runup (but might punish after by defecting)--> if yes, i, if not, 0
    # [:,:,12] is it like c_i but allows for punishmend without leaving in either phase
    # [:,:,13] is it like 12 but allows also for more non-D handshakes
    # [:,:,14] is it like 13 but also allows for non-pure D terminal cycles (but the terminal cycle needs a C somewhere)

    terminal_cycle = zeros(Int64,size(df,1),nstrat)
    for i = 1:size(df,1)
        for j = 1:nstrat
            #println(i)
            s = collect(df[i,(3*(j-1)+1):3*j]);
            s = convert(Array{UInt128}, s)
            stra = number_to_strat_3(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
#            stra = number_to_strat_5(s,max_n, s5tn_v1_UInt,s5tn_v2_UInt)
#            stra = convert_strat_5_to_3(stra)
            if stra[1,1]!=0 
                #display(stra)
                leave,c_tc, statelist, statelist_pre = detect_terminal_cycle(stra)
                if leave !=1
                    terminal_cycle[i,j] = 1 
                    stayspace[i,j,1] = length(statelist_pre)
                    stayspace[i,j,2] = length(statelist)
                    if length(statelist_pre)>0
                        stayspace[i,j,3] = 2 - mean(stra[1,statelist_pre]); #cooperation level in the terminal cycle
                        stayspace[i,j,5] = mean(stra[1,stra[2,statelist_pre]].==0)
                        stayspace[i,j,7] = mean(stra[1,stra[2,statelist_pre]].==2)
                    end
                    stayspace[i,j,4] = c_tc
                    stayspace[i,j,6] = mean(stra[1,stra[3,statelist]].==0)
                    stayspace[i,j,8] = mean(stra[1,stra[3,statelist]].==2)
                    if (1-stayspace[i,j,3])*prod(stayspace[i,j,4:6])==1
                        stayspace[i,j,9] = length(statelist_pre)
                    end
                    if (1-stayspace[i,j,3])*prod(stayspace[i,j,4:5])*sum(stayspace[i,j,[6,8]])==1
                        stayspace[i,j,10] = length(statelist_pre)
                    end
                    if (1-stayspace[i,j,3])*prod(stayspace[i,j,[4,6]])*sum(stayspace[i,j,[5,7]])==1
                        stayspace[i,j,11] = length(statelist_pre)
                    end
                    if ((1-stayspace[i,j,3])*stayspace[i,j,[4]]*sum(stayspace[i,j,[6,8]])*sum(stayspace[i,j,[5,7]]))[1]==1
                        stayspace[i,j,12] = length(statelist_pre)
                    end
                    if (stayspace[i,j,[4]].*sum(stayspace[i,j,[6,8]])*sum(stayspace[i,j,[5,7]]))[1]==1
                        stayspace[i,j,13] = length(statelist_pre)
                    end
                    if ((stayspace[i,j,[4]].>0)*sum(stayspace[i,j,[6,8]])*sum(stayspace[i,j,[5,7]]))[1]==1
                        stayspace[i,j,14] = length(statelist_pre)
                    end

                end
            end
        end
    end
    return terminal_cycle, stayspace
end 