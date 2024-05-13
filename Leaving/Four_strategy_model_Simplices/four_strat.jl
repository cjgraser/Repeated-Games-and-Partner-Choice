#using Pkg
#Pkg.add("Optim")
using Optim

function objective(x)
    x11, x12, x13, x14, x22, x23, x24, x33, x34, x44 = x
    b11, b12, b13, b14, b22, b23, b24, b33, b34, b44 = B

    constraints = [        (1-x11 -x12 -x13 -x14 -x22 -x23 -x24 -x33 -x34 -x44),
        (x1 - (x11 + 0.5*x12 + 0.5*x13 + 0.5*x14)),
        (x2 - (x22 + 0.5*x12 + 0.5*x23 + 0.5*x24)),
        (x3 - (x33 + 0.5*x13 + 0.5*x23 + 0.5*x34)),
        (x4 - (x44 + 0.5*x14 + 0.5*x24 + 0.5*x34)),
        (b11*x11 - (((x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 + 0.5*x14*b14)*(x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 + 0.5*x14*b14))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b12*x12 - (((x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 + 0.5*x14*b14)*(x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 + 0.5*x24*b24))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b13*x13 - (((x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 + 0.5*x14*b14)*(x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 + 0.5*x34*b34))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b14*x14 - (((x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 + 0.5*x14*b14)*(x44*b44 + 0.5*x14*b14 + 0.5*x24*b24 + 0.5*x34*b34))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b22*x22 - (((x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 + 0.5*x24*b24)*(x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 + 0.5*x24*b24))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b23*x23 - (((x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 + 0.5*x24*b24)*(x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 + 0.5*x34*b34))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b24*x24 - (((x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 + 0.5*x24*b24)*(x44*b44 + 0.5*x14*b14 + 0.5*x24*b24 + 0.5*x34*b34))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b33*x33 - (((x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 + 0.5*x34*b34)*(x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 + 0.5*x34*b34))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b34*x34 - (((x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 + 0.5*x34*b34)*(x44*b44 + 0.5*x14*b14 + 0.5*x24*b24 + 0.5*x34*b34))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44)))),
        (b44*x44 - (((x44*b44 + 0.5*x14*b14 + 0.5*x24*b24 + 0.5*x34*b34)*(x44*b44 + 0.5*x14*b14 + 0.5*x24*b24 + 0.5*x34*b34))/((x11*b11 + x12*b12 + x13*b13 + x14*b14 + x22*b22 + x23*b23 + x24*b24 + x33*b33 + x34*b34 + x44*b44))))
    ]

    return sum(constraint^2 for constraint in constraints)
end



function objective_three(x,X)
    # the general function has 4 --> but here we use that two strategies have the same properties
    x1, x2, x3  = X
    x11, x12, x13,  x22, x23,  x33 = x

    constraints = [        (1-x11 -x12 -x13 -x22 -x23 -x33 ),
        (x1 - (x11 + 0.5*x12 + 0.5*x13)),
        (x2 - (x22 + 0.5*x12 + 0.5*x23)),
        (x3 - (x33 + 0.5*x13 + 0.5*x23)),
        (b11*x11 - (((x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 )*(x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 ))/((x11*b11 + x12*b12 + x13*b13 +  x22*b22 + x23*b23  + x33*b33 )))),
        (0.5*b12*x12 - (((x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 )*(x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 ))/((x11*b11 + x12*b12 + x13*b13 +  x22*b22 + x23*b23  + x33*b33 )))),
        (0.5*b13*x13 - (((x11*b11 + 0.5*x12*b12 + 0.5*x13*b13 )*(x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 ))/((x11*b11 + x12*b12 + x13*b13 +  x22*b22 + x23*b23  + x33*b33 )))),
        (b22*x22 - (((x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 )*(x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 ))/((x11*b11 + x12*b12 + x13*b13 +  x22*b22 + x23*b23  + x33*b33 )))),
        (0.5*b23*x23 - (((x22*b22 + 0.5*x12*b12 + 0.5*x23*b23 )*(x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 ))/((x11*b11 + x12*b12 + x13*b13 +  x22*b22 + x23*b23  + x33*b33 )))),
        (b33*x33 - (((x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 )*(x33*b33 + 0.5*x13*b13 + 0.5*x23*b23 ))/((x11*b11 + x12*b12 + x13*b13 +  x22*b22 + x23*b23  + x33*b33 )))),
    ]

    return sum(constraint^2 for constraint in constraints)
end


using Combinatorics
multisets(n, k) = map(A -> [sum(A .== i) for i in 1:n],
                      with_replacement_combinations(1:n, k))

delta = 0.9;

strategy1 = ones(Int64,3,1); strategy1[1]=2
strategy2 = ones(Int64,3,3); strategy2[:,1] = [2 ;2;2]; strategy2[:,2] = [1 ; 2; 3] ; strategy2[:,3] = [0 ; 1; 1]; strategy2
strategy3 = ones(Int64,3,3); strategy3[:,1] = [2 ;2;2]; strategy3[:,2] = [1 ; 2; 3] ;strategy3[:,3] = [2 ; 3; 3] ; strategy3
strategy4 = ones(Int64,3,2); strategy4[:,1] = [2 ;2;2]; strategy4[:,2] = [1 ; 2; 2] ; strategy4

compute_stationary_pair_distributions = function(state)

    x1 = state[1]
    x2 = state[2]
    x3 = state[3]+state[4]
    lower = zeros(6)
    upper = ones(6)

    initial_x = [x1^2, 2*x1*x2, 2*x1*x3, x2*x2, 2*x2*x3, x3*x3]
    objective_three(initial_x)

    opt = optimize(objective_three, lower, upper, initial_x, Fminbox(LBFGS()))#,options)
    res = Optim.minimizer(opt)
    # res gives us     x11, x12, x13,  x22, x23,  x33 
    #now map this to a 4-type matrix 
    #     x11, x12, x13, x14, x22, x23, x24, x33, x34, x44 = x

    proportion3 = state[3]/(state[3]+ state[4])

    res4 = [res[1] res[2] proportion3*res[3] (1-proportion3)*res[3] res[4] proportion3*res[5] (1-proportion3)*res[5] (proportion3^2)*res[6] 2*(proportion3*(1-proportion3))*res[6] ((1-proportion3)^2)*res[6]]

    res4 = res4./(sum(res4)) #because of numerical errors   
    return res4
end



compute_stationary_pair_distributions_3 = function(state)

    x1 = state[1]
    x2 = state[2]
    x3 = state[3]
    lower = zeros(6)
    upper = ones(6)


    initial_x = [x1^2, 2*x1*x2, 2*x1*x3, x2*x2, 2*x2*x3, x3*x3]
    objective_function(x) = objective_three(x, [x1,x2,x3])

    opt = optimize(objective_function, lower, upper, initial_x, Fminbox(LBFGS()))#,options)
    res = Optim.minimizer(opt)
    res3 = res./(sum(res))    
    return res3
end


# function that extends Pair_Freqs3 to Pair_Freqs4

extend_to_4 = function(States3, States3_integer, Pair_Freqs3)
    States4 = [zeros(4)]; States4[1][1]=1; 
    States4_integer  = [zeros(Int64,4)]; States4_integer[1][1]=States3_integer[1][1]; 
    Pair_Freqs4 = [zeros(10)]; Pair_Freqs4[1][1] = 1 


    for i = 2:length(States3_integer)
        println(i)
        n3 = States3_integer[i][3]
        for j = n3:-1:0
            append!(States4_integer, [[States3_integer[i][1], States3_integer[i][2], j, n3-j]])
            append!(States4, [[States3_integer[i][1], States3_integer[i][2], j, n3-j]./N])
            res = Pair_Freqs3[i]
            
            if n3>0
                proportion3 = j/n3
            else 
                proportion3 = 0
            end

            res4 = [res[1], res[2], proportion3*res[3], (1-proportion3)*res[3], res[4], proportion3*res[5], (1-proportion3)*res[5], (proportion3^2)*res[6], 2*(proportion3*(1-proportion3))*res[6], ((1-proportion3)^2)*res[6]]
            append!(Pair_Freqs4, [res4])
        end
    end


    return States4, States4_integer, Pair_Freqs4
end


include("C:/Users/chris/OneDrive - UvA/Julian_Code/Full_Project/Leaving/Explicit_Markov_Chain/functions_to_fill_markov_matrix.jl")
#PM = [2 0; 3 1]
#delta = 0.8;
PM = [3 0; 4 1]
delta = 0.9;
tolerance = 0.0001;
PAYOFFS = zeros(4,4)
EBS = zeros(10)
EBS[1],PAYOFFS[1,1],a =get_ebs_and_ps(strategy1,strategy1,delta,PM,tolerance)
EBS[2],PAYOFFS[1,2],PAYOFFS[2,1] =get_ebs_and_ps(strategy1,strategy2,delta,PM,tolerance)
EBS[3],PAYOFFS[1,3],PAYOFFS[3,1] =get_ebs_and_ps(strategy1,strategy3,delta,PM,tolerance)
EBS[4],PAYOFFS[1,4],PAYOFFS[4,1] =get_ebs_and_ps(strategy1,strategy4,delta,PM,tolerance)
EBS[5],PAYOFFS[2,2],a =get_ebs_and_ps(strategy2,strategy2,delta,PM,tolerance)
EBS[6],PAYOFFS[2,3],PAYOFFS[3,2] =get_ebs_and_ps(strategy2,strategy3,delta,PM,tolerance)
EBS[7],PAYOFFS[2,4],PAYOFFS[4,2] =get_ebs_and_ps(strategy2,strategy4,delta,PM,tolerance)
EBS[8],PAYOFFS[3,3],a =get_ebs_and_ps(strategy3,strategy3,delta,PM,tolerance)
EBS[9],PAYOFFS[3,4],PAYOFFS[4,3] =get_ebs_and_ps(strategy3,strategy4,delta,PM,tolerance)
EBS[10],PAYOFFS[4,4],a =get_ebs_and_ps(strategy4,strategy4,delta,PM,tolerance)

B4 = (-delta*(-EBS.+1)).+1
B = B4[[1, 2, 3, 5, 6, 8]]
b11,b12,b13,b22,b23,b33 = B



N = 40
States3 = multisets(3,N)./N
States3_integer = multisets(3,N)


Pair_Freqs3 = compute_stationary_pair_distributions_3.(States3)
States4, States4_integer, Pair_Freqs4 =extend_to_4(States3, States3_integer, Pair_Freqs3)


per_state_payoffs = [[PAYOFFS[1,1],0 ,0,0]] 
for i = 2:length(States4_integer)
    pf = Pair_Freqs4[i]
    s = States4[i]

    p1 = (pf[1]*PAYOFFS[1,1] + 0.5*pf[2]*PAYOFFS[1,2] + 0.5*pf[3]*PAYOFFS[1,3]+ 0.5*pf[4]*PAYOFFS[1,4])/(s[1])
    p2 = (pf[5]*PAYOFFS[2,2] + 0.5*pf[2]*PAYOFFS[2,1] + 0.5*pf[6]*PAYOFFS[2,3]+ 0.5*pf[7]*PAYOFFS[2,4])/(s[2])
    p3 = (pf[8]*PAYOFFS[3,3] + 0.5*pf[6]*PAYOFFS[3,2] + 0.5*pf[3]*PAYOFFS[3,1]+ 0.5*pf[9]*PAYOFFS[3,4])/(s[3])
    p4 = (pf[10]*PAYOFFS[4,4] + 0.5*pf[4]*PAYOFFS[4,1] + 0.5*pf[7]*PAYOFFS[4,2]+ 0.5*pf[9]*PAYOFFS[4,3])/(s[4])
    payoffs_here = [p1, p2 , p3 , p4]
    append!(per_state_payoffs, [payoffs_here])
end


# now this generates the selection markov matrix

Markov = zeros(length(States4_integer),length(States4_integer))

for i = 1:size(Markov,1)
    println(i)
    for j = 1:size(Markov,2)
        change = States4_integer[i].- States4_integer[j]
        if sum(change.==1) ==1 && sum(change.==-1) ==1 && sum(change.==0)==2#so exactly one goes up, and one goes down by one
            up_index = findfirst(change.==-1)
            down_index = findfirst(change.==1)
            if States4_integer[i][up_index]>=1 #only things that are present can win in selection
                
                total_payoff = States4_integer[i][States4_integer[i].>0]'*per_state_payoffs[i][States4_integer[i].>0]
                winner_payoff = States4_integer[i][up_index]'*per_state_payoffs[i][up_index]
                Markov[i,j] = (States4_integer[i][down_index]/N)*(winner_payoff/total_payoff) #probability that down_index dies times sum of payoff of up_index individuals divided by total payoffs
            end
        end
    end
end

for i = 1:size(Markov,1)
    Markov[i,i] = 1-sum(Markov[i,setdiff(1:size(Markov,1),i)])
end







# now this generates the mutation markov matrix

mut_Markov = zeros(length(States4_integer),length(States4_integer))

for i = 1:size(mut_Markov,1)
    println(i)
    for j = 1:size(mut_Markov,2)
        change = States4_integer[i].- States4_integer[j]
        if sum(change.==1) ==1 && sum(change.==-1) ==1 && sum(change.==0)==2 #so exactly one goes up, and one goes down by one
            down_index = findfirst(change.==1)
            mut_Markov[i,j] = (1/3)*States4[i][down_index] #probability that down_index mutates, and then times 1/3, because there are three possibilities
        end
    end
end


for i = 1:size(Markov,1)
    mut_Markov[i,:] = mut_Markov[i,:]./sum(mut_Markov[i,:])
end


######
###### lower dimensional mutation matrix --> to explore facets of the simplex (if needed)

mut123_Markov = deepcopy(mut_Markov)
mut124_Markov = deepcopy(mut_Markov)
mut134_Markov = deepcopy(mut_Markov)
mut234_Markov = deepcopy(mut_Markov)


for i = 1:size(mut_Markov,1)
    println(i)
    for j = 1:size(mut_Markov,2)
        change = States4_integer[i].- States4_integer[j]
        if sum(change.==1) ==1 && sum(change.==-1) ==1 && sum(change.==0)==2#so exactly one goes up, and one goes down by one
            #
            up_index = findfirst(change.==-1)
            down_index = findfirst(change.==1)
            if up_index ==4
                mut123_Markov[i,j]=0
            elseif up_index ==3
                mut124_Markov[i,j]=0
            elseif up_index ==2
                mut134_Markov[i,j]=0
            elseif up_index ==1
                mut234_Markov[i,j]=0
            end        
        end
    end
end


for i = 1:size(Markov,1)
    mut123_Markov[i,:] = mut123_Markov[i,:]./sum(mut123_Markov[i,:])
    mut124_Markov[i,:] = mut124_Markov[i,:]./sum(mut124_Markov[i,:])
    mut134_Markov[i,:] = mut134_Markov[i,:]./sum(mut134_Markov[i,:])
    mut234_Markov[i,:] = mut234_Markov[i,:]./sum(mut234_Markov[i,:])
end



create_selection_mutation_matrix = function(Markov, mut_Markov,mu)
    M = (Markov.*(1-mu)).+ (mut_Markov.*mu)
    for i = 1:size(M,1)
        M[i,i] = 1-sum(M[i,setdiff(1:size(M,1),i)])
    end
    return M 
end


using DiscreteMarkovChains


mu_vector = [0.00001, 0.0001, 0.001, 0.0025 , 0.005 ,0.007 ,0.009, 0.01, 0.013, 0.016, 0.02, 0.025, 0.03,0.04, 0.05,0.06,0.08, 0.1]
SDs = Any[]
for i = 1:length(mu_vector)
    println(i)
    MM = create_selection_mutation_matrix(Markov, mut_Markov,mu_vector[i])
    chain = DiscreteMarkovChain(MM)
    SDs = append!(SDs,[stationary_distribution(chain)])
end


find_pure_indices = function(States4_integer, N)
    index1 = 0
    index2 = 0
    index3 = 0
    index4= 0
    for i = 1:length(States4_integer)
        if States4_integer[i][1]==N
            index1 = i 
        end
        if States4_integer[i][2]==N
            index2 = i 
        end
        if States4_integer[i][3]==N
            index3 = i 
        end
        if States4_integer[i][4]==N
            index4 = i 
        end
    end
    return index1,index2,index3,index4
end

index1, index2, index3, index4 = find_pure_indices(States4_integer, N)
States4_integer[index1]
States4_integer[index2]
States4_integer[index3]
States4_integer[index4]

SD_pure_states = zeros(4,length(mu_vector))

for i = 1:length(mu_vector)
    SD_pure_states[:,i] = SDs[i][[index1,index2,index3,index4]]
end

find_state_in_list = function(States4_integer,state)
    index = 0
    for i = 1:length(States4_integer)
        if States4_integer[i]==state
            index = i
        end
    end
    return index
end


mv = zeros(size(Markov,2))
for i = 1:length(mv)
    mv[i] = Markov[i,i]
end



##############3
##############
############### And we can also use the functions above to create plots of assortment as a function of frequency

N = 400
D_proportion = collect(0:1/N:1)
Assortment = zeros(length(D_proportion))
Payoffs_D = zeros(length(D_proportion))
Payoffs_C1 = zeros(length(D_proportion))

for i = 2:length(D_proportion)-1
    s = [D_proportion[i], 1-D_proportion[i] ,0]
    pf = compute_stationary_pair_distributions_3(s)
    Assortment[i] = (1/(1-D_proportion[i]))*(pf[4]) - (1/(D_proportion[i]))*(1/2)*pf[2]
    Payoffs_D[i] = (pf[1]*PAYOFFS[1,1] + 0.5*pf[2]*PAYOFFS[1,2])/(s[1])
    Payoffs_C1[i] = (pf[4]*PAYOFFS[2,2] + 0.5*pf[2]*PAYOFFS[2,1]) /(s[2])
    
end

p11 = PAYOFFS[1,1]
p13 = PAYOFFS[1,3]
p31 = PAYOFFS[3,1]
p33 = PAYOFFS[3,3]

Payoffs_D_no_assortment= zeros(N-1)
Payoffs_D_PUN_no_assortment = zeros(N-1)

for i = 1:(N-1)
    Payoffs_D_no_assortment[i] = (i/N)*p11 + ((N-i)/N)*p13
    Payoffs_D_PUN_no_assortment[i] = (i/N)*p31 + ((N-i)/N)*p33
end


delta = 0.9
compute_stationary_pair_distributions_3([0.4,0.6,0])

0.3275415957310722+ 0.14491680642030874+ 0.5275415942785167

b22