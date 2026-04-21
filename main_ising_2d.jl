using Random
using Plots
using LaTeXStrings
using GLM
using DataFrames
using CSV
using StatsBase
using LsqFit


current_dir = @__DIR__

J = 1

#I want to use Wolff algorithm to simulate a 2D Ising model


"""
For all this cases I apply metropolis for T>2.6 and then the rest wolff

"""
#To do the results
"""
#This one will be for doing calculations
T_array_1 = 5.0:-0.001:3.000
T_array_2 = 2.999:-0.001:2.600
T_array_3 = 2.599:-0.001:2.000
T_array_4 = 1.999:-0.001:0.1
"""


#To do the graphs
"""
#This one is for representation
T_array_1 = 5.0:-0.1:3.0
T_array_2 = 2.95:-0.05:2.65
T_array_3 = 2.6:-0.05:2.0
T_array_4 = 1.9:-0.1:0.1
"""



#For the time correlation

T_array_1 = 3.00:-0.01:2.8
T_array_2 = 2.79:-0.01:2.61
T_array_3 = 2.60:-0.01:2.30
T_array_4 = 2.29:-0.01:2.00


"""
Here I apply only wolff
"""

"""
#For the critical tempertaure
T_array_1 = 2.269:-0.0001:2.250
T_array_2 = 2.249:-0.0001:2.200
T_array_3 = 2.199:-0.0001:2.100
T_array_4 = 2.099:-0.0001:1.900
"""




T_array = vcat(T_array_1, T_array_2, T_array_3, T_array_4)

MCS_thermalization = 2000

#Estos valores los he escogido tras hacer el análisis de los tiempos de correlación. He dado margen en los tiempos de correlación.
MCS_1 = 40 #20 #5  #quizás poner aquí otra cosa
MCS_2 = 40 #10 #40 
MCS_3 = 200 #5 #200
MCS_4 = 200 #2 #2

#We want to take 1000 measurements for each temperature


#We first create the arrays with the neighbours
function neighbours_creator(L)
    n_r = zeros(L*L) #n1
    n_t = zeros(L*L) #n2
    n_l = zeros(L*L) #n3
    n_b = zeros(L*L) #n4
    

    for ix in 1:L
        for iy in 1:L
            i = (iy - 1) * L + ix
            ix1 = ix + 1
            if ix1 == L + 1
                ix1 = 1
            end
            n_r[i] = (iy - 1) * L + ix1
            
            iy2 = iy + 1
            if iy2 == L + 1
                iy2 = 1
            end
            n_t[i] = (iy2 - 1) * L + ix

            ix3 = ix - 1
            if ix3 == 0
                ix3 = L
            end
            n_l[i] = (iy - 1) * L + ix3

            iy4 = iy - 1
            if iy4 == 0
                iy4 = L
            end
            n_b[i] = (iy4 - 1) * L + ix
        end
    end
    
    return n_r, n_t, n_l, n_b
end


# Here is where is the simulation where all the data is stored for each size of the lattice. Now I should do a program were 
# I call this program and I store the matrix for each size in a dataframe so then I can do all the calculation properly.
#I should get 5 or 6 different dataframes, depending if I include the simulation with L = 128.

"""

Just try it for L = 4 because I need to include the rest of things, like the correlation time or the correlation function.

IMPORTANT STUFF I HAVEN'T ADDED YET:

-AS I HAVE SAID BEFORE, I HAVEN'T ADDED THE CORRELATION TIME NEITHER THE CORRELATION FUNCTION 

-WHEN I'M DOING THE MEAN IS THE STANDAR MEAN OR IS THE MEAN FOLLOWING THE BOLTZMANN DISTRIBUTION??? -----> LOOK AT PAGES 100 TO 103

-I HAVE TO CREATE ANOTHER PROGRAM TO CALL THIS ONE AND THEN DO ALL THE SIMULTIONS.


"""


"""
I have used here the limit for counts as 1. Is a Monte Carlo Step just one flip of a cluster or is until you flip N spins?
"""



function simulation(L,N,s_array,MCS_1, MCS_2 ,MCS_3 , MCS_4, MCS_thermalization,T_array)

    correlation_matrix = zeros(length(T_array), 2)
    #Here we creat different matrixes to store the results
    n_r, n_t, n_l, n_b = neighbours_creator(L)
    neigh = Int.(hcat(n_r, n_t, n_l, n_b))

    # Ensure the shape is (4, 16)
    neigh = transpose(neigh)
    
    t_index = 1

    mean_std_H_M_temperature_matrix = zeros(length(T_array),11) #Matrix to store the mean and variance of H and M for each temperature. The fifth value is the temperature.

    r = r1k(L)

    g_r = zeros(length(T_array),N)

    """
    Here I only want to do τ_M and τ_H the first time to estimate the number of monte carlo steps
    """

    τ_H_array = zeros(length(T_array))
    τ_M_array = zeros(length(T_array))
    τ_H_2_array = zeros(length(T_array))
    τ_M_2_array = zeros(length(T_array))
    τ_M_4_array = zeros(length(T_array))
    


    #First we start the system in a random state
    for i in 1:N
        u = rand()
        if u>=0.5
            s_array[i] = +1
        else
            s_array[i] = -1
        end
    end

    #High temperatures part
    for T in T_array_1

        H_M_matrix = zeros(1000,2) #matrix we will store each measurment of H and M for each value of the temperature
        # Por el momento lo dejo así
        β = 1/T
        pa = 1-exp(-2*J*β) #We need to add here the temperature


        #We also want to measure the correlation length. The correlation time will be measured apart, using just 1 MC step.
        s1_sj = zeros(1000, N)


        for _ in 1:MCS_thermalization
            i = rand(1:N) #Chooses a random spin to flip it.
            s_array=metropolis(s_array, neigh, T) #
        end

        # We save the results of the measurment
        H_M_matrix[1,:] .= hamiltonian_and_M(s_array, neigh, N) #Does this correctly

        
        s1_sj[1, :] .=  [s_array[1] * s_array[j] for j in 1:N]

        
        



        #Now I will do this 999 times more

        for n in 2:1000

            #Here we differenciate between the region of temperatures near the critical point and the rest of the range of temperatures.
            #here we use wolff

            #here we use metropolis
            for _ in 1:MCS_1
                i = rand(1:N)
                s_array= metropolis(s_array, neigh, T)
            end

            #We save the results of the measurment
            H_M_matrix[n,:] .= hamiltonian_and_M(s_array, neigh, N)

            
            s1_sj[n, :] .=  [s_array[1] * s_array[j] for j in 1:N]


        end


        
        #Once we have taken 1000 measurments we measuere the mean and the variance of the energy and the magnetization
        mean_std_H_M_temperature_matrix[t_index, :] .= mean(H_M_matrix[:,1]), std(H_M_matrix[:,1]), mean(H_M_matrix[:,2]), std(H_M_matrix[:,2]), mean(H_M_matrix[:,1].^2), std(H_M_matrix[:,1].^2) ,mean(H_M_matrix[:,2].^2), std(H_M_matrix[:,2].^2) ,T_array[t_index], mean(H_M_matrix[:,2].^4), std(H_M_matrix[:,2].^4)

        

        
        #We also want to compute the correlation length so we compute the correlation function g(r)

        m_avg = mean(H_M_matrix[:,2])/N
        

        s1_sj_avg = vec(mean(s1_sj, dims = 1))

        g_avg = s1_sj_avg .- m_avg^2 #We have now a vector of the value of g_r for j∈[1,N]

        # Create a DataFrame with r and g_avg
        df = DataFrame(r = r, g_avg = g_avg)

        # Group by r and calculate the average of g_avg for each group. This means creat a DataFrame with all the g_avg that correspond to the spins at distance r.
        df_grouped = combine(groupby(df, :r), :g_avg => mean)

        # Extract the unique r values and the corresponding average g_avg values
        r_unique = df_grouped.r
        g_avg_mean = abs.(df_grouped.g_avg_mean) # because I get negative numbers due to numerical errors

        p = sortperm(r_unique)
        r_unique = r_unique[p]
        g_avg_mean = g_avg_mean[p]

        # Define the exponential fitting function
        exp_func(x, p) = exp.(-x ./ p[1])
        fit_result = curve_fit(exp_func, r_unique, g_avg_mean, [1.0])

        # Extract the fitted parameter ξ with error
        ξ = fit_result.param[1]; σ_ξ = estimate_errors(fit_result)[1]

        correlation_matrix[t_index,1] = ξ
        correlation_matrix[t_index,2] = σ_ξ

        

        
        

        """
        Look at main luis lines 920 to 971

        HOW DO I GET THE ERROR OF G_R????????
        
        
        
        Here I will compute the correlation time for each value of T given L. The idea is to do first a simulation with just MCS, compute the correlation lenght and put it in a csv file.
        """

        τ_H_array[t_index] = Correlation_time(H_M_matrix[:,1])
        τ_M_array[t_index] = Correlation_time(H_M_matrix[:,2])
        τ_H_2_array[t_index] = Correlation_time(H_M_matrix[:,1].^2)
        τ_M_2_array[t_index] = Correlation_time(H_M_matrix[:,2].^2)
        τ_M_4_array[t_index] = Correlation_time(H_M_matrix[:,2].^4)


        t_index = t_index + 1
    end 

    #Temperatures near-after the critical point
    for T in T_array_2


        H_M_matrix = zeros(1000,2) #matrix we will store each measurment of H and M for each value of the temperature
        # Por el momento lo dejo así
        β = 1/T
        pa = 1-exp(-2*J*β) #We need to add here the temperature


        #We also want to measure the correlation length. The correlation time will be measured apart, using just 1 MC step.
        s1_sj = zeros(1000, N)


        for _ in 1:MCS_thermalization
            i = rand(1:N) #Chooses a random spin to flip it.
            s_array=metropolis(s_array, neigh, T) 
        end

        # We save the results of the measurment
        H_M_matrix[1,:] .= hamiltonian_and_M(s_array, neigh, N) #Does this correctly

        s1_sj[1, :] .=  [s_array[1] * s_array[j] for j in 1:N]


        #Now I will do this 999 times more

        for n in 2:1000

            #Here we differenciate between the region of temperatures near the critical point and the rest of the range of temperatures.
            #here we use wolff
            for _ in 1:MCS_2
                i = rand(1:N) #Chooses a random spin to flip it.
                s_array=metropolis(s_array, neigh, T) #
            end

            #We save the results of the measurment
            H_M_matrix[n,:] .= hamiltonian_and_M(s_array, neigh, N)

            s1_sj[n, :] .=  [s_array[1] * s_array[j] for j in 1:N]
            


        end

        
        #Once we have taken 1000 measurments we measuere the mean and the variance of the energy and the magnetization
        mean_std_H_M_temperature_matrix[t_index, :] .= mean(H_M_matrix[:,1]), std(H_M_matrix[:,1]), mean(H_M_matrix[:,2]), std(H_M_matrix[:,2]), mean(H_M_matrix[:,1].^2), std(H_M_matrix[:,1].^2) ,mean(H_M_matrix[:,2].^2), std(H_M_matrix[:,2].^2), T_array[t_index], mean(H_M_matrix[:,2].^4), std(H_M_matrix[:,2].^4)

        

        #We also want to compute the correlation length so we compute the correlation function g(r)

        m_avg = mean(H_M_matrix[:,2])/N
        

        s1_sj_avg = vec(mean(s1_sj, dims = 1))

        g_avg = s1_sj_avg .- m_avg^2 #We have now a vector of the value of g_r for j∈[1,N]

        # Create a DataFrame with r and g_avg
        df = DataFrame(r = r, g_avg = g_avg)

        # Group by r and calculate the average of g_avg for each group. This means creat a DataFrame with all the g_avg that correspond to the spins at distance r.
        df_grouped = combine(groupby(df, :r), :g_avg => mean)

        # Extract the unique r values and the corresponding average g_avg values
        r_unique = df_grouped.r
        g_avg_mean = abs.(df_grouped.g_avg_mean) # because I get negative numbers due to numerical errors

        p = sortperm(r_unique)
        r_unique = r_unique[p]
        g_avg_mean = g_avg_mean[p]

        # Define the exponential fitting function
        exp_func(x, p) = exp.(-x / p[1])
        fit_result = curve_fit(exp_func, r_unique, g_avg_mean, [1.0])

        # Extract the fitted parameter ξ with error
        ξ = fit_result.param[1]; σ_ξ = estimate_errors(fit_result)[1]

        correlation_matrix[t_index,1] = ξ
        correlation_matrix[t_index,2] = σ_ξ


        """
        Look at main luis lines 920 to 971

        HOW DO I GET THE ERROR OF G_R????????
        
        
        
        Here I will compute the correlation time for each value of T given L. The idea is to do first a simulation with just MCS, compute the correlation lenght and put it in a csv file.
        """

        τ_H_array[t_index] = Correlation_time(H_M_matrix[:,1])
        τ_M_array[t_index] = Correlation_time(H_M_matrix[:,2])
        τ_H_2_array[t_index] = Correlation_time(H_M_matrix[:,1].^2)
        τ_M_2_array[t_index] = Correlation_time(H_M_matrix[:,2].^2)
        τ_M_4_array[t_index] = Correlation_time(H_M_matrix[:,2].^4)

        

        t_index = t_index + 1
    end 


    #Temperatures around the critical point
    for T in T_array_3


        H_M_matrix = zeros(1000,2) #matrix we will store each measurment of H and M for each value of the temperature
        # Por el momento lo dejo así
        β = 1/T
        pa = 1-exp(-2*J*β) #We need to add here the temperature


        #We also want to measure the correlation length. The correlation time will be measured apart, using just 1 MC step.
        s1_sj = zeros(1000, N)


        for _ in 1:MCS_thermalization
            #check = falses(N)
            #i = rand(1:N) #Chooses a random spin to flip it.
            i = rand(1:N)
            s_array = wolff(i, s_array, neigh, pa)
        end

        # We save the results of the measurment
        H_M_matrix[1,:] .= hamiltonian_and_M(s_array, neigh, N) #Does this correctly

        s1_sj[1, :] .=  [s_array[1] * s_array[j] for j in 1:N]


        #Now I will do this 999 times more

        for n in 2:1000

            #Here we differenciate between the region of temperatures near the critical point and the rest of the range of temperatures.
            #here we use wolff
            for _ in 1:MCS_3
                #check = falses(N)
                #i = rand(1:N)
                i = rand(1:N)
                s_array = wolff(i, s_array, neigh, pa)
            end

            #We save the results of the measurment
            H_M_matrix[n,:] .= hamiltonian_and_M(s_array, neigh, N)

            s1_sj[n, :] .=  [s_array[1] * s_array[j] for j in 1:N]
           


        end


        
        #Once we have taken 1000 measurments we measuere the mean and the variance of the energy and the magnetization
        mean_std_H_M_temperature_matrix[t_index, :] .= mean(H_M_matrix[:,1]), std(H_M_matrix[:,1]), mean(H_M_matrix[:,2]), std(H_M_matrix[:,2]), mean(H_M_matrix[:,1].^2), std(H_M_matrix[:,1].^2) ,mean(H_M_matrix[:,2].^2), std(H_M_matrix[:,2].^2), T_array[t_index], mean(H_M_matrix[:,2].^4), std(H_M_matrix[:,2].^4)
        #We also want to compute the correlation length so we compute the correlation function g(r)

        
        
        m_avg = mean(H_M_matrix[:,2])/N
        

        s1_sj_avg = vec(mean(s1_sj, dims = 1))

        g_avg = s1_sj_avg .- m_avg^2 #We have now a vector of the value of g_r for j∈[1,N]

        # Create a DataFrame with r and g_avg
        df = DataFrame(r = r, g_avg = g_avg)

        # Group by r and calculate the average of g_avg for each group. This means creat a DataFrame with all the g_avg that correspond to the spins at distance r.
        df_grouped = combine(groupby(df, :r), :g_avg => mean)

        # Extract the unique r values and the corresponding average g_avg values
        r_unique = df_grouped.r
        g_avg_mean = abs.(df_grouped.g_avg_mean) # because I get negative numbers due to numerical errors

        p = sortperm(r_unique)
        r_unique = r_unique[p]
        g_avg_mean = g_avg_mean[p]

        # Define the exponential fitting function
        exp_func(x, p) = exp.(-x / p[1])
        fit_result = curve_fit(exp_func, r_unique, g_avg_mean, [1.0])

        # Extract the fitted parameter ξ with error
        ξ = fit_result.param[1]; σ_ξ = estimate_errors(fit_result)[1]

        correlation_matrix[t_index,1] = ξ
        correlation_matrix[t_index,2] = σ_ξ




        """
        Look at main luis lines 920 to 971

        HOW DO I GET THE ERROR OF G_R????????
        
        
        
        Here I will compute the correlation time for each value of T given L. The idea is to do first a simulation with just MCS, compute the correlation lenght and put it in a csv file.
        """

        τ_H_array[t_index] = Correlation_time(H_M_matrix[:,1])
        τ_M_array[t_index] = Correlation_time(H_M_matrix[:,2])
        τ_H_2_array[t_index] = Correlation_time(H_M_matrix[:,1].^2)
        τ_M_2_array[t_index] = Correlation_time(H_M_matrix[:,2].^2)
        τ_M_4_array[t_index] = Correlation_time(H_M_matrix[:,2].^4)



        t_index = t_index + 1
    end 


    #Low temperatures part
    for T in T_array_4


        H_M_matrix = zeros(1000,2) #matrix we will store each measurment of H and M for each value of the temperature
        # Por el momento lo dejo así
        β = 1/T
        pa = 1-exp(-2*J*β) #We need to add here the temperature


        #We also want to measure the correlation length. The correlation time will be measured apart, using just 1 MC step.
        s1_sj = zeros(1000, N)


        for _ in 1:MCS_thermalization
            #check = falses(N)
            #i = rand(1:N) #Chooses a random spin to flip it.
            i = rand(1:N)
            s_array = wolff(i, s_array, neigh, pa)
        end

        # We save the results of the measurment
        H_M_matrix[1,:] .= hamiltonian_and_M(s_array, neigh, N) #Does this correctly

        s1_sj[1, :] .=  [s_array[1] * s_array[j] for j in 1:N]

        #Now I will do this 999 times more

        for n in 2:1000

            #Here we differenciate between the region of temperatures near the critical point and the rest of the range of temperatures.
            #here we use wolff
            for _ in 1:MCS_4
                #check = falses(N)
                #i = rand(1:N)
                i = rand(1:N)
                s_array = wolff(i, s_array, neigh, pa)
            end

            #We save the results of the measurment
            H_M_matrix[n,:] .= hamiltonian_and_M(s_array, neigh, N)

            

            s1_sj[n, :] .=  [s_array[1] * s_array[j] for j in 1:N]

        end

        #Once we have taken 1000 measurments we measuere the mean and the variance of the energy and the magnetization
        mean_std_H_M_temperature_matrix[t_index, :] .= mean(H_M_matrix[:,1]), std(H_M_matrix[:,1]), mean(H_M_matrix[:,2]), std(H_M_matrix[:,2]), mean(H_M_matrix[:,1].^2), std(H_M_matrix[:,1].^2) ,mean(H_M_matrix[:,2].^2), std(H_M_matrix[:,2].^2), T_array[t_index], mean(H_M_matrix[:,2].^4), std(H_M_matrix[:,2].^4)
        #We also want to compute the correlation length so we compute the correlation function g(r)

        m_avg = mean(H_M_matrix[:,2])/N
        

        s1_sj_avg = vec(mean(s1_sj, dims = 1))

        g_avg = s1_sj_avg .- m_avg^2 #We have now a vector of the value of g_r for j∈[1,N]

        # Create a DataFrame with r and g_avg
        df = DataFrame(r = r, g_avg = g_avg)

        # Group by r and calculate the average of g_avg for each group. This means creat a DataFrame with all the g_avg that correspond to the spins at distance r.
        df_grouped = combine(groupby(df, :r), :g_avg => mean)

        # Extract the unique r values and the corresponding average g_avg values
        r_unique = df_grouped.r
        g_avg_mean = abs.(df_grouped.g_avg_mean) # because I get negative numbers due to numerical errors

        p = sortperm(r_unique)
        r_unique = r_unique[p]
        g_avg_mean = g_avg_mean[p]

        # Define the exponential fitting function
        exp_func(x, p) = exp.(-x / p[1])
        fit_result = curve_fit(exp_func, r_unique, g_avg_mean, [1.0])

        # Extract the fitted parameter ξ with error
        ξ = fit_result.param[1]; σ_ξ = estimate_errors(fit_result)[1]

        correlation_matrix[t_index,1] = ξ
        correlation_matrix[t_index,2] = σ_ξ



        """
        Look at main luis lines 920 to 971

        HOW DO I GET THE ERROR OF G_R????????
        
        
        
        Here I will compute the correlation time for each value of T given L. The idea is to do first a simulation with just MCS, compute the correlation lenght and put it in a csv file.
        """

        τ_H_array[t_index] = Correlation_time(H_M_matrix[:,1])
        τ_M_array[t_index] = Correlation_time(H_M_matrix[:,2])
        τ_H_2_array[t_index] = Correlation_time(H_M_matrix[:,1].^2)
        τ_M_2_array[t_index] = Correlation_time(H_M_matrix[:,2].^2)
        τ_M_4_array[t_index] = Correlation_time(H_M_matrix[:,2].^4)



        t_index = t_index + 1
    end
    #Now I save all the data in a dataframe

    
    df_results = DataFrame(:H => mean_std_H_M_temperature_matrix[:,1], :std_H => mean_std_H_M_temperature_matrix[:,2], :M => mean_std_H_M_temperature_matrix[:,3], :std_M => mean_std_H_M_temperature_matrix[:,4], :H_2 => mean_std_H_M_temperature_matrix[:,5], :std_H_2 => mean_std_H_M_temperature_matrix[:,6], :M_2 => mean_std_H_M_temperature_matrix[:,7], :std_M_2 => mean_std_H_M_temperature_matrix[:,8], :T => mean_std_H_M_temperature_matrix[:,9], :M_4 => mean_std_H_M_temperature_matrix[:,10], :std_M_4 => mean_std_H_M_temperature_matrix[:,11])

    
    
    file_path =joinpath(current_dir, "simulation_data_$(L)_ising_2d.csv")
    
    CSV.write(file_path, df_results)

    df_g_r = DataFrame(correlation_matrix, :auto)

    file_path = joinpath(current_dir, "correlation_function_$(L)_ising_2d.csv")

    CSV.write(file_path, df_g_r)

    df_r = DataFrame(r = r)

    file_path = joinpath(current_dir, "distance_r_not_unique_$(L)_ising_2d.csv")

    CSV.write(file_path, df_r)


    
    

    """
    Here I will create a df to save the correlation time given τ_H and τ_M for each value of the temperature given the size of the lattice.
    """
    
    df_τ = DataFrame(:τ_H => τ_H_array, :τ_M => τ_M_array, :τ_H_2 => τ_H_2_array, :τ_M_2 => τ_M_2_array, :T => T_array, :τ_M_4 => τ_M_4_array)


    file_path = joinpath(current_dir, "correlation_time_L_$(L)_ising_2d.csv")

    CSV.write(file_path, df_τ)
    


end

# This is a recursive algorithm that follows the wolff method. We start on a random site on the lattice and flip the spin.
# Then we look at its neighbours, and the neighbours of the neighbours and so on. The probability of going to neighbout spin is pa.

#I should a counter to know how many spins have been flipped.

#He escondido count

"""
Mirar si poner un visited en wolff

function wolff(i, s_array, neigh, pa, check)
    check[i] = true
    s_array[i] = -s_array[i]
    for k in 1:4
        j = neigh[k,i]
        if s_array[j] == -s_array[i] && !check[j] && rand() < pa
            s_array= wolff(j, s_array, neigh, pa, check) 
        end
    end
    return s_array 
end
"""
function wolff(i, s_array, neigh, pa)
    N = length(s_array)
    stack = [i]  # Initialize stack with the starting spin
    visited = falses(N)  # Track visited spins
    visited[i] = true  # Mark the starting spin as visited
    s_array[i] = -s_array[i]  # Flip the starting spin

    while !isempty(stack)
        current = pop!(stack)
        for k in 1:4  # Check all 4 neighbors
            neighbor = neigh[k, current]
            if s_array[neighbor] == -s_array[current] && !visited[neighbor] && rand() < pa
                s_array[neighbor] = -s_array[neighbor]  # Flip the neighbor spin
                push!(stack, neighbor)  # Add neighbor to the stack
                visited[neighbor] = true  # Mark neighbor as visited
            end
        end
    end

    return s_array
end





#To compute the hamiltonian 
function hamiltonian_and_M(s_array , neigh, N)
    energy = 0
    M = 0
    for i in 1:N
        M = M + s_array[i] #computing the magnetization
        for j in 1:4
            energy = energy - s_array[i]*s_array[neigh[j,i]] #computing the energy
        end
    end

    energy = (J/2)*energy
    M = abs(M)
    return energy, M
end

#Here we compute the correlation time. We will use to know how many Monte Carlo steps are needed to make the measuremnts uncorrelated, if you could say so. 
function Correlation_time(Observable)
    # Correlation function 
    ρg = autocor(Observable, 1:1)
    τeq = ρg[1]/(1-ρg[1])
    return τeq
end


function metropolis(s_array, neigh, T)
    N = length(s_array) # Number of spins

    # we compute the acceptance probability for the possible values of bi and store them in a list
    h = [exp(-2 * j / T) for j in [2,4]] # only for j = 2, 4 because for -4, -2, 0 we always accept the flip

    for _ in 1:(N) # Note: We are doing just one MCS
        i = rand(1:N)
        bi = 0
        for j in 1:4
            bi = bi + s_array[i] * s_array[neigh[j,i]]
        end
        if bi <= 0 # if bi <= 0, then exp(-2 * bi / T) >= 1, so min(1, exp(-2 * bi / T)) = 1; we always accept the flip
            s_array[i] = -s_array[i]
        else # h = min(1, exp(-2 * bi / T)) = exp(-2 * bi / T)
            if rand() < h[Int(div(bi, 2))] # avoiding to calculate exp(-2 * bi / T) each time
                s_array[i] = -s_array[i]
            end
        end
    end
    return s_array
end



#A function to compute the distance between each site of the lattice with the site (1,1), having into account that we are with periodic boundary conditions in a square lattice
function r1k(L)
    ds = []
    for j in 1:L
        for i in 1:L
            if i != 1 && j != 1 # If not in horizontal or vertical border, diagonal distance
                d = sqrt((i-1)^2 + (j-1)^2)
            else # If in horizontal or vertical border, minimum distance
                d_non_BC = max(abs(i-1),abs(j-1))
                d = min(d_non_BC, L-d_non_BC)
            end
            push!(ds, d)
        end
    end
    return ds
end



"""

First I should do a simulation with 1 MCS to compute the value of the correlation time of M and H for each size L so I can do a better approximation of what montecarlo I should take.

"""

L_list = [64,128]
for L in L_list 
    """
    Me falta 64. Ahora como lo tengo las gráficas no quedan aesthetic, así que quizás lo mejor sea reducir el número de intervalos
    de timpo para la representación de la gráfica pero los cálculos hacerlos con los csv con más presción.
    Simplemente guardarlos aparte y ya.

    """


    N = Int(L*L)



    s_array = zeros(L*L) #array storing the value of the spin for each site

    simulation(L,N,s_array,MCS_1, MCS_2, MCS_3 , MCS_4, MCS_thermalization,T_array)
end


