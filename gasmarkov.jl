#using  Gadfly
#using Cairo

#my_theme = Theme(
 #   background_color=colorant"white",
#    default_color=colorant"orange"
#)


function initialize(N,L)
    i = ∛N
    positions = zeros(N,3)
    index = 1
    for j ∈ 0:L/(i):L-L/(i)
        for k ∈ 0:L/(i):L-L/(i)
            for m ∈ 0:L/(i):L-L/(i)
                positions[index,1]= j - L/2
                positions[index,2] = k - L/2
                positions[index,3] = m - L/2
                index+=1
            end
        end
    end
    return positions
end
random(from,to) = rand(Float64)*(to-from) + from 

normal_rand(std) = √(-2(std^2)*log(1-rand(Float64)))*cospi(2rand(Float64))

lennard_jones(r) =  4(r^-12 - r^-6)

function calculate_potential(positions,N,L)
    potential = 0
    for i ∈ 1:N
        for j ∈ i+1:N
            r = positions[i,:]-positions[j,:]-round((positions[i,:]-positions[j,:])/L)
            dist = (r[1]^2+r[2]^2 +r[3]^2)^2
           
            if dist<L/2 && dist !=0
                potential += 2lennard_jones(dist)
            end

        end
    end
    return potential
end

function find_delta(temperature,num_particles,L)
    delta = 0.1
    steps = 0
    acceptance = 0
    acceptance_ratio = 0
    β = 1/temperature

    positions = initialize(num_particles,L)
    potential = calculate_potential(positions,num_particles,L)

    while acceptance_ratio<45 || acceptance_ratio > 55
        new_positions = zeros(num_particles,3)
        for i ∈ 1:3
            for j ∈ 1:num_particles
                new_positions[j,i] = positions[j,i]+random(-delta/2,delta/2)
                new_positions[j,i] = new_positions[j,i]-round(new_positions[j,i]/L)
            end
        end

        new_potential = calculate_potential(new_positions,num_particles,L)
        p = min(1,e^(-β*(new_potential - potential)))
        xi = rand(Float64)
        if xi < p
            positions = new_positions
            potential = new_potential
            acceptance+=1
        end
        steps+=1
        if steps%10000==0
            acceptance_ratio = acceptance/100
            acceptance = 0
            delta = delta*(acceptance_ratio/50)
        end 
    end
    return delta
end


function simulation(num_particles, density, M)

    #########################
    #temp = linspace(1.5,4,5)
    #Vt = zeros(length(temp))
    #Vt2 = zeros(length(temp))
    #βs = 1./temp
    #normalization = zeros(length(temp))
    #######################

    temp = 1.5
    V = num_particles/density
    L = ∛V
    println("Attendi, trovo il migliore delta....")
    delta = find_delta(temp[1],num_particles,L)
    println("Delta migliore trovato: $delta, adesso inizio la simulazione")
    acceptance = 0

    β = 1/temp
    #inizializzazione
    positions = initialize(num_particles,L)
   
    potential = calculate_potential(positions,num_particles,L)
    
    #Vt[1] =  potential
    #Vt2[1] =  potential^2

    Vprova = potential
    
   # for t in 2:length(temp)
    #    Vt[t] =  potential*e^((β - βs[t])Vt[1])
    #    Vt2[t] =  (potential^2)*e^((β - βs[t])Vt[1])
    #    normalization = e^((β - βs[t])Vt[1])
    #end

    #loop principale della catena
    for k ∈ 1:M

        #nuove posizioni
        new_positions = zeros(num_particles,3)
        for i ∈ 1:3
            for j ∈ 1:num_particles
                new_positions[j,i] = positions[j,i]+random(-delta/2,delta/2)
                new_positions[j,i] = new_positions[j,i]-round(new_positions[j,i]/L)
            end
        end

        new_potential = calculate_potential(new_positions,num_particles,L)
        p = min(1,e^(-β*(new_potential - potential)))
        xi = rand(Float64)
        if xi < p
            positions = new_positions
            potential = new_potential
            acceptance+=1
        end
       # Vt[1] += new_potential
        #Vt2[1]+= new_potential^2
        Vprova += new_potential
        #for t in 2:length(temp)
        #    Vt[t] +=  potential*e^((β - βs[t])new_potential)
        #    Vt2[t] +=  (potential^2)*e^((β - βs[t])new_potential)
        #    normalization += e^((β - βs[t])new_potential)
       # end

    end

    #for t in 2:length(temp)
    #        Vt[t] /=  normalization 
    #        Vt2[t] /=  normalization 
    #       
    #end
    #Vt[1]/=M
    #Vt2[1]/=M
    println("Potenziale medio: $(Vprova/M) acceptance ratio: $(acceptance*100/M)")

end

simulation(125,0.1,100000)