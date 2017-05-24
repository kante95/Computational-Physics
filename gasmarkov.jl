using Gadfly
using Cairo

my_theme = Theme(
   background_color=colorant"white",
    default_color=colorant"orange"
)


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
            dist = (r[1]^2+r[2]^2 +r[3]^2)^(1/2)
           
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


    while acceptance_ratio<35 || acceptance_ratio > 65
        positions = initialize(num_particles,L)
        potential = calculate_potential(positions,num_particles,L)
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
        if steps%1000==0
            acceptance_ratio = acceptance/10
            #println("Nuovo delta: $delta, $acceptance_ratio")
            if(acceptance==0)
                delta = delta/10
                acceptance = 0
            else
                acceptance = 0
                delta = delta*(acceptance_ratio/50)
            end
            #println("Nuovo delta: $delta, $acceptance_ratio")
        end 
    end
    return delta
end


function simulation(num_particles, density, M,nrank)

    temp_simul = 2
    #########################
    temp = linspace(0.5,3,50)
    Vt = zeros(length(temp))
    Vt2 = zeros(length(temp))
    βs = 1./temp
    normalization = zeros(length(temp))
    #######################

    V = num_particles/density
    L = ∛V
    println("Processore $nrank: Attendi, trovo il migliore delta....")
    #delta = -0.23density + 0.0315
    delta = 0.05/density^(1/3)#find_delta(temp_simul,num_particles,L)
    #delta = 0.03
    println("Processore $nrank: Delta migliore trovato: $delta, adesso inizio la simulazione")
    acceptance = 0

    β = 1/temp_simul
    #inizializzazione
    positions = initialize(num_particles,L)   
    potential = calculate_potential(positions,num_particles,L)
    
    Vt_simul =  potential
    Vt2_simul =  potential^2

    #Vprova = potential
    
    for t in 1:length(temp)
        Vt[t] =  potential*e^((β - βs[t])Vt_simul)
        Vt2[t] =  (potential^2)*e^((β - βs[t])Vt_simul)
        normalization[t] = e^((β - βs[t])Vt_simul)
    end

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
        Vt_simul += potential
        Vt2_simul+= potential^2
        #Vprova += potential
        for t in 1:length(temp)
            Vt[t] +=  potential*e^((β - βs[t])potential)
            Vt2[t] +=  (potential^2)*e^((β - βs[t])potential)
            normalization[t] += e^((β - βs[t])potential)
        end
        println(k)

    end

    for t in 1:length(temp)
            Vt[t] /=  normalization[t]
            Vt2[t] /=  normalization[t]       
    end
    Vt_simul/=M
    Vt2_simul/=M

    println("Poteniziali simulazione: $Vt_simul, $Vt2_simul, Potenziale medio: $(Vt), potenziale medio quadro: $(Vt2), acceptance ratio: $(acceptance*100/M)")

    Cv = zeros(length(temp))
    for i ∈ 1: length(temp)
        Cv[i] = (Vt2[i]-Vt[i]^2)/(temp[i]^2)
    end
    #println("Cv: $Cv, Cv simulato: $((Vt2_simul-Vt_simul^2)/(temp_simul^2))")

    #plot
    #p = plot(x=temp,y=Cv,Geom.point, Geom.line,Guide.ylabel("Cv"),Guide.xlabel("Temperatura"),my_theme)
    #draw(PNG("graphs/Cv$(density).png", 17inch, 10inch), p)


    cv_max = findmax(Cv)
    temp_max = temp[cv_max[2]]

    println("Processore $nrank: Temperatura massima: $temp_max, $cv_max")

    return temp_max
end

simulation(125,0.3,50000,1)