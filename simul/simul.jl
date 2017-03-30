#Box-Muller formula
normal_rand(std) = √(-2(std^2)*log(1-rand(Float64)))*cospi(2rand(Float64))

function initialize(N,L,temperature)
    i = ∛N
    positions = zeros(N,3)
    velocities = zeros(N,3)
    index = 1
    for j ∈ 0:L/(i-1):i+1
        for k ∈ 0:L/(i-1):i+1
            for m ∈ 0:L/(i-1):i+1
                positions[index,1]= j - L/2
                positions[index,2] = k - L/2
                positions[index,3] = m - L/2
                velocities[index,1]= normal_rand(√temperature)
                velocities[index,2] = normal_rand(√temperature)
                velocities[index,3] = normal_rand(√temperature)
                index+=1
            end
        end
    end
    #RANDOM PARTICLES! to fix
    #for j ∈ 0:N-i^3
    #    append!(position,rand(-L/2:L/2,3))
    #end

    return positions,velocities
end

lennard_jones(r) =  4(r^-12 - r^-6)

function lennard_jones_derivative!(r,returning)
    den = r[1]^2+r[2]^2 +r[3]^2
    returning[1] = -4(-12r[1]*den^-7 + 6r[1]*den^-4)
    returning[2] = -4(-12r[2]*den^-7 + 6r[2]*den^-4)
    returning[3] = -4(-12r[3]*den^-7 + 6r[3]*den^-4)
end

function calculate_accellerations_and_potential(positions,N,L)
    accelerations = zeros(positions)
    potentials = zeros(125)
    force = Array{Float64}(3)
    for i ∈ 1:N
        for j ∈ i+1:N
            r = positions[i,:]-positions[j,:]-round((positions[i,:]-positions[j,:])/L)
            dist = √(r[1]^2+r[2]^2 +r[3]^2)
            if dist<L/2
                lennard_jones_derivative!(r,force)
                accelerations[i,:] += force
                accelerations[j,:] -= force
                potentials[i] += lennard_jones(dist)
                potentials[j] += lennard_jones(dist)
            end
        end
    end
    return accelerations,potentials
end


function simulation(num_particles, density, temperature, time, step)
    V = num_particles/density
    L = ∛V
    print("Inizializzo lo stato iniziale....\n")
    positions,velocities = initialize(num_particles, L,temperature)
    #velocity verlet algorithm
    current_accelerations,current_potentials = calculate_accellerations_and_potential(positions,num_particles,L)
    for t ∈ step:step:time
        #print("$positions")
        #print("Simulando secondo "+str(t)+" s")
        #calcolo le nuove posizioni
        #print("$shape(positions) $shape(velocities)")
        positions = positions + velocities*step + 0.5current_accelerations*step^2
        positions = positions - L*round(positions/L)
        #accelerazioni dovute alle nuove posizioni
        new_accelerations,new_potentials = calculate_accellerations_and_potential(positions,num_particles,L)
        #nuove velocita che si calcolano con le vecchie e le nuove accelerazioni
        velocities = velocities + 0.5(current_accelerations + new_accelerations)step
        #aggiorno le accelerazioni per il prossimo ciclo
        current_accelerations = new_accelerations
        current_potentials = new_potentials
        energy = 0.5sum(velocities[:][1]^2+velocities[:][2]^2+velocities[:][3]^2) + sum(current_potentials)
        print("Energy at second $t: $energy \n")
    end
        #average =  what da fuck have i to put here??
    #time average
    #P = 
    #return P
end

#for ρ ∈ [0.1,0.01,0.001]
    simulation(125,0.01,273,5.0,0.001)
#end
