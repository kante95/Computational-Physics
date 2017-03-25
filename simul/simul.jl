function initialize(N,L)
    i = ∛N
    positions = zeros(N,3)
    index = 1
    for j ∈ 0:L/(i-1):i+1
        for k ∈ 0:L/(i-1):i+1
            for m ∈ 0:L/(i-1):i+1
                positions[index,1]= j - L/2
                positions[index,2] = k - L/2
                positions[index,3] = m - L/2
                index+=1
            end
        end
    end
    #RANDOM PARTICLES! to fix
    #for j ∈ 0:N-i^3
    #    append!(position,rand(-L/2:L/2,3))
    #end
    velocities = zeros(positions)
    return positions,velocities
end

function lennard_jones(r)
    return 4(r^-12 - r^-6)
end

function lennard_jones_derivative(r)
    den = r[1]^2+r[2]^2 +r[3]^2
    return -[4(-14r[1]*den^-7 + 8r[1]*den^-4),4(-14r[2]*den^-7 + 8r[2]*den^-4),4(-14r[3]*den^-7 + 8r[3]*den^-4)]
end

function calculate_accellerations_and_potential(positions,N,L)
    accelerations = zeros(positions)
    potentials = zeros(125)
    for i ∈ 1:N
        for j ∈ i+1:N
            r = positions[i,:]-positions[j,:]-round((positions[i,:]-positions[j,:])/L)
            dist = √(r[1]^2+r[2]^2 +r[3]^2)
            if dist<L/2
                accelerations[i,:] += lennard_jones_derivative(r)
                accelerations[j,:] -= lennard_jones_derivative(r)
                potentials[i] += lennard_jones(dist)
                potentials[j] += lennard_jones(dist)
            end
        end
    end
    return accelerations,potentials
end

function simulation(num_particles,L, density, temperature, time, step)
    print("Inizializzo lo stato iniziale....\n")
    positions,velocities = initialize(num_particles, L)
    #mesh = collect(step,time,step)
    #energy = zeros(lenght(mesh))
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
    #plt.plot(mesh,energy)
    #plt.show()
    #time average
    #P = 
    #return P
end

simulation(125,5,1,273,5,0.01)
