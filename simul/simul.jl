using ArgParse

#Box-Muller formula
normal_rand(std) = √(-2(std^2)*log(1-rand(Float64)))*cospi(2rand(Float64))

function initialize(N,L,temperature)
    i = ∛N
    positions = zeros(N,3)
    velocities = zeros(N,3)
    index = 1
    for j ∈ 0:L/(i):L-L/(i)
        for k ∈ 0:L/(i):L-L/(i)
            for m ∈ 0:L/(i):L-L/(i)
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
lennard_jones_derivative_scalar(r) = -48r^-13 +24r^-7

function lennard_jones_derivative!(r,returning)
    den = r[1]^2+r[2]^2 +r[3]^2
    returning[1] = -4(-12r[1]*den^-7 + 6r[1]*den^-4)
    returning[2] = -4(-12r[2]*den^-7 + 6r[2]*den^-4)
    returning[3] = -4(-12r[3]*den^-7 + 6r[3]*den^-4)
end

function calculate_accellerations_and_potential_and_partial_pressure(positions,N,L)
    accelerations = zeros(positions)
    potentials = zeros(125)
    force = Array{Float64}(3)
    partial_pressure = 0    
    for i ∈ 1:N
        for j ∈ i+1:N
            r = positions[i,:]-positions[j,:]-round((positions[i,:]-positions[j,:])/L)
            dist = (r[1]^2+r[2]^2 +r[3]^2)^2
           
            if dist<L/2 && dist !=0
                lennard_jones_derivative!(r,force)
                accelerations[i,:] += force
                accelerations[j,:] -= force
                potentials[i] += lennard_jones(dist)
                potentials[j] += lennard_jones(dist)
            end
            partial_pressure += lennard_jones_derivative_scalar(dist)*dist

        end
    end
    return accelerations,potentials, partial_pressure
end


function simulation(num_particles, density, temperature, time, step)
    V = num_particles/density
    L = ∛V
    pressure = 0
    #print("Inizializzo lo stato iniziale....\n")
    positions,velocities = initialize(num_particles, L,temperature)
    #velocity verlet algorithm
    current_accelerations,current_potentials,partial_pressure = calculate_accellerations_and_potential_and_partial_pressure(positions,num_particles,L)
    pressure += partial_pressure
    for t ∈ step:step:time
        #print("$positions\n")
        #print("Simulando secondo "+str(t)+" s")
        #calcolo le nuove posizioni
        #print("$shape(positions) $shape(velocities)")
        positions = positions + velocities*step + 0.5current_accelerations*step^2
        positions = positions - L*round(positions/L)
        #accelerazioni dovute alle nuove posizioni
        new_accelerations,new_potentials,partial_pressure = calculate_accellerations_and_potential_and_partial_pressure(positions,num_particles,L)
        pressure += partial_pressure
        #nuove velocita che si calcolano con le vecchie e le nuove accelerazioni
        velocities = velocities + 0.5(current_accelerations + new_accelerations)step
        #aggiorno le accelerazioni per il prossimo ciclo
        current_accelerations = new_accelerations
        current_potentials = new_potentials
        energy = 0.5sum(velocities[:][1]^2+velocities[:][2]^2+velocities[:][3]^2) + sum(current_potentials)
        print("Energy at second $t: $energy, Partial pressure: $pressure\n")
    end
    return density*temperature - (1/3)*pressure/length(0:step:time)    
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--temp"
            help = "Temperature to simulate, if not provided a range of temperature will be simulated"
            arg_type = Int
            default = 0
        "--density"
            help = "Numbers of density to simulate"
            arg_type = Int
            default = 20
        "--timestep"
            help = "Time step"
            arg_type = Float64
            default = 0.00001
        "--step"
            help = "Number of step to simulate"
            arg_type = Int
            default = 100000

    end
    return parse_args(s)
end


parsed_args = parse_commandline()

if parsed_args["temp"]!=0
    temp = parsed_args["temp"] 
else
    temp = logspace(0,3,8)
end

ρ = logspace(-2,0,parsed_args["density"])
pression = zeros(Float64,length(ρ))
volume = zeros(Float64,length(ρ))

for t∈temp
    for i=1:length(ρ)
        pression[i] = simulation(125,ρ[i],t,parsed_args["step"]*parsed_args["timestep"],parsed_args["timestep"])
        volume[i] = 125/ρ[i]
        c = ρ[i]
        print("Sono arrivato ad temp: $t, desità: $c\n")
        writedlm("/tmp/simulation$t.txt",hcat(pression,volume), ",")
    end
end

