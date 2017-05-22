import MPI

include("gasmarkov.jl")

ρ = linspace(0.01,0.3,12)
temp = zeros(length(ρ))

MPI.Init()
comm = MPI.COMM_WORLD

ncore = MPI.Comm_size(comm)
nrank = MPI.Comm_rank(comm)

nstart = nrank*length(ρ)/ncore+1
nend = nstart + length(ρ)/ncore-1

for i=nstart:nend
    intero = convert(Int,i)
    temp[intero] = simulation(125,ρ[intero],50000,nrank)
    println("Processore: $nrank, densità: $(ρ[intero]),temperatura: $(temp[intero])\n")
end

MPI.Barrier(comm)


p = plot(x=ρ,y=temp,Geom.point, Geom.line,Guide.ylabel("Temperatura"),Guide.xlabel("Densità"),my_theme)
    draw(PNG("graphs/Finale2.png", 17inch, 10inch), p)


println(temp)

MPI.Finalize()
 
