#Primo esercizio, calcolo di π

points = 10000000
random(from,to) = rand(Float64)*(to-from) + from

N = 0
for i=1:points
    x = random(-0.5,0.5)
    y = random(-0.5,0.5)
    d = √(x^2+y^2)
    #print("$x $y $d\n")
    if d<0.5
        N +=1
    end
end
π = 4N/(points)
print("π calcolato: $π, vero: $pi\n")


#Secondo esercizio, integrale di 1/2π x^2y^2 e^{-x^2+y^2}/2

#normal_rand(square_std) = √(-2(square_std)*log(1-rand(Float64)))*cospi(2rand(Float64))
#f(x) = x^2

#I = 0
#dI = 0
#for i=1:points
#    x = normal_rand(1)
#    I += f(x)
#    dI += f(x)^2
#end
#I = (I/points)^2
#dI = sqrt((dI/points-I^2)/points)

#print("Valore di 1/(2π)∫x^2y^2e^{-(x^2+y^2)/2}: $I ± $dI, vero: 1\n")


#################################
normal_rand(square_std) = √(-2(square_std)*log(1-rand(Float64)))*cospi(2rand(Float64))
f(x,y) = x^2*y^2
I = 0
dI = 0
for i=1:points
    x = normal_rand(1)
    y = normal_rand(1)
    I += f(x,y)
    dI += f(x,y)^2
end
I = (I/points)^2
dI = sqrt((dI/points-I^2)/points)

print("Valore di 1/(2π)∫x^2y^2e^{-(x^2+y^2)/2}: $I ± $dI, vero: 1\n")





#Terzo esercizo, integrale e^x-1 tra 0 e 1 con densità uniforme e lineare

g(x) = e^x-1
I = 0
dI = 0
for i=1:points
    x = rand(Float64)
    I += g(x)
    dI += g(x)^2
end
I = (I/points)
dI = sqrt((dI/points-I^2)/points)
tru = e-2
print("Valore di ∫e^x-1: $I ± $dI, vero: $tru densità uniforme\n")


uniform_rand() = √rand(Float64)
h(x) = (e^x-1)/(2x)
I = 0
dI = 0
for i=1:points
    x = uniform_rand()
    I += h(x)
    dI += h(x)^2
end
I = (I/points)
dI = sqrt((dI/points-I^2)/points)
tru = e-2
print("Valore di ∫e^x-1: $I ± $dI, vero: $tru densità lineare\n")

