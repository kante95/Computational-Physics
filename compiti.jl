#Primo esercizio, calcolo di π

points = 10000000
function random(from,to)
    return rand(Float64)*(to-from) + from
end
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
π = N/(points*0.5*0.5)
print("π calcolato: $π, vero: $pi\n")


#Secondo esercizio, integrale di 1/2π x^2y^2 e^{-x^2+y^2}/2

normal_rand(square_std) = √(-2(square_std)*log(1-rand(Float64)))*cospi(2rand(Float64))
f(x) = x^2

I = 0
for i=1:points
    I += f(normal_rand(1))
end
I = (I/points)^2
print("Valore di 1/(2π)∫x^2y^2e^{-(x^2+y^2)/2}: $I, vero: 1\n")

#Terzo esercizo, integrale e^x-1 tra 0 e 1 con densità costante e uniforme

g(x) = e^x-1
I = 0
for i=1:points
    I += g(rand(Float64))
end
I = (I/points)
tru = e-2
print("Valore di ∫e^x-1: $I, vero: $tru densità uniforme\n")


uniform_rand() = √rand(Float64)
h(x) = (e^x-1)/(2x)
I = 0
for i=1:points
    I += h(uniform_rand())
end
I = (I/points)
tru = e-2
print("Valore di ∫e^x-1: $I, vero: $tru densità uniforme\n")

