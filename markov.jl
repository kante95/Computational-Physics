using  Gadfly
using Cairo

my_theme = Theme(
    background_color=colorant"white",
    default_color=colorant"orange"
)

Gadfly.push_theme(:dark)

random(from,to) = rand(Float64)*(to-from) + from

function find_delta()
	delta = 2
	i = 0
	acceptance = 0
	acceptance_ratio = 0
	while acceptance_ratio<45 || acceptance_ratio >55
		x = random(-2,2)
		d = random(-delta/2,delta/2) 
		p = min(1,e^(x^2/2 -(x+d)^2/2))
		xi = rand(Float64)
		if xi < p
			x += d
			acceptance+=1
		end
		i+=1
		if i%100000==0
			acceptance_ratio = acceptance/1000
			acceptance = 0
			delta = delta*(acceptance_ratio/50)
		end
	end
	return delta
end

acceptance = 0
N = 1000000
x = zeros(N)

delta = find_delta()
x[1] = random(-delta,delta)
V = x[1]^2/2
xmedio = x[1]
x2medio = x[1]^2
dxmedio = x[1]^2
dx2medio = x[1]^4

for i in 2:N 
	d = random(-delta/2,delta/2) 
	p = min(1,e^(x[i-1]^2/2 -(x[i-1]+d)^2/2))
	xi = rand(Float64)
	if xi < p
		x[i] = x[i-1]+d
		acceptance+=1
	else
		x[i] = x[i-1]
	end	
	V+=x[i]^2/2
	xmedio += x[i]
	x2medio += x[i]^2
	dxmedio += x[i]^2
	dx2medio += x[i]^4
end
V = V/N
xmedio = xmedio/N
x2medio = x2medio/N
dxmedio = sqrt((dxmedio/N-xmedio^2)/N)
dx2medio = sqrt((dx2medio/N-x2medio^2)/N)
acceptance_ratio = acceptance/N
println("Valore medio x: $xmedio +/- $dxmedio, valore medio x^2: $x2medio +/- $dxmedio acceptance ratio: $(acceptance_ratio*100)%")

p = plot(x=x, Geom.histogram,my_theme)
draw(PNG("histogram1.png", 17inch, 10inch), p)

p = plot(x=x.^2, Geom.histogram,my_theme)
draw(PNG("histogram2.png", 17inch, 10inch), p)

l = 50
c = zeros(l)

D = x2medio - xmedio^2

for i in 1:l
	j = 1
	while j+i<N
		c[i]+= 	x[j]*x[j+i]
		j+=1
	end
	c[i]/=j-1	
	c[i] = (c[i] - xmedio^2)/D
end

p = plot(x=collect(1:l),y=c, Geom.point, Geom.line,Guide.ylabel("Ci"),Guide.xlabel("i"),my_theme)
draw(PNG("autocorrelation.png", 17inch, 10inch), p)
