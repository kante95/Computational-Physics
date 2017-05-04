random(from,to) = rand(Float64)*(to-from) + from

println("Dimmi un delta:")
delta = float(readline(STDIN))
acceptance = 0
#delta = 0.5
N = 1000000
x = zeros(N)
x[1] = random(-7,7)
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
	dxmedio = x[i]^2
	dx2medio = x[i]^4
end
V = V/N
xmedio = xmedio/N
x2medio = x2medio/N
dxmedio = sqrt(abs(dxmedio/N-xmedio^2)/N)
dx2medio = sqrt(abs(dx2medio/N-x2medio^2)/N)
acceptance_ratio = acceptance/N
println("Valore medio x: $xmedio +/- $dxmedio, valore medio x^2: $x2medio +/- $dxmedio acceptance radio: $(acceptance_ratio*100)%")
