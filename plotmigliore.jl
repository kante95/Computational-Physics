using Gadfly
using Cairo
using DataFrames

my_theme = Theme(
   background_color=colorant"white",
    default_color=colorant"orange"
)


#filename = ["0.01.csv","0.04.csv","0.07.csv","0.025.csv","0.055.csv","0.085.csv","0.0175.csv","0.0325.csv","0.0475.csv","0.0625.csv","0.0775.csv","0.0925.csv","0.1.csv"]
cd("datifinali")
filename = readdir()
#print(filename)
density = [0.505000,0.752500,0.257500,0.010000,0.787857,0.540357,0.292857,0.045357,0.823214,0.575714,0.328214,0.080714,0.858571,0.611071,0.363571,0.116071,0.893929,0.646429,0.398929,0.151429,0.929286,0.681786,0.434286,0.186786, 0.964643,0.717143,0.469643,0.222143]
temp =    [2.100000,1.710000,2.100000,1.500000,2.100000,2.190000,2.130000,2.280000,2.070000,2.100000,2.130000,2.190000,2.040000,1.680000,2.160000,2.130000,2.040000,1.680000,2.130000,2.100000,2.040000,1.740000,2.130000,2.130000,2.040000,1.740000,2.070000,2.130000 ]

                                                                                                                                                                                                                                            

df = DataFrame()

for i in filename
    df1= readtable("$i",header = false)
    df2 = DataFrame(Densità = repeat([i[8:end-4]], inner = [50]))
    df1 = hcat(df1,df2)
    #print(names(df1))
    #df2 = DataFrame(x=1:10, y=g(1:10), label="g(x)")
    df = vcat(df, df1)
end

cd("..")
p = plot(df, x=:x2, y=:x1, Geom.line, Geom.point,Guide.ylabel("Cv"),Guide.xlabel("Temeperatura"), color=:Densità,my_theme)
draw(PNG("Plotcv.png", 17inch, 10inch), p) 

p = plot(x=density,y=temp,Geom.point, Geom.line,Guide.ylabel("Temperatura"),Guide.xlabel("Densità"),my_theme)
draw(PNG("Tempcrit.png", 17inch, 10inch), p) 
