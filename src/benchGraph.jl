using Plots

pythonSeries = [0.04839897155761719, 0.21028518676757812, 1.6918182373046875, 11.152029037475586, 94.91419792175293]
juliaSeries = [0.000172202, 0.000905400, 0.010241, 0.110558, 1.183]
cSeries = [0.012, 0.054, 0.556, 5.449, 25.346]
n = [10^x for x=2:6]


plot(n, pythonSeries, seriestype = :scatter, labels = "Python", markersize = 6) 
plot!(n, cSeries, seriestype = :scatter, labels = "C", markersize = 6)
plot!(n, juliaSeries,  seriestype = :scatter, labels = "Julia", markersize = 6)
plot!(legend = :bottomright, xscale = :log , yscale = :log, ytitle = "Tempo de execução(ms)", xtitle = "Número de entradas")
savefig("bench")
