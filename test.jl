using Plots


include("functions.jl")


ρ = 0.0
σ = 2
Y = 200
pdfY, y = Tauchen(ρ,σ,Y)


sim = simMC(y,pdfY,500,y[100])

plot(sim)


mean(sim)
