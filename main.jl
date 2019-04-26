
include("functions.jl")
#Uncomment the section below to test this code

#Defining parameters
α = 0.36
β = 0.9
δ = .025
σ = .2
μ = 1.0
r = 0.02
w =  1.0
O = 3.0#-Inf #Outside option

#Cash in hand grid:
nX = 1500 #number of Cash in Hand gridpoints
X = range(-0.0,stop = 30, length = nX)


policy_k,policy_b,policy_c,policy_def,V,Vgrid,policygrid,xbar,q = monetary(X,O,w,σ)




using Plots
plot(X,policy_k.(X)) #log utility, exogenous labor case
plot(X,policy_b.(X))
plot(X,policy_c.(X))
plot(X,x.(policy_k.(X),policy_b.(X),1.0))
plot(X,V.(X))
#plot(A[1]:0.05:A[end],[policy_a.(A[1]:0.05:A[end],1) policy_a.(A[1]:0.05:A[end],0) A[1]:0.05:A[end]],label =["Employed", "Unemployed","45"],legend = :bottomright)

#plot(A,[policy_a.(A) Z*α*β.*(A).^α A])
|
