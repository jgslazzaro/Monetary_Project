
include("functions.jl")
#Uncomment the section below to test this code

#Defining parameters
α = 0.36
β = 0.92
δ = .1
σ = .5
μ = 1.0
r = 0.1
w =  .450
O = .01#-Inf #Outside option

1/(1+r)
#Cash in hand grid:
nX = 300 #number of Cash in Hand gridpoints
X = range(-10.0,stop = 50, length = nX)

policy_k,policy_b,policy_c,policy_def,V,Vgrid,policygrid,xbar,q = monetary(X,O,w,σ)


nstar(k,z;α=α,w=w) = (z*(1-α)/w)^(1/α)*k
X2(k::Float64,b::Float64,z;w::Float64=w,δ::Float64=δ,α::Float64 = α) = z*k^α * nstar(k,z)^(1-α)-w*nstar(k,z) + (1-δ)*k - b



using Plots
plot(X,policy_k.(X),label = "Capital policy") #log utility, exogenous labor case
plot(X,policy_b.(X),label = "Debt policy")
plot(X,policy_c.(X),label = "Consumption policy")

plot(X,V.(X),label = "Value Function")
plot(X,q.(policy_k.(X),policy_b.(X),xbar))
plot(X,[x.(policy_k.(X),policy_b.(X),6.) x.(policy_k.(X),policy_b.(X),.5)],label=["X, z=1" "X, z=0.5"])
#plot(A[1]:0.05:A[end],[policy_a.(A[1]:0.05:A[end],1) policy_a.(A[1]:0.05:A[end],0) A[1]:0.05:A[end]],label =["Employed", "Unemployed","45"],legend = :bottomright)


pdef = 1 .-q.(policy_k.(X),policy_b.(X),xbar)*(1+r)



#plot(A,[policy_a.(A) Z*α*β.*(A).^α A])
