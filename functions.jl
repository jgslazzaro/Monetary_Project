include("CRRA_utility.jl")
using Distributions,Optim, Interpolations,  FastGaussQuadrature

function Tauchen(ρ::Float64,σ::Float64,Y::Int64,μ::Float64 = 0.0,m::Float64 = 3.0)
    #This function is to discretize an AR(1) process following Tauchen(1986) method
    # y_{t+1} = μ + ρy_t + ϵ
    #ϵ~N(0,σ^2)
    #Y is the number of y states
if Y>1
    ybar = μ/(1-ρ)
    ymax= ybar + m*(σ^2/(1-ρ^2))^(1/2) #maximum y
    ymin= ybar - m*(σ^2/(1-ρ^2))^(1/2) #minimum y

    Δ = (ymax-ymin)/(Y-1)# #distance between each y
    y=ymin:Δ:ymax #vector of possible states of p

    d=Normal()

    pdfY=ones(Y,Y) #preallocate memory and create the transition matrix in the following loop
    for i in 1:Y
        pdfY[i,1]=cdf(d,(y[1] + Δ/2 -ρ*y[i]) / σ^0.5);
        pdfY[i,Y]=1-cdf(d,(y[Y] -Δ/2 - ρ*y[i]) / σ^0.5);
        for j in 2:Y-1
            pdfY[i,j]=cdf(d,(y[j] + Δ/2 - ρ*y[i])/σ^0.5) - cdf(d,(y[j] - Δ/2 - ρ*y[i]) / σ^0.5);
        end
    end
else
    y=μ
    pdfY=1.0
end

    return pdfY, y
end


function simMC(S,pdf,T,s0)
    #This function simulates a Markov chain
    #S possible states
    #pdF transition matrix of the states
    #T simulation length
    #s0 initial state
    nS = length(S)
    ssim = fill(s0, T)
    r = rand(T)
    s=1
    #Simulating the economy
    for t=2:T
        s = findfirst(S.==ssim[t-1])

        ps = pdf[s,1]
        for i=1:nS
            if r[t]<=ps
                ssim[t]=S[i]
                break
            else
                ps+=pdf[s,i+1]
            end
        end
    end
    return ssim
end


function VFI(X::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
    O::Float64,w::Float64,σ::Float64,xbar::Float64;α::Float64 = α,β::Float64 = β,
    δ::Float64=δ, μ::Float64 =μ,integranodes::Int64=50,inner_optimizer = BFGS(), tol::Float64 = 1e-6)

    nX = length(X)
    #Defining consumption function:
    c(x::Float64,k1::Float64,b1::Float64,xbar::Float64=xbar) = x-k1+q(k1,b1,xbar)*b1

    #nodes, weights = gausshermite(integranodes)#Gauss Legendre nodes and weights,this function is just a Quadrature table
    #Expected value function given states and a value function V
    #=#function EV(k::Float64,b::Float64,σ::Float64 ;nodes = nodes,weights=weights)
        #find expected value of V using Gauss-Hermite method
        expected = 0.0
        for i=1:length(nodes)
            expected += weights[i]*V(x(k,b,exp(sqrt(2.0)*σ*nodes[i])))
        end
        expected = π^(-1/2)*expected
        return expected
    end =#

    weights,nodes = Tauchen(0.0,σ,integranodes)
    nodes = exp.(nodes[1,:])
    function EV(k::Float64,b::Float64,σ::Float64 ;nodes = nodes,weights=weights)
        #find expected value of V using Tauchen discretized process
        expected = 0.0
        for i=1:length(nodes)
            expected += weights[i]*V(x(k,b,nodes[i]))
        end

        return expected
    end


    #Function to be maximized by the solver
    function Vf(S::Array{Float64,1},X::Float64,σ::Float64;β::Float64 = β,xbar::Float64=xbar)
        k1,b1 = S
        value = u(c(X,k1,b1)) + β * EV(k1,b1,σ)
        return -value
    end

    #solver stuff
    initial = [.5,.5]
    lower = [0.001,-Inf]
    upper = [Inf,Inf]
    #predefining variables and loop stuff
    policy = ones(nX,3) #last dimension indicates if it is the policy for k,b,default
    distance = 1
    Vgrid = ones(nX).*u.(X)
    Vgridf = copy(Vgrid)
    #Guess for Value function
    itp = LinearInterpolation((X),Vgridf, extrapolation_bc=Line())
    V(x::Float64) = itp(x)
    iteration =0
    while distance > tol
        #global policy, Vgrid,Vgridf,initial,iteration,distance,itp
        Vgridf = copy(Vgrid)

        for x0 = 1:nX #Parallel for!!
            maxV = optimize( s-> Vf(s,X[x0],σ), lower,upper,initial,Fminbox(inner_optimizer))
            #autodiff=:forward)
            policy[x0,1:2] = maxV.minimizer
            if -maxV.minimum >= O
                Vgrid[x0] = -maxV.minimum
                policy[x0,3] =0.0
            else
                Vgrid[x0] = O
                policy[x0,3] = 1.0
            end
            initial = policy[x0,1:2]
        end
        distance = maximum(abs.(Vgrid-Vgridf))
        itp = LinearInterpolation((X),Vgrid, extrapolation_bc=Line())
        V(x::Float64) = itp(x)
        iteration += 1
        println("In iteration $(iteration), distance is $(distance)")
    end
    #Finally, find labor, capital and consumption policies:
    itpb = LinearInterpolation((X),policy[:,2], extrapolation_bc=Line())
    policy_b(x) = itpb(x)
    itpk = LinearInterpolation((X),policy[:,1], extrapolation_bc=Line())
    policy_k(x) = itpk(x)
    policy_c(x) = c(x,policy_k(x),policy_b(x))
    function policy_def(x;O=O)
        if V(x)>=O
            return 0.0
        else
            return 1.0
        end
    end

    return  policy_k,policy_b,policy_c,policy_def,V,Vgrid,policy
end

#Cash in Hands function
x(k::Float64,b::Float64,z;w::Float64=w,δ::Float64=δ,α::Float64 = α) = (z[1]*(1-α)/w)^(1/α)*k*(z[1]*(z[1]*(1-α)/w)^(1-α)-w)-b+(1-δ)*k
#Defining the Default probability function
DefProb(k1,b1,xbar,σ) = cdf(LogNormal(0.0,σ),
    optimize( z-> abs(xbar - x(k1,b1,z)),eps(),30.0).minimizer)
#finds the probability of z such that X tommorrow is below the threshold
q(k1::Float64,b1::Float64,xbar::Float64;r::Float64=r,σ::Float64 = σ) = (1-DefProb(k1,b1,xbar,σ)) / (1+r)

function monetary(X::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
    O::Float64, w::Float64, σ::Float64; α::Float64 = α,β::Float64 = β, δ::Float64=δ,
    μ::Float64 =μ, r::Float64=r, tol::Float64 = 1e-6)


    nX=length(X)
    xbar = X[1]
    qgrid = 1/(1+r)*ones(nX)
    policy_k,policy_b,policy_c,policy_def,V,Vgrid,policygrid= VFI(X,O,w,σ,xbar;α = α,
    β = β, δ=δ, μ =μ,inner_optimizer = BFGS(), tol = 1e-6)
    xbar = X[findxbar(policygrid[:,3])] #find the threshold
    distance = maximum(abs.(q.(policy_k.(X),policy_b.(X),xbar) - qgrid))
    println("Distance is $(distance)")
    while distance >tol
        @time policy_k,policy_b,policy_c,policy_def,V,Vgrid,policygrid= VFI(X,O,w,σ,xbar;α = α,
        β = β, δ=δ, μ =μ,inner_optimizer = BFGS(), tol = 1e-6)
        qgrid = q.(policy_k.(X),policy_b.(X),xbar)
        xbar = X[findxbar(policygrid[:,3])] #find the threshold
        distance = maximum(abs.(q.(policy_k.(X),policy_b.(X),xbar) - qgrid))
        println("Distance is $(distance)")
    end
    return policy_k,policy_b,policy_c,policy_def,V,Vgrid,policygrid,xbar,q

end

function findxbar(policygrid)
    xbar = length(policygrid)
    for i=1:length(policygrid)
        if policygrid[i] == 0.0
            xbar = i
            break
        end
    end
    return xbar
end
