#Author: João Lazzaro
#This code defines typical CRRA utility functions


function u(c::Float64,l::Float64;η::Real = η,μ::Real = μ)
    if (c<=0) || (η<1 && l<=0)
        if μ ==1
            u=log(eps(0.)^η * l^(1-η))+ 1e5*min(c,l) - eps(0.)
        else
            u=((eps(0.)^η * l^(1-η))^(1-μ) )/ (1-μ)-eps(0.) + 1e5*1e5*min(c,l)
        end
    elseif μ == 1 && η!=1.0
        u = log(c^η * l^(1-η))

    elseif μ != 1 && η!=1.0
        u = ((c^η * l^(1-η))^(1-μ) )/ (1-μ)
    else
        u = u(c)
    end
    return u
end

function u(c::Float64;μ::Real = μ)
    if (c<=0 && μ !=0.0)
        if μ ==1
            u=log(eps(0.))+ 1e10*c - eps(0.)
        else
            u=(eps(0.)^(1-μ) )/(1-μ)-eps(0.) + 1e5*c
        end
    elseif μ == 1
        u = log(c)

    elseif μ != 1
        u = (c^(1-μ) )/ (1-μ)

    end
    return u
end

#derivative of u with respect to c
function uc(c,l;η = η,μ = μ)
    if η!=1 && l <= 0
        return 0.0
    elseif c<=0
        return 1e10*(-c) +1e5
    else
        return (η * c^(η-1) * l^(1-η)) * (c^η * l^(1-η))^(-μ)
    end
end

#derivative of u with respect to l
function ul(c,l;η = η,μ = μ)
    if η!=1
        if c<=0
            return 0.0
        elseif l<=0
            return 1e10 * (-l) +1e5
        else
            return ((1-η) * c^η * l^(-η))*(c^η * l^(1-η))^(-μ)
        end
    else
        return 0.0
    end

end
