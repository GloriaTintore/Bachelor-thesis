using Plots

# constant Cp

struct Ideal
    R:: Float64
end

potential_temperature(params::Ideal, p, T) = T.*(p./params.Pr).^(-params.R/params.Cp) #computing potential temperature, p and T can be either vectors (?) or numbers
exner(params, p) = (p./params.Pr).^(params.R/params.Cp) #computing the coeficients of exner with c = (p/pr)^(R/Cp)
inverse_PT(params, p, theta) = theta*(exner(params,p)) #computing T from theta, need numbers not vectors
pressure(params) = [100000 - 100000/params.N*i for i in(0:params.N)]
temperature(params, P) = params.T0*(P./params.Pr).^(params.y*params.R/params.g)
Cp(gas::Ideal, T) = gas.Cp
# variable Cp

struct VarCp
    R:: Float64
end

potential_temperature(params::VarCp, p, T) = T.*(p./params.Pr).^(-params.R/params.Cp) #computing potential temperature, p and T can be either vectors (?) or numbers
# Cp(gas::Ideal, T) = gas.Cp0*(...)

function h(params, p, theta) 
    Cp0, v, T0, Pr, R = params.Cp0, params.v, params.T0, params.Pr, params.R 
    return Cp0/((v + 1)*T0^v)*(theta^v - v*T0^v*log(Pr/p)^(R/Cp0))^((v+1)/v)
end

function dh(params, p, theta) 
    Cp0, v, T0, Pr, R = params.Cp0, params.v, params.T0, params.Pr, params.R 
    return Cp0/(T0^v)*theta^(v-1)*(theta^v-v*T0^v*log(Pr/p)^(R/Cp0))^(1/v)
end

function warm!(params, m, T, dt)
    T[1] = T[1] + params.Q*dt/(m[1]*Cp(gas,T[1]))  # get new T1 with formula Q = CpdT in layer 1
end

enthalpy_T_Cpvar(params, T) = (T/params.T0)^params.v*T

function energy_2layers(gas, m, T)
#    energy = params.Cp0*(((T[1]/params.T0)^params.v)*T[1]*m[1]+((T[2]/params.T0)^params.v)*T[2]*m[2])/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
    energy = ( m[1]*enthalphy_T(gas, T[1]) + m[2]*enthalphy_T(gas, T[2]) ) / sum(m)
#    ((T[2]/params.T0)^params.v)*T[2]*m[2])/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
end

function energy_Nlayers(params, m, T)
    energy = params.Cp0*sum((T./params.T0).^params.v.*T.*m)/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
end

function simplecolumn(params, dt, p, T, m, f, df)
    t = 0
    init_theta = potential_temperature.(Ref(params), p, T)
    while t <= params.Tf 
        warm!(params, m, T, dt) #we warm only the first layer
        theta = potential_temperature.(Ref(params), p, T)
        before = energy_2layers(params, m, T) #getting the energy of the two layers before the adjustment
        Newton_adjust!(params, f, df, p, T, theta) #adjusting the two layers if theta 1 > theta 2
        after = energy_2layers(params, m, T) # energy after
        @info "energy" before after after-before
        t += dt
    end
    theta = potential_temperature.(Ref(params), p, T)
    plot([init_theta, theta], p, label=["initial PT" "final PT"], yflip= true)
    title!("Potential temperature profile")
    xlabel!("potential temperature")
    ylabel!("pressure (Pa)")
end

function simplecolumnN(params, dt, p, T, m, f, df)
    t = 0
    init_theta = potential_temperature.(Ref(params), p, T)
    thetas = []
    while t <= params.Tf 
        warm!(params, m, T, dt) #we warm only the first layer
        push!(thetas, potential_temperature.(Ref(params), p, T))
        theta = potential_temperature.(Ref(params), p, T)
        before = energy_Nlayers(params, m, T) #getting the energy of the two layers before the adjustment
        for n in (2:params.N-1)
            Newton_adjustN!(params, f, df, p, T, theta, n) #adjusting the two layers if theta 1 > theta 2
        end
        after = energy_Nlayers(params, m, T) # energy after
        @info "energy" before after after-before
        t += dt
    end
    theta = potential_temperature.(Ref(params), p, T)
    plot([init_theta, theta], p, label=["initial PT" "final PT"], yflip= true)
    plot!(thetas, p, label = "warmed PT")
    title!("Potential temperature profile")
    xlabel!("potential temperature")
    ylabel!("pressure (Pa)")
end

function Newton_adjust!(params, f, df, p, T, theta)
    if theta[1] > theta[2]
        thetas = zeros(params.k_max)
        thetas[1] = potential_temperature(params, params.Pr, params.T0)
        k = 1
        while abs(f(params, p, theta, thetas[k])) > params.eps &&  k < params.k_max
            thetas[k+1] = thetas[k] - f(params, p, theta, thetas[k])/df(params, p, theta, thetas[k])
            k += 1
        end
        theta[1:2] .= last(thetas[1:k])
        T[1:2] = inverse_PT(params, p , theta[1])[1:2]
    end
end

function Newton_adjustN!(params, f, df, p, T, theta, n)
    if theta[1] > theta[n]
        thetas = zeros(params.k_max)
        thetas[1] = theta[1]
        k = 1
        while abs(f(params, p, theta, thetas[k], n)) > params.eps &&  k < params.k_max
            thetas[k+1] = thetas[k] - f(params, p, theta, thetas[k], n)/df(params, p, theta, thetas[k], n)
            k += 1
        end
        theta[1:n] .= last(thetas[1:k])
        T[1:n] = inverse_PT(params, p, theta[1])[1:n]
    end
end

function main()
#    params = (Cp = 1000.0, R = 287.0, N = 10, Pr = 100000.0, Q = 1000.0, Tf = 30000.0, g = 9.8, y = 0.005, T0 = 293.0, v = 0.35, Cp0 = 1000.0, k_max = 20, eps = 0.0001)
    gas = Ideal(1000.0, 287.0, Pr = 100000.0)
    params = N = 10, Q = 1000.0, Tf = 30000.0, g = 9.8, y = 0.005, T0 = 293.0, v = 0.35, Cp0 = 1000.0, k_max = 20, eps = 0.0001)
    P_int = pressure(params)
    T_int = temperature(params, P_int)

    p = zeros(params.N-1)
    m = zeros(params.N-1)
    T = zeros(params.N-1)

    for i in (1:params.N-1)
        p[i] = (P_int[i] + P_int[i+1])/2 #computing p in layers, only 9 values
        m[i] = (P_int[i] - P_int[i+1])/params.g #computing the mass using dp/g 
        T[i] = (T_int[i] + T_int[i+1])/2
    end

    dest_deriv(params, p, T)

    f(params, p, theta, x) = h(params, p[1], x) + h(params, p[2], x) - h(params, p[1], theta[1]) - h(params, p[2], theta[2])
    df(params, p, theta, x) = dh(params, p[1], x) + dh(params, p[2], x)

    fN(params, p, theta, x, n) = sum(h(params, p[i], x) for i in (1:n)) - sum(h(params, p[i], theta[i]) for i in (1:n))
    dfN(params, p, theta, x, n) = sum(dh(params, p[i], x) for i in (1:n))
    simplecolumnN(params, 10000, p, T, m, fN, dfN) #warming and adjusting the layers
end

using ForwardDiff: derivative

function dest_deriv(params, p, T)
    theta = potential_temperature.(Ref(params), p, T)
    f(x) = h(params, p[1], x) + h(params, p[2], x) - h(params, p[1], theta[1]) - h(params, p[2], theta[2])
    df(x) = dh(params, p[1], x) + dh(params, p[2], x)
    x=theta[1]
    @info "f" f(x) df(x) derivative(f,x)
end

main()
