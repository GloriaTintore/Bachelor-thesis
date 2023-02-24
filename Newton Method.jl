using Plots
#common functions
pressure(params) = [100000 - 100000/params.N*i for i in(0:params.N)]
temperature(gas, params, P) = gas.T0*(P./gas.Pr).^(params.y*gas.R/params.g)
exner(gas, p) = (p./gas.Pr).^(gas.R/gas.Cp0) #computing the coeficients of exner with c = (p/pr)^(R/Cp)

#for constant Cp
struct Ideal
    R::Float64
    Cp0::Float64
    Pr::Float64
    T0::Float64
end

potential_temperature(gas::Ideal, p, T) = T.*(p./gas.Pr).^(-gas.R/gas.Cp0) #computing potential temperature, p and T can be either vectors (?) or numbers
inverse_PT(gas::Ideal, p, theta) = theta*(exner(gas,p)) #computing T from theta, need numbers not vectors
Cp(gas::Ideal, T) = gas.Cp0

function adjust_2layers!(gas::Ideal, params, f, df, p, T, theta)
    if theta[1] > theta[2]
        coef = exner(gas, p)
        theta[1:2] .= (T[1] + T[2])/(coef[1] + coef[2])
        T[1:2] = inverse_PT(gas, p , theta[1])[1:2]
        print(theta)
    end
end

function adjust_Nlayers!(gas::Ideal, params, f, df, p, T, theta, n)
    if theta[1] > theta[n]
        coef = exner(gas, p)
        theta[1:n] .= sum(T[1:n])/sum(coef[1:n])
        T[1:n] = inverse_PT(gas, p, theta[1])[1:n]
    end
end

#for varying Cp
struct NoIdeal
    R :: Float64
    Cp0 :: Float64
    Pr :: Float64
    T0 :: Float64
    nu :: Float64
end

potential_temperature(gas::NoIdeal, p, T) = (T^gas.nu + gas.nu*gas.T0^gas.nu*log.((gas.Pr/p)^(gas.R/gas.Cp0)))^(1/gas.nu)
inverse_PT(gas::NoIdeal, p, theta) = (theta^gas.nu .- gas.nu*gas.T0^gas.nu*log.((gas.Pr./p).^(gas.R/gas.Cp0))).^(1/gas.nu)
Cp(gas::NoIdeal, T) = gas.Cp0*(T./gas.T0).^gas.nu

function adjust_2layers!(gas::NoIdeal, params, f, df, p, T, theta)
    if theta[1] > theta[2]
        thetas = zeros(params.k_max)
        thetas[1] = theta[1]
        k = 1
        while abs(f(gas, p, theta, thetas[k])) > params.eps &&  k < params.k_max
            thetas[k+1] = thetas[k] - f(gas, p, theta, thetas[k])/df(gas, p, theta, thetas[k])
            k += 1
        end
        theta[1:2] .= last(thetas[1:k])
        T[1:2] = inverse_PT(gas, p, theta[1])[1:2]
    end
end

function adjust_Nlayers!(gas::NoIdeal, params, f, df, p, T, theta, n)
    if theta[1] > theta[n]
        thetas = zeros(params.k_max)
        thetas[1] = theta[1]
        k = 1
        while abs(f(gas, p, theta, thetas[k], n)) > params.eps &&  k < params.k_max
            thetas[k+1] = thetas[k] - f(gas, p, theta, thetas[k], n)/df(gas, p, theta, thetas[k], n)
            k += 1
        end
        theta[1:n] .= last(thetas[1:k])
        T[1:n] = inverse_PT(gas, p, theta[1])[1:n]
    end
end


#function warm for both ideal and no ideal using Cp()
function warm!(gas, params, m, T)
    T[1] = T[1] + params.Q*params.dt/(m[1]*Cp(gas,T[1]))
end

enthalpy_T(gas::Ideal, T) = Cp(gas, T)*T #Cp0 cst
enthalpy_T(gas::NoIdeal, T) = Cp(gas, T)/(gas.nu+1).*T

function h(gas::NoIdeal, p, theta) 
    Cp0, nu, T0, Pr, R = gas.Cp0, gas.nu, gas.T0, gas.Pr, gas.R 
    return Cp0/((nu + 1)*T0^nu)*(theta^nu - nu*T0^nu*log((Pr/p)^(R/Cp0)))^((nu+1)/nu)
end

function dh(gas::NoIdeal, p, theta) 
    Cp0, nu, T0, Pr, R = gas.Cp0, gas.nu, gas.T0, gas.Pr, gas.R 
    return Cp0/(T0^nu)*theta^(nu-1)*(theta^nu-nu*T0^nu*log((Pr/p)^(R/Cp0)))^(1/nu)
end

function energy_2layers(gas, m, T)
#    energy = params.Cp0*(((T[1]/params.T0)^params.v)*T[1]*m[1]+((T[2]/params.T0)^params.v)*T[2]*m[2])/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
    return (m[1]*enthalpy_T(gas, T[1]) + m[2]*enthalpy_T(gas, T[2])) / sum(m)
#    ((T[2]/params.T0)^params.v)*T[2]*m[2])/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
end

function energy_Nlayers(gas, m, T)
    return sum(m.*enthalpy_T(gas, T))/sum(m)
    #energy = params.Cp0*sum((T./params.T0).^params.v.*T.*m)/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
end

function simplecolumn(gas, params, p, T, m, f, df)
    t = 0
    init_theta = potential_temperature.(Ref(gas), p, T)
    while t <= params.Tf 
        warm!(gas, params, m, T) #we warm only the first layer
        theta = potential_temperature.(Ref(gas), p, T)
        #print(T, theta)
        before = energy_2layers(gas, m, T) #getting the energy of the two layers before the adjustment
        adjust_2layers!(gas, params, f, df, p, T, theta) #adjusting the two layers if theta 1 > theta 2
        after = energy_2layers(gas, m, T) # energy after
        @info "energy" before after round(after-before)
        t += params.dt
    end
    theta = potential_temperature.(Ref(gas), p, T)
    plot([init_theta, theta], p, label=["initial PT" "final PT"], yflip= true)
    title!("Potential temperature profile")
    xlabel!("potential temperature")
    ylabel!("pressure (Pa)")
end

function simplecolumnN(gas, params, p, T, m, f, df)
    t = 0
    init_theta = potential_temperature.(Ref(gas), p, T)
    thetas = []
    while t <= params.Tf 
        warm!(gas, params, m, T) #we warm only the first layer
        push!(thetas, potential_temperature.(Ref(gas), p, T))
        theta = potential_temperature.(Ref(gas), p, T)
        before = energy_Nlayers(gas, m, T) #getting the energy of the two layers before the adjustment
        for n in (2:params.N-1)
            adjust_Nlayers!(gas, params, f, df, p, T, theta, n) #adjusting the two layers if theta 1 > theta 2
        end
        after = energy_Nlayers(gas, m, T) # energy after
        @info "energy" before after round(after-before)
        t += params.dt
    end
    theta = potential_temperature.(Ref(gas), p, T)
    plot([init_theta, theta], p, label=["initial PT" "final PT"], yflip= true)
    plot!(thetas, p, label = "warmed PT")
    title!("Potential temperature profile")
    xlabel!("potential temperature")
    ylabel!("pressure (Pa)")
end



function main()
#   params = (Cp = 1000.0, R = 287.0, N = 10, Pr = 100000.0, Q = 1000.0, Tf = 30000.0, g = 9.8, y = 0.005, T0 = 293.0, v = 0.35, Cp0 = 1000.0, k_max = 20, eps = 0.0001)
    #gas = Ideal(287.0, 1000.0, 100000.0, 293.0)
    gas = NoIdeal(287.0, 1000.0, 100000.0, 293.0, 0.35)
    params = (N = 10, Q = 1000.0, Tf = 30000.0, g = 9.8, k_max = 20, eps = 0.0001, dt = 10000, y = 0.007)
    P_int = pressure(params)
    T_int = temperature(gas, params, P_int)

    p = zeros(params.N-1)
    m = zeros(params.N-1)
    T = zeros(params.N-1)

    for i in (1:params.N-1)
        p[i] = (P_int[i] + P_int[i+1])/2 #computing p in layers, only 9 values
        m[i] = (P_int[i] - P_int[i+1])/params.g #computing the mass using dp/g 
        T[i] = (T_int[i] + T_int[i+1])/2
    end

    f(gas, p, theta, x) = h(gas, p[1], x) + h(gas, p[2], x) - h(gas, p[1], theta[1]) - h(gas, p[2], theta[2])
    df(gas, p, theta, x) = dh(gas, p[1], x) + dh(gas, p[2], x)

    fN(gas, p, theta, x, n) = sum(h(gas, p[i], x) for i in (1:n)) - sum(h(gas, p[i], theta[i]) for i in (1:n))
    dfN(gas, p, theta, x, n) = sum(dh(gas, p[i], x) for i in (1:n))
    simplecolumnN(gas, params, p, T, m, fN, dfN) #warming and adjusting the layers
end

using ForwardDiff: derivative

function dest_deriv(gas, p, T)
    theta = potential_temperature.(Ref(gas), p, T)
    f(x) = h(gas, p[1], x) + h(gas, p[2], x) - h(gas, p[1], theta[1]) - h(gas, p[2], theta[2])
    df(x) = dh(gas, p[1], x) + dh(gas, p[2], x)
    x=theta[1]
    @info "f" f(x) df(x) derivative(f,x)
end

main()
