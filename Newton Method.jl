using Plots
potential_temperature(params, p, T) = T.*(p./params.Pr).^(-params.R/params.Cp) #computing potential temperature, p and T can be either vectors (?) or numbers

exner(params, p) = (p./params.Pr).^(params.R/params.Cp) #computing the coeficients of exner with c = (p/pr)^(R/Cp)

inverse_PT(params, p, theta) = theta*(exner(params,p)) #computing T from theta, need numbers not vectors

pressure(params) = [100000 - 8000*i for i in(0:params.N)]
temperature(params, P) = params.T0*(P./params.Pr).^(params.y*params.R/params.g)

h(params, p, theta) = params.Cp0/((params.v + 1)*params.T0^params.v)*(theta^params.v - params.v*params.T0^params.v*log(params.Pr/p)^(params.R/params.Cp0))^((params.v+1)/params.v)
dh(params, p, theta) = params.Cp0/(params.T0^params.v)*theta^(params.v-1)*(theta^params.v-params.v*params.T0^params.v*log(params.Pr/p)^params.R/params.Cp0)^(1/params.v)

function warm!(params, m, T, dt)
    T[1] = T[1] + params.Q*dt/(m[1]*params.Cp)  # get new T1 with formula Q = CpdT in layer 1
end

function energy_2layers(params, m, T)
    energy = params.Cp0*(((T[1]/params.T0)^params.v)*T[1]*m[1]+((T[2]/params.T0)^params.v)*T[2]*m[2])/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
end

function energy_Nlayers(params, m, T)
    energy = params.Cp*(sum(T.*m))/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
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
    params = (Cp = 1000.0, R = 287.0, N = 10, Pr = 100000.0, Q = 1000.0, Tf = 30000.0, g = 9.8, y = 0.005, T0 = 293.0, v = 0.35, Cp0 = 1000.0, k_max = 40, eps = 0.001)
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
    f(params, p, theta, x) = h(params, p[1], x) + h(params, p[2], x) - h(params, p[1], theta[1]) - h(params, p[2], theta[2])
    df(params, p, theta, x) = dh(params, p[1], x) + dh(params, p[2], x)

    fN(params, p, theta, x, n) = sum(h(params, p[i], x) for i in (1:n)) - sum(h(params, p[i], theta[i]) for i in (1:n))
    dfN(params, p, theta, x, n) = sum(dh(params, p[i], x) for i in (1:n))
    simplecolumnN(params, 10000, p, T, m, fN, dfN) #warming and adjusting the layers
end

main()