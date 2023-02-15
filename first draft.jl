using Plots
potential_temperature(params, p, T) = T.*(p./params.Pr).^(-params.R/params.Cp) #computing potential temperature, p and T can be either vectors (?) or numbers

exner(params, p) = (p./params.Pr).^(params.R/params.Cp) #computing the coeficients of exner with c = (p/pr)^(R/Cp)

inverse_PT(params, p, theta) = theta*(exner(params,p)) #computing T from theta, need numbers not vectors

function warm!(params, m, T, dt)
    T[1] = T[1] + params.Q*dt/(m[1]*params.Cp)  # get new T1 with formula Q = CpdT in layer 1
end

function adjust_2layers!(params, p, T, theta)
    if theta[1] > theta[2]
        print("adjust") 
        coef = exner(params, p)
        theta[1:2] .= (T[1] + T[2])/(coef[1] + coef[2])
        T[1:2] = inverse_PT(params, p , theta[1])[1:2]
        print(theta)
    end
end

function adjust_3layers!(params, p, T, theta)
    if theta[1] > theta[3]
        print("adjust")
        coef = exner(params, p)
        theta[1:3] .= (T[1] + T[2] + T[3])/(coef[1] + coef[2] + coef[3])
        T[1:3] = inverse_PT(params, p , PT13)[1:3]
    end
end

function adjust_Nlayers!(params, p, T, theta, n)
    if theta[1] > theta[n]
        print("adjust")
        coef = exner(params, p)
        theta[1:n] .= sum(T[1:n])/sum(coef[1:n])
        T[1:n] = inverse_PT(params, p, theta[1])[1:n]
        #print(theta)
    end
end

function energy_2layers(params, m, T)
    energy = params.Cp*(T[1]*m[1]+T[2]*m[2])/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
end

function energy_Nlayers(params, m, T)
    energy = params.Cp*(sum(T.*m))/sum(m) # computing the energy with the formula Cp/g(T1+T2)dp
end

function simplecolumn(params, dt, p, T, m)
    t = 0
    init_theta = potential_temperature.(Ref(params), p, T)
    while t <= params.Tf 
        warm!(params, m, T, dt) #we warm only the first layer
        theta = potential_temperature.(Ref(params), p, T)
        before = energy_2layers(params, m, T) #getting the energy of the two layers before the adjustment
        adjust_2layers!(params, p, T, theta) #adjusting the two layers if theta 1 > theta 2
        after = energy_2layers(params, m, T) # energy after
        @info "energy" before after after-before
        t += dt
    end
    theta = potential_temperature.(Ref(params), p, T)
    plot([theta, init_theta], p)
end

function simplecolumnN(params, dt, p, T, m)
    t = 0
    init_theta = potential_temperature.(Ref(params), p, T)
    thetas = []
    while t <= params.Tf 
        warm!(params, m, T, dt) #we warm only the first layer
        push!(thetas, potential_temperature.(Ref(params), p, T))
        theta = potential_temperature.(Ref(params), p, T)
        before = energy_Nlayers(params, m, T) #getting the energy of the two layers before the adjustment
        for i in (2:params.N-1)
            adjust_Nlayers!(params, p, T, theta, i) #adjusting the two layers if theta 1 > theta 2
        end
        after = energy_Nlayers(params, m, T) # energy after
        @info "energy" before after after-before
        t += dt
    end
    theta = potential_temperature.(Ref(params), p, T)
    plot([init_theta, theta], p, label=["initial PT" "final PT"], yflip= true)
    plot!(thetas, p, label ="warmed PT")
    title!("Potential temperature profile")
    xlabel!("potential temperature")
    ylabel!("pressure (Pa)")
end

pressure(params) = [100000 - 10000*i for i in(0:params.N)]
temperature(params, P) = params.T0*(P./params.P0).^(params.y*params.R/params.g)

function main()
    params = (Cp = 1000.0, R = 287.0, N = 10, Pr = 100000.0, Q = 1000.0, Tf = 30000.0, g = 9.8, y = 0.007, T0 = 293.0, P0 = 100000)
    P_int = pressure(params)
    T_int = temperature(params, P_int)

    #T_int = [293 - params.y*i for i in (0:params.N)]
    #P_int = [100000 - 10000*i for i in (0:params.N)]

    p = zeros(9)
    m = zeros(9)
    T = zeros(9)

    for i in (1:9)
        p[i] = (P_int[i] + P_int[i+1])/2 #computing p in layers, only 9 values
        m[i] = (P_int[i] - P_int[i+1])/params.g #computing the mass using dp/g 
        T[i] = (T_int[i] + T_int[i+1])/2
    end
    #theta = potential_temperature.(Ref(params), p, T) #computing the potential temperature with p in layers
    simplecolumnN(params, 10000, p, T, m) #warming and adjusting the layers
    #theta = potential_temperature(params, p, T)
    #print(T, p, theta)
end

main()


