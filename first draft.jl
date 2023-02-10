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
        theta[1:3] = (T[1] + T[2] + T[3])/(coef[1] + coef[2] + coef[3])
        T[1:3] = inverse_PT(params, p , PT13)[1:3]
    end
end

function adjust_Nlayers!(params, p, T, theta, n)
    if theta[n-1] > theta[n]
        print("adjust")
        coef = exner(params, p)
        theta[1:n] = sum(T[1:n])/sum(coef[1:n])
        T[1:n] = inverse_PT(params, p, PT)[1:n]
    end
end

function energy_2layers(params, m, T)
    energy = params.Cp*(T[1]*m[1]+T[2]*m[2]) # computing the energy with the formula Cp/g(T1+T2)dp
end

function simplecolumn(params, dt, p, T, m)
    t = 0
    while t <= params.Tf 
        warm!(params, m, T, dt) #we warm only the first layer
        theta = potential_temperature.(Ref(params), p, T)
        before = energy_2layers(params, m, T) #getting the energy of the two layers before the adjustment
        adjust_2layers!(params, p, T, theta) #adjusting the two layers if theta 1 > theta 2
        after = energy_2layers(params, m, T) # energy after
        @info "energy" before after after-before
        t += dt
    end
end

function main()
    params = (Cp = 1000.0, R = 287.0, N = 10, Pr = 100000.0, Q = 1000.0, Tf = 6000.0, g = 9.8)
    p_int = [100000.0, 90000.0, 80000.0, 70000.0, 60000.0, 50000.0, 40000.0, 30000.0, 20000.0, 10000.0]
    T = [293.0, 289.0, 283.0, 277.0, 271.0, 265.0, 259.0, 253.0, 247.0]

    p = zeros(9)
    m = zeros(9)
    for i in (1:9)
        p[i] = (p_int[i] + p_int[i+1])/2 #computing p in layers, only 9 values
        m[i] = (p_int[i] - p_int[i+1])/params.g #computing the mass using dp/g 
    end
    #theta = potential_temperature.(Ref(params), p, T) #computing the potential temperature with p in layers
    simplecolumn(params, 1000, p, T, m) #warming and adjusting the layers
end

main()


