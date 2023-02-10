using Plots

include("PerfectGas.jl")

struct Params
    cp :: Float64
    r :: Float64
    N :: Int
end

function ff(x::Float64)
    return x*x
end

ff(x::Int) = x+x

function ff(x)
    return x+1
end

function f(x)
    for i in (1,3,7)
        x[i] = i
    end
    return nothing
end


if false
    params = (Cp=1000.0, R=287.0, N=10)

    @info "X" params

    print("toto")

    params = (params... , N=20)

    params = Params(1.5f0, 2.5, 3)
end

#theend