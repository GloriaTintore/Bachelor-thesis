function Newton(f, df, theta0, k_max, eps, p, theta)
    x = zeros(k_max)
    thetas = zeros(k_max)
    thetas[1] = theta0
    k = 1
    while abs(f(p, theta, thetas[k])) > eps &&  k < k_max
        thetas[k+1] = thetas[k] - f(p, theta, thetas[k])/df(p, theta, thetas[k])
        k += 1
    end
    xks = thetas[:k+1]
end

f(p, theta, x) = h(p[1], x) + h(p[2], x) - h(p[1], theta[1]) - h(p[2], theta[2])

df(p, theta, x) = dh(p[1], x) + dh(p[2], x)


function main()
    p = [100000 - 10000*i for i in(0:9)]
    theta = [300 - 5*i for i in (0:9)]

    xks = Newton(f, df, 340, 15, 1^(-5), p, theta)
    print(xks)
end

main()