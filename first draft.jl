params = (Cp = 1000.0, R = 287.0, N = 10, Pr = 1000.0, Q = 100.0, Tf = 10.0)


pressure = [1000.0, 900.0, 800.0, 700.0, 600.0, 500.0, 400.0, 300.0, 200.0, 100.0]
#temperature = [293.0, 287.0, 281.0, 275.0, 269.0, 263.0, 257.0, 251.0, 245.0, 239.0]

temperature = [293.0 for i in (1:10)]

theta = temperature.*(pressure/params.Pr).^(-params.R/params.Cp)
print(theta)

function simplecolumn(temperature, Tf)
    dt = 1
    while dt <= Tf
        temperature[1] = params.Q/params.Cp + temperature[1]
        theta[1] = temperature[1]*(pressure[1]/params.Pr)^(-params.R/params.Cp) 
        dt += 1
        
    print(theta)

end

simple(temperature, paramas.Tf)
