function Example1()
    function f(dx,x,p,t)
        λ = p
        dx[1] = x[1] + x[3]*x[2] 
        dx[2] = x[4]*x[1] + x[2]
        dx[3] = 0.0
        dx[4] = 0.0
        dx[5] = 0.0
    end

    evr = Array{Float64}(undef,0,4)
    s1(u, t, integrator) = integrator.p[1]*u[1] - u[2]  
    s2(u, t, integrator) = -0.36*u[1] + u[2] 

    function h1(integrator) 
        integrator.u[3] = integrator.u[4] 
        integrator.u[4] = integrator.u[3]
        integrator.u[5] = -integrator.u[5]
        evr = vcat(evr, [integrator.t 1 1 1])
    end
    function h2(integrator) 
        integrator.u[3] = integrator.u[4] 
        integrator.u[4] = integrator.u[3]
        integrator.u[5] = -integrator.u[5]
        evr = vcat(evr, [integrator.t 1 2 1])
    end
    cb1 = ContinuousCallback(s1,nothing,affect_neg! =h1)
    cb2 = ContinuousCallback(s2,nothing,affect_neg! =h2)

    u0 = [0.0,1,-100,10,1]
    p = [2.75]
    ode_fun = ODEFunction(f)
    prob = ODEProblem(ode_fun, u0, (0.0, 0.2), p)
    sol = solve(prob, Rodas4(), callback = CallbackSet(cb1,cb2),dtmax=1e-4)
    return sol, evr
end

sol,evr = Example1();
plot(sol.t,sol[1,:])
plot!(sol.t,sol[2,:])
plot(sol.t,sol[5,:])

plot(sol,idxs=(1,2))