using DifferentialEquations
using Plots
using LinearAlgebra
using Roots

begin 
    function f(du,u,p,t)
        Ud,w,Pdc,Cd = p

        du[1] = u[1] - Ud*sin(w*t)#*exp(-t)
        du[2] = 0.0 # u_off
        du[3] = 0.0 # t_sum
        du[4] = 0.0 # t_on
        du[5] = 0.0 # t_off
        du[6] = 0.0 # half_cycle
        du[7] = 0.0 # a
        du[8] = 0.0 # b
        du[9] = u[9] - Ud#*exp(-t)
        du[10] = 0.0 # hilf_t_on
        du[11] = 0.0 # hilf_t_off
    end

    u0 = [0.0, 230*sqrt(2),0.0,0.0034,0.00511210289264389,0.0,0.0,0.0,230*sqrt(2),0.0,0.0]
    p =  [230*sqrt(2),100*pi,1000.0,700e-6]

    M = Diagonal([0, 1, 1, 1, 1, 1, 1, 1,0,1,1])

    function fault_start(integrator)
        integrator.p[1] = 0.5*integrator.p[1]
    end

    function fault_end(integrator)
        integrator.p[1] = 2*integrator.p[1]
    end

    s0(u, t, integrator) = true
    function h0(integrator) 
        integrator.u[3] = integrator.uprev[3] + integrator.t -integrator.tprev   
    end

    s1(u, t, integrator) = u[3] >= 0.01
    function h1(integrator) 
        Ud,w,Pdc,Cd = integrator.p
        integrator.u[6] = 1.0
        integrator.u[3] = mod(integrator.t,0.01)
        if integrator.u[4] > 0
            u_off = integrator.u[9]*sin(w*integrator.u[5])
            integrator.u[2] = sqrt(u_off^2 - 2*Pdc*(0.01-integrator.u[5])/Cd)
        else
            u_off = integrator.u[9]*sin(w*integrator.u[5])
            integrator.u[2] = maximum(real(sqrt(Complex(integrator.u[2]^2 - 2*Pdc*(0.01)/Cd))))
        end
    end    

    s2(u, t, integrator) = u[10] == 1.0 && u[11] == 1.0
    function h2(integrator)
        integrator.u[6] = 0.0
        integrator.u[10] = 0.0
        integrator.u[11] = 0.0
    end 

    s3(u, t, integrator) = u[6] == 1.0 || mod(t,0.01) < u[5]
    function h3(integrator) 
        Ud,w,Pdc,Cd = integrator.p
        integrator.u[5] = (pi + asin(2*Pdc/(w*Cd*integrator.u[9]^2)))/(2*w)
        integrator.u[11] = 1.0
    end 

    s4(u, t, integrator) = u[6] == 1.0 || mod(t,0.01) < u[4]
    function h4(integrator)
        Ud,w,Pdc,Cd = integrator.p
        uofft2 = integrator.u[2]
        u0 = integrator.u[9]
        x1 = 0.004517042542168
        x2 = -0.084973441092720
        x3 =  1.367157627046003
        x4 = -0.580239811988436

        A = u0^2*w^3*x4
        B = u0^2*w^2*x3
        C = u0^2*w*x2 + 2*Pdc/Cd
        D = u0^2*x1 - uofft2^2  

        integrator.u[4] = CardanosFormular(A,B,C,D)
        integrator.u[10] = 1.0
    end 

    s5(u, t, integrator) = u[4] <= 0.0
    function h5(integrator)
        integrator.u[7] = 0.0
        integrator.u[8] = 0.0
    end

    s6(u, t, integrator) = u[4] > 0.0
    function h6(integrator)
        Ud,w,Pdc,Cd = p
        T = 0.02
        Ud = integrator.u[9]
        te = integrator.u[4]
        ta = integrator.u[5]

        Teil1 = Ud*(1/2)*w*Cd*(ta-te+1/(2*w)*(sin(2*w*ta)-sin(2*w*te)))
        Teil2 = Pdc/(Ud*w)*(log(abs(sin(w*ta)))-log(abs(sin(w*te))))
        an = 4/T*(Teil1+Teil2)
    
        Teil1b = Ud*(1/4)*Cd*(cos(2*w*te)-cos(2*w*ta))
        Teil2b = Pdc/Ud*(ta-te)
        bn = 4/T*(Teil1b+Teil2b)

        integrator.u[7] = an*Ud/2/Pdc
        integrator.u[8] = bn*Ud/2/Pdc
    end

    cb_fault = PresetTimeCallback([0.04],fault_start)
    cb_clear = PresetTimeCallback([0.1],fault_end)

    cb0 = DiscreteCallback(s0,h0)
    cb1 = DiscreteCallback(s1,h1)
    cb2 = DiscreteCallback(s2,h2)
    cb3 = DiscreteCallback(s3,h3)
    cb4 = DiscreteCallback(s4,h4)
    cb5 = DiscreteCallback(s5,h5)
    cb6 = DiscreteCallback(s6,h6)

    ode_fun = ODEFunction(f,mass_matrix = M)
    prob = ODEProblem(ode_fun, u0, (0.0, 0.14), p)
    sol = solve(prob, Rodas4P(), callback = CallbackSet(cb_fault,cb_clear,cb0,cb1,cb2,cb3,cb4,cb5,cb6),dtmax=1e-3);
    nothing
end

plot(sol,idxs=[9])

plot(sol.t.-sol[9,:])



function CardanosFormular(A::Float64,B::Float64,C::Float64,D::Float64)
    p = (9*A*C-3*B^2)./(9*A^2);
    q = (2*B^3 -9*A*B*C + 27*D*A^2)./(27*A^3)
    delta = (q/2)^2 +  (p/3)^3

    roots = Float64
    if delta >0
        u = Complex(-q/2+sqrt(delta))^(1/3)
        v = Complex(-q/2-sqrt(delta))^(1/3)
        x1 = u+v- B/(3*A) 
        x2 = -(u+v)/2 - B/(3*A) + 1im*(u-v)/2*sqrt(3)
        x3 = -(u+v)/2 - B/(3*A) - 1im*(u-v)/2*sqrt(3)
        roots = minimum(real([x1 x2 x3]))
    elseif delta == 0 && p == 0
        x2 = - B/(3*A)
        roots = real(x2)
    elseif delta == 0 && p != 0
        x1 = 3*q/p- B/(3*A)
        x23 = -3*q/(2*p)- B/(3*A)
        roots = real(x23)
    elseif delta < 0
        x1 = -sqrt(-4/3*p).*cos(1/3.0*acos(-q/2.0*sqrt(-27.0/(p.^3)))+pi/3)-B/(3*A)
        x2 =  sqrt(-4/3*p).*cos(1/3.0*acos(-q/2.0*sqrt(-27.0/(p.^3))))     -B./(3*A)
        x3 = -sqrt(-4/3*p).*cos(1/3.0*acos(-q/2.0*sqrt(-27.0/(p.^3)))-pi/3)-B./(3*A)
        roots =  real(x1)
    end
    roots
end

a = -3.598215230287720e+07
b = 2.698660986575224e+05
c = 0.619874292343603
d = -1.724611446133485
CardanosFormular(a,b,c,d)

function fton(ton)
    ton^3*(u0^2*w^3*x4) + ton^2*(u0^2*w^2*x3) + ton*(u0^2*w*x2 + 2*Pdc/Cd) + u0^2*x1 - uofft2^2  
end