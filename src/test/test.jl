using ModelingToolkit, Plots, DifferentialEquations

@variables t
@connector function Pin(;name)
    sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name=name)
end

function Ground(;name)
    @named g = Pin()
    eqs = [g.v ~ 0]
    compose(ODESystem(eqs, t, [], []; name=name), g)
end

function OnePort(;name)
    @named p = Pin()
    @named n = Pin()
    sts = @variables v(t)=1.0 i(t)=1.0
    eqs = [
           v ~ p.v - n.v
           0 ~ p.i + n.i
           i ~ p.i
          ]
    compose(ODESystem(eqs, t, sts, []; name=name), p, n)
end

function Resistor(;name, R = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
           v ~ i * R
          ]
    extend(ODESystem(eqs, t, [], ps; name=name  ), oneport)
end

function Capacitor(;name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
           D(v) ~ i / C
          ]
    extend(ODESystem(eqs, t, [], ps; name=name  ), oneport)
end

function ACVoltage(;name, V = 1.0, freq = 1.0, phase = 0.0)
    @named oneport = OnePort()
    @unpack v,i = oneport
    ps = @parameters V=V freq=freq phase=phase
    eqs = [
           V*sin(2*pi*freq*t+phase) ~ v
          ]
    extend(ODESystem(eqs, t, [], ps; name=name ), oneport)
end

@named R = Resistor(R= 1)
@named C = Capacitor(C= 1e-3)
@named source = ACVoltage(V=1*sqrt(2),freq = 50.0)
@named ground = Ground()

 rc_eqs = [
            connect(source.p,R.p)
            connect(R.n,C.p)
            connect(C.n,source.n)
            connect(source.n,ground.g)
          ]


@named rc_model = ODESystem(rc_eqs, t,[R.v],[];
                            systems = [R,C,source,ground],
                            continuous_events = [C.v ~  1.0] => [source.freq ~ 1.0])

sys = structural_simplify(rc_model)

u0 = [
    C.v => 0.0
  ]
prob = ODEProblem(sys, u0, (0,2.0))

sol = solve(prob, Rodas4())
plot(sol,vars = [R.i], layout=(1,1))

ModelingToolkit.continuous_events(sys)


function dq(a,b,c,θ)
    θ += pi/2
    faktor = sqrt(2/3)
    d = faktor * (sin(θ)*a + sin(θ-2*pi/3)*b + sin(θ-4*pi/3)*c)
    q = faktor * (cos(θ)*a + cos(θ-2*pi/3)*b + cos(θ-4*pi/3)*c)
    z = faktor * (sqrt(0.5)*  (a+b+c))

    return d,q,z
end

theta = 0:1e-3:6*pi
ua = sqrt(1)*cos.(theta)
ub = sqrt(1)*cos.(theta.-2*pi/3)
uc = sqrt(1)*cos.(theta.+2*pi/3)

ia = 0.5*sqrt(2)*sin.(theta)
ib = 0.5*sqrt(2)*sin.(theta.-2*pi/3)
ic = 0.5*sqrt(2)*sin.(theta.+2*pi/3)


dqz = dq.(ua,ub,uc,theta)
idqz = dq.(ia,ib,ic,theta)

pdqz= dqz.*idqz

dqz = [collect(i) for i in dqz]
for i in enumerate(dqz)
    dqz
end
function tup2mat(tup)
    u = zeros(length(tup))
    v = zeros(length(tup))
    w = zeros(length(tup))

    for (ind,val) in enumerate(tup)
        u[ind] = val[1]
        v[ind] = val[2]
        w[ind] = val[3]
    end
    return [u v w]
end
function invdq(d,q,z,θ)
    θ += pi/2
    faktor = sqrt(2/3)
    a = faktor * (sin(θ)*d + cos(θ)*q + sqrt(0.5)*z)
    b = faktor * (sin(θ-2*pi/3)*d + cos(θ-2*pi/3)*q + sqrt(0.5)*z)
    c = faktor * (sin(θ-4*pi/3)*d + cos(θ-4*pi/3)*q + sqrt(0.5)*z)

    return a,b,c
end

abc = invdq.(d,q,z,theta)

tmp = tup2mat(abc)

plot(tmp)
plot!(ua)

t = 0:1e-4:0.04
ua = sqrt(1)*cos.(2*pi*50 .*t)
ub = sqrt(1)*cos.(2*pi*50 .*t .-2*pi/3)
uc = sqrt(1)*cos.(2*pi*50 .*t .+2*pi/3)

plot(t,ua)
plot!(t,ub)
plot!(t,uc)


function CalcParallel(R)
    tmp = R[1]
    for ind in 2:length(R)
        tmp = parallel(R[ind],tmp)
    end
    return tmp
end


function parallel(R1::Union{ComplexF64,Float64},R2::Union{ComplexF64,Float64})
    return (R1*R2)/(R1+R2)
end

r1 = [9.6+1im*64,9.6+1im*64,16+1im*64,16+1im*96]
a = [1,1]
CalcParallel(r1)

cd(joinpath(DEPOT_PATH[1], "registries", "General")) do
           deps = Pkg.dependencies()
           registry = Pkg.TOML.parse(read("Registry.toml", String))
           general_pkgs = registry["packages"]

           constrained = Dict{String, Tuple{VersionNumber,VersionNumber}}()
           for (uuid, dep) in deps
               suuid = string(uuid)
               dep.is_direct_dep || continue
               dep.version === nothing && continue
               haskey(general_pkgs, suuid) || continue
               pkg_meta = general_pkgs[suuid]
               pkg_path = joinpath(pkg_meta["path"], "Versions.toml")
               versions = Pkg.TOML.parse(read(pkg_path, String))
               newest = maximum(VersionNumber.(keys(versions)))
               if newest > dep.version
                   constrained[dep.name] = (dep.version, newest)
               end
           end

           return constrained
       end


function find_holdback(held_back_package, recursive = 0, installed = collect(keys(Pkg.installed())), inset=0)
   @assert held_back_package ∈ installed "The package you are querying is not installed"
   for package in installed
       REQUIRE = Pkg.dir(package, "REQUIRE")
       isfile(REQUIRE) || continue
       for line in eachline(REQUIRE)
           splitline = split(line, ' ')
           splitline[1] == held_back_package || continue
           for i = 1:inset print("    ") end
           if length(splitline) > 2
               print_with_color(:red, @display("Package %-35s requires %s\n", package, line))
           else
               @display("Package %-35s requires %s\n", package, line)
           end
           recursive > 0 && find_holdback(package, recursive-1, installed, inset+1)
       end
   end
end


using DifferentialEquations
using LinearAlgebra

function f(du,u,p,t)
    du[1] = u[2]
    du[2] = u[2] - (u[1] + 1/3*u[2]^3  - u[2])
end

mm = Diagonal([true, false])

fun  = ODEFunction(f,mass_matrix=mm)

u0 = [2,-2/3]
prob = ODEProblem(fun, u0,(1.5,5),nothing)
sol = solve(prob,Rodas4(),dtmax=1e-5, initializealg = BrownFullBasicInit())

plot([scatter(x=sol.t,y=sol[1,:]),scatter(x=sol.t,y=sol[2,:])])