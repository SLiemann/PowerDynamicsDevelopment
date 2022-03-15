using PowerDynamics
#using OrderedCollections: OrderedDict
using Plots
#import PowerDynamics: PiModel
using DifferentialEquations
#using CSV #read PF DataFrames
#using DataFrames #for CSV
#using Distributed
using ForwardDiff

#Ubase = 380e3
#Sbase = 100e6
#Zbase = (Ubase^2)/Sbase

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/LTVS_Test_System_TapParam.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/sensitivity_analyses/Local_Sensitivity.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/utility/utility_functions.jl")
end


begin
    pg = LTVS_Test_System()
    Qmax   = [Inf, Inf, Inf,Inf, Inf,Inf*sqrt(1-0.95^2)]
    Qmin   = -Qmax
    U,δ,ic0 = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2,max_tol = 1e-6)
    Uc = U.*exp.(1im*δ/180*pi)
    Ykk = NodalAdmittanceMatrice(pg)
    Ic = abs.(Ykk*Uc./5.5)
    S  = round.(Uc.*(conj.(Ykk)*conj.(Uc)),digits=3)
    pg, ic0 = InitializeInternalDynamics(pg,ic0)
end

prob = ODEProblem{true}(rhs(pg),ic0,(0.0,1.0))
integrator = init(prob,Rodas4(),reltol=1e-6,abstol=1e-6)

function sum_of_solution(x)
    #_prob = remake(prob,u0=x)
    step!(integrator)
    return integrator.sol[end]
    #solve(_prob,Rodas4(),reltol=1e-6,abstol=1e-6)[end] #,reltol=1e-6,abstol=1e-6,saveat=0.1
end


using DiffResults
result = DiffResults.JacobianResult(ic0)

@time dx = ForwardDiff.jacobian!(result,sum_of_solution,ic0)
@time dx = ForwardDiff.jacobian(sum_of_solution,ic0)
DiffResults.value(result)

_prob = remake(prob,tspan=(0,1e-6),p=[])
solve(_prob,Rodas4(),reltol=1e-6,abstol=1e-6)
