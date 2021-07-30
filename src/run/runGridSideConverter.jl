using OrderedCollections: OrderedDict
using Plots
using DifferentialEquations
using DataFrames #for CSV
using Distributed
using PowerDynamics
using LightGraphs
using LinearAlgebra
using Roots
using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate, guess


@everywhere using IfElse

include("C:/Users/Micha/Documents/Master-Studium/HiWi/PowerDynamicsDevelopment/PowerDynamicsDevelopment/src/include_custom_nodes_lines_utilities.jl")
include("C:/Users/Micha/Documents/Master-Studium/HiWi/PowerDynamicsDevelopment/PowerDynamicsDevelopment/src/grids/Inverter_Test_System.jl")

@everywhere Ubase = 380e3
@everywhere Sbase = 100e6
@everywhere Zbase = (Ubase^2)/Sbase
powergrid = Inverter_Test_System()

powergrid2 = Inverter_Test_System2()
Qmax   = [Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf,]
Qmin   = -Qmax

U,δ1,ic = PowerFlowClassic(powergrid2, iwamoto = false, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2,max_tol = 1e-8)

Ykk = NodalAdmittanceMatrice(powergrid)
Uc = U.*exp.(1im*δ1/180*pi)
I_c = Ykk*Uc
S = conj(Ykk*Uc).*Uc
pg, ic0 = InitializeInternalDynamics(powergrid,I_c,ic)


problem = ODEProblem{true}(rhs(pg),ic1,tspan)

timespan= (0.0,10.0)
# simulating a tripped line between node 1 and 2
fault2 = LineFailure(line_name="Line_1-2", tspan_fault=(1.,5.))
solution2 = simulate(fault2, pg, ic0, timespan)
plot2 = create_plot(solution2)
savefig(plot2, "/Dokumente/GridSideConverter.pdf")
display(plot2)
