using PowerDynamics: SlackAlgebraic, FourthOrderEq, VoltageDependentLoad, PiModelLine, StaticLine, Transformer, PowerGrid#, write_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using PowerDynamics: symbolsof, initial_guess, guess, RLLine
using PowerDynamics: StaticPowerTransformer, DynamicPowerTransformer, SixOrderMarconatoMachine,SixOrderMarconatoMachineSin,SixOrderMarconatoMachineAVROEL
using PowerDynamics
using PowerDynamics: rhs, State
using OrderedCollections: OrderedDict
using Plots
import PowerDynamics: PiModel
using DifferentialEquations
using CSV #read PF DataFrames
using DataFrames #for CSV
using Distributed
@everywhere using IfElse

include("operationpoint/PowerFlow.jl")
include("grids/LTVS_Test_System.jl")

pg = LTVS_Test_System()

#Load Flow
Qmax   = [Inf, Inf, Inf,Inf, Inf,sqrt(1-0.9^2)*Inf]
Qmin   = -Qmax
U,δ1,ic = PowerFlowClassic(pg,iwamoto = true, Qmax = Qmax, Qmin = Qmin, Qlimit_iter_check = 2)
Ykk = NodalAdmittanceMatrice(pg)
Uc  = U.*exp.(1im*δ1/180*pi)
S = conj(Ykk*Uc).*Uc
