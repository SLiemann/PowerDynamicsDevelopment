using PowerDynamics
using DifferentialEquations
using Plots
#using ModelingToolkit
#using CSV
#using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
end
pg = GFC_Test_Grid()
U,δ,ic0 = PowerFlowClassic(pg)
Ykk = NodalAdmittanceMatrice(pg)
Uc = U.*exp.(1im*δ/180*pi)
I_c = Ykk*Uc
S = conj(Ykk*Uc).*Uc

pg1 ,ic = InitializeInternalDynamics(pg,I_c,ic0)
