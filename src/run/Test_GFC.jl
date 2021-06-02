using PowerDynamics
using DifferentialEquations
using Plots
#using ModelingToolkit
#using CSV
#using DataFrames

begin
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/include_costum_nodes_lines_utilities.jl")
    include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/grids/GFC_Test_Grid.jl")
    #include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/nodes/GridFormingConverter.jl")
end
include("C:/Users/liemann/github/PowerDynamicsDevelopment/src/nodes/GridFormingConverter.jl")
pg = GFC_Test_Grid()

function test()
    using(PowerDynamics)
end
