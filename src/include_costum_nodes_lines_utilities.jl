import Base: @__doc__
import PowerDynamics: AbstractLine
## those are needed for including custom nodes directly instead copying them into PowerDynamics.jö
#import PowerDynamics: dimension, symbolsof, construct_vertex, construct_edge
#import NetworkDynamics: ODEVertex, StaticEdge
##############
#Linesd
#include("lines/DynamicPowerTransformer.jl")
#include("lines/PiModelLineParam.jl")
#include("lines/StaticPowerTransformer.jl")
#include("lines/StaticPowerTransformerTapParam.jl")

#Nodes
#include("nodes/OLTC.jl")
#include("nodes/SlackAlgebraicParam.jl")
#include("nodes/SimpleRecoveryLoad.jl")
#include("nodes/SimpleRecoveryLoadParam.jl")
#include("nodes/SixOrderMarconatoMachine.jl")
#include("nodes/SixOrderMarconatoMachineAVROEL.jl")
#include("nodes/GridFormingConverter.jl")
#include("nodes/GridFormingConverterCSA.jl")
#include("nodes/GridFormingConverterCSAAntiWindup.jl")
#include("nodes/GridFormingConverterParam.jl")
#include("nodes/GridSideConverter.jl")

#Operation Point
include("operationpoint/PowerFlow.jl")
include("operationpoint/InitializeInternals.jl")


#Utility functions
include("utility/utility_functions.jl")
