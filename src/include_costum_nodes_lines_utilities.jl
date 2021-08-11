import Base: @__doc__
import PowerDynamics: AbstractLine
## those are need for including custom nodes directly instead copying them into PowerDynamics.jö
import PowerDynamics: dimension, symbolsof, construct_vertex
import NetworkDynamics: ODEVertex
##############
#Lines
#include("lines/DynamicPowerTransformer.jl")
#include("lines/PiModelLineParam.jl")
#include("lines/StaticPowerTransformer.jl")
#include("lines/StaticPowerTransformerTapParam.jl")

#Nodes
#include("nodes/OLTC.jl")
#include("nodes/SimpleRecoveryLoad.jl")
#include("nodes/SimpleRecoveryLoadParam.jl")
#include("nodes/SixOrderMarconatoMachine.jl")
#include("nodes/SixOrderMarconatoMachineAVROEL.jl")
include("nodes/GridFormingConverterCSA.jl")
#include("nodes/GridFormingConverterParam.jl")
#Operation Point
include("operationpoint/InitializeInternals.jl")
include("operationpoint/PowerFlow.jl")

#Utility functions
include("utility/utility_functions.jl")
