import Base: @__doc__
## those are need for including custom nodes directly instead copying them into PowerDynamics.jl
import PowerDynamics: dimension, symbolsof, construct_vertex
import PowerDynamics: AbstractLine, AbstractNode
import NetworkDynamics: ODEVertex
##############
#Lines
include("lines/DynamicPowerTransformer.jl")
include("lines/PiModelLineParam.jl")
include("lines/StaticPowerTransformer.jl")
include("lines/StaticPowerTransformerTapParam.jl")

#Nodes
include("nodes/GridFormingConverter.jl")
include("nodes/GridFormingConverterCSA.jl")
include("nodes/GridFormingConverterParam.jl")
include("nodes/SimpleRecoveryLoad.jl")
include("nodes/SimpleRecoveryLoadParam.jl")
include("nodes/SixOrderMarconatoMachine.jl")
include("nodes/SixOrderMarconatoMachineAVROEL.jl")
include("nodes/GridFormingConverterCSA.jl")
include("nodes/GridSideConverter.jl")
include("nodes/SlackAlgebraicParam.jl")

#Utility functions
include("utility/utility_functions.jl")

#Operation Point
include("operationpoint/InitializeInternals.jl")
include("operationpoint/PowerFlow.jl")
