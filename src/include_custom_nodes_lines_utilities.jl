import Base: @__doc__
using PowerDynamics
import PowerDynamics: AbstractLine, AbstractNode
## those are needed for including custom nodes directly instead copying them into PowerDynamics.jö
import PowerDynamics: dimension, symbolsof, construct_vertex, construct_edge
import NetworkDynamics: ODEVertex, StaticEdge
##############
#Linesd
include("lines/DynamicPowerTransformer.jl")
include("lines/PiModelLineParam.jl")
include("lines/StaticPowerTransformer.jl")
include("lines/StaticPowerTransformerTapParam.jl")

#Nodes
#include("nodes/OLTC.jl")
include("nodes/SlackAlgebraicParam.jl")
include("nodes/SimpleRecoveryLoad.jl")
include("nodes/SimpleRecoveryLoadParam.jl")
include("nodes/SixOrderMarconatoMachine.jl")
include("nodes/SixOrderMarconatoMachineAVROEL.jl")
include("nodes/GridFormingConverter.jl")
include("nodes/GridFormingConverterCSA.jl")
include("nodes/GridFormingConverterCSAAntiWindup.jl")
include("nodes/GridFormingConverterParam.jl")
include("nodes/GridSideConverter.jl")
include("nodes/GFMCurrentPrio.jl")
include("nodes/SynchronousMachineGENSAL.jl")
include("nodes/oPFC.jl")
include("nodes/MatchingControl.jl")
include("nodes/MatchingControlRed.jl")
include("nodes/dVOC.jl")
include("nodes/droop.jl")
include("nodes/VSM.jl")
include("nodes/GeneralVoltageDependentLoad.jl")
include("nodes/gentpj.jl")
include("nodes/gentpjAVROEL.jl")
include("nodes/ThreePhaseFault.jl")
include("nodes/ThreePhaseFaultContinouos.jl")


#Operation Point
include("operationpoint/PowerFlow.jl")
include("operationpoint/InitializeInternals.jl")


#Utility functions
include("utility/utility_functions.jl")
