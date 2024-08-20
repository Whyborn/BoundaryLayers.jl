module BoundaryLayers

using DifferentialEquations, DiffEqCallbacks
import NaNMath

include("FluidModels.jl")
include("BVPMethod.jl")
include("Interface.jl")

end
