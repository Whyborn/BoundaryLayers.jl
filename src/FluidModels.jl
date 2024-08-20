# Author: Lachlan Whyborn
# Last Modified: Sun 18 Aug 2024 11:15:01 AM AEST

# Define the components of the fluid models
struct GasModel{T}
    γ::T    # Ratio of specific heats
    Pr::T   # Prandtl number
    S::T    # Sutherland temperature
    μ::T    # Reference viscosity
    R::T    # Specific gas constant
end

InbuiltModels = Dict(
    :Air => GasModel(1.4, 0.71, 110.4, 1.715e-5, 287.15)
   )

function GasModel(Gas::Symbol)
    """
    Retrieve a gas model from the list of pre-defined gas models.
    """
    return InbuiltModels[Gas]
end

function GasModel(GasDict::Dict)
    return GasModel(GasDict[:γ], GasDict[:Pr], GasDict[:S], GasDict[:μ], GasDict[:R])
end

export GasModel

struct CompressibleFreestream{T}
    T::T
    M::T
    ρ::T
    Gas::GasModel
end

export CompressibleFreestream

struct AdiabaticBC{T}
    Tw::T
end

function AdiabaticBC() AdiabaticBC(300.0) end

export AdiabaticBC
