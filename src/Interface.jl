# Author: Lachlan Whyborn
# Last Modified: Tue 20 Aug 2024 09:38:49 PM AEST

function ComputeBoundaryLayer(Freestream::CompressibleFreestream, BC::AdiabaticBC, InitCond::NTuple{2, Float64} = (0.1, 3.0), η_span = (eps(), 5.0))
    """
    Compute the boundary layer profile for the given freestream and boundary conditions.
    """

    # The parameters are T∞, M∞, Pr, Tw, γ, S
    parameters = [Freestream.T,
                  Freestream.M,
                  Freestream.Gas.Pr,
                  BC.Tw,
                  Freestream.Gas.γ,
                  Freestream.Gas.S]

    # Initial conditions
    u₀ = [0.0, 0.0, 0.1, 3.0, 0.0]

    # Positive domain call back
    pos_callback = PositiveDomain(u₀)
    # SS_callback = TerminateSteadyState()

    # Assume a η span for the height- can't start at 0, so use an eps()
    η_span = (eps(), 5.0)
    bvp = DifferentialEquations.TwoPointBVProblem(
                            CompressibleODESystem!,
                            (CompAdiabaticWallBC!, CompAdiabaticFarBC!),
                            u₀,
                            η_span,
                            parameters;
                            bcresid_prototype = (zeros(3), zeros(2)),
                            callback = CallbackSet(pos_callback) 
                            )

    
end

export ComputeBoundaryLayer
