# Author: Lachlan Whyborn
# Last Modified: Tue 20 Aug 2024 09:38:46 PM AEST

# This is in the form provided in
# "A CFD Tutorial in Julia: Introduction to Compressible Laminar Boundary-Layer Flows"
function CompressibleODESystem!(du, u, p, η)
    # Destructure the inputs to match equations
    T∞, M∞, Pr, Tw, γ, S = p
    y1, y2, y3, y4, y5 = u

    # Equation 71
    du[1] = y2
    du[2] = y3
    du[3] = -y3 * ((y5 / (2 * y4)) - (y5 / (y4 + S / T∞))) - y1 * y3 * ((y4 + S / T∞) / (NaNMath.sqrt(y4) * (1 + S / T∞)))
    du[4] = y5
    du[5] = -y5^2 * (1 / (2 * y4) - (1 / (y4 + S / T∞))) - Pr * y1 * y5 * (y4 + S / T∞) / (NaNMath.sqrt(y4) * (1 + S / T∞)) - (γ - 1) * Pr * (M∞ * y3)^2
end

# Possible boundary conditions
# Adiabatic BCs, Eqs 72-76
function CompAdiabaticWallBC!(residual, u, p)
    residual[1] = u[1]
    residual[2] = u[2]
    residual[3] = u[5]
end

function CompAdiabaticFarBC!(residual, u, p)
    residual[1] = u[2] - 1
    residual[2] = u[4] - 1
end

# Isothermal BCs, Eqs 77-81
function CompIsothermalWallBC!(residual, u, p)
    T∞, _, _, Tw, _, _ = p

    residual[1] = u[1]
    residual[2] = u[2]
    residual[3] = u[4] - Tw / T∞
end

function CompIsothermalFarBC!(residual, u, p)
    residual[1] = u[2] - 1
    residual[2] = u[4] - 1
end
