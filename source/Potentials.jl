########################################
# DEFINING THE POTENTIALS AND THEIR DERIVATIVES.
########################################
#
#   Isochrone 
#
"""
    psiIsochrone(r::Float64)

Isochrone potential.

# Arguments
- `r`: spherical radius.

# Output
- value of the potential.
"""
function psiIsochrone(r::Float64)
    return -(G * MTot) / (BIsochrone + sqrt(BIsochrone^(2) + r^(2)))
end
"""
    dPsiIsochroneDR(r::Float64)

Derivative of the isochrone potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the derivative.
"""
function dPsiIsochroneDR(r::Float64)
	a = sqrt(BIsochrone^(2) + r^(2))
    return (G * MTot * r) / (a * (BIsochrone + a)^(2))
end
"""
    d2PsiIsochroneDR2(r::Float64)

Second derivative of the isochrone potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the second derivative.
"""
function d2PsiIsochroneDR2(r::Float64)
	b2 = BIsochrone^(2)
	r2 = r^2
	a = sqrt(b2 + r2)
    (G * MTot * (BIsochrone * b2 + b2 * a - 2.0 * r2 * a)) / (a^(3) * (BIsochrone + a)^(3))
end


#
#   Plummer 
#
"""
    psiPlummer(r::Float64)

Plummer potential.

# Arguments
- `r`: spherical radius.

# Output
- value of the potential.
"""
function psiPlummer(r::Float64)
    return -(G * MTot) / sqrt(BPlummer^(2) + r^(2))
end
"""
    dPsiPlummerDR(r::Float64)

Derivative of the Plummer potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the derivative.
"""
function dPsiPlummerDR(r::Float64)
    return (G * MTot * r) / ((BPlummer^(2) + r^(2))^(1.5))
end
"""
    d2PsiPlummerDR2(r::Float64)

Second derivative of the plummer potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the second derivative.
"""
function d2PsiPlummerDR2(r::Float64)
	b2 = BPlummer^(2)
	r2 = r^2
    return G * MTot * (b2 - 2.0 * r2)/((b2 + r2)^(2.5))
end


#
#   Hernquist
#
"""
    psiHernquist(r::Float64)

Hernquist potential.

# Arguments
- `r`: spherical radius.

# Output
- value of the potential.
"""
function psiHernquist(r::Float64)
    return -(G * MTot) / (BHernquist + r)
end
"""
    dPsiHernquistDR(r::Float64)

Derivative of the Hernquist potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the derivative.
"""
function dPsiHernquistDR(r::Float64)
    return (G * MTot) / ((BHernquist + r)^2)
end
"""
    d2PsiHernquistDR2(r::Float64)

Second derivative of the Hernquist potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the second derivative.
"""
function d2PsiHernquistDR2(r::Float64)
    return - (2.0 * G * MTot) / ((BHernquist + r)^3)
end
"""
    rhoHernquist(r::Float64)

Hernquist density.

# Arguments
- `r`: spherical radius.

# Output
- value of the density.
"""
function rhoHernquist(r::Float64)
    return MTot / (2.0 * pi * r * (1.0 + r)^3)
end


#
#   NFW
#
"""
    psiNFW(r::Float64)

NFW potential. The profile is cutoff at a maximum radius RMaxNFW. The total mass within is fixed at Mtot. These conditions fix the value of Rho0NFW for given Mtot, RMaxNFW, ANFW.

# Arguments
- `r`: spherical radius.

# Output
- value of the potential.
"""
const Rho0NFW = MTot / (4.0 * pi * ANFW^3 * (log(1.0 + RMaxNFW / ANFW) - RMaxNFW / ANFW / (1.0 + RMaxNFW / ANFW)))

function psiNFW(r::Float64)
	rescaledR = r / ANFW
    return -(4.0 * pi * G * Rho0NFW * ANFW^2) * log(1.0 + rescaledR) / rescaledR
end
"""
    dPsiNFWDR(r::Float64)

Derivative of the NFW potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the derivative.
"""
function dPsiNFWDR(r::Float64)
	rescaledR = r / ANFW
    return -(4.0 * pi * G * Rho0NFW * ANFW) * (rescaledR / (1.0 + rescaledR) - log(1.0 + r)) / rescaledR^2 
end
"""
    d2PsiNFWDR2(r::Float64)

Second derivative of the NFW potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the second derivative.
"""
function d2PsiNFWDR2(r::Float64)
	rescaledR = r / ANFW
    return -(4.0 * pi * G * Rho0NFW) * (- rescaledR * (2.0 + 3.0 * rescaledR) / (1.0 + rescaledR)^2 + 2.0 * log(1.0 + rescaledR)) / rescaledR^3
end


#
#   Truncated cusp
#
"""
    psiTruncatedCusp(r::Float64)

Truncated cusp potential: psi(r) = - (1 - exp(-r)) / r.

# Arguments
- `r`: spherical radius.

# Output
- value of the potential.
"""
const MCut = 1.0
const RCut = 1.0
function psiTruncatedCusp(r::Float64)
    rescaledR = r / RCut
    return - G * MCut / RCut * (1.0 - exp(-rescaledR)) / rescaledR
end
"""
    dPsiTruncatedCuspDR(r::Float64)

Derivative of the truncated cusp potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the derivative.
"""
function dPsiTruncatedCuspDR(r::Float64)
    rescaledR = r / RCut
    return - G * MCut / RCut^2 * ((1.0 + rescaledR) * exp(-rescaledR) - 1.0) / rescaledR^2
end
"""
    d2PsiTruncatedCuspDR2(r::Float64)

Second derivative of the truncated cusp potential w.r.t. r.

# Arguments
- `r`: spherical radius.

# Output
- value of the second derivative.
"""
function d2PsiTruncatedCuspDR2(r::Float64)
    rescaledR = r / RCut
    return - G * MCut / RCut^3 * (2.0 / rescaledR^3 - exp(-rescaledR) * (1.0 / rescaledR + 2.0 / rescaledR^2 + 2.0 / rescaledR^3)) 
end


##################################################
# DICTIONARIES CONTAINING THE POTENTIAL AND ITS DERIVATIVES
##################################################
potential = Dict("Isochrone" => psiIsochrone, "Plummer" => psiPlummer, "Hernquist" => psiHernquist, "NFW" => psiNFW, "TruncatedCusp" => psiTruncatedCusp)
dPotentialDR = Dict("Isochrone" => dPsiIsochroneDR, "Plummer" => dPsiPlummerDR, "Hernquist" => dPsiHernquistDR, "NFW" => dPsiNFWDR, "TruncatedCusp" => dPsiTruncatedCuspDR)
d2PotentialDR2 = Dict("Isochrone" => d2PsiIsochroneDR2, "Plummer" => d2PsiPlummerDR2, "Hernquist" => d2PsiHernquistDR2, "NFW" => d2PsiNFWDR2, "TruncatedCusp" => d2PsiTruncatedCuspDR2)


##################################################
# DEFINITION OF THE POTENTIAL
##################################################
const Pot = potential[PotentialType]
const DPotDR = dPotentialDR[PotentialType]
const D2PotDR2 = d2PotentialDR2[PotentialType]

