########################################
# DEFINING THE CHANGE OF COORDINATES rp,ra -> u,v.
########################################


#
# Constants characterising the transformation
# 
"""
	uslope: slope of the linear part of rp(u). Also equal to rp(u0).
	u0: value of u at the transition between linear and logarithmic.
	umin, umax: minimum and maximum values of u, corresponding resp. to the maximum and minimum radii rmax, rmin.
"""
const USlope = RBasis / (DeltaUV * NMax * Q)
const U0 = DeltaUV * NMax * Q - 1.0
const UMin = (U0 + 1.0) - RMax / USlope
const UMax = U0 + log(USlope / RMin)


#
# Change of variables (rp,ra)->(u,v)
#
"""
	rpFromU(u::Float64)

Function rp(u).

# Arguments
- `u`: first dimension of the new variables.

# Output
- value of rp(u).
"""
function rpFromU(u::Float64)
	uMinU0 = u - U0
	U0OverSigma = U0 / Sigma
	uMinU0OverSigma = uMinU0 / Sigma
	return USlope / 2.0 * (exp(- uMinU0) * (erf(U0OverSigma) + erf(uMinU0OverSigma)) + (1.0 - uMinU0) * (erfc(U0OverSigma) + erfc(uMinU0OverSigma)))
end
"""
	raFromUV(u::Float64, v::Float64)

Function ra(u,v).

# Arguments
- `u,v`: new variables u,v.

# Output
- value of ra(u,v).
"""
function raFromUV(u::Float64, v::Float64)
	return rpFromU(u) + rpFromU(v)
end


#
# Gradients of the change of variables.
#
"""
	dRpDU(u::Float64)

First derivative of rp w.r.t. u, dRp/du.

# Arguments
- `u`: first dimension of the new variables.

# Output
- value of the derivative drp/du(u).
"""
function dRpDU(u::Float64)
	uMinU0 = u - U0
	U0OverSigma = U0 / Sigma
	uMinU0OverSigma = uMinU0 / Sigma
	return USlope / 2.0 * (- exp(- uMinU0) * (erf(U0OverSigma) + erf(uMinU0OverSigma)) - (erfc(U0OverSigma) + erfc(uMinU0OverSigma)) + 2.0 / (sqrt(pi) * Sigma) * exp(-uMinU0OverSigma^2) * (uMinU0 - 1.0 + exp(- uMinU0)))
end
"""
	d2RpDU2(u)

Second derivative of rp w.r.t. u, d2rp/du2.

# Arguments
- `u`: first dimension of the new variables.

# Output
- value of the second derivative d2rp/du(u).
"""
function d2RpDU2(u::Float64)
	uMinU0 = u - U0
	U0OverSigma = U0 / Sigma
	uMinU0OverSigma = uMinU0 / Sigma
	return USlope / (sqrt(pi) * Sigma^3) * exp(- uMinU0OverSigma^2) * 2.0 * (- exp(- uMinU0) * (uMinU0 + Sigma^2) + Sigma^2 - uMinU0 * (uMinU0 - 1.0)) + USlope / 2.0 * exp(- uMinU0) * (erf(uMinU0OverSigma) + erf(U0OverSigma))
end



