########################################
# DEFINING THE RADIAL BASES
########################################

#
# Radial basis from Weinberg 1989.
# The zeros of the Bessel functions are allocated, so that they are not re-computed.
#
"""
	besselZeros!()

Function filling the table besselZeros with the zeros of the Bessel function J_{ell - 0.5}

# Arguments
None

# Output
None
"""
const besselZeros = zeros(Float64, EllMax + 1, NMax)
function besselZeros!()
	for ell=0:EllMax, n=1:NMax
		besselZeros[ell + 1, n] = besselj_zero(ell - 0.5, n)
	end
end
# Some bug here with EllMax = 10: besselj evaluated at a negative point?
if BasisType == "Bessel"
	besselZeros!()
end
"""
	amplitudeBesselBasis(ell::Int64, n::Int64)

Prefactor of the basis elements. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.

# Output
- value of the amplitude.
"""
function amplitudeBesselBasis(ell::Int64, n::Int64)
	zero = besselZeros[ell + 1, n]
	return - sqrt(8.0 * pi / RBasis) / (zero * sqrt(pi / (2.0 * zero)) * abs(besselj(ell + 0.5, zero)))
end
"""
	uBasisBessel(ell::Int64, n::Int64, r::Float64)

Radial part of the potential basis function. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the potential basis function.
"""
function uBasisBessel(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = besselZeros[ell + 1, n] * r / RBasis
	return amplitudeBesselBasis(ell, n) * sqrt(pi / (2.0 * rescaledRadius)) * besselj(ell + 0.5, rescaledRadius)
end
"""
	dBasisBessel(ell::Int64, n::Int64, r::Float64)

Radial part of the density basis function. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the density basis function.
"""
function dBasisBessel(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = besselZeros[ell + 1, n] * r / RBasis
	return - amplitudeBesselBasis(ell, n) * (besselZeros[ell + 1, n] / RBasis)^2 / (4 * pi) * sqrt(pi / (2.0 * rescaledRadius)) * besselj(ell + 0.5, rescaledRadius)
end
"""
	dUBasisBesselDR(ell::Int64, n::Int64, r::Float64)

First derivative of the potential basis function w.r.t. r, dU/dr. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the derivative of the potential basis function.
"""
function dUBasisBesselDR(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = besselZeros[ell + 1, n] * r / RBasis
	return amplitudeBesselBasis(ell, n) / (2.0 * r) * sqrt(pi / (2.0 * rescaledRadius)) * (- besselj(ell + 0.5, rescaledRadius) + rescaledRadius * (besselj(ell - 0.5, rescaledRadius) - besselj(ell + 1.5, rescaledRadius)))
end


#
# Radial basis from Hernquist & Ostriker, 1992.
#
"""
	amplitudeHernquistBasis(ell::Int64, n::Int64)

Prefactor of the basis elements. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.

# Output
- value of the amplitude.
"""
function amplitudeHernquistBasis(ell::Int64, n::Int64)
	fourEllPlusThree = 4 * ell + 3
	halfFourEllPlusThree = fourEllPlusThree / 2.0
	return 2.0^(fourEllPlusThree) * gamma(halfFourEllPlusThree) * sqrt(G / RBasis * gamma(n + 1) * (n + halfFourEllPlusThree) / (gamma(n + fourEllPlusThree) * (n / 2.0 * (n + fourEllPlusThree) + (ell + 1.0) * (halfFourEllPlusThree - 0.5))))
end
"""
	uBasisHernquist(ell::Int64, n::Int64, r::Float64)

Radial part of the potential basis function. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the potential basis function.
"""
function uBasisHernquist(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = r / RBasis
	rescaledRadiusPlusOne = 1.0 + rescaledRadius
	chi = (rescaledRadius - 1.0) / rescaledRadiusPlusOne
	nEllCombination = n + 4 * ell + 3
	ellPlusOne = ell + 1
	twoEllPlusOne = 2 * ell + 1
	return amplitudeHernquistBasis(ell, n) * rescaledRadius^ell / rescaledRadiusPlusOne^(twoEllPlusOne) * sqrt(pi) / 2.0^(2 * twoEllPlusOne) * gamma(nEllCombination - n) / (gamma(twoEllPlusOne + 0.5) * gamma(2 * ellPlusOne)) * gegenbauer(n, twoEllPlusOne + 0.5, chi)
end
"""
	dBasisHernquist(ell::Int64, n::Int64, r::Float64)

Radial part of the density basis function. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the density basis function.
"""
function dBasisHernquist(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = r / RBasis
	rescaledRadiusPlusOne = 1.0 + rescaledRadius
	return - uBasisHernquist(ell, n, r) / (sqrt(G) * 2.0 * pi * RBasis^2) * (n / 2.0 * (n + 4 * ell + 3) + (ell + 1) * (2 * ell + 1)) / (rescaledRadius * rescaledRadiusPlusOne^2)
end
"""
	dUBasisHernquistDR(ell::Int64, n::Int64, r::Float64)

First derivative of the potential basis function w.r.t. r, dU/dr. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the derivative of the potential basis function.
"""
function dUBasisHernquistDR(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = r / RBasis
	rescaledRadiusPlusOne = 1.0 + rescaledRadius
	chi = (rescaledRadius - 1.0) / rescaledRadiusPlusOne
	nEllCombination = n + 4 * ell + 3
	ellPlusOne = ell + 1
	twoEllPlusOne = 2 * ell + 1
	twoEllPlusThree = twoEllPlusOne + 2
	return amplitudeHernquistBasis(ell, n) / RBasis * rescaledRadius^(ell - 1) / rescaledRadiusPlusOne^(twoEllPlusThree) * sqrt(pi) / 2.0^(2 * twoEllPlusOne) * gamma(nEllCombination - n) / (gamma(twoEllPlusOne + 1.5) * gamma(twoEllPlusThree)) * ((nEllCombination - n)^2 / 2.0 * rescaledRadius * (nEllCombination - n + 1) * gegenbauer(n - 1, twoEllPlusThree - 0.5, chi) - rescaledRadiusPlusOne * (rescaledRadius + ell * (rescaledRadius - 1.0)) * (twoEllPlusOne + 0.5) * 2 * ellPlusOne * gegenbauer(n, twoEllPlusOne + 0.5, chi))
end


##################################################
# Definition of the radial basis from Clutton-Brock, 1973
##################################################
#
# Radial basis from Clutton-Brock, 1973.
#
"""
	gegenbauer(n::Int64, l::Float64, x::Float64)

Definition of the Gegenbauer polynomials, G_n^l(x). 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `x`: variable of the polynomial.

# Output
- value of the polynomial.
"""
function gegenbauer(n::Int64, l::Float64, x::Float64)
	return gamma(2.0 * l + n) / gamma(2.0 * l) * gamma(l + 0.5) / gamma(l + 0.5 + n) * jacobi(x, n, l - 0.5, l - 0.5)
end
"""
	gegenbauerLargeEllN(n::Int64, l::Float64, x::Float64)

Approximation of the Gegenbauer polynomials, used when the arguments are too large for the gamma function.

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `x`: variable of the polynomial.

# Output
- value of the approximation.
"""
function gegenbauerLargeEllN(n::Int64, l::Float64, x::Float64)
	return (2.0 * l + n - 1.0)^(l - 0.5) * ((2.0 * l + n - 1.0) / (l + n - 0.5))^(l + n) * exp(- l + 0.5) * gamma(l + 0.5) / gamma(2.0 * l) * jacobi(x, n, l - 0.5, l - 0.5)
end
"""
	amplitudeCluttonBasis(ell::Int64, n::Int64)

Prefactor of the basis elements. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.

# Output
- value of the amplitude.
"""
function amplitudeCluttonBasis(ell::Int64, n::Int64)
	KCB = 4 * n * (n + 2 * ell + 2) + (2 * ell + 1) * (2 * ell + 3)
	ICB = pi * 2.0^(- 2 * ell - 1) * gamma(n + 2 * ell + 2) / (gamma(n + 1) * (n + ell + 1) * gamma(ell + 1)^2)
	return - sqrt(4.0 * pi * 2.0^(2 * ell + 3) / (KCB * ICB * RBasis))
end
"""
	amplitudeCluttonBasisLargeEllN(ell::Int64, n::Int64)

Approximate amplitude for large ell and n. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.

# Output
- value of the approximate amplitude.
"""
function amplitudeCluttonBasisLargeEllN(ell::Int64, n::Int64)
	KCB = 4 * n * (n + 2 * ell + 2) + (2 * ell + 1) * (2 * ell + 3)
	ICB = pi * 2.0^(- 2 * ell - 1) * (2.0 * ell + n + 1.0)^(2.0 * ell + 1.0) * ((2.0 * ell + n + 1.0) / n)^(n + 0.5) * exp(- 2.0 * ell - 1.0) / ((n + ell + 1) * gamma(ell + 1)^2)
	return - sqrt(4.0 * pi * 2.0^(2 * ell + 3) / (KCB * ICB * RBasis))
end
"""
	uBasisClutton(ell::Int64, n::Int64, r::Float64)

Radial part of the potential basis function. 
When 2.0 * ell + n >= 165.0 (experimental criterion), the approximate versions of the Gegenbauer polynomials and of the amplitude are used.

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the potential basis function.
"""
function uBasisClutton(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = r / RBasis
	r2PlusOne = rescaledRadius^2 + 1.0
	chi = (rescaledRadius^2 - 1.0) / r2PlusOne
	if 2.0 * ell + n >= 165.0
		amplitude = amplitudeCluttonBasisLargeEllN(ell, n)
		gegen = gegenbauerLargeEllN(n, ell + 1.0, chi)
	else
		amplitude = amplitudeCluttonBasis(ell, n)
		gegen = gegenbauer(n, ell + 1.0, chi)
	end
	return amplitude * rescaledRadius^ell / r2PlusOne^(ell + 0.5) * gegen
end
"""
	dBasisClutton(ell::Int64, n::Int64, r::Float64)

Radial part of the density basis function. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the density basis function.
"""
function dBasisClutton(ell::Int64, n::Int64, r::Float64)
	KCB = 4 * n * (n + 2 * ell + 2) + (2 * ell + 1) * (2 * ell + 3)
	rescaledRadius = r / RBasis
	r2PlusOne = rescaledRadius^2 + 1.0
	chi = (rescaledRadius^2 - 1.0) / r2PlusOne
	if 2.0 * ell + n >= 165.0
		amplitude = amplitudeCluttonBasisLargeEllN(ell, n)
		gegen = gegenbauerLargeEllN(n, ell + 1.0, chi)
	else
		amplitude = amplitudeCluttonBasis(ell, n)
		gegen = gegenbauer(n, ell + 1.0, chi)
	end
	return - KCB / (4.0 * pi * RBasis^2) * amplitude * rescaledRadius^ell / r2PlusOne^(ell + 2.5) * gegen
end
"""
	dUBasisCluttonDR(ell::Int64, n::Int64, r::Float64)

First derivative of the potential basis function w.r.t. r, dU/dr. 

# Arguments
- `ell`: harmonic number ell.
- `n`: order of the radial basis function.
- `r`: spherical radius.

# Output
- value of the derivative of the potential basis function.
"""
function dUBasisCluttonDR(ell::Int64, n::Int64, r::Float64)
	rescaledRadius = r / RBasis
	r2PlusOne = rescaledRadius^2 + 1.0
	chi = (rescaledRadius^2 - 1.0) / r2PlusOne
	if 2.0 * ell + n >= 165.0
		amplitude = amplitudeCluttonBasisLargeEllN(ell, n)
		gegen1 = gegenbauerLargeEllN(n - 1, ell + 2.0, chi)
		gegen2 = gegenbauerLargeEllN(n, ell + 1.0, chi)
	else
		amplitude = amplitudeCluttonBasis(ell, n)
		gegen1 = gegenbauer(n - 1, ell + 2.0, chi)
		gegen2 = gegenbauer(n, ell + 1.0, chi)
	end
	return amplitude / RBasis * rescaledRadius^(ell - 1) / r2PlusOne^(ell + 2.5) * (8.0 * (ell + 1) * rescaledRadius^2 * gegen1 - ((1 + ell) * rescaledRadius^4 + rescaledRadius^2 - ell) * gegen2)
end


##################################################
# DICTIONARIES CONTAINING THE RADIAL BASIS AND ITS DERIVATIVE.
##################################################
uBasis = Dict("Bessel" => uBasisBessel, "Hernquist" => uBasisHernquist, "Clutton" => uBasisClutton)
dBasis = Dict("Bessel" => dBasisBessel, "Clutton" => dBasisClutton, "Hernquist" => dBasisHernquist)
dUBasisDR = Dict("Bessel" => dUBasisBesselDR, "Hernquist" => dUBasisHernquistDR, "Clutton" => dUBasisCluttonDR)


##################################################
# DEFINITION OF THE RADIAL BASIS.
##################################################
const UBasis = uBasis[BasisType]
const DBasis = dBasis[BasisType]
const DUBasisDR = dUBasisDR[BasisType]


#
# Two constants which depend on the choice of basis: the initial radial order (NStart) and the total number of elements (NBasisElements).
#
const NStart = begin
	if BasisType == "Bessel"
		1
	else
		0
	end
end
const NBasisElements = NMax + 1 - NStart
