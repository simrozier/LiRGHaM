########################################
# DEFINITION OF THE FULL BASIS, INCLUDING THE ANGULAR HARMONICS.
########################################

#
# Type definition: named tuple for the harmonic numbers.
#
Harmonics = NamedTuple{(:m, :ell, :n),Tuple{Int64, Int64, Int64}}


#
# Defining spherical harmonics.
#
"""
	PlmGeneral(ell,m,x)

Associated Legendre polynomial, P_ell^m(x). This definition extends Plm to m<0.

# Arguments
- `ell`: harmonic number ell.
- `m`: harmonic number m.
- `x`: variable of the polynomial.

# Output
- value of the polynomial.
"""
function PlmGeneral(ell, m, x)
	if m >= 0
		return Plm(ell, m, x)
	else
		return (-1)^m * gamma(ell + m + 1) / gamma(ell - m + 1) * Plm(ell, -m, x)
	end
end
"""
	sphericalY(ell::Int64, m::Int64, theta::Float64, phi::Float64)

Spherical harmonic, Y_ell^m(theta,phi).

# Arguments
- `ell`: harmonic number ell.
- `m`: harmonic number m.
- `theta`: polar angle.
- `phi`: azimuthal angle.

# Output
- value of the spherical harmonic.
"""
function sphericalY(ell::Int64, m::Int64, theta::Float64, phi::Float64)
	return sqrt((2 * ell + 1) * gamma(ell - m + 1) / (4.0 * pi * gamma(ell + m + 1))) * PlmGeneral(ell, m, cos(theta)) * exp(im * m * phi)
end


#
# Defining the full spherical basis.
#
"""
	phiBasis(r::Float64,phi::Float64,theta::Float64,args::Harmonics)

Spherical potential basis.  
WARNING: the order of the variables is r, phi, theta.

# Arguments
- `r`: spherical radius.
- `phi`: azimuthal angle.
- `theta`: polar angle.
- `args`: variable of type Harmonics containing the harmonic numbers, (m,ell,n).

# Output
- value of the potential basis function.
"""
function phiBasis(r::Float64,phi::Float64,theta::Float64,args::Harmonics)
	return UBasis(args.ell, args.n,r) * sphericalY(args.ell,args.m,theta,phi)
end
"""
	rhoBasis(r,phi,theta,args)

Spherical density basis.  
WARNING: the order of the variables is r, phi, theta.

# Arguments
- `r`: spherical radius.
- `phi`: azimuthal angle.
- `theta`: polar angle.
- `args`: variable of type Harmonics containing the harmonic numbers, (m,ell,n).

# Output
- value of the density basis function.
"""
function rhoBasis(r::Float64,phi::Float64,theta::Float64,args::Harmonics)
	return DBasis(args.ell, args.n,r) * sphericalY(args.ell,args.m,theta,phi)
end
