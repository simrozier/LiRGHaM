########################################
# COMPUTING THE RESPONSE VECTOR AND ITS CORRESPONDING DENSITY.
########################################


"""
	inverseIMinusM!(MMat::Array{Array{Complex{Float64},2},1}, inverse::Array{Array{Complex{Float64},2},1})

Function filling the argument `inverse` with the inverse of the argument `matrix`. `matrix` is given as a list of square matrices, corresponding to the values of the response matrix at different times, and is considered as a single triangular block-Toeplitz matrix to compute the inverse. Since the inverse will also be a triangular block-Toeplitz matrix, the `inverse` also has the structure of a list of square matrices.

# Arguments
- `MMat`: list of square matrices, corresponding to the response matrix at different times.
- `inverse`: list of square matrices to fill.

# Output
None
"""
function inverseIMinusM!(MMat::Array{Array{Complex{Float64},2},1}, inverse::Array{Array{Complex{Float64},2},1})
	inverse[1] = inv(IMat - DeltaT * MMat[1])
	for k=2:GridTimesSize
		matrixSum = zeros(Complex{Float64},NBasisElements,NBasisElements)
		for i=1:k-1
			matrixSum += - DeltaT * MMat[i + 1] * inverse[k - i]
		end
		inverse[k] = - inverse[1] * matrixSum
	end
end
"""
	response!(inverse::Array{Array{Complex{Float64},2},1}, perturber::Array{Array{Complex{Float64},1},1},response::Array{Array{Complex{Float64},1},1})

Function filling the argument `response` with the product of the argument `inverse` minus the identity matrix with `perturber`. `inverse` is given as a list of square matrices, considered as a triangular block-Toeplitz matrix. 

# Arguments
- `inverse`: list of square matrices, corresponding to the matrix [(I-M)^-1] at different times.
- `perturber`: list of vectors corresponding to the perturber. The outer dimension concerns different time steps, while the inner dimension concerns the order of radial basis functions.
- `response`: list of vectors to fill.

# Output
None
"""
function response!(inverse::Array{Array{Complex{Float64},2},1}, perturber::Array{Array{Complex{Float64},1},1},response::Array{Array{Complex{Float64},1},1})
	for k=1:GridTimesSize
		matrixSum = (inverse[1] - IMat) * perturber[k]
		for i=2:k
			matrixSum += inverse[i] * perturber[k - i + 1]
		end
		response[k] = matrixSum
	end
end
"""
	bareResponse!(matrix::Array{Array{Complex{Float64},2},1}, perturber::Array{Array{Complex{Float64},1},1},response::Array{Array{Complex{Float64},1},1})

Function filling the argument `response` with the product of the argument `matrix` with `perturber`. `matrix` is given as a list of square matrices, considered as a triangular block-Toeplitz matrix. 

# Arguments
- `matrix`: list of square matrices, corresponding to the matrix [(I-M)^-1] at different times.
- `perturber`: list of vectors corresponding to the perturber. The outer dimension concerns different time steps, while the inner dimension concerns the order of radial basis functions.
- `response`: list of vectors to fill.

# Output
None
"""
function bareResponse!(matrix::Array{Array{Complex{Float64},2},1}, perturber::Array{Array{Complex{Float64},1},1},response::Array{Array{Complex{Float64},1},1})
	for k=1:GridTimesSize
		matrixSum = matrix[1] * perturber[k]
		for i=2:k
			matrixSum += matrix[i] * perturber[k - i + 1]
		end
		response[k] = matrixSum
	end
end


"""
	deprojectedDensitySingleLM(r::Float64, theta::Float64, phi::Float64, ell::Int64, m::Int64, responseTable::Array{Complex{Float64},1})

Deprojected density corresponding to a vector with a single value of the harmonic numbers ell,m. For m!=0, the deprojection is counted twice to account for the +-m symmetry.

# Arguments
- `r,theta,phi`: spherical radius, polar angle, azimuthal angle in spherical coordinates.
- `ell,m`: hamonic numbers ell,m.
- `responseTable`: vector of values in the space of basis functions, corresponding to a single ell,m harmonic and all radial orders n.

# Output
value of the density at spherical coordinates (r,theta,phi).
"""
function deprojectedDensitySingleLM(r::Float64, theta::Float64, phi::Float64, ell::Int64, m::Int64, responseTable::Array{Complex{Float64},1})
	Ylm = computeYlm(theta, phi, lmax=ell)
	sumdens = 0.0
	for n=1:NBasisElements
		sumdens += responseTable[n] * DBasis(ell, n - 1 + NStart, r)
	end
	sumdens *= sphericalY(ell, m, theta, phi)
	return ifelse(m==0, sumdens, 2.0 * real(sumdens))
end
"""
	deprojectedDensityRadialSingleLM(r::Float64, ell::Int64, m::Int64, responseTable::Array{Complex{Float64},1})

Deprojected radial part of the density, corresponding to a vector with a single value of the harmonic numbers ell,m. This function should be multiplied with spherical harmonics to recover the full density.

# Arguments
- `r`: spherical radius.
- `ell,m`: hamonic numbers ell,m.
- `responseTable`: vector of values in the space of basis functions, corresponding to a single ell,m harmonic and all radial orders n.

# Output
value of the radial part of the density at spherical radius r.
"""
function deprojectedDensityRadialSingleLM(r::Float64, ell::Int64, m::Int64, responseTable::Array{Complex{Float64},1})
	sumdens = 0.0
	for n=1:NBasisElements
		sumdens += responseTable[n] * DBasis(ell, n - 1 + NStart, r)
	end
	return sumdens
end