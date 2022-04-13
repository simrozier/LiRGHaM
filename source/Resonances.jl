########################################
# DEFINING LISTS OF RESONANCE VECTORS AND THE FULL RESONANCE MATRIX.
########################################


#
# Functions for the number of resonances and the table of all resonances for given values of ell and n1max. Note that resonances with opposite signs are considered as equivalent.
#
"""
	numberOfResonances(ell::Int64, n1Bound::Int64)

Number of different resonance vectors for a given value of ell and n1Bound.
Resonances with opposite signs are considered as a single resonance.

# Arguments
- `ell`: harmonic number ell.
- `n1Bound`: maximum value of n1, the Fourier number of the radial angle.

# Output
- number of different resonance vectors.
"""
function numberOfResonances(ell::Int64, n1Bound::Int64)
	return (ell + 1) * n1Bound + ceil(Int, ell / 2)
end
"""
	resonanceTable(ell::Int64, n1Bound::Int64)

List of different resonance vectors for a given value of ell and n1Bound.
Resonances with opposite signs are considered as a single resonance.
The vectors are ordered starting with n2=0 if ell is even, then n2=ell downwards.

# Arguments
- `ell`: harmonic number ell.
- `n1Bound`: maximum value of n1, the Fourier number of the radial angle.

# Output
- list of different resonance vectors.
"""
function resonanceTable(ell::Int64, n1Bound::Int64)
	numberOfRes = numberOfResonances(ell, n1Bound)
	resTable = [(0,0) for i=1:numberOfRes]
	accu = 1
	# Starting with n2 = 0 if ell is even.
	if iseven(ell)
		for n1=1:n1Bound
			resTable[accu] = (n1, 0)
			accu+=1
		end
	end
	# Then n2 != 0
	for n1=-n1Bound:n1Bound
		for n2=ell:-2:1
			resTable[accu] = (n1, n2)
			accu+=1
		end
	end
	return resTable
end
"""
	fullResonances(ell::Int64, n1Bound::Int64)

List of different resonance vectors for a given value of ell and n1Bound.
Resonances with opposite signs are considered as a single resonance.
The vectors are ordered starting with n2=0 upwards.

# Arguments
- `ell`: harmonic number ell.
- `n1Bound`: maximum value of n1, the Fourier number of the radial angle.

# Output
- list of different resonance vectors.
"""
function fullResonances(ell::Int64, n1Bound::Int64)
	res = []
	for n2=ell:-2:1
		for n1=-N1Max:N1Max
			append!(res,[(n1,n2)])
		end
	end
	if iseven(ell)
		for n1=1:N1Max
			append!(res,[(n1,0)])
		end
	end
	return res
end


#
# Computing the response matrix for every single resonance.
#
"""
	fillResonanceMatrix!(resonanceMat::Array{Array{Array{Array{Complex{Float64},2},1},1},1})

Function filling the list resonanceMat with the values of the response matrix for all single resonance vector. The dimensions of the list correspond, from the inside out, to: the radial orders (np,nq) [2 dimensions]; the time steps t [1 dimension]; the difference resonance vectors (n1,n2) [1 dimension]; the harmonic number ell [1 dimension].

# Arguments
- `resonanceMat`: list of matrices to fill.

# Output
None
"""
function fillResonanceMatrix!(resonanceMat::Array{Array{Array{Array{Complex{Float64},2},1},1},1})
	for ell=0:EllMax
		resList = resonanceTable(ell, N1Max)
		numRes = numberOfResonances(ell, N1Max)
		for i=1:numRes
			resonance = resList[i]
			tabMMat = [zeros(Complex{Float64},NBasisElements,NBasisElements) for j=1:GridTimesSize]
			println("ell = ", ell, "; resonance = ", resonance, "; Response matrix computation")
			flush(stdout)
			@time tabMMatResonances!(ell, [resonance], tabMMat)
			flush(stdout)
			resonanceMat[ell + 1][i] = tabMMat
		end
	end
end


