########################################
# FUNCTIONS TO COMPUTE THE RESPONSE AND SAVE IT IN DIFFERENT FORMATS.
########################################


"""
	computeDumpResponseRadial!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)

Function taking a list of resonances, a response matrix for each of the single resonances and a perturber, computing the self-gravitating response and saving the radial response density in a file for each pair of harmonics ell, m.
Format: the density is computed and saved on a grid of r in 1D. The full density can be recovered by multiplying with the right spherical harmonic.
The list of resonances should be consistent: all n2 must have the same parity. We use a convention for the list of resonances where n2>=0 (knowing that each resonance is identified with its opposite).

# Arguments
- `resonanceList`: list of resonance vectors (n_1, n_2) on which to compute the response.
- `resonanceMatrix`: response matrix for each ell and each of the individual resonances, ordered as in resonanceTable(ell, N1Max).
- `ProjectedLMCTrajectory`: array of values of the LMC's projection onto the bi-orthogonal basis for each of its positions along its motion.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function computeDumpResponseRadial!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)
	len = length(resonanceList)
	ellmin = maximum([resonanceList[i][2] for i=1:len])
	for ell=ellmin:2:EllMax
		fullResonanceList = resonanceTable(ell, N1Max)
		tabMMat = [zeros(Complex{Float64},NBasisElements,NBasisElements) for i=1:GridTimesSize]
		for i=1:len
			indexInFullList = findall(x->x==resonanceList[i], fullResonanceList)[1]
			tabMMat += resonanceMatrix[ell + 1][indexInFullList]
		end
		inverseIMinusM!(tabMMat, inverseIMinusM)
		for m=ell:-2:0
			perturberLM = [ProjectedLMCTrajectory[i][ell + 1, m + 1, :] for i=1:GridTimesSize]
			response!(inverseIMinusM, perturberLM, responseSingleLMResonances)
			writeResponseDensityRadialSingleLMResonances!(ell, m, resonanceList, responseSingleLMResonances, Acc)
		end
	end
end
"""
	computeDumpBareResponseRadial!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)

Function taking a list of resonances, a response matrix for each of the single resonances and a perturber, computing the bare response and saving the radial response density in a file for each pair of harmonics ell, m.
Format: the density is computed and saved on a grid of r in 1D. The full density can be recovered by multiplying with the right spherical harmonic.
The list of resonances should be consistent: all n2 must have the same parity. We use a convention for the list of resonances where n2>=0 (knowing that each resonance is identified with its opposite).

# Arguments
- `resonanceList`: list of resonance vectors (n_1, n_2) on which to compute the response.
- `resonanceMatrix`: response matrix for each ell and each of the individual resonances, ordered as in resonanceTable(ell, N1Max).
- `ProjectedLMCTrajectory`: array of values of the LMC's projection onto the bi-orthogonal basis for each of its positions along its motion.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function computeDumpBareResponseRadial!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)
	len = length(resonanceList)
	ellmin = maximum([resonanceList[i][2] for i=1:len])
	for ell=ellmin:2:EllMax
		fullResonanceList = resonanceTable(ell, N1Max)
		tabMMat = [zeros(Complex{Float64},NBasisElements,NBasisElements) for i=1:GridTimesSize]
		for i=1:len
			indexInFullList = findall(x->x==resonanceList[i], fullResonanceList)[1]
			tabMMat += resonanceMatrix[ell + 1][indexInFullList]
		end
		for m=ell:-2:0
			perturberLM = [ProjectedLMCTrajectory[i][ell + 1, m + 1, :] for i=1:GridTimesSize]
			bareResponse!(tabMMat, perturberLM, responseSingleLMResonances)
			writeBareResponseDensityRadialSingleLMResonances!(ell, m, resonanceList, responseSingleLMResonances, Acc)
		end
	end
end
"""
	computeDumpResponseCoeff!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)

Function taking a list of resonances, a response matrix for each of the single resonances and a perturber, computing the self-gravitating response and saving the list of vector coefficients of the response for each pair of harmonics ell, m.
The list of resonances should be consistent: all n2 must have the same parity. We use a convention for the list of resonances where n2>=0 (knowing that each resonance is identified with its opposite).

# Arguments
- `resonanceList`: list of resonance vectors (n_1, n_2) on which to compute the response.
- `resonanceMatrix`: response matrix for each ell and each of the individual resonances, ordered as in resonanceTable(ell, N1Max).
- `ProjectedLMCTrajectory`: array of values of the LMC's projection onto the bi-orthogonal basis for each of its positions along its motion.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function computeDumpResponseCoeff!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)
	len = length(resonanceList)
	ellmin = maximum([resonanceList[i][2] for i=1:len])
	for ell=ellmin:2:EllMax
		fullResonanceList = resonanceTable(ell, N1Max)
		tabMMat = [zeros(Complex{Float64},NBasisElements,NBasisElements) for i=1:GridTimesSize]
		for i=1:len
			indexInFullList = findall(x->x==resonanceList[i], fullResonanceList)[1]
			tabMMat += resonanceMatrix[ell + 1][indexInFullList]
		end
		inverseIMinusM!(tabMMat, inverseIMinusM)
		for m=ell:-2:0
			perturberLM = [ProjectedLMCTrajectory[i][ell + 1, m + 1, :] for i=1:GridTimesSize]
			response!(inverseIMinusM, perturberLM, responseSingleLMResonances)
			writeResponseCoeffSingleLMResonances!(ell, m, resonanceList, responseSingleLMResonances, Acc)
		end
	end
end
"""
	computeDumpBareResponseCoeff!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)

Function taking a list of resonances, a response matrix for each of the single resonances and a perturber, computing the bare response and saving the list of vector coefficients of the response for each pair of harmonics ell, m.
The list of resonances should be consistent: all n2 must have the same parity. We use a convention for the list of resonances where n2>=0 (knowing that each resonance is identified with its opposite).

# Arguments
- `resonanceList`: list of resonance vectors (n_1, n_2) on which to compute the response.
- `resonanceMatrix`: response matrix for each ell and each of the individual resonances, ordered as in resonanceTable(ell, N1Max).
- `ProjectedLMCTrajectory`: array of values of the LMC's projection onto the bi-orthogonal basis for each of its positions along its motion.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function computeDumpBareResponseCoeff!(resonanceList, resonanceMatrix::Array{Array{Array{Array{Complex{Float64},2},1},1},1}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)
	len = length(resonanceList)
	ellmin = maximum([resonanceList[i][2] for i=1:len])
	for ell=ellmin:2:EllMax
		fullResonanceList = resonanceTable(ell, N1Max)
		tabMMat = [zeros(Complex{Float64},NBasisElements,NBasisElements) for i=1:GridTimesSize]
		for i=1:len
			indexInFullList = findall(x->x==resonanceList[i], fullResonanceList)[1]
			tabMMat += resonanceMatrix[ell + 1][indexInFullList]
		end
		for m=ell:-2:0
			perturberLM = [ProjectedLMCTrajectory[i][ell + 1, m + 1, :] for i=1:GridTimesSize]
			bareResponse!(tabMMat, perturberLM, responseSingleLMResonances)
			writeBareResponseCoeffSingleLMResonances!(ell, m, resonanceList, responseSingleLMResonances, Acc)
		end
	end
end
