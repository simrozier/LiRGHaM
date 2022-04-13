########################################
# DEFINING THE PROJECT PATH; LOADING ALL FILES FOR THE PROJECT.
########################################Data

#
# ACTION REQUIRED: replace the path in "CodeDirectory" with the local path to the parent folder.
#
CodeDirectory = "/LiRGHaM/source/" 
# "./" 


#
# Parameters in the case where the code is run in parallel.
#
const OnCluster = false			# Whether the code is run in parallel on a CPU cluster.
const NThrds = Threads.nthreads()	# Number of parallel threads of the current julia session.


#
# Subsequent paths to the folder specific to the LMC response and corresponding data drops.
#
LMCDirectory = CodeDirectory * "MW_LMC_Response/"
DataDirectory = LMCDirectory * "Data/"


#
# Include all relevant files.
#
println("Loading packages")
flush(stdout)
@time include(CodeDirectory * "Packages.jl")
flush(stdout)
include(CodeDirectory * "utils.jl")
include(LMCDirectory * "MW_LMC_Parameters.jl")
include(LMCDirectory * "LMC_Parameters.jl")
include(CodeDirectory * "Potentials.jl")
include(CodeDirectory * "DistributionFunctions.jl")
include(CodeDirectory * "uvCoordinates.jl")
include(CodeDirectory * "OrbitalFunctions.jl")
include(CodeDirectory * "Grid.jl")
include(CodeDirectory * "RadialBasis.jl")
include(CodeDirectory * "HarmonicBasis.jl")
include(CodeDirectory * "WMat.jl")
include(CodeDirectory * "ResponseMatrix.jl")
include(CodeDirectory * "Response.jl")
include(CodeDirectory * "Resonances.jl")
include(LMCDirectory * "LMC_Perturber.jl")
include(LMCDirectory * "LMC_Projection.jl")
include(LMCDirectory * "LMC_Response.jl")
include(LMCDirectory * "Read_Drop_Data.jl")


########################################
# COMPUTING THE PROJECTION, THE MATRIX AND THE RESPONSE
########################################
const ProjectedLMCTrajectory = [zeros(Complex{Float64}, EllMax + 1, EllMax + 1, NBasisElements) for i=1:GridTimesSize]
#
# First step: project the LMC onto the bi-orthogonal basis. The names of the lists of positions are LMCTrajectory for the fiducial trajectory and LMCAccTrajectory for the accelerated reference frame.
#
if projectionQ == "computed"
	println("Projection of the LMC")
	flush(stdout)
	@time projectLMCTrajectory!(LMCTrajectory, ProjectedLMCTrajectory)
	flush(stdout)
	#
	# Dropping the perturber, in vector and radial forms.
	#
	writePerturberLMCAllEll!(ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
	writePerturberLMCRadial!(ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
elseif projectionQ == "read"
	println("Reading the LMC projection")
	flush(stdout)
	@time readPerturberLMC!(ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
	flush(stdout)
else
	println("No projection computed.")
	flush(stdout)
end
#
# Second step: compute the resonance matrix.
#
const tabMMat = [zeros(Complex{Float64},NBasisElements,NBasisElements) for j=1:GridTimesSize]
const resonanceMatrix = [[[zeros(Complex{Float64},NBasisElements,NBasisElements) for j=1:GridTimesSize] for i=1:numberOfResonances(ell, N1Max)] for ell=0:EllMax]
if resonanceMatrixQ == "computed"
	println("Computing the matrix: ")
	flush(stdout)
	@time fillResonanceMatrix!(resonanceMatrix)
	flush(stdout)
	writeResponseMatrixResonances!(resonanceMatrix)
elseif resonanceMatrixQ == "read"
	println("Reading the matrix: ")
	flush(stdout)
	@time readResponseMatrixResonances!(resonanceMatrix)
	flush(stdout)
else
	println("No matrix computed.")
	flush(stdout)
end


const IMat = zeros(Complex{Float64},NBasisElements,NBasisElements)
for np=1:NBasisElements
    IMat[np,np] = 1.0 + 0.0im 
end
const inverseIMinusM = [zeros(Complex{Float64},NBasisElements,NBasisElements) for i=1:GridTimesSize]
const responseSingleLMResonances = [zeros(Complex{Float64},NBasisElements) for i=1:GridTimesSize]
#
# Computing and dumping the response for some sets of resonances.
#
if responseQ == "computed"
	# Full resonances of each ell.
	for ell=0:EllMax
		computeDumpResponseRadial!(fullResonances(ell, N1Max), resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
		computeDumpBareResponseRadial!(fullResonances(ell,N1Max), resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
		computeDumpResponseCoeff!(fullResonances(ell,N1Max), resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
		computeDumpBareResponseCoeff!(fullResonances(ell,N1Max), resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
	end
	# Single resonances. Considering EllMax and EllMax - 1 allows for all possible single resonances.
	for resonance in resonanceTable(EllMax, N1Max)
		computeDumpBareResponseRadial!([resonance], resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
		computeDumpBareResponseCoeff!([resonance], resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
	end
	for resonance in resonanceTable(EllMax - 1, N1Max)
		computeDumpBareResponseRadial!([resonance], resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
		computeDumpBareResponseCoeff!([resonance], resonanceMatrix, ProjectedLMCTrajectory, (version == "MW_LMC_Acc"))
	end
else
	println("No response computed.")
	flush(stdout)
end


