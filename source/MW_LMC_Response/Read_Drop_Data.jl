########################################
# FUNCTIONS USED TO READ AND WRITE DATA INTO FILES IN DIFFERENT FORMATS.
########################################

#
# Unique grid for the radial sampling, when radial data have to be saved.
#
nPointsR = 200
nPointsPhi = 50
rMinPlot = 0.001
rMaxPlot = 10.401


"""
	readPerturberLMC!(perturberTable::Array{Array{Complex{Float64},3},1}, Acc=false)

Function reading the right file "Perturber_LMC_..." and storing the data in the array perturberTable.

# Arguments
- `perturberTable`: array where the perturber vector must be stored.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function readPerturberLMC!(perturberTable::Array{Array{Complex{Float64},3},1}, Acc=false)
	for ell=0:EllMax
		println("Rank of current ell: ", ell, " / ", EllMax)
		flush(stdout)
		@time for m=0:ell
			stringPerturber = readdlm(DataDirectory * "Perturber_LMC" * ifelse(Acc, "_Acc", "") * "_Ell" * string(ell) * "_m" * string(m) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt")
			for t=1:GridTimesSize
				for n=1:NBasisElements
					perturberTable[t][ell + 1, m + 1, n] = parse(Complex{Float64},string(stringPerturber[t, 3*(n-1) + 1])*stringPerturber[t,3*(n-1) + 2]*stringPerturber[t,3*(n-1) + 3])
				end
			end
		end
		flush(stdout)
	end
end


"""
	writePerturberLMC!(ell::Int64, perturberTable::Array{Array{Complex{Float64},3},1}, Acc=false)

Function writing the values of the perturbing vector to a file for one value of ell and all possible m.

# Arguments
- `ell`: harmonic number ell.
- `perturberTable`: array of values of the perturber vector.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writePerturberLMC!(ell::Int64, perturberTable::Array{Array{Complex{Float64},3},1}, Acc=false)
	for m=0:ell
		writedlm(DataDirectory * "Perturber_LMC" * ifelse(Acc, "_Acc", "") * "_Ell" * string(ell) * "_m" * string(m) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt", [perturberTable[i][ell + 1, m + 1, :] for i=1:GridTimesSize])
	end
end
"""
	writePerturberLMCAllEll!(perturberTable::Array{Array{Complex{Float64},3},1}, Acc=false)

Function writing the values of the perturbing vector to a file for all values of ell.

# Arguments
- `perturberTable`: array of values of the perturber vector.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writePerturberLMCAllEll!(perturberTable::Array{Array{Complex{Float64},3},1}, Acc=false)
	for ell=0:EllMax
		writePerturberLMC!(ell, perturberTable, Acc)
	end
end
"""
	writePerturberLMCRadialEllM!(ell::Int64, m::Int64, perturberTable::Array{Array{Complex{Float64},1},1}, Acc=false)

Function computing the radial density from the values of the perturbing vector for a single value of ell, m and writing the radial density values to a file.

# Arguments
- `ell,m`: harmonic numbers ell, m.
- `perturberTable`: array of values of the perturbing vector for a single value of ell,m, all values of the radial number n and the time t.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writePerturberLMCRadialEllM!(ell::Int64, m::Int64, perturberTable::Array{Array{Complex{Float64},1},1}, Acc=false)
	densityRadial = zeros(Complex{Float64},nPointsR + 1, GridTimesSize)
	for time=1:GridTimesSize
		perturberAtTimeT = perturberTable[time]
		for ir=1:nPointsR + 1
			r = rMinPlot + (ir - 1) * (rMaxPlot-rMinPlot)/nPointsR
			densityRadial[ir, time] = deprojectedDensityRadialSingleLM(r, ell, m, perturberAtTimeT)
		end
	end
	writedlm(DataDirectory * "Perturber_LMC_" * ifelse(Acc, "Acc_", "") * "Radial_Ell" * string(ell) * "_m" * string(m) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt", densityRadial)
end
"""
	writePerturberLMCRadial!(projectedTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)

Function computing the radial density from the values of the perturbing vector for all values of ell, m and writing the radial density values to corresponding files.

# Arguments
- `perturberTable`: array of values of the perturber vector.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writePerturberLMCRadial!(projectedTrajectory::Array{Array{Complex{Float64},3},1}, Acc=false)
	for ell=0:EllMax
		for m=0:ell
			perturberTable = [projectedTrajectory[time][ell + 1, m + 1, :] for time=1:GridTimesSize]
			writePerturberLMCRadialEllM!(ell, m, perturberTable, Acc)
		end
	end
end


#
# Loading and dropping the resonance matrix. Some strategy must be put in place for loading, because the matrix can be very large.
#
"""
	parseList!(listOfArrays,listOfArraysParsed)

Function (fine-tuned) taking as input a list of strings containing the list of responce matrices for all resonances, all times and a single value of ell, and parsing that string in the right way.
# Fine tuning 
- each element of the listOfArrays argument is a string with the response matrix for one resonance and all times. In these strings, the 26 first characters are "Array{Complex{Float64},2}[", and the last is "]". We remove both these chunks. 
- then, the string is split by ",". We end up with a list of strings, each carrying a response matrix at a given time. Inside the string, elements are separated by a blank space, and lines are separated by a ";". The matrix begins with "[" and ends with "]".
- therefore, it just remains to eval(Meta.parse()) each of these strings to parse each matrix, and store the value at the right place in listOfArraysParsed.

# Arguments
- `listOfArrays`: list of strings containing the matrix values for a single ell and all resonances.
- `listOfArraysParsed`: list to store the parsed matrices for a single ell and all resonances.

# Output
None
"""
function parseList!(listOfArrays,listOfArraysParsed)
	for j=1:length(listOfArrays)
		arrayNoHead = listOfArrays[j][27:length(listOfArrays[j]) - 1]
		listOfMatricesStrings = split(arrayNoHead, ",")
		println("Resonance ", j, "/", length(listOfArrays))
		@time for k=1:length(listOfMatricesStrings)
			listOfArraysParsed[j][k] = eval(Meta.parse(listOfMatricesStrings[k]))
		end
		flush(stdout)
	end
end
"""
	readResponseMatrixResonances!(resonanceMatrixLoc::Array{Array{Array{Array{Complex{Float64},2},1},1},1})

Function reading matrix values from the "ResponseMatrixAllEllAllResonances..." file.

# Arguments
- `resonanceMatrixLoc`: array to store the matrix values.

# Output
None
"""
function readResponseMatrixResonances!(resonanceMatrixLoc::Array{Array{Array{Array{Complex{Float64},2},1},1},1})
	op = open(DataDirectory * "ResponseMatrixAllEllAllResonances" * "_beta" * string(Beta) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt")
	i = 0
	for line in eachline(op)
		i += 1
		println("Rank of current ell: ", i, " / ", EllMax + 1)
		flush(stdout)
		listOfArrays = split(line, "\t")
		listOfArraysParsed = [[zeros(Complex{Float64},NBasisElements,NBasisElements) for j=1:GridTimesSize] for r=1:length(listOfArrays)]
		@time parseList!(listOfArrays,listOfArraysParsed)
		flush(stdout)
		resonanceMatrix[i] = listOfArraysParsed
		# println("Dimensions: Resonances - ", length(resonanceMatrix[i]), "; Times - ", length(resonanceMatrix[i][1]), "; Elements - ", length(resonanceMatrix[i][1][1]))
		# println("Matrix value: ",resonanceMatrix[i][1][2][1,1])
		println("Reading succeeded for line ", i)
		flush(stdout)
		# listOfArrays = split(responseMatrixString[i], "\t")
		# listOfArraysParsed = [Meta.parse(listOfArrays[j]) |> eval for j=1:length(listOfArrays)]
		# resonanceMatrix[i] = listOfArraysParsed
	end
end
"""
	writeResponseMatrixResonances!(resonanceMatrixLoc::Array{Array{Array{Array{Complex{Float64},2},1},1},1})

Function writing matrix values to the "ResponseMatrixAllEllAllResonances..." file.

# Arguments
- `resonanceMatrixLoc`: array to carrying the matrix values for each ell, each resonance, each time.

# Output
None
"""
function writeResponseMatrixResonances!(resonanceMatrixLoc::Array{Array{Array{Array{Complex{Float64},2},1},1},1})
	writedlm(DataDirectory * "ResponseMatrixAllEllAllResonances" * "_beta" * string(Beta) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt", resonanceMatrixLoc)
end


#
# Dropping the response values. 
#
"""
	writePerturberLMC!(ell::Int64, perturberTable::Array{Array{Complex{Float64},3},1}, Acc=false)

Function which takes an array of tuples and turns it into a string readable by Mathematica, with l instead of "[" and "(", u instead of "," and r instead of "]" and ")".

# Arguments
- `tupleArray`: list of 2-tuples.

# Output
- `stringResult`: string resulting from replacing each "[" and "(" by "l", each "," by "u" and each "]" and ")" by "r".
"""
function stringFromTupleArray(tupleArray)
	len = length(tupleArray)
	stringResult = "l"
	for i=1:len
		n1, n2 = tupleArray[i]
		stringResult *= "l" * string(n1) * "u" * string(n2) * "r"
		if i < len
			stringResult *= "u"
		else
			stringResult *= "r"
		end
	end
	return stringResult
end
"""
	writeResponseDensityRadialSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1})

Function writing the values of the self-gravitating radial response density to a file for one value of ell, m and a given list of resonances. The response density is computed in a 1D grid of radii.

# Arguments
- `ell`: harmonic number ell.
- `m`: 
- `resonanceList`: list of resonances on which to compute the self-gravitating response.
- `responseTable`: array of values of the response vector for each time and all radial numbers n.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writeResponseDensityRadialSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1}, Acc=false)
	densityRadial = zeros(Complex{Float64},nPointsR + 1, GridTimesSize)
	for time=1:GridTimesSize
		responseAtTimeT = responseTable[time]
		for ir=1:nPointsR + 1
			r = rMinPlot + (ir - 1) * (rMaxPlot-rMinPlot)/nPointsR
			densityRadial[ir, time] = deprojectedDensityRadialSingleLM(r, ell, m, responseAtTimeT)
		end
	end
	writedlm(DataDirectory * "ResponseDensityRadial" * ifelse(Acc, "_Acc", "") * "_To_LMC_Ell" * string(ell) * "_m" * string(m) * "_resonances" * stringFromTupleArray(resonanceList) * "_beta" * string(Beta) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt", densityRadial)
end
"""
	writeBareResponseDensityRadialSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1})

Function writing the values of the bare radial response density to a file for one value of ell, m and a given list of resonances. The response density is computed in a 1D grid of radii.

# Arguments
- `ell`: harmonic number ell.
- `m`: 
- `resonanceList`: list of resonances on which to compute the bare response.
- `responseTable`: array of values of the response vector for each time and all radial numbers n.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writeBareResponseDensityRadialSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1}, Acc=false)
	densityRadial = zeros(Complex{Float64},nPointsR + 1, GridTimesSize)
	for time=1:GridTimesSize
		responseAtTimeT = responseTable[time]
		for ir=1:nPointsR + 1
			r = rMinPlot + (ir - 1) * (rMaxPlot-rMinPlot)/nPointsR
			densityRadial[ir, time] = deprojectedDensityRadialSingleLM(r, ell, m, responseAtTimeT)
		end
	end
	writedlm(DataDirectory * "BareResponseDensityRadial" * ifelse(Acc, "_Acc", "") * "_To_LMC_Ell" * string(ell) * "_m" * string(m) * "_resonances" * stringFromTupleArray(resonanceList) * "_beta" * string(Beta) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt", densityRadial)
end
"""
	writeResponseCoeffSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1})

Function writing the values of the self-gravitating response vector coefficients to a file for one value of ell, m and a given list of resonances. 

# Arguments
- `ell`: harmonic number ell.
- `m`: 
- `resonanceList`: list of resonances on which to compute the self-gravitating response.
- `responseTable`: array of values of the response vector for each time and all radial numbers n.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writeResponseCoeffSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1}, Acc=false)
	writedlm(DataDirectory * "ResponseCoeff" * ifelse(Acc, "_Acc", "") * "_To_LMC_Ell" * string(ell) * "_m" * string(m) * "_resonances" * stringFromTupleArray(resonanceList) * "_beta" * string(Beta) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt", permutedims(reshape(hcat(responseTable...), (NBasisElements, GridTimesSize))))
end
"""
	writeBareResponseCoeffSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1})

Function writing the values of the bare response vector coefficients to a file for one value of ell, m and a given list of resonances. 

# Arguments
- `ell`: harmonic number ell.
- `m`: 
- `resonanceList`: list of resonances on which to compute the bare response.
- `responseTable`: array of values of the response vector for each time and all radial numbers n.
- `Acc`: whether the LMC's trajectory is computed in a static (false) or in an accelerated (true) MW reference frame.

# Output
None
"""
function writeBareResponseCoeffSingleLMResonances!(ell::Int64, m::Int64, resonanceList, responseTable::Array{Array{Complex{Float64},1},1}, Acc=false)
	writedlm(DataDirectory * "BareResponseCoeff" * ifelse(Acc, "_Acc", "") * "_To_LMC_Ell" * string(ell) * "_m" * string(m) * "_resonances" * stringFromTupleArray(resonanceList) * "_beta" * string(Beta) * "_nmax" * string(NMax) * "_Ubasis" * string(BasisType) * "_Rbasis" * string(RBasis) * ".txt", permutedims(reshape(hcat(responseTable...), (NBasisElements, GridTimesSize))))
end