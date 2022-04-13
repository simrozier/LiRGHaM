##################################################
# PROJECTING THE LMC POTENTIAL ONTO THE BI-ORTHOGONAL BASIS TO BUILD THE PERTURBING VECTOR.
##################################################


"""
	projectionOneBasisElement(func, argsf::Tuple{Float64,Float64,Float64}, argsg::Harmonics, atol::Float64=1e-5, rtol::Float64=1e-5)

Function taking a 3D scalar field and projecting it onto an element of the bi-orthogonal basis.

# Arguments
- `func`:  Scalar field to be projected. This should represent the density of a perturber, with the arguments: r, phi, theta, argsCen. Here, (r,phi,theta) are the spherical coordinates (with phi the azimuth and theta the polar angle), and argsCen is a 3-tuple of Float64, indicating the 3D cartesian coordinates of the LMC's centre.
- `argsf`: Value of argsCen, the 3-tuple of coordinates of the LMC's centre.
- `argsg`: Harmonic numbers for the chosen element of the basis.
- `atol, rtol`: accuracy and precision parameters for the integration. 

# Output
- projection of the field func onto the element of the basis given by argsg. This projection is equal to - INT[d^3x func^* phiBasis]
"""
function projectionOneBasisElement(func, argsf::Tuple{Float64,Float64,Float64}, argsg::Harmonics, atol::Float64=1e-5, rtol::Float64=1e-5)
	result, error, probability, neval, fail, nregions = dotprod_spherical(func, phiBasis, 0.0 + 1.0 * Sanitizer, 2.0 * RVirMW / aMW, argsf, argsg, atol, rtol)
	return - (result[1] + im * result[2])
end
"""
	projectionFull!(func, argsf::Tuple{Float64,Float64,Float64}, projection::Array{Complex{Float64}})

Function taking a 3D scalar field and projecting it onto the full bi-orthogonal basis.

# Arguments
- `func`:  Scalar field to be projected. This should represent the density of a perturber, with the arguments: r, phi, theta, argsCen. Here, (r,phi,theta) are the spherical coordinates (with phi the azimuth and theta the polar angle), and argsCen is a 3-tuple of Float64, indicating the 3D cartesian coordinates of the LMC's centre.
- `argsf`: Value of argsCen, the 3-tuple of coordinates of the LMC's centre.
- `projection`: array to fill with the full set of values of the projection. The dimensions correspond to the harmonic numbers ell and m as the first two dimensions, and the radial order n as the third.

# Output
None
"""
function projectionFull!(func, argsf::Tuple{Float64,Float64,Float64}, projection::Array{Complex{Float64}})
	for ell=0:EllMax, m=ell:-2:0, n=0:NMax 
		projection[ell + 1, m + 1, n + 1] = projectionOneBasisElement(func, argsf, Harmonics((m, ell, n)))
	end
end


"""
	LMCdensitySpherical(r::Float64, phi::Float64, theta::Float64, argsCen::Tuple{Float64,Float64,Float64})

Function computing the LMC's density in spherical coordinates, depending on its current position in cartesian coordinates.

# Arguments
- `r,phi,theta`: the spherical coordinates (with phi the azimuth and theta the polar angle).
- `argsCen`: the 3-tuple of coordinates of the LMC's centre.

# Output
- LMC's density in theoretical units.
"""
function LMCdensitySpherical(r::Float64, phi::Float64, theta::Float64, argsCen::Tuple{Float64,Float64,Float64})
		x, y, z = sph2cart(r, phi, theta)
		return rhoLMC(x, y, z, argsCen[1], argsCen[2], argsCen[3])
	end


"""
	projectionLMC!(xCen::Float64, yCen::Float64, zCen::Float64, projectionLMC::Array{Complex{Float64}})

Function taking the position of the LMC's center in cartesian coordinates and projecting it onto a BOB.

# Arguments
- `xCen,yCen,zCen`: cartesian coordinates of the LMC's centre.
- `projectionLMC`: vector to fill with the full set of values of the projection.

# Output
None
"""
function projectionLMC!(xCen::Float64, yCen::Float64, zCen::Float64, projectionLMC::Array{Complex{Float64}})
	projectionFull!(LMCdensitySpherical, (xCen, yCen, zCen), projectionLMC)
end


"""
	projectLMCTrajectory!(LMCPositionTable::Array{Float64,2}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1})

Function taking the trajectory of the LMC as a 2D list of (x,y) coordinates in its orbital plane and projecting it onto a BOB at each time step.

# Arguments
- `LMCPositionTable`: list of 2D cartesian coordinates of the LMC's centre along the time evolution.
- `ProjectedLMCTrajectory`: array to fill with the full set of values of the projection. The outer dimension corresponds to time steps, and the three inner to the ell, m harmonic numbers and the n radial number.

# Output
None
"""
function projectLMCTrajectory!(LMCPositionTable::Array{Float64,2}, ProjectedLMCTrajectory::Array{Array{Complex{Float64},3},1})
	 Threads.@threads for i=1:GridTimesSize
		xLMC, yLMC = LMCPositionTable[i,:]
		if OnCluster == false
			println("Time step = ", i, "; Projection")
		end
		@time projectionLMC!(xLMC, yLMC, 0.0, ProjectedLMCTrajectory[i])
		flush(stdout)
	end
end

