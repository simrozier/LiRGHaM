########################################
# COMPUTING THE LMC'S ORBIT.
########################################


"""
	rhoLMC(x::Float64, y::Float64, z::Float64, xCen::Float64, yCen::Float64, zCen::Float64)

Density of the LMC in theoretical units, off-centered at a cartesian position xCen, yCen, zCen given in theoretical units.
A prefactor is required to rescale the density to that of the LMC (a fraction of 1). A radial rescaling is require to get to the LMC radius (a fraction of 1).

# Arguments
- `x,y,z`: cartesian coordinates.
- `xCen,yCen,zCen`: position of the LMC's centre in cartesian coordinates.

# Output
- value of the density.
"""
function rhoLMC(x::Float64, y::Float64, z::Float64, xCen::Float64, yCen::Float64, zCen::Float64)
	return (MLMC / aLMC^3) / RhoUnit * rhoHernquist(sqrt((x - xCen)^2 + (y - yCen)^2 + (z - zCen)^2) / (aLMC / aMW))
end


#
# Orbit of the LMC in a static MW potential.
#
"""
	leapFrogSpherical!(tMin::Float64, tMax::Float64, deltaT::Float64, rPeri::Float64, vPeri::Float64, xvTrajectory::Array{Float64,2})

Function filling the list xvTrajectory with the trajectory of the LMC in the MW potential, computed using a leap-frog integration scheme, considering rPeri and vPeri as the initial conditions.
Here, the MW potential is static.

# Arguments
- `tMin,tMax,deltaT`: definition of the integration time grid. Minimal time, maximal time and time step.
- `rPeri`: LMC's pericentric radius.
- `vPeri`: LMC's tangential velocity at pericentre.
- `xvTrajectory`: list of trajectory points to fill.

# Output
None
"""
function leapFrogSpherical!(tMin::Float64, tMax::Float64, deltaT::Float64, rPeri::Float64, vPeri::Float64, xvTrajectory::Array{Float64,2})
	NumberOfTimeSteps = floor(Int64,(tMax - tMin) / deltaT) + 1
	x0, y0 = rPeri, 0.0
	vx0, vy0 = 0.0, vPeri
	# Initialization
	x, y = x0, y0
	vx, vy = vx0, vy0
	r, v = sqrt(x^2 + y^2), sqrt(vx^2 + vy^2)
	ax, ay = - x / r * DPotDR(r), - y / r * DPotDR(r)
	xvTrajectory[1,:] = [x0 y0 vx0 vy0]
	for i=1:NumberOfTimeSteps-1
		vxh = vx + ax * deltaT / 2.0
		vyh = vy + ay * deltaT / 2.0
		x = x + vxh * deltaT
		y = y + vyh * deltaT
		r = sqrt(x^2 + y^2)
		ax, ay = - x / r * DPotDR(r), - y / r * DPotDR(r)
		vx = vxh + ax * deltaT / 2.0
		vy = vyh + ay * deltaT / 2.0
		xvTrajectory[i + 1,:] = [x y vx vy]
	end
end


#
#Building a reference trajectory for the LMC. Initial time: 41.7742. Final: 61.7141, theoretical units.
# Realistic orbit achieved between time steps 408-608.
#
const LMCReferenceTrajectory = zeros(Float64,NTStepsIntegration,4)
if (version == "MW_LMC_notAcc") && (projectionQ == "computed")
	println("Computing the trajectory.")
	flush(stdout)
	@time leapFrogSpherical!(TMin, TMaxIntegration, DeltaTIntegration, RpLMC / aMW, VpLMC / VUnit, LMCReferenceTrajectory)
	flush(stdout)
	const LMCTrajectory = LMCReferenceTrajectory[408:10:608,1:2]
	LMCTrajectory[:,1] = - LMCTrajectory[:,1]
end


#
# Building an orbit of the LMC while letting the MW undergo a reflex motion.
# At each time step, the MW cusp has the dynamics of a fictitious point mass in  the potential of the LMC. The dynamics of the LMC is that of a small Hernquist sphere in the potential of a larger Hernquist sphere centered at the cusp (two inter-penetrating spheres).
#
"""
	forceSpheres(rLMC::Float64)

Force applied by the LMC on the MW.

# Arguments
- `rLMC`: distance between the LMC and the MW.

# Output
- value of the force.
"""
function forceSpheres(rLMC::Float64)
	integralForce, errorForce = hcubature(Rz -> DPotDR(norm(Rz)) * Rz[2] / norm(Rz) * Rz[1] * rhoLMC(Rz[1],0.0,Rz[2],0.0,0.0,rLMC), (0.0, -10000.0), (10000.0, 10000.0))
	return - 2.0 * pi * integralForce
end
"""
	dPsiLMCDR(r::Float64)

Rescaled version of dPsi/dr, corresponding to the LMC with the same potential as the MW.

# Arguments
- `r`: spherical radius.

# Output
- value of the derivative of the potential.
"""
function dPsiLMCDR(r::Float64)
	return MLMC / MMW * (aMW / aLMC)^2 * DPotDR(r * aMW / aLMC)
end

"""
	leapFrogAccelerated!(tMin::Float64, tMax::Float64, deltaT::Float64, rPeri::Float64, vPeri::Float64, xvTrajectory::Array{Float64,2})

Function filling the list xvTrajectory with the trajectory of the LMC and the MW, computed using a leap-frog integration scheme, considering rPeri and vPeri as the initial conditions for the LMC.
Here, the MW potential is moving.

# Arguments
- `tMin,tMax,deltaT`: definition of the integration time grid. Minimal time, maximal time and time step.
- `rPeri`: LMC's pericentric radius.
- `vPeri`: LMC's tangential velocity at pericentre.
- `xvTrajectory`: list of trajectory points to fill. Each element corresponds to a single time step and is composed of: the 2D phase space of the MW (x, y, vx, vy); the 2D phase space of the LMC (x, y, vx, vy).

# Output
None
"""
function leapFrogAccelerated!(tMin::Float64, tMax::Float64, deltaT::Float64, rPeri::Float64, vPeri::Float64, xvTrajectory::Array{Float64,2})
	NumberOfTimeSteps = floor(Int64,(tMax - tMin) / deltaT) + 1
	x10, y10 = 0.0, 0.0
	vx10, vy10 = 0.0, 0.0
	x20, y20 = rPeri, 0.0
	vx20, vy20 = 0.0, vPeri
	# Initialization
	x1, y1 = x10, y10
	x2, y2 = x20, y20
	x, y = x2 - x1, y2 - y1
	vx1, vy1 = vx10, vy10
	vx2, vy2 = vx20, vy20
	r = sqrt(x^2 + y^2)
	ax1, ay1 = x / r * dPsiLMCDR(r), y / r * dPsiLMCDR(r)
	ax2, ay2 = x / r * MMW / MLMC * forceSpheres(r), y / r * MMW / MLMC * forceSpheres(r)
	xvTrajectory[1,:] = [x10 y10 vx10 vy10 x20 y20 vx20 vy20]
	for i=1:NumberOfTimeSteps-1
		vx1h = vx1 - ax1 * deltaT / 2.0
		vy1h = vy1 - ay1 * deltaT / 2.0
		vx2h = vx2 - ax2 * deltaT / 2.0
		vy2h = vy2 - ay2 * deltaT / 2.0
		x1 = x1 - vx1h * deltaT
		y1 = y1 - vy1h * deltaT
		x2 = x2 - vx2h * deltaT
		y2 = y2 - vy2h * deltaT
		x, y = x2 - x1, y2 - y1
		r = sqrt(x^2 + y^2)
		ax1, ay1 = x / r * dPsiLMCDR(r), y / r * dPsiLMCDR(r)
		ax2, ay2 = x / r * MMW / MLMC * forceSpheres(r), y / r * MMW / MLMC * forceSpheres(r)
		vx1 = vx1h - ax1 * deltaT / 2.0
		vy1 = vy1h - ay1 * deltaT / 2.0
		vx2 = vx2h - ax2 * deltaT / 2.0
		vy2 = vy2h - ay2 * deltaT / 2.0
		xvTrajectory[i + 1,:] = [x1 y1 vx1 vy1 x2 y2 vx2 vy2]
	end
end


#
# Building a new trajectory for the LMC. 
# Realistic orbit achieved between time steps 1604-1804.
#
const LMCAcceleratedTrajectory = zeros(Float64,NTStepsIntegrationAcc,8)
if (version == "MW_LMC_Acc") && (projectionQ == "computed")
	println("Computing the accelerated trajectory.")
	flush(stdout)
	@time leapFrogAccelerated!(TMin, TMaxIntegrationAcc, DeltaTIntegration, RpLMC / aMW, VpLMC / VUnit, LMCAcceleratedTrajectory)
	flush(stdout)
	const LMCAccTrajectory = LMCAcceleratedTrajectory[1604:10:1804,5:6] - LMCAcceleratedTrajectory[1604:10:1804,1:2]
	LMCAccTrajectory[:,1] = - LMCAccTrajectory[:,1]
	const LMCTrajectory = LMCAccTrajectory
end
