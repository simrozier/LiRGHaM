########################################
# BUILDING THE U-V GRID CORRESPONDING TO A MIXED LINEAR-LOGARITHMIC GRID IN RP-RA SPACE.
########################################

const NUGridLarge = floor(Int, (UMax - UMin) / DeltaUV + 1.0)
const gridUVInit = [[u,v] for u=UMin:DeltaUV:UMax, v=UMin:DeltaUV:UMax]

"""
	isInUVGrid(uv::Array{Float64})

Function taking a (u,v) pair and determining wether or not they should belong to the grid.
Excluded points are:
- too circular orbits (ra - rp < 10^-3 rp)
- too large orbits (ra > RMax)
- orbits with invalid Q parameters (case of Osipkov-Merritt distribution functions).

# Arguments
- `uv` : pair of coordinates in the u,v plane.

# Output
- boolean value, true if the point should be included in the grid.
"""
function isInUVGrid(uv::Array{Float64})
	u, v = uv
	rp, ra = rpFromU(u), raFromUV(u,v)
	isNotTooCircular = (ra - rp > 1.0e-3 * rp)
	isWithinRMax = (ra < RMax)
	e, l = eLFromRpRa(rp,ra)
	hasValidSahaQ = ((DFType != "IsochroneOM") && (DFType != "PlummerOM")) || (begin
		0.0 <= sahaQ(e, l) <= 1.0
	end)
	return isNotTooCircular && hasValidSahaQ && isWithinRMax
end
"""
	makeGridUV!()

Function filling the array gridUVToFill with the allowed (u,v) values, and storing the total number of points in the grid as the last element of gridUVToFill.

# Arguments
None

# Output
None
"""
const gridUVToFill = zeros(Float64, 2, NUGridLarge * NUGridLarge + 1)
function makeGridUV!()
	gridUVToFill[1, NUGridLarge * NUGridLarge + 1] = 0.0
	for i=1:NUGridLarge * NUGridLarge
		uv = gridUVInit[i]
		if isInUVGrid(uv)
			gridUVToFill[1, NUGridLarge * NUGridLarge + 1] += 1.0
			nUVGridLoc = convert(Int64, gridUVToFill[1, NUGridLarge * NUGridLarge + 1])
			gridUVToFill[:,nUVGridLoc] = uv
		end
	end
end


#
# Modifying the gridUVToFill array and retrieving the number of points.
#
makeGridUV!()
const NUVGrid = convert(Int64, gridUVToFill[1, NUGridLarge * NUGridLarge + 1])
const gridUV = gridUVToFill[:,1:NUVGrid]


#
# Building the grid of (rp,ra).
#
"""
	gridRpRa!()

Function filling the array gridRpRa with the (rp,ra) values.

# Arguments
None

# Output
None
"""
const gridRpRa = zeros(Float64, 2, NUVGrid)
function gridRpRa!()
	for i=1:NUVGrid
		u, v = gridUV[:,i]
		gridRpRa[:, i] = [rpFromU(u), raFromUV(u,v)]
	end
end
gridRpRa!()


#
# Grid of values of VMax for each value of u. This allows us to know where the edge of the integration plane is located.
#
"""
	makeGridVMax!()

Function filling the array GridVMax with the maximum v value present in gridUVToFill for each value of u.

# Arguments
None

# Output
None
"""
const GridVMax = [u for u=UMin:DeltaUV:UMax]
function makeGridVMax!()
	for i=1:NUVGrid
		u, v = gridUV[:,i]
		index = floor(Int64, (u - UMin) / DeltaUV + 1.0)
		GridVMax[index] = maximum(gridUV[:,gridUV[1,:].==u][2,:])
	end
end
makeGridVMax!()



