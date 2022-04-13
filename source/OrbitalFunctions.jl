########################################
# DEFINITION OF SOME CONSERVED QUANTITIES OF THE ORBITS
# e.g. energy and angular momentum, frequencies, and their derivatives.
########################################


#
# Energy and angular momentum.
#
"""
	eLFromRpRa(rp::Float64,ra::Float64)

Energy and angular momentum as a function of the peri and apocentre.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `e,l`: orbital energy and angular momentum.
"""
function eLFromRpRa(rp::Float64,ra::Float64)
    squareRpRa = (rp / ra)^2
    oneMinusSquareRpRa = 1.0 - squareRpRa
    prp = Pot(rp)
    pra = Pot(ra)
    #####
    e = (pra - squareRpRa * prp) / oneMinusSquareRpRa
    l = rp * sqrt(2.0 * (pra - prp) / oneMinusSquareRpRa)
    #####
    return e, l
end
"""
	dEDRpRa(rp::Float64,ra::Float64)

First derivatives of the energy w.r.t. rp, ra.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `dEDRp,dEDRa`: first derivatives of the energy w.r.t. rp, ra.
"""
function dEDRpRa(rp::Float64,ra::Float64)
    squareRpRa = (rp / ra)^2
    oneMinusSquareRpRa = 1.0 - squareRpRa
    oneMinusSquareRpRaSquared = oneMinusSquareRpRa^2
    diffPotRpRa = Pot(ra) - Pot(rp)
    #####
    dEDRp = squareRpRa * (2.0 / rp * diffPotRpRa - oneMinusSquareRpRa * DPotDR(rp)) / oneMinusSquareRpRaSquared
    dEDRa = squareRpRa * (- 2.0 / ra * diffPotRpRa + oneMinusSquareRpRa / squareRpRa * DPotDR(ra)) / oneMinusSquareRpRaSquared
    #####
    return dEDRp, dEDRa
end
"""
	dLDRpRa(rp::Float64,ra::Float64)

First derivatives of the angular momentum w.r.t. rp, ra.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `dLDRp,dLDRa`: first derivatives of the angular momentum w.r.t. rp, ra.
"""
function dLDRpRa(rp::Float64,ra::Float64)
    squareRpRa = (rp / ra)^2
    oneMinusSquareRpRa = 1.0 - squareRpRa
    oneMinusSquareRpRaSquared = oneMinusSquareRpRa^2
    diffPotRpRa = Pot(ra) - Pot(rp)
    prefactor = 1.0 * sqrt(oneMinusSquareRpRa / (2.0 * diffPotRpRa))
    #####
    dLDRp = prefactor * (2.0 * diffPotRpRa - rp * oneMinusSquareRpRa * DPotDR(rp)) / oneMinusSquareRpRaSquared
    dLDRa = prefactor * rp / ra * squareRpRa * (-2.0 * diffPotRpRa + ra * oneMinusSquareRpRa / squareRpRa * DPotDR(ra)) / oneMinusSquareRpRaSquared
    #####
    return dLDRp, dLDRa
end
"""
	d2EDRpDRa(rp::Float64,ra::Float64)

Second derivatives of the energy w.r.t. rp, ra. The output is (d2E/drp2, d2E/drpdra, d2E/dra2).

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `d2EDRp2`: second derivative of the energy w.r.t. rp, d2E/drp2.
- `d2EDRpDRaloc`: second derivative of the energy w.r.t. rp and ra, d2E/drpdra.
- `d2EDRa2`: second derivative of the energy w.r.t. ra, d2E/dra2.
"""
function d2EDRpDRa(rp::Float64,ra::Float64)
    squareRpRa = (rp / ra)^2
    oneMinusSquareRpRa = 1.0 - squareRpRa
    oneMinusSquareRpRaCubed = oneMinusSquareRpRa^3
    diffPotRpRa = Pot(ra) - Pot(rp)
    dPotDRp = DPotDR(rp)
    dPotDRa = DPotDR(ra)
    prefactor = 1.0 / (ra^2 * oneMinusSquareRpRaCubed)
    #####
    d2EDRp2 = prefactor * (2.0 * (1.0 + 3.0 * squareRpRa) * diffPotRpRa - rp * oneMinusSquareRpRa * (4.0 * dPotDRp + rp * oneMinusSquareRpRa * D2PotDR2(rp)))
    d2EDRpDRaloc = prefactor * 2.0 * rp / ra * (-2.0 * (1.0 + squareRpRa) * diffPotRpRa + oneMinusSquareRpRa * (ra * dPotDRa + rp * dPotDRp))
    d2EDRa2 = prefactor * (2.0 * squareRpRa * (3.0 + squareRpRa) * diffPotRpRa + ra * oneMinusSquareRpRa * (-4.0 * squareRpRa * dPotDRa + ra * oneMinusSquareRpRa * D2PotDR2(ra)))
    #####
    return d2EDRp2, d2EDRpDRaloc, d2EDRa2
end
"""
	d2LDRpDRa(rp::Float64,ra::Float64)

Second derivatives of the angular momentum w.r.t. rp, ra.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `d2LDRp2`: second derivative of the angular momentum w.r.t. rp, d2E/drp2.
- `d2LDRpDRaloc`: second derivative of the angular momentum w.r.t. rp and ra, d2E/drpdra.
- `d2LDRa2`: second derivative of the angular momentum w.r.t. ra, d2E/dra2.
"""
function d2LDRpDRa(rp::Float64,ra::Float64)
    squareRpRa = (rp / ra)^2
    oneMinusSquareRpRa = 1.0 - squareRpRa
    oneMinusSquareRpRaSquared = oneMinusSquareRpRa^2
    diffPotRpRa = Pot(ra) - Pot(rp)
    dPotDRp = DPotDR(rp)
    dPotDRa = DPotDR(ra)
    prefactor = 1.0 / (2.0 * sqrt(2.0)) * sqrt(oneMinusSquareRpRa^(-5) / diffPotRpRa)
    #####
    d2LDRp2 = - prefactor * (-12.0 * rp / ra^2 * diffPotRpRa + 2.0 * oneMinusSquareRpRa * (2.0 * dPotDRp + rp * oneMinusSquareRpRa * D2PotDR2(rp)) + rp * oneMinusSquareRpRaSquared * dPotDRp^2 / diffPotRpRa)
    d2LDRpDRaloc = prefactor * (-12.0 * squareRpRa / ra * diffPotRpRa + 2.0 * oneMinusSquareRpRa * (dPotDRa + rp / ra * squareRpRa * dPotDRp) + rp * oneMinusSquareRpRaSquared * dPotDRp * dPotDRa / diffPotRpRa)
    d2LDRa2 = prefactor * rp / ra * (12.0 * squareRpRa / ra * diffPotRpRa + 2.0 * oneMinusSquareRpRa * (-2.0 * squareRpRa * dPotDRa + ra * oneMinusSquareRpRa * D2PotDR2(ra)) - ra * oneMinusSquareRpRaSquared * dPotDRa^2 / diffPotRpRa)
    #####
    return d2LDRp2, d2LDRpDRaloc, d2LDRa2
end


#
# Radial action, orbital frequencies and their derivatives.
#
"""
	j1FromRpRa(rp::Float64,ra::Float64)

Radial action as a function of rp, ra. 
This function uses adaptive Gauss-Kronrod integration from QuadGK.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- radial action of the orbit.
"""
function j1FromRpRa(rp::Float64,ra::Float64)
	e, l = eLFromRpRa(rp,ra)
	####
	integral, error = quadgk(r -> sqrt(2.0 * (e - Pot(r)) - (l / r)^2), rp, ra, rtol=1.0e-6)
	return 1.0 / pi * integral
end
"""
	anomaly(u::Float64)

Fake anomaly used in the radial integration of the frequencies and their derivatives.
anomaly(u) = u*((3/2) - u^(2)/(2))
Warning: this u is a dummy integration variable, not to be mixed up with the (u,v) grid.

# Arguments
- `u`: dummy variable between -1 and 1.

# Output
- value of the anomaly.
"""
function anomaly(u::Float64)
    u * (1.5 - 0.5 * u^(2))
end
"""
	dAnomalyDU(u::Float64)

Derivative of the anomaly w.r.t. u.

# Arguments
- `u`: dummy variable between -1 and 1.

# Output
- value of the derivative.
"""
function dAnomalyDU(u::Float64)
    1.5 * (1.0 - u^(2))
end
"""
	omega12FromRpRa(rp::Float64,ra::Float64)

Orbital frequencies as a function of rp,ra. 
This function uses adaptive Gauss-Kronrod integration from QuadGK. 
The boundaries are excluded in order to avoid numerical issues, which are sanitised by the change of variables. The Sanitizer is defined as a const.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `omega1,omega2`: orbital frequencies Omega_1 and Omega_2.
"""
const Sanitizer = 1.0e-3
function omega12FromRpRa(rp::Float64,ra::Float64)
	e, l = eLFromRpRa(rp,ra)
	sigmaRpRa = (rp + ra) / 2.0
	deltaRpRa = (ra - rp) / 2.0
	####
	function o1Integrand(u::Float64)
		r = sigmaRpRa + deltaRpRa * anomaly(u)
		return deltaRpRa * dAnomalyDU(u) / sqrt(abs(2.0 * (e - Pot(r)) - (l/r)^2))
	end
	function o2Integrand(u::Float64)
		r = sigmaRpRa + deltaRpRa * anomaly(u)
		return deltaRpRa * dAnomalyDU(u) / (r^2 * sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2)))
	end
	####
	integral1, error1 = quadgk(o1Integrand, -1 + Sanitizer, 1 - Sanitizer, rtol=Sanitizer)
	integral2, error2 = quadgk(o2Integrand, -1 + Sanitizer, 1 - Sanitizer, rtol=Sanitizer)
	####
	omega1 = pi / integral1
	omega2 = omega1 * l * integral2 / pi
	return omega1, omega2
end
"""
	dJ1DRpRa(rp::Float64,ra::Float64)

First derivatives of the radial action w.r.t. rp,ra. 

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `dJ1DRp, dJ1DRa`: first derivatives of the radial action w.r.t. rp and ra.
"""
function dJ1DRpRa(rp::Float64,ra::Float64)
	omega1, omega2 = omega12FromRpRa(rp,ra)
	dEDRp, dEDRa = dEDRpRa(rp,ra)
	dLDRp, dLDRa = dLDRpRa(rp,ra)
	####
	dJ1DRp = (dEDRp - omega2 * dLDRp) / omega1
	dJ1DRa = (dEDRa - omega2 * dLDRa) / omega1
	return dJ1DRp, dJ1DRa
end
"""
	dOmega12DRpRa(rp::Float64,ra::Float64)

First derivatives of the orbital frequencies w.r.t. rp,ra. 

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `dOmega1DRp, dOmega1DRa`: first derivatives of the radial frequency w.r.t. rp and ra.
- `dOmega2DRp, dOmega2DRa`: first derivatives of the angular frequency w.r.t. rp and ra.
"""
function dOmega12DRpRa(rp::Float64,ra::Float64)
	e, l = eLFromRpRa(rp,ra)
	omega1, omega2 = omega12FromRpRa(rp,ra)
	dEDRp, dEDRa = dEDRpRa(rp,ra)
	dLDRp, dLDRa = dLDRpRa(rp,ra)
	sigmaRpRa = (rp + ra) / 2.0
	deltaRpRa = (ra - rp) / 2.0
	####
	function dOmega1DRpIntegrand(u::Float64)
		anom = anomaly(u)
		oneMinusAnom = 1.0 - anom
		r = sigmaRpRa + deltaRpRa * anom
		vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
		dAnomdu = dAnomalyDU(u)
		return 0.5 * (-dAnomdu / vr - deltaRpRa * dAnomdu / vr^(3) * (2.0 * dEDRp - DPotDR(r) * oneMinusAnom - 2.0 * l * dLDRp / r^2 + l^2 / r^3 * oneMinusAnom))
	end
	function dOmega1DRaIntegrand(u::Float64)
		anom = anomaly(u)
		onePlusAnom = 1.0 + anom
		r = sigmaRpRa + deltaRpRa * anom
		vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
		dAnomdu = dAnomalyDU(u)
		return 0.5 * (dAnomdu / vr - deltaRpRa * dAnomdu / vr^(3) * (2.0 * dEDRa - DPotDR(r) * onePlusAnom - 2.0 * l * dLDRa / r^2 + l^2 / r^3 * onePlusAnom))
	end
	function dOmega2DRpIntegrand(u::Float64)
		anom = anomaly(u)
		oneMinusAnom = 1.0 - anom
		r = sigmaRpRa + deltaRpRa * anom
		vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
		r2 = r^2
		r3 = r * r2
		dAnomdu = dAnomalyDU(u)
		return dAnomdu * (- 0.5 / (r2 * vr) - deltaRpRa * oneMinusAnom / (r3 * vr) - deltaRpRa / (r2 * vr^(3)) * (dEDRp - 0.5 * DPotDR(r) * oneMinusAnom - l * dLDRp / r2 + l^2 / (2.0 * r3) * oneMinusAnom))
	end
	function dOmega2DRaIntegrand(u::Float64)
		anom = anomaly(u)
		onePlusAnom = 1.0 + anom
		r = sigmaRpRa + deltaRpRa * anom
		vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
		r2 = r^2
		r3 = r * r2
		dAnomdu = dAnomalyDU(u)
		return dAnomdu * (0.5 / (r2 * vr) - deltaRpRa * onePlusAnom / (r3 * vr) - deltaRpRa / (r2 * vr^(3)) * (dEDRa - 0.5 * DPotDR(r) * onePlusAnom - l * dLDRa / r2 + l^2 / (2.0 * r3) * onePlusAnom))
	end
	####
	integralOmega1Rp, errorOmega1Rp = quadgk(dOmega1DRpIntegrand, -1 + Sanitizer, 1 - Sanitizer, rtol=Sanitizer)
	integralOmega1Ra, errorOmega1Ra = quadgk(dOmega1DRaIntegrand, -1 + Sanitizer, 1 - Sanitizer, rtol=Sanitizer)
	integralOmega2Rp, errorOmega2Rp = quadgk(dOmega2DRpIntegrand, -1 + Sanitizer, 1 - Sanitizer, rtol=Sanitizer)
	integralOmega2Ra, errorOmega2Ra = quadgk(dOmega2DRaIntegrand, -1 + Sanitizer, 1 - Sanitizer, rtol=Sanitizer)
	####
	dOmega1DRp = -omega1^2 / pi * integralOmega1Rp
	dOmega1DRa = -omega1^2 / pi * integralOmega1Ra
	dOmega2DRp = omega2 / omega1 * dOmega1DRp + omega2 / l * dLDRp + omega1 * l / pi * integralOmega2Rp
	dOmega2DRa = omega2 / omega1 * dOmega1DRa + omega2 / l * dLDRa + omega1 * l / pi * integralOmega2Ra
	return dOmega1DRp, dOmega1DRa, dOmega2DRp, dOmega2DRa
end


#
# Jacobians of the successive changes of coordinates.
#
"""
	jacobianELToRpRa(rp::Float64,ra::Float64)

Jacobian of the (E, L)->(rp, ra) change of coordinates.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- Jacobian of the transformation.
"""
function jacobianELToRpRa(rp::Float64,ra::Float64)
	dEDRp, dEDRa = dEDRpRa(rp,ra)
	dLDRp, dLDRa = dLDRpRa(rp,ra)
	return abs(dEDRp * dLDRa - dEDRa * dLDRp)
end
"""
	dJacobianELToRpRaDRpRa(rp::Float64,ra::Float64)

First derivatives of the Jacobian (E, L)->(rp, ra) w.r.t. rp,ra.

# Arguments
- `rp,ra`: orbital peri- and apocentre.

# Output
- `dJacDRp, dJacDRa`: first derivatives of the Jacobian w.r.t. rp and ra.
"""
function dJacobianELToRpRaDRpRa(rp::Float64,ra::Float64)
	dEDRp, dEDRa = dEDRpRa(rp,ra)
	dLDRp, dLDRa = dLDRpRa(rp,ra)
	jacEL = dEDRp * dLDRa - dEDRa * dLDRp
	d2EDRp2, d2EDRpDRaloc, d2EDRa2 = d2EDRpDRa(rp,ra)
	d2LDRp2, d2LDRpDRaloc, d2LDRa2 = d2LDRpDRa(rp,ra)
	#####
	dJacDRp = sign(jacEL) * (d2EDRp2 * dLDRa + dEDRp * d2LDRpDRaloc - d2EDRpDRaloc * dLDRp - dEDRa * d2LDRpDRaloc)
	dJacDRa = sign(jacEL) * (d2EDRpDRaloc * dLDRa + dEDRp * d2LDRa2 - d2EDRa2 * dLDRp - dEDRa * d2LDRpDRaloc)
	return dJacDRp, dJacDRa
end
"""
	jacobianRpRaToUV(u::Float64,v::Float64)

Jacobian of the (rp, ra)->(u, v) change of coordinates.

# Arguments
- `u,v`: new coordinates.

# Output
- Jacobian of the transformation.
"""
function jacobianRpRaToUV(u::Float64,v::Float64)
	return abs(dRpDU(u) * dRpDU(v))
end
"""
	dJacobianRpRaToUVDUV(u::Float64,v::Float64)

First derivatives of the Jacobian (rp, ra)->(u, v) w.r.t. u,v. 

# Arguments
- `u,v`: new coordinates.

# Output
- `dJacDU, dJacDV`: first derivatives of the Jacobian w.r.t. u and v.
"""
function dJacobianRpRaToUVDUV(u::Float64,v::Float64)
	dRpDUloc = dRpDU(u)
	dRpDV = dRpDU(v)
	jacRpRa = dRpDUloc * dRpDV
	#####
	dJacDU = sign(jacRpRa) * d2RpDU2(u) * dRpDV
	dJacDV = sign(jacRpRa) * d2RpDU2(v) * dRpDUloc
	return dJacDU, dJacDV
end

