########################################
# DEFINITION OF THE PHASE SPACE DISTRIBUTION FUNCTIONS.
########################################

#
# Definition of the Osipkov-Merritt Q function
#
"""
	sahaQ(e::Float64, l::Float64)

This function is used in the Osipkov-Merrit distribution functions of all potentials.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- value of Q(e,l).
"""
function sahaQ(e::Float64, l::Float64)
	return -(e + l^2 / (2.0 * Ra^2))
end


#
# Isochrone potential, Osipkov-Merritt DF.
#
"""
	dFIsochroneOsipkovMerritt(e::Float64, l::Float64)

Osipkov-Merritt DF for the isochrone sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- value of the distribution function.
"""
function dFIsochroneOsipkovMerritt(e::Float64, l::Float64)
	qSah = sahaQ(e, l)
	alpha = 1.0 / Ra^2
	sqSah = sqrt(qSah)
	qSah2 = qSah^2
	return ((3.0 * asin(sqSah) * (-9.0 + 17.0 * alpha + (28.0 - 44.0 * alpha) * qSah + (16.0 - 8.0 * alpha) * qSah2)) / sqrt(1.0 - qSah) + sqSah * (27.0 + 77.0 * alpha - (66.0 + 286.0 * alpha) * qSah + (320.0 + 136.0 * alpha) * qSah2 - (240.0 + 32.0 * alpha) * qSah2 * qSah + 64.0 * qSah2 * qSah2)) / (128.0 * sqrt(2.0) * pi^3 * (1.0 - qSah)^4)
end
"""
	dDFIsochroneOsipkovMerrittDEL(e::Float64, l::Float64)

First derivatives of the Osipkov-Merritt DF for the isochrone sphere. 

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `dDFDE, dDFDL`: first derivatives of the DF w.r.t. E and L.
"""
function dDFIsochroneOsipkovMerrittDEL(e::Float64, l::Float64)
	qSah = sahaQ(e, l)
	alpha = 1.0 / Ra^2
	sqrtQ = sqrt(qSah)
	oneMinusQ = 1 - qSah
	oneMinusQSqrt = sqrt(oneMinusQ)
	q2 = qSah^2
	q3 = qSah^3
	q4 = qSah^4
	arcSinSqrtQ = asin(sqrtQ)
	dDFDQ = - 1.0 / (256.0 * sqrt(2) * pi^3 * oneMinusQ^5 * sqrtQ * oneMinusQSqrt) * (oneMinusQSqrt * qSah * (-75.0 - 1318.0 * qSah + 720.0 * q2 - 336.0 * q3 + 64.0 * q4) - 15.0 * sqrtQ * (- 5.0 + 52.0 * qSah + 16.0 * q2) * arcSinSqrtQ + alpha * (oneMinusQSqrt * (- 128.0 + 451.0 * qSah + 774.0 * q2 - 184.0 * q3 + 32.0 * q4) + 15.0 * sqrtQ * (- 13.0 + 68.0 * qSah + 8.0 * q2) * arcSinSqrtQ))
	dDFDE = - dDFDQ
	dDFDL = - l * alpha * dDFDQ
	return dDFDE, dDFDL
end
"""
	d2DFIsochroneOsipkovMerrittDEDL(e::Float64, l::Float64)

Second derivatives of the Osipkov-Merritt DF for the isochrone sphere. 

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `d2DFDE2`: second derivative of the DF w.r.t. E, d2DF/dE2.
- `d2DFDEDL`: second derivative of the DF w.r.t. E and L, d2DF/dEdL.
- `d2DFDL2`: second derivative of the DF w.r.t. L, d2DF/dL2.
"""
function d2DFIsochroneOsipkovMerrittDEDL(e::Float64, l::Float64)
	qSah = sahaQ(e, l)
	alpha = 1.0 / Ra^2
	sqrtQ = sqrt(qSah)
	oneMinusQ = 1 - qSah
	oneMinusQSqrt = sqrt(oneMinusQ)
	q2 = qSah^2
	q3 = qSah^3
	q4 = qSah^4
	arcSinSqrtQ = asin(sqrtQ)
	qSqrtQ = qSah * sqrtQ
	dDFDQ = - 1.0 / (256.0 * sqrt(2) * pi^3 * oneMinusQ^5 * sqrtQ * oneMinusQSqrt) * (oneMinusQSqrt * qSah * (-75.0 - 1318.0 * qSah + 720.0 * q2 - 336.0 * q3 + 64.0 * q4) - 15.0 * sqrtQ * (- 5.0 + 52.0 * qSah + 16.0 * q2) * arcSinSqrtQ + alpha * (oneMinusQSqrt * (- 128.0 + 451.0 * qSah + 774.0 * q2 - 184.0 * q3 + 32.0 * q4) + 15.0 * sqrtQ * (- 13.0 + 68.0 * qSah + 8.0 * q2) * arcSinSqrtQ))
	d2DFDQ2 = 1.0 / (512.0 * pi^3 * sqrt(2) * oneMinusQSqrt * oneMinusQ^6 * qSqrtQ) * (oneMinusQSqrt * q2 * (5409.0 + 5866.0 * qSah - 1248.0 * q2 + 432.0 * q3 - 64.0 * q4) + 105.0 * qSqrtQ * (7.0 + 76.0 * qSah + 16.0 * q2) * arcSinSqrtQ - alpha * (oneMinusQSqrt * (128.0 - 1152.0 * qSah + 7401.0 * q2 + 4618.0 * q3 - 696.0 * q4 + 96.0 * q4 * qSah) + 105.0 * qSqrtQ * (- 1.0 + 92.0 * qSah + 8.0 * q2) * arcSinSqrtQ))
	d2DFDE2 = d2DFDQ2
	d2DFDEDL = l * alpha * d2DFDQ2
	d2DFDL2 = - alpha * dDFDQ + l^2 * alpha^2 * d2DFDQ2
	return d2DFDE2, d2DFDEDL, d2DFDL2
end


#
# Plummer potential, Osipkov-Merritt DF.
#
"""
	dFPlummerOsipkovMerritt(e::Float64, l::Float64)

Osipkov-Merritt DF for the Plummer sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- value of the distribution function.
"""
function dFPlummerOsipkovMerritt(e::Float64, l::Float64)
	qSah = sahaQ(e, l)
	ra2 = Ra^2
	alpha = 1.0 / ra2
	return (24.0 * sqrt(2.0) * (1.0 - alpha + 7.0 / (16.0 * ra2 * qSah^2)) * qSah^3.5) / (7.0 * pi^3)
end
"""
	dDFPlummerOsipkovMerrittDEL(e::Float64, l::Float64)

First derivatives of the Osipkov-Merritt DF for the Plummer sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `dDFDE, dDFDL`: first derivatives of the DF w.r.t. E and L.
"""
function dDFPlummerOsipkovMerrittDEL(e::Float64, l::Float64)
	qSah = sahaQ(e, l)
	ra2 = Ra^2
	alpha = 1.0 / ra2
	dDFDQ = 3.0 * sqrt(qSah) / (2.0 * sqrt(2.0) * pi^3 * ra2) * (3.0 + 16.0 * qSah^2 * (ra2 - 1.0))
	dDFDE = - dDFDQ
	dDFDL = - l * alpha * dDFDQ
	return dDFDE, dDFDL
end
"""
	d2DFPlummerOsipkovMerrittDEDL(e::Float64, l::Float64)

Second derivatives of the Osipkov-Merritt DF for the Plummer sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `d2DFDE2`: second derivative of the DF w.r.t. E, d2DF/dE2.
- `d2DFDEDL`: second derivative of the DF w.r.t. E and L, d2DF/dEdL.
- `d2DFDL2`: second derivative of the DF w.r.t. L, d2DF/dL2.
"""
function d2DFPlummerOsipkovMerrittDEDL(e::Float64, l::Float64)
	qSah = sahaQ(e, l)
	ra2 = Ra^2
	alpha = 1.0 / ra2
	dDFDQ = 3.0 * sqrt(qSah) / (2.0 * sqrt(2.0) * pi^3 * ra2) * (3.0 + 16.0 * qSah^2 * (ra2 - 1.0))
	d2DFDQ2 = (9.0 + 240.0 * qSah^2 * (ra2 - 1.0)) / (4.0 * sqrt(2.0) * pi^3 * sqrt(qSah) * ra2)
	d2DFDE2 = d2DFDQ2
	d2DFDEDL = l * alpha * d2DFDQ2
	d2DFDL2 = - alpha * dDFDQ + l^2 * alpha^2 * d2DFDQ2
	return d2DFDE2, d2DFDEDL, d2DFDL2
end


#
# Plummer potential, Dejonghe 1987 DF.
#
"""
	prefactorPlummerDejonghe(sign::Int64)

Prefactor of the DF, depending on the sign of L^2 / (2.0 * E) + 1.

# Arguments
- `sign`: +1 or -1, depending on the sign of L^2 / (2.0 * E) + 1.

# Output
- value of the prefactor.
"""
function prefactorPlummerDejonghe(sign::Int64)
	if sign == 1
		res = 1.0 / gamma(4.5 - QPlummer)
	else
		res = 1.0 / (gamma(1.0 - QPlummer / 2.0) * gamma(4.5 - QPlummer / 2.0))
	end
	return (3.0 * gamma(6.0 - QPlummer)) / (8.0 * sqrt(2.0) * pi^2.5) * res
end
"""
	dFPlummerDejonghe(e::Float64, l::Float64)

Dejonghe DF for the Plummer sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- value of the distribution function.
"""
function dFPlummerDejonghe(e::Float64, l::Float64)
	x = -l^2 / (2.0 * e)
	if x <= 1
		res = prefactorPlummerDejonghe(1) * _₂F₁(QPlummer / 2.0, QPlummer - 3.5, 1.0, x)
	else
		res = prefactorPlummerDejonghe(-1) * x^(-QPlummer / 2.0) * _₂F₁(QPlummer / 2.0, QPlummer / 2.0, 4.5 - QPlummer / 2.0, 1.0 / x)
	end
	return (-e)^(3.5 - QPlummer) * res
end
"""
	dDFPlummerDejongheDEL(e::Float64, l::Float64)

First derivatives of the Dejonghe DF for the Plummer sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `dDFDE, dDFDL`: first derivatives of the DF w.r.t. E and L.
"""
function dDFPlummerDejongheDEL(e::Float64, l::Float64)
	x = -l^2 / (2.0 * e)
	if x <= 1
		dDFDE = prefactorPlummerDejonghe(1) * (2.0 * QPlummer - 7.0) / 2.0 * (- e)^(2.5 - QPlummer) * (_₂F₁(QPlummer - 3.5, QPlummer / 2.0, 1.0, x) + x * QPlummer / 2.0 * _₂F₁(QPlummer / 2.0 + 1.0, QPlummer - 2.5, 2.0, x))
		dDFDLFactor = prefactorPlummerDejonghe(1) * x / (2.0 * l) * QPlummer * (2.0 * QPlummer - 7.0) * _₂F₁(1.0 + QPlummer / 2.0, QPlummer - 2.5, 2.0, x)
	else
		dDFDE = prefactorPlummerDejonghe(-1) / (2.0 * (QPlummer - 9.0)) * (-e)^(2.5 - QPlummer) * x^(-QPlummer / 2.0 - 1.0) * (QPlummer^2 * _₂F₁(1.0 + QPlummer / 2.0, 1.0 + QPlummer / 2.0, 5.5 - QPlummer / 2.0, 1.0 / x) + x * (63.0 - 16.0 * QPlummer + QPlummer^2) * _₂F₁(QPlummer / 2.0, QPlummer / 2.0, 4.5 - QPlummer / 2.0, 1.0 / x))
		dDFDLFactor = prefactorPlummerDejonghe(-1) * QPlummer / (l * (QPlummer - 9.0)) * x^(- 1.0 - QPlummer / 2.0) * (QPlummer * _₂F₁(1.0 + QPlummer / 2.0, 1.0 + QPlummer / 2.0, 5.5 - QPlummer / 2.0, 1.0 / x) - x * (QPlummer - 9.0) * _₂F₁(QPlummer / 2.0, QPlummer / 2.0, 4.5 - QPlummer / 2.0, 1.0 / x))
	end
	return dDFDE, dDFDLFactor * (-e)^(3.5 - QPlummer)
end
"""
	d2DFPlummerDejongheDEDL(e::Float64, l::Float64)

Second derivatives of the Dejonghe DF for the Plummer sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `d2DFDE2`: second derivative of the DF w.r.t. E, d2DF/dE2.
- `d2DFDEDL`: second derivative of the DF w.r.t. E and L, d2DF/dEdL.
- `d2DFDL2`: second derivative of the DF w.r.t. L, d2DF/dL2.
"""
function d2DFPlummerDejongheDEDL(e::Float64, l::Float64)
	x = -l^2 / (2.0 * e)
	if x <= 1
		d2DFDE2 = prefactorPlummerDejonghe(1) / 4.0 * (2.0 * QPlummer - 7.0) * (2.0 * QPlummer - 5.0) * (- e)^(1.5 - QPlummer) * (_₂F₁(QPlummer - 3.5, QPlummer / 2.0, 1.0, x) + x * QPlummer * _₂F₁(QPlummer / 2.0 + 1.0, QPlummer - 2.5, 2.0, x) + x^2 * QPlummer / 4.0 * (1.0 + QPlummer / 2.0) * _₂F₁(QPlummer / 2.0 + 2.0, QPlummer - 1.5, 3.0, x))
		d2DFDEDL = prefactorPlummerDejonghe(1) / 32.0 * (- e)^(1.5 - QPlummer) * l * QPlummer * (2.0 * QPlummer - 7.0) * (2.0 * QPlummer - 5.0) * (4.0 * _₂F₁(1.0 + QPlummer / 2.0, QPlummer - 2.5, 2.0, x) + x * (QPlummer + 2.0) * _₂F₁(QPlummer / 2.0 + 2.0, QPlummer - 1.5, 3.0, x))
		d2DFDL2Factor = prefactorPlummerDejonghe(1) / (16.0 * (-e)) * QPlummer * (2.0 * QPlummer - 7.0) * (4.0 * _₂F₁(1.0 + QPlummer / 2.0, QPlummer - 2.5, 2.0, x) - x * (10.0 + QPlummer - 2.0 * QPlummer^2) * _₂F₁(2.0 + QPlummer / 2.0, QPlummer - 1.5, 3.0, x))
	else
		d2DFDE2 = prefactorPlummerDejonghe(-1) / (4.0 * (QPlummer - 9.0) * (QPlummer - 11.0)) * (-e)^(1.5 - QPlummer) * x^(-QPlummer / 2.0 - 2.0) * (QPlummer^2 * (2.0 + QPlummer)^2 * _₂F₁(2.0 + QPlummer / 2.0, 2.0 + QPlummer / 2.0, 6.5 - QPlummer / 2.0, 1.0 / x) + 2.0 * x * QPlummer^2 * (77.0 - 18.0 * QPlummer + QPlummer^2) * _₂F₁(1.0 + QPlummer / 2.0, 1.0 + QPlummer / 2.0, 5.5 - QPlummer / 2.0, 1.0 / x) + x^2 * (3465.0 - 1888.0 * QPlummer + 374.0 * QPlummer^2 - 32.0 * QPlummer^3 + QPlummer^4) * _₂F₁(QPlummer / 2.0, QPlummer / 2.0, 4.5 - QPlummer / 2.0, 1.0 / x))
		d2DFDEDL = prefactorPlummerDejonghe(-1) / (2.0 * l * (QPlummer - 9.0) * (QPlummer - 11.0)) * (-e)^(2.5 - QPlummer) * x^(-QPlummer / 2.0 - 2.0) * QPlummer * (QPlummer * (2.0 + QPlummer)^2 * _₂F₁(2.0 + QPlummer / 2.0, 2.0 + QPlummer / 2.0, 6.5 - QPlummer / 2.0, 1.0 / x) - 9.0 * x * (QPlummer - 11.0) * _₂F₁(1.0 + QPlummer / 2.0, 1.0 + QPlummer / 2.0, 5.5 - QPlummer / 2.0, 1.0 / x) - x^2 * (- 693.0 + 239.0 * QPlummer - 27.0 * QPlummer^2 + QPlummer^3) * _₂F₁(QPlummer / 2.0, QPlummer / 2.0, 4.5 - QPlummer / 2.0, 1.0 / x))
		d2DFDL2Factor = prefactorPlummerDejonghe(-1) * QPlummer / ((- 2.0 * e) * (QPlummer - 9.0) * (QPlummer - 11.0)) * x^(- 3.0 - QPlummer / 2.0) * (QPlummer * (2.0 + QPlummer)^2 * _₂F₁(2.0 + QPlummer / 2.0, 2.0 + QPlummer / 2.0, 6.5 - QPlummer / 2.0, 1.0 / x) - x * QPlummer * ( - 33.0 - 19.0 * QPlummer + 2.0 * QPlummer^2) * _₂F₁(1.0 + QPlummer / 2.0, 1.0 + QPlummer / 2.0, 5.5 - QPlummer / 2.0, 1.0 / x) + x^2 * (99.0 + 79.0 * QPlummer - 19.0 * QPlummer^2 + QPlummer^3) * _₂F₁(QPlummer / 2.0, QPlummer / 2.0, 4.5 - QPlummer / 2.0, 1.0 / x))
	end
	return d2DFDE2, d2DFDEDL, d2DFDL2Factor  * (-e)^(3.5 - QPlummer)
end


#
# Hernquist potential, Baes & Van Hese 2007 DF with beta_0 = beta_\infty = Beta.
#
"""
	prefactorHernquistBaes()

Prefactor of the DF.

# Arguments
None

# Output
- value of the prefactor.
"""
function prefactorHernquistBaes()
	return (2.0^(-2.5 + Beta) * BHernquist * (G * MTot)^(-4.0 + 2.0 * Beta) * gamma(5.0 - 2.0 * Beta)) / (pi^2.5 * gamma(1.0 - Beta) * gamma(3.5 - Beta))
end
"""
	dFHernquistBaes(e::Float64, l::Float64)

Baes & Van Hese DF for the Hernquist sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- value of the distribution function.
"""
function dFHernquistBaes(e::Float64, l::Float64)
	return ((-e)^(2.5 - Beta) * prefactorHernquistBaes() * _₂F₁(1.0 - 2.0 * Beta, 5.0 - 2.0 * Beta, 3.5 - Beta, -((BHernquist * e) / (G * MTot))))/l^(2.0 * Beta)
end
"""
	dDFHernquistBaesDEL(e::Float64, l::Float64)

First derivatives of the Baes & Van Hese DF for the Hernquist sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `dDFDE, dDFDL`: first derivatives of the DF w.r.t. E and L.
"""
function dDFHernquistBaesDEL(e::Float64, l::Float64)
	rescaledEnergy = - BHernquist * e / (G * MTot)
	dDFDE = - prefactorHernquistBaes() * (-e)^(1.5 - Beta) * l^(- 2.0 * Beta) * (2.5 - Beta) * (_₂F₁(1.0 - 2.0 * Beta, 5.0 - 2.0 * Beta, 3.5 - Beta, rescaledEnergy) + rescaledEnergy * 2.0 * (1.0 - 2.0 * Beta) / (3.5 - Beta) * _₂F₁(2.0 - 2.0 * Beta, 6.0 - 2.0 * Beta, 4.5 - Beta, rescaledEnergy))
	dDFDL = - prefactorHernquistBaes() * 2.0 * Beta * (-e)^(2.5 - Beta) * l^(-1.0 - 2.0 * Beta) * _₂F₁(1.0 - 2.0 * Beta, 5.0 - 2.0 * Beta, 3.5 - Beta, rescaledEnergy)
	return dDFDE, dDFDL
end
"""
	d2DFHernquistBaesDEDL(e::Float64, l::Float64)

Second derivatives of the Baes & Van Hese DF for the Hernquist sphere.

# Arguments
- `e, l`: orbital energy and angular momentum.

# Output
- `d2DFDE2`: second derivative of the DF w.r.t. E, d2DF/dE2.
- `d2DFDEDL`: second derivative of the DF w.r.t. E and L, d2DF/dEdL.
- `d2DFDL2`: second derivative of the DF w.r.t. L, d2DF/dL2.
"""
function d2DFHernquistBaesDEDL(e::Float64, l::Float64)
	rescaledEnergy = - BHernquist * e / (G * MTot)
	d2DFDE2 = prefactorHernquistBaes() * (-e)^(0.5 - Beta) * l^(- 2.0 * Beta) * (5 - 2.0 * Beta) * ((3.0 - 2.0 * Beta) / 4.0 * _₂F₁(1.0 - 2.0 * Beta, 5.0 - 2.0 * Beta, 3.5 - Beta, rescaledEnergy) + 2.0 * rescaledEnergy * (5 - 2.0 * Beta) * (1.0 - 2.0 * Beta) / (7.0 - 2.0 * Beta) * _₂F₁(2.0 - 2.0 * Beta, 6.0 - 2.0 * Beta, 4.5 - Beta, rescaledEnergy) + 16.0 * rescaledEnergy^2 * (3.0 - Beta) * (1.0 - Beta) * (1.0 - 2.0 * Beta) / ((9.0 - 2.0 * Beta) * (7.0 - 2.0 * Beta)) * _₂F₁(3.0 - 2.0 * Beta, 7.0 - 2.0 * Beta, 5.5 - Beta, rescaledEnergy))
	d2DFDEDL = prefactorHernquistBaes() * (-e)^(1.5 - Beta) * l^(- 1.0 - 2.0 * Beta) * Beta * (5 - 2.0 * Beta) * (_₂F₁(1.0 - 2.0 * Beta, 5.0 - 2.0 * Beta, 3.5 - Beta, rescaledEnergy) + 4.0 * rescaledEnergy * (1.0 - 2.0 * Beta) * _₂F₁(2.0 - 2.0 * Beta, 6.0 - 2.0 * Beta, 4.5 - Beta, rescaledEnergy))
	d2DFDL2 = prefactorHernquistBaes() * 2.0 * Beta * (-e)^(2.5 - Beta) * l^(- 2.0 - 2.0 * Beta) * (1.0 + 2.0 * Beta)* _₂F₁(1.0 - 2.0 * Beta, 5.0 - 2.0 * Beta, 3.5 - Beta, rescaledEnergy)
	return d2DFDE2, d2DFDEDL, d2DFDL2
end


##################################################
# DICTIONARIES CONTAINING THE DFS AND THEIR DERIVATIVES.
##################################################
distributionFunction = Dict("IsochroneOM" => dFIsochroneOsipkovMerritt, "PlummerDejonghe" => dFPlummerDejonghe, "PlummerOM" => dFPlummerOsipkovMerritt, "HERNQUISTBaes" => dFHernquistBaes)
dDFDEL = Dict("IsochroneOM" => dDFIsochroneOsipkovMerrittDEL, "PlummerDejonghe" => dDFPlummerDejongheDEL, "PlummerOM" => dDFPlummerOsipkovMerrittDEL, "HERNQUISTBaes" => dDFHernquistBaesDEL)
d2DFDEDL = Dict("IsochroneOM" => d2DFIsochroneOsipkovMerrittDEDL, "PlummerDejonghe" => d2DFPlummerDejongheDEDL, "PlummerOM" => d2DFPlummerOsipkovMerrittDEDL, "HERNQUISTBaes" => d2DFHernquistBaesDEDL)


##################################################
# DEFINITION OF THE DF.
##################################################
const DF = distributionFunction[DFType]
const DDFDEL = dDFDEL[DFType]
const D2DFDEDL = d2DFDEDL[DFType]

