########################################
# DEFINITION OF THE RUNGE-KUTTA SCHEMES FOR THE COMPUTATION OF WMAT, XMAT AND THEIR DERIVATIVES.
########################################

"""
    WMat(ell::Int64, np::Int64,n1::Int64,n2::Int64,rp::Float64,ra::Float64)

Function WMat, the in-plane angular Fourier transform of the radial basis elements. 

# Arguments
- `ell`: harmonic number ell.
- `np`: order of the radial basis element.
- `n1,n2`: angular Fourier numbers.
- `rp,ra`: orbital peri- and apocentre.

# Output
- `s1`: value of WMat, computed from a Runge-Kutta integration.
"""
function WMat(ell::Int64, np::Int64,n1::Int64,n2::Int64,rp::Float64,ra::Float64)
    e, l = eLFromRpRa(rp,ra)
    omega1, omega2 = omega12FromRpRa(rp,ra)
    sigma, delta = (ra + rp) * 0.5, (ra - rp) * 0.5
    ########################################
    # Step 0 -- Initialisation
    s1, s2, s3 = 0.0, 0.0, 0.0 # State vectors: s1=W; s2=theta1; s3=theta2-psi
    u = - 1.0 + EpsilonWMat # Time for the integral. Avoiding the edges, where vr is ill-defined.
    r = sigma + delta * anomaly(u)
    dTheta1DU = delta * dAnomalyDU(u) * omega1 / sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
    dTheta2DU = (omega2 - l / (r^2)) * dTheta1DU / omega1
    uBasisValue = UBasis(ell, np, r)
    ########################################
    for i=1:NStepsWMat
        ########################################
        # Step 1
        # Velocity vectors for s1, s2, s3.
        k1_1 = HWMat * (1.0 / pi) * dTheta1DU * uBasisValue * cos(n1 * s2 + n2 * s3)
        k2_1 = HWMat * dTheta1DU
        k3_1 = HWMat * dTheta2DU
        ########################################
        # Step 2
        u += 0.5 * HWMat
        r = sigma + delta * anomaly(u) 
        dTheta1DU = delta * dAnomalyDU(u) * omega1 / sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dTheta2DU = (omega2 - l / (r^2)) * dTheta1DU / omega1
        uBasisValue = UBasis(ell, np, r)
        ########################################
        k1_2 = HWMat * (1.0 / pi) * dTheta1DU * uBasisValue * cos(n1 * (s2 + 0.5 * k2_1) + n2 * (s3 + 0.5 * k3_1))
        k2_2 = HWMat * dTheta1DU
        k3_2 = HWMat * dTheta2DU 
        ########################################
        # Step 3
        k1_3 = HWMat * (1.0 / pi) * dTheta1DU * uBasisValue * cos(n1 * (s2 + 0.5 * k2_2) + n2 * (s3 + 0.5 * k3_2)) 
        k2_3 = k2_2 
        k3_3 = k3_2 
        ########################################
        # Step 4
        u += 0.5 * HWMat 
        r = sigma + delta * anomaly(u) 
        dTheta1DU = delta * dAnomalyDU(u) * omega1 / sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dTheta2DU = (omega2 - l / (r^2)) * dTheta1DU / omega1
        uBasisValue = UBasis(ell, np, r)
        ########################################
        k1_4 = HWMat * (1.0 / pi) * dTheta1DU * uBasisValue * cos(n1 * (s2 + k2_3) + n2 * (s3 + k3_3)) 
        k2_4 = HWMat * dTheta1DU 
        k3_4 = HWMat * dTheta2DU 
        ########################################
        # Update
        s1 += (k1_1 + 2.0 * k1_2 + 2.0 * k1_3 + k1_4) / 6.0
        s2 += (k2_1 + 2.0 * k2_2 + 2.0 * k2_3 + k2_4) / 6.0 
        s3 += (k3_1 + 2.0 * k3_2 + 2.0 * k3_3 + k3_4) / 6.0 
        # if s1 == Inf
        #     println("s1 = Inf; value of i, u, r, dtheta1du, dtheta2du: ", i, " ", u, " ", r, " ", dTheta1DU, " ", dTheta2DU)
        # end
    end
    return s1
end


"""
    dWMatDRpRa(ell::Int64, np::Int64,n1::Int64,n2::Int64,rp::Float64,ra::Float64)

First derivatives of WMat w.r.t. rp and ra.

# Arguments
- `ell`: harmonic number ell.
- `np`: order of the radial basis element.
- `n1,n2`: angular Fourier numbers.
- `rp,ra`: orbital peri- and apocentre.

# Output
- `s1P,s1A`: derivatives of WMat w.r.t. rp and ra, dW/drp,dW/dra.
"""
function dWMatDRpRa(ell::Int64, np::Int64,n1::Int64,n2::Int64,rp::Float64,ra::Float64)
    e, l = eLFromRpRa(rp,ra)
    omega1, omega2 = omega12FromRpRa(rp,ra)
    dEDRp, dEDRa = dEDRpRa(rp,ra)
    dLDRp, dLDRa = dLDRpRa(rp,ra)
    dOmega1DRp, dOmega1DRa, dOmega2DRp, dOmega2DRa = dOmega12DRpRa(rp,ra)
    sigma, delta = (ra + rp) * 0.5, (ra - rp) * 0.5
    ########################################
    # Step 0 -- Initialisation
    s1P, s1A, s2, s3, s4P, s4A, s5P, s5A = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 # State vectors: s1P=dW/drp; s1A=dW/dra; s2=theta1; s3=theta2-psi; s4P=dtheta1/drp; s4A=dtheta1/dra; s5P=d(theta2-psi)/drp; ds5A=d(theta-psi)/dra
    u = - 1.0 + EpsilonWMat # Time for the integral. Avoiding the edges, where vr is ill-defined.
    r = sigma + delta * anomaly(u)
    vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
    dRDRp, dRDRa, dRDU = 0.5 * (1.0 - anomaly(u)), 0.5 * (1.0 + anomaly(u)), delta * dAnomalyDU(u)
    d2RDRpDU, d2RDRaDU = - 0.5 * dAnomalyDU(u), 0.5 * dAnomalyDU(u)
    dTheta1DR = omega1 / vr
    dTheta2DR = (omega2 - l / (r^2)) / vr
    dTheta1DU = dRDU * dTheta1DR
    dTheta2DU = dRDU * dTheta2DR
    d2Theta1DRpDR, d2Theta1DRaDR = dOmega1DRp / vr - omega1 / vr^3 * (dEDRp - l / r^2 * dLDRp + (l^2 / r^3 - DPotDR(r)) * dRDRp), dOmega1DRa / vr - omega1 / vr^3 * (dEDRa - l / r^2 * dLDRa + (l^2 / r^3 - DPotDR(r)) * dRDRa)
    d2Theta2DRpDR, d2Theta2DRaDR = (dOmega2DRp - dLDRp / r^2 + 2.0 * l / r^3 * dRDRp) / vr - (omega2 - l / r^2) / vr^3 * (dEDRp - dLDRp * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRp), (dOmega2DRa - dLDRa / r^2 + 2.0 * l / r^3 * dRDRa) / vr - (omega2 - l / r^2) / vr^3 * (dEDRa - dLDRa * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRa)
    d2Theta1DRpDU, d2Theta1DRaDU = d2Theta1DRpDR * dRDU + dTheta1DR * d2RDRpDU, d2Theta1DRaDR * dRDU + dTheta1DR * d2RDRaDU 
    d2Theta2DRpDU, d2Theta2DRaDU = d2Theta2DRpDR * dRDU + dTheta2DR * d2RDRpDU, d2Theta2DRaDR * dRDU + dTheta2DR * d2RDRaDU 
    uBasisValue = UBasis(ell, np, r)
    dUBasisDRValue = DUBasisDR(ell, np, r)
    ########################################
    for istep=1:NStepsWMat 
        ########################################
        # Step 1
        valSin, valCos = sincos(n1 * s2 + n2 * s3)
        prefactorP, prefactorA = n1 * s4P + n2 * s5P, n1 * s4A + n2 * s5A 
        #####
        # Velocity vectors for s1P, s1A, s2, s3, s4P, s4A, s5P, s5A
        k1P_1 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * uBasisValue + dTheta1DU * dRDRp * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorP * valSin)
        k1A_1 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * uBasisValue + dTheta1DU * dRDRa * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorA * valSin)
        k2_1  = HWMat * dTheta1DU 
        k3_1  = HWMat * dTheta2DU 
        k4P_1 = HWMat * d2Theta1DRpDU 
        k4A_1 = HWMat * d2Theta1DRaDU
        k5P_1 = HWMat * d2Theta2DRpDU 
        k5A_1 = HWMat * d2Theta2DRaDU 
        ########################################
        # Step 2
        u += 0.5 * HWMat 
        r = sigma + delta * anomaly(u)
        vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dRDRp, dRDRa, dRDU = 0.5 * (1.0 - anomaly(u)), 0.5 * (1.0 + anomaly(u)), delta * dAnomalyDU(u)
        d2RDRpDU, d2RDRaDU = - 0.5 * dAnomalyDU(u), 0.5 * dAnomalyDU(u)
        dTheta1DR = omega1 / vr
        dTheta2DR = (omega2 - l / (r^2)) / vr
        dTheta1DU = dRDU * dTheta1DR
        dTheta2DU = dRDU * dTheta2DR
        d2Theta1DRpDR, d2Theta1DRaDR = dOmega1DRp / vr - omega1 / vr^3 * (dEDRp - l / r^2 * dLDRp + (l^2 / r^3 - DPotDR(r)) * dRDRp), dOmega1DRa / vr - omega1 / vr^3 * (dEDRa - l / r^2 * dLDRa + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta2DRpDR, d2Theta2DRaDR = (dOmega2DRp - dLDRp / r^2 + 2.0 * l / r^3 * dRDRp) / vr - (omega2 - l / r^2) / vr^3 * (dEDRp - dLDRp * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRp), (dOmega2DRa - dLDRa / r^2 + 2.0 * l / r^3 * dRDRa) / vr - (omega2 - l / r^2) / vr^3 * (dEDRa - dLDRa * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta1DRpDU, d2Theta1DRaDU = d2Theta1DRpDR * dRDU + dTheta1DR * d2RDRpDU, d2Theta1DRaDR * dRDU + dTheta1DR * d2RDRaDU 
        d2Theta2DRpDU, d2Theta2DRaDU = d2Theta2DRpDR * dRDU + dTheta2DR * d2RDRpDU, d2Theta2DRaDR * dRDU + dTheta2DR * d2RDRaDU 
        uBasisValue = UBasis(ell, np, r)
        dUBasisDRValue = DUBasisDR(ell, np, r)
        #####
        valSin, valCos = sincos(n1 * (s2 + 0.5 * k2_1) + n2 * (s3 + 0.5 * k3_1)) 
        prefactorP, prefactorA   = n1 * (s4P + 0.5 * k4P_1) + n2 * (s5P + 0.5 * k5P_1), n1 * (s4A + 0.5 * k4A_1) + n2 * (s5A + 0.5 * k5A_1) 
        #####
        k1P_2 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * uBasisValue + dTheta1DU * dRDRp * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorP * valSin)
        k1A_2 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * uBasisValue + dTheta1DU * dRDRa * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorA * valSin)
        k2_2  = HWMat * dTheta1DU 
        k3_2  = HWMat * dTheta2DU 
        k4P_2 = HWMat * d2Theta1DRpDU 
        k4A_2 = HWMat * d2Theta1DRaDU 
        k5P_2 = HWMat * d2Theta2DRpDU 
        k5A_2 = HWMat * d2Theta2DRaDU 
        ########################################
        # Step 3
        valSin, valCos = sincos(n1 * (s2 + 0.5 * k2_2) + n2 * (s3 + 0.5 * k3_2))
        prefactorP, prefactorA = n1 * (s4P + 0.5 * k4P_2) + n2 * (s5P + 0.5 * k5P_2), n1 * (s4A + 0.5 * k4A_2) + n2 * (s5A + 0.5 * k5A_2)
        #####
        k1P_3 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * uBasisValue + dTheta1DU * dRDRp * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorP * valSin)
        k1A_3 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * uBasisValue + dTheta1DU * dRDRa * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorA * valSin)
        k2_3  = HWMat * dTheta1DU 
        k3_3  = HWMat * dTheta2DU 
        k4P_3 = HWMat * d2Theta1DRpDU 
        k4A_3 = HWMat * d2Theta1DRaDU 
        k5P_3 = HWMat * d2Theta2DRpDU 
        k5A_3 = HWMat * d2Theta2DRaDU 
        ########################################
        # Step 4
        u += 0.5 * HWMat 
        r = sigma + delta * anomaly(u)
        vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dRDRp, dRDRa, dRDU = 0.5 * (1.0 - anomaly(u)), 0.5 * (1.0 + anomaly(u)), delta * dAnomalyDU(u)
        d2RDRpDU, d2RDRaDU = - 0.5 * dAnomalyDU(u), 0.5 * dAnomalyDU(u)
        dTheta1DR = omega1 / vr
        dTheta2DR = (omega2 - l / (r^2)) / vr
        dTheta1DU = dRDU * dTheta1DR
        dTheta2DU = dRDU * dTheta2DR
        d2Theta1DRpDR, d2Theta1DRaDR = dOmega1DRp / vr - omega1 / vr^3 * (dEDRp - l / r^2 * dLDRp + (l^2 / r^3 - DPotDR(r)) * dRDRp), dOmega1DRa / vr - omega1 / vr^3 * (dEDRa - l / r^2 * dLDRa + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta2DRpDR, d2Theta2DRaDR = (dOmega2DRp - dLDRp / r^2 + 2.0 * l / r^3 * dRDRp) / vr - (omega2 - l / r^2) / vr^3 * (dEDRp - dLDRp * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRp), (dOmega2DRa - dLDRa / r^2 + 2.0 * l / r^3 * dRDRa) / vr - (omega2 - l / r^2) / vr^3 * (dEDRa - dLDRa * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta1DRpDU, d2Theta1DRaDU = d2Theta1DRpDR * dRDU + dTheta1DR * d2RDRpDU, d2Theta1DRaDR * dRDU + dTheta1DR * d2RDRaDU 
        d2Theta2DRpDU, d2Theta2DRaDU = d2Theta2DRpDR * dRDU + dTheta2DR * d2RDRpDU, d2Theta2DRaDR * dRDU + dTheta2DR * d2RDRaDU 
        uBasisValue = UBasis(ell, np, r)
        dUBasisDRValue = DUBasisDR(ell, np, r)
        #####
        valSin, valCos = sincos(n1 * (s2 + k2_3) + n2 * (s3 + k3_3))
        prefactorP, prefactorA = n1 * (s4P + k4P_3) + n2 * (s5P + k5P_3), n1 * (s4A + k4A_3) + n2 * (s5A + k5A_3)
        #####
        k1P_4 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * uBasisValue + dTheta1DU * dRDRp * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorP * valSin) 
        k1A_4 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * uBasisValue + dTheta1DU * dRDRa * dUBasisDRValue) * valCos - dTheta1DU * uBasisValue * prefactorA * valSin) 
        k2_4  = HWMat * dTheta1DU 
        k3_4  = HWMat * dTheta2DU 
        k4P_4 = HWMat * d2Theta1DRpDU 
        k4A_4 = HWMat * d2Theta1DRaDU 
        k5P_4 = HWMat * d2Theta2DRpDU 
        k5A_4 = HWMat * d2Theta2DRaDU 
        ########################################
        # Update of the state vectors
        s1P += (k1P_1 + 2.0 * k1P_2 + 2.0 * k1P_3 + k1P_4) / 6.0 
        s1A += (k1A_1 + 2.0 * k1A_2 + 2.0 * k1A_3 + k1A_4) / 6.0 
        s2  += (k2_1  + 2.0 * k2_2  + 2.0 * k2_3  + k2_4 ) / 6.0 
        s3  += (k3_1  + 2.0 * k3_2  + 2.0 * k3_3  + k3_4 ) / 6.0 
        s4P += (k4P_1 + 2.0 * k4P_2 + 2.0 * k4P_3 + k4P_4) / 6.0 
        s4A += (k4A_1 + 2.0 * k4A_2 + 2.0 * k4A_3 + k4A_4) / 6.0 
        s5P += (k5P_1 + 2.0 * k5P_2 + 2.0 * k5P_3 + k5P_4) / 6.0 
        s5A += (k5A_1 + 2.0 * k5A_2 + 2.0 * k5A_3 + k5A_4) / 6.0 
    end
    return s1P, s1A 
end


"""
    XMat(n1::Int64,n2::Int64,rp::Float64,ra::Float64)

Function XMat, the in-plane angular Fourier transform of the spherical radius. 

# Arguments
- `n1,n2`: angular Fourier numbers.
- `rp,ra`: orbital peri- and apocentre.

# Output
- `s1`: value of XMat, computed from a Runge-Kutta integration.
"""
function XMat(n1::Int64,n2::Int64,rp::Float64,ra::Float64)
    e, l = eLFromRpRa(rp,ra)
    omega1, omega2 = omega12FromRpRa(rp,ra)
    sigma, delta = (ra + rp) * 0.5, (ra - rp) * 0.5
    ########################################
    # Step 0 -- Initialisation
    s1, s2, s3 = 0.0, 0.0, 0.0 # State vectors: s1=W; s2=theta1; s3=theta2-psi
    u = - 1.0 + EpsilonWMat # Time for the integral. Avoiding the edges, where vr is ill-defined.
    r = sigma + delta * anomaly(u)
    dTheta1DU = delta * dAnomalyDU(u) * omega1 / sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
    dTheta2DU = (omega2 - l / (r^2)) * dTheta1DU / omega1
    ########################################
    for i=1:NStepsWMat
        ########################################
        # Step 1
        # Velocity vectors for s1, s2, s3.
        k1_1 = HWMat * (1.0 / pi) * dTheta1DU * r * cos(n1 * s2 + n2 * s3)
        k2_1 = HWMat * dTheta1DU
        k3_1 = HWMat * dTheta2DU
        ########################################
        # Step 2
        u += 0.5 * HWMat
        r = sigma + delta * anomaly(u) 
        dTheta1DU = delta * dAnomalyDU(u) * omega1 / sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dTheta2DU = (omega2 - l / (r^2)) * dTheta1DU / omega1        ########################################
        k1_2 = HWMat * (1.0 / pi) * dTheta1DU * r * cos(n1 * (s2 + 0.5 * k2_1) + n2 * (s3 + 0.5 * k3_1))
        k2_2 = HWMat * dTheta1DU
        k3_2 = HWMat * dTheta2DU 
        ########################################
        # Step 3
        k1_3 = HWMat * (1.0 / pi) * dTheta1DU * r * cos(n1 * (s2 + 0.5 * k2_2) + n2 * (s3 + 0.5 * k3_2)) 
        k2_3 = k2_2 
        k3_3 = k3_2 
        ########################################
        # Step 4
        u += 0.5 * HWMat 
        r = sigma + delta * anomaly(u) 
        dTheta1DU = delta * dAnomalyDU(u) * omega1 / sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dTheta2DU = (omega2 - l / (r^2)) * dTheta1DU / omega1
        ########################################
        k1_4 = HWMat * (1.0 / pi) * dTheta1DU * r * cos(n1 * (s2 + k2_3) + n2 * (s3 + k3_3)) 
        k2_4 = HWMat * dTheta1DU 
        k3_4 = HWMat * dTheta2DU 
        ########################################
        # Update
        s1 += (k1_1 + 2.0 * k1_2 + 2.0 * k1_3 + k1_4) / 6.0
        s2 += (k2_1 + 2.0 * k2_2 + 2.0 * k2_3 + k2_4) / 6.0 
        s3 += (k3_1 + 2.0 * k3_2 + 2.0 * k3_3 + k3_4) / 6.0 
    end
    return s1
end


"""
    dXMatDRpRa(n1::Int64,n2::Int64,rp::Float64,ra::Float64)

First derivatives of XMat w.r.t. rp and ra.

# Arguments
- `n1,n2`: angular Fourier numbers.
- `rp,ra`: orbital peri- and apocentre.

# Output
- `s1P,s1A`: derivatives of XMat w.r.t. rp and ra, dX/drp,dX/dra.
"""
function dXMatDRpRa(n1::Int64,n2::Int64,rp::Float64,ra::Float64)
    e, l = eLFromRpRa(rp,ra)
    omega1, omega2 = omega12FromRpRa(rp,ra)
    dEDRp, dEDRa = dEDRpRa(rp,ra)
    dLDRp, dLDRa = dLDRpRa(rp,ra)
    dOmega1DRp, dOmega1DRa, dOmega2DRp, dOmega2DRa = dOmega12DRpRa(rp,ra)
    sigma, delta = (ra + rp) * 0.5, (ra - rp) * 0.5
    ########################################
    # Step 0 -- Initialisation
    s1P, s1A, s2, s3, s4P, s4A, s5P, s5A = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 # State vectors: s1P=dW/drp; s1A=dW/dra; s2=theta1; s3=theta2-psi; s4P=dtheta1/drp; s4A=dtheta1/dra; s5P=d(theta2-psi)/drp; ds5A=d(theta-psi)/dra
    u = - 1.0 + EpsilonWMat # Time for the integral. Avoiding the edges, where vr is ill-defined.
    r = sigma + delta * anomaly(u)
    vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
    dRDRp, dRDRa, dRDU = 0.5 * (1.0 - anomaly(u)), 0.5 * (1.0 + anomaly(u)), delta * dAnomalyDU(u)
    d2RDRpDU, d2RDRaDU = - 0.5 * dAnomalyDU(u), 0.5 * dAnomalyDU(u)
    dTheta1DR = omega1 / vr
    dTheta2DR = (omega2 - l / (r^2)) / vr
    dTheta1DU = dRDU * dTheta1DR
    dTheta2DU = dRDU * dTheta2DR
    d2Theta1DRpDR, d2Theta1DRaDR = dOmega1DRp / vr - omega1 / vr^3 * (dEDRp - l / r^2 * dLDRp + (l^2 / r^3 - DPotDR(r)) * dRDRp), dOmega1DRa / vr - omega1 / vr^3 * (dEDRa - l / r^2 * dLDRa + (l^2 / r^3 - DPotDR(r)) * dRDRa)
    d2Theta2DRpDR, d2Theta2DRaDR = (dOmega2DRp - dLDRp / r^2 + 2.0 * l / r^3 * dRDRp) / vr - (omega2 - l / r^2) / vr^3 * (dEDRp - dLDRp * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRp), (dOmega2DRa - dLDRa / r^2 + 2.0 * l / r^3 * dRDRa) / vr - (omega2 - l / r^2) / vr^3 * (dEDRa - dLDRa * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRa)
    d2Theta1DRpDU, d2Theta1DRaDU = d2Theta1DRpDR * dRDU + dTheta1DR * d2RDRpDU, d2Theta1DRaDR * dRDU + dTheta1DR * d2RDRaDU 
    d2Theta2DRpDU, d2Theta2DRaDU = d2Theta2DRpDR * dRDU + dTheta2DR * d2RDRpDU, d2Theta2DRaDR * dRDU + dTheta2DR * d2RDRaDU 
    ########################################
    for istep=1:NStepsWMat 
        ########################################
        # Step 1
        valSin, valCos = sincos(n1 * s2 + n2 * s3)
        prefactorP, prefactorA = n1 * s4P + n2 * s5P, n1 * s4A + n2 * s5A 
        #####
        # Velocity vectors for s1P, s1A, s2, s3, s4P, s4A, s5P, s5A
        k1P_1 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * r + dTheta1DU * dRDRp) * valCos - dTheta1DU * r * prefactorP * valSin)
        k1A_1 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * r + dTheta1DU * dRDRa) * valCos - dTheta1DU * r * prefactorA * valSin)
        k2_1  = HWMat * dTheta1DU 
        k3_1  = HWMat * dTheta2DU 
        k4P_1 = HWMat * d2Theta1DRpDU 
        k4A_1 = HWMat * d2Theta1DRaDU
        k5P_1 = HWMat * d2Theta2DRpDU 
        k5A_1 = HWMat * d2Theta2DRaDU 
        ########################################
        # Step 2
        u += 0.5 * HWMat 
        r = sigma + delta * anomaly(u)
        vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dRDRp, dRDRa, dRDU = 0.5 * (1.0 - anomaly(u)), 0.5 * (1.0 + anomaly(u)), delta * dAnomalyDU(u)
        d2RDRpDU, d2RDRaDU = - 0.5 * dAnomalyDU(u), 0.5 * dAnomalyDU(u)
        dTheta1DR = omega1 / vr
        dTheta2DR = (omega2 - l / (r^2)) / vr
        dTheta1DU = dRDU * dTheta1DR
        dTheta2DU = dRDU * dTheta2DR
        d2Theta1DRpDR, d2Theta1DRaDR = dOmega1DRp / vr - omega1 / vr^3 * (dEDRp - l / r^2 * dLDRp + (l^2 / r^3 - DPotDR(r)) * dRDRp), dOmega1DRa / vr - omega1 / vr^3 * (dEDRa - l / r^2 * dLDRa + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta2DRpDR, d2Theta2DRaDR = (dOmega2DRp - dLDRp / r^2 + 2.0 * l / r^3 * dRDRp) / vr - (omega2 - l / r^2) / vr^3 * (dEDRp - dLDRp * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRp), (dOmega2DRa - dLDRa / r^2 + 2.0 * l / r^3 * dRDRa) / vr - (omega2 - l / r^2) / vr^3 * (dEDRa - dLDRa * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta1DRpDU, d2Theta1DRaDU = d2Theta1DRpDR * dRDU + dTheta1DR * d2RDRpDU, d2Theta1DRaDR * dRDU + dTheta1DR * d2RDRaDU 
        d2Theta2DRpDU, d2Theta2DRaDU = d2Theta2DRpDR * dRDU + dTheta2DR * d2RDRpDU, d2Theta2DRaDR * dRDU + dTheta2DR * d2RDRaDU 
        #####
        valSin, valCos = sincos(n1 * (s2 + 0.5 * k2_1) + n2 * (s3 + 0.5 * k3_1)) 
        prefactorP, prefactorA   = n1 * (s4P + 0.5 * k4P_1) + n2 * (s5P + 0.5 * k5P_1), n1 * (s4A + 0.5 * k4A_1) + n2 * (s5A + 0.5 * k5A_1) 
        #####
        k1P_2 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * r + dTheta1DU * dRDRp) * valCos - dTheta1DU * r * prefactorP * valSin)
        k1A_2 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * r + dTheta1DU * dRDRa) * valCos - dTheta1DU * r * prefactorA * valSin)
        k2_2  = HWMat * dTheta1DU 
        k3_2  = HWMat * dTheta2DU 
        k4P_2 = HWMat * d2Theta1DRpDU 
        k4A_2 = HWMat * d2Theta1DRaDU 
        k5P_2 = HWMat * d2Theta2DRpDU 
        k5A_2 = HWMat * d2Theta2DRaDU 
        ########################################
        # Step 3
        valSin, valCos = sincos(n1 * (s2 + 0.5 * k2_2) + n2 * (s3 + 0.5 * k3_2))
        prefactorP, prefactorA = n1 * (s4P + 0.5 * k4P_2) + n2 * (s5P + 0.5 * k5P_2), n1 * (s4A + 0.5 * k4A_2) + n2 * (s5A + 0.5 * k5A_2)
        #####
        k1P_3 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * r + dTheta1DU * dRDRp) * valCos - dTheta1DU * r * prefactorP * valSin)
        k1A_3 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * r + dTheta1DU * dRDRa) * valCos - dTheta1DU * r * prefactorA * valSin)
        k2_3  = HWMat * dTheta1DU 
        k3_3  = HWMat * dTheta2DU 
        k4P_3 = HWMat * d2Theta1DRpDU 
        k4A_3 = HWMat * d2Theta1DRaDU 
        k5P_3 = HWMat * d2Theta2DRpDU 
        k5A_3 = HWMat * d2Theta2DRaDU 
        ########################################
        # Step 4
        u += 0.5 * HWMat 
        r = sigma + delta * anomaly(u)
        vr = sqrt(abs(2.0 * (e - Pot(r)) - (l / r)^2))
        dRDRp, dRDRa, dRDU = 0.5 * (1.0 - anomaly(u)), 0.5 * (1.0 + anomaly(u)), delta * dAnomalyDU(u)
        d2RDRpDU, d2RDRaDU = - 0.5 * dAnomalyDU(u), 0.5 * dAnomalyDU(u)
        dTheta1DR = omega1 / vr
        dTheta2DR = (omega2 - l / (r^2)) / vr
        dTheta1DU = dRDU * dTheta1DR
        dTheta2DU = dRDU * dTheta2DR
        d2Theta1DRpDR, d2Theta1DRaDR = dOmega1DRp / vr - omega1 / vr^3 * (dEDRp - l / r^2 * dLDRp + (l^2 / r^3 - DPotDR(r)) * dRDRp), dOmega1DRa / vr - omega1 / vr^3 * (dEDRa - l / r^2 * dLDRa + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta2DRpDR, d2Theta2DRaDR = (dOmega2DRp - dLDRp / r^2 + 2.0 * l / r^3 * dRDRp) / vr - (omega2 - l / r^2) / vr^3 * (dEDRp - dLDRp * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRp), (dOmega2DRa - dLDRa / r^2 + 2.0 * l / r^3 * dRDRa) / vr - (omega2 - l / r^2) / vr^3 * (dEDRa - dLDRa * l / r^2 + (l^2 / r^3 - DPotDR(r)) * dRDRa)
        d2Theta1DRpDU, d2Theta1DRaDU = d2Theta1DRpDR * dRDU + dTheta1DR * d2RDRpDU, d2Theta1DRaDR * dRDU + dTheta1DR * d2RDRaDU 
        d2Theta2DRpDU, d2Theta2DRaDU = d2Theta2DRpDR * dRDU + dTheta2DR * d2RDRpDU, d2Theta2DRaDR * dRDU + dTheta2DR * d2RDRaDU 
        #####
        valSin, valCos = sincos(n1 * (s2 + k2_3) + n2 * (s3 + k3_3))
        prefactorP, prefactorA = n1 * (s4P + k4P_3) + n2 * (s5P + k5P_3), n1 * (s4A + k4A_3) + n2 * (s5A + k5A_3)
        #####
        k1P_4 = HWMat * (1.0 / pi) * ((d2Theta1DRpDU * r + dTheta1DU * dRDRp) * valCos - dTheta1DU * r * prefactorP * valSin) 
        k1A_4 = HWMat * (1.0 / pi) * ((d2Theta1DRaDU * r + dTheta1DU * dRDRa) * valCos - dTheta1DU * r * prefactorA * valSin) 
        k2_4  = HWMat * dTheta1DU 
        k3_4  = HWMat * dTheta2DU 
        k4P_4 = HWMat * d2Theta1DRpDU 
        k4A_4 = HWMat * d2Theta1DRaDU 
        k5P_4 = HWMat * d2Theta2DRpDU 
        k5A_4 = HWMat * d2Theta2DRaDU 
        ########################################
        # Update of the state vectors
        s1P += (k1P_1 + 2.0 * k1P_2 + 2.0 * k1P_3 + k1P_4) / 6.0 
        s1A += (k1A_1 + 2.0 * k1A_2 + 2.0 * k1A_3 + k1A_4) / 6.0 
        s2  += (k2_1  + 2.0 * k2_2  + 2.0 * k2_3  + k2_4 ) / 6.0 
        s3  += (k3_1  + 2.0 * k3_2  + 2.0 * k3_3  + k3_4 ) / 6.0 
        s4P += (k4P_1 + 2.0 * k4P_2 + 2.0 * k4P_3 + k4P_4) / 6.0 
        s4A += (k4A_1 + 2.0 * k4A_2 + 2.0 * k4A_3 + k4A_4) / 6.0 
        s5P += (k5P_1 + 2.0 * k5P_2 + 2.0 * k5P_3 + k5P_4) / 6.0 
        s5A += (k5A_1 + 2.0 * k5A_2 + 2.0 * k5A_3 + k5A_4) / 6.0 
    end
    return s1P, s1A 
end


"""
    pFactor(ell::Int64, np::Int64)

Prefactor of the barycentric drift: radial integral of the density basis elements.

# Arguments
- `ell`: harmonic number ell.
- `np`: order of the radial basis element.

# Output
- `integralPFactor`: value of p = int_0^{8 RMax} D_ell^np(r) dr.
"""
function pFactor(ell::Int64, np::Int64)
    integralPFactor, errorPFactor = quadgk(r -> DBasis(ell, np,r), 0.0 + Sanitizer, 8.0 * RMax, rtol=Sanitizer)
    return integralPFactor
end
const tabPFactor = zeros(Float64,NBasisElements)


"""
    tabPFactor!(ell::Int64)

Function filling the table tabPFactor with the values of the prefactors.

# Arguments
- `ell`: harmonic number ell.

# Output
None
"""
function tabPFactor!(ell::Int64)
    for np=1:NBasisElements
        tabPFactor[np] = pFactor(ell, np - 1 + NStart)
    end
end 


#
# Filling the table of prefactors. Only the value ell=1 is required, because the barycentric drift only acts on this harmonic.
#
if EllMax>0
    println("Acceleration p factor, timing:")
    flush(stdout)
    @time tabPFactor!(1)
    flush(stdout)
else
    println("No acceleration p factor computed.")
end


