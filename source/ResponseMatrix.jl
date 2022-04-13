##################################################
# COMPUTATION OF THE RESPONSE MATRIX.
##################################################

#
# Tables containing the times.
#
const GridTimes = [t for t=TMin:DeltaT:TMax]
const GridTimesSize = length(GridTimes)


#
# Computing Ylm(pi/2,0)
# 
"""
    CMatrix!()

Function filling the CMatrix table with the Ylm(pi/2,0) values.

# Arguments
None

# Output
None
"""
const CMatrix = zeros(Float64,EllMax + 1,EllMax + 1) # Container for the C-coefficients !! HARMONIC 
function CMatrix!()
    for m=1:(EllMax + 1), l=m:(EllMax + 1)
        CMatrix[l,m] = sphericalY(l - 1, m - 1, pi / 2.0, 0.0)
    end
end
CMatrix!() # Initialising CMatrix with the pre-computed coefficients
"""
    CYlm(ell::Int64, m::Int64)

Wrapped function returning the spherical harmonic prefactors.

# Arguments
- `ell`: harmonic number ell.
- `m`: harmonic number m.

# Output
- value of Ylm(pi/2,0). 
"""
function CYlm(ell::Int64, m::Int64)
    if iseven(ell)
        return CMatrix[ell + 1, abs(m) + 1]
    else # Yl,m(pi/2) = - Yl,-m(pi/2,0) for odd m
        return CMatrix[ell + 1, abs(m) + 1] * sign(m)
    end
end


#
# Pre-computing chunks of values of WMat and its derivatives.
# 
"""
    tabWMat!(ell::Int64, n1::Int64,n2::Int64,rp::Float64,ra::Float64,tabWMat::Array{Float64,1})

Function filling the list tabWMat with all values of WMat for given (n1,n2,rp,ra), i.e. over all the radial basis indices.

# Arguments
- `ell`: harmonic number ell.
- `n1,n2`: angular Fourier numbers.
- `rp,ra`: orbital peri- and apocentre.
- `tabWMat`: table to fill with the values of WMat.

# Output
None
"""
function tabWMat!(ell::Int64, n1::Int64,n2::Int64,rp::Float64,ra::Float64,tabWMat::Array{Float64,1})
    for np=1:NBasisElements
        tabWMat[np] = WMat(ell, np - 1 + NStart,n1,n2,rp,ra)
    end
end
"""
    tabDWMatDRpRa!(ell::Int64, n1::Int64,n2::Int64,rp::Float64,ra::Float64,tabDWMatDRpRa::Array{Float64,2})

Function filling the list tabDWMatDRpRa with all values of the derivatives of WMat w.r.t. rp and ra for given (n1,n2,rp,ra), i.e. over all the radial basis indices.

# Arguments
- `ell`: harmonic number ell.
- `n1,n2`: angular Fourier numbers.
- `rp,ra`: orbital peri- and apocentre.
- `tabDWMatDRpRa`: table to fill with the values of the derivatives, dW/drp,dW/dra. Filled so that tabDWMatDRpRa[1,np] is dW/drp for the given np, and tabDWMatDRpRa[2,np] is dW/dra for the given np.

# Output
None
"""
function tabDWMatDRpRa!(ell::Int64, n1::Int64,n2::Int64,rp::Float64,ra::Float64,tabDWMatDRpRa::Array{Float64,2})
    for np=1:NBasisElements
        tabDWMatDRpRa[1,np], tabDWMatDRpRa[2,np] = dWMatDRpRa(ell, np - 1 + NStart,n1,n2,rp,ra)
    end
end
const tabWMat = [zeros(Float64,NBasisElements) for i=1:NThrds]
const tabDWMatDRpRa = [zeros(Float64,2,NBasisElements) for i=1:NThrds]


#
# Computing the response matrix.
# 
"""
    tabMMat!(ell::Int64, tabMMat::Array{Array{Complex{Float64},2},1})

Function filling the list of matrices tabMMat with all values of the response matrix, for a given value of ell. The outer dimension of tabMMat corresponds to the different time steps, the inner 2 dimensions to the basis radial orders, np and nq.
This function is parallelised on multiple threads, over the scanning of the action space through the index i.
This function considers all possible resonances by default.

# Arguments
- `ell`: harmonic number ell.
- `tabMMat`: list of matrices to fill.

# Output
None
"""
function tabMMat!(ell::Int64, tabMMat::Array{Array{Complex{Float64},2},1})
    tabMMatLarge = [[zeros(Complex{Float64},NBasisElements,NBasisElements) for i=1:GridTimesSize] for j=1:NThrds]
    println("Grid size: ", NUVGrid)
    flush(stdout)
    Threads.@threads for i=1:NUVGrid
        ThrID = Threads.threadid()
        if OnCluster == false
            print("Integration domain: ", i, " / ", NUVGrid, '\r')
        elseif (floor(10 * (i + 1) / NUVGrid) != floor(10 * i / NUVGrid))
            println("Integration domain: ", floor(10 * (i + 1) / NUVGrid) * 10, "%")
            flush(stdout)
        end
        u, v = gridUV[:,i]
        indexInGridU = floor(Int64, (u - UMin) / DeltaUV + 1.0)
        rp, ra = gridRpRa[:,i]
        e, l = eLFromRpRa(rp,ra)
        dRpDUloc = dRpDU(u)
        dRaDV = dRpDU(v)
        dEDRp, dEDRa = dEDRpRa(rp,ra)
        dLDRp, dLDRa = dLDRpRa(rp,ra)
        dRpDE = 1.0 / (dEDRp - dEDRa * dLDRp / dLDRa)
        dRpDL = 1.0 / (dLDRp - dLDRa * dEDRp / dEDRa)
        dRaDE = 1.0 / (dEDRa - dEDRp * dLDRa / dLDRp)
        dRaDL = 1.0 / (dLDRa - dLDRp * dEDRa / dEDRp)
        d2EDRp2, d2EDRpDRaloc, d2EDRa2 = d2EDRpDRa(rp,ra)
        d2LDRp2, d2LDRpDRaloc, d2LDRa2 = d2LDRpDRa(rp,ra)
        dDRpDEDRa = - dRpDE^2 * (d2EDRpDRaloc - d2EDRa2 * dLDRp / dLDRa - dEDRa * d2LDRpDRaloc / dLDRa + dEDRa * dLDRp * d2LDRa2 / dLDRa^2)
        dDRpDLDRa = - dRpDL^2 * (d2LDRpDRaloc - d2LDRa2 * dEDRp / dEDRa - dLDRa * d2EDRpDRaloc / dEDRa + dLDRa * dEDRp * d2EDRa2 / dEDRa^2)
        dDRaDEDRp = - dRaDE^2 * (d2EDRpDRaloc - d2EDRp2 * dLDRa / dLDRp - dEDRp * d2LDRpDRaloc / dLDRp + dEDRp * dLDRa * d2LDRp2 / dLDRp^2)
        dDRaDLDRp = - dRaDL^2 * (d2LDRpDRaloc - d2LDRp2 * dEDRa / dEDRp - dLDRp * d2EDRpDRaloc / dEDRp + dLDRp * dEDRa * d2EDRp2 / dEDRp^2)
        dDRaDEDRa = - dRaDE^2 * (d2EDRa2 - d2EDRpDRaloc * dLDRa / dLDRp - dEDRp * d2LDRa2 / dLDRp + dEDRp * dLDRa * d2LDRpDRaloc / dLDRp^2)
        dDRaDLDRa = - dRaDL^2 * (d2LDRa2 - d2LDRpDRaloc * dEDRa / dEDRp - dLDRp * d2EDRa2 / dEDRp + dLDRp * dEDRa * d2EDRpDRaloc / dEDRp^2)
        dDRpDEDRp = - dRpDE^2 * (d2EDRp2 - d2EDRpDRaloc * dLDRp / dLDRa - dEDRa * d2LDRp2 / dLDRa + dEDRa * dLDRp * d2LDRpDRaloc / dLDRa^2)
        dDRpDLDRp = - dRpDL^2 * (d2LDRp2 - d2LDRpDRaloc * dEDRp / dEDRa - dLDRa * d2EDRp2 / dEDRa + dLDRa * dEDRp * d2EDRpDRaloc / dEDRa^2)
        omega1, omega2 = omega12FromRpRa(rp,ra)
        dOmega1DRp, dOmega1DRa, dOmega2DRp, dOmega2DRa = dOmega12DRpRa(rp,ra)
        jacELToRpRa = jacobianELToRpRa(rp,ra)
        dJacELDRp, dJacELDRa = dJacobianELToRpRaDRpRa(rp,ra)
        jacRpRaToUV = jacobianRpRaToUV(u,v)
        dJacRpRaDU, dJacRpRaDV = dJacobianRpRaToUVDUV(u,v)
        FEL = DF(e,l)
        dDFDE, dDFDL = DDFDEL(e,l)
        d2DFDE2, d2DFDEDL, d2DFDL2 = D2DFDEDL(e,l)
        ##########
        # We benefit from the symmetry W[-n1,-n2]=W[n1,n2] to reduce the number of evaluations of WMat
        # We also note that the resonance (n1,n2)=(0,0) does not contribute to the response matrix, and we make sure not to consider it
        for n2=ell:-2:0 # Loop over the resonant index n2 !! We use the fact that it goes by steps of two
            pref = 2.0 * (2.0 * pi)^(3) * CYlm(ell, n2)^(2) / (2.0 * ell + 1.0) # Value of the prefactor
            n1Bound = - N1Max # Value of the bound for the n1 summation if n2>0
            ##########
            if (n2==0)
                n1Bound = 1
            end
            ##########
            for n1=n1Bound:N1Max
                nDotOmega = n1 * omega1 + n2 * omega2
                dNDotOmegaDRp = n1 * dOmega1DRp + n2 * dOmega2DRp
                dNDotOmegaDRa = n1 * dOmega1DRa + n2 * dOmega2DRa
                nDotDDFDJ = n1 * omega1 * dDFDE + n2 * (dDFDL + omega2 * dDFDE)
                dnDotDDFDJDRp = n1 * (dOmega1DRp * dDFDE + omega1 * (d2DFDE2 * dEDRp + d2DFDEDL * dLDRp)) + n2 * (d2DFDL2 * dLDRp + d2DFDEDL * dEDRp + dOmega2DRp * dDFDE + omega2 * (d2DFDE2 * dEDRp + d2DFDEDL * dLDRp)) 
                dnDotDDFDJDRa = n1 * (dOmega1DRa * dDFDE + omega1 * (d2DFDE2 * dEDRa + d2DFDEDL * dLDRa)) + n2 * (d2DFDL2 * dLDRa + d2DFDEDL * dEDRa + dOmega2DRa * dDFDE + omega2 * (d2DFDE2 * dEDRa + d2DFDEDL * dLDRa)) 
                nDotDRpDJ = n1 * omega1 * dRpDE + n2 * (dRpDL + omega2 * dRpDE)
                dnDotDRpDJDRa = n1 * (dOmega1DRa * dRpDE + omega1 * dDRpDEDRa + n2 * (dDRpDLDRa + dOmega2DRa * dRpDE + omega2 * dDRpDEDRa)) 
                dnDotDRpDJDRp = n1 * (dOmega1DRp * dRpDE + omega1 * dDRpDEDRp + n2 * (dDRpDLDRp + dOmega2DRp * dRpDE + omega2 * dDRpDEDRp))
                nDotDRaDJ = n1 * omega1 * dRaDE + n2 * (dRaDL + omega2 * dRaDE)
                dnDotDRaDJDRp = n1 * (dOmega1DRp * dRaDE + omega1 * dDRaDEDRp + n2 * (dDRaDLDRp + dOmega2DRp * dRaDE + omega2 * dDRaDEDRp)) 
                dnDotDRaDJDRa = n1 * (dOmega1DRa * dRaDE + omega1 * dDRaDEDRa + n2 * (dDRaDLDRa + dOmega2DRa * dRaDE + omega2 * dDRaDEDRa))
                g0 = pref * jacRpRaToUV * jacELToRpRa * (l / omega1) * nDotDDFDJ # Value of the numerator without Wnp and Wnq !! 
                dG0DU = pref * (dJacRpRaDU * jacELToRpRa * (l / omega1) * nDotDDFDJ + jacRpRaToUV * dRpDUloc * ((dJacELDRp + dJacELDRa) * l / omega1 * nDotDDFDJ + jacELToRpRa * (dLDRp + dLDRa) / omega1 * nDotDDFDJ - jacELToRpRa * l * (dOmega1DRp + dOmega1DRa) / omega1^2 * nDotDDFDJ + jacELToRpRa * l / omega1 * (dnDotDDFDJDRp + dnDotDDFDJDRa)))
                dG0DV = pref * (dJacRpRaDV * jacELToRpRa * (l / omega1) * nDotDDFDJ + jacRpRaToUV * dRaDV * (dJacELDRa * l / omega1 * nDotDDFDJ + jacELToRpRa * dLDRa / omega1 * nDotDDFDJ - jacELToRpRa * l * dOmega1DRa / omega1^2 * nDotDDFDJ + jacELToRpRa * l / omega1 * dnDotDDFDJDRa))
                g0EdgeU = pref * abs(dRaDV) * jacELToRpRa * (l / omega1) * FEL * nDotDRpDJ
                dG0EdgeUDV = pref * (sign(dRaDV) * d2RpDU2(v) * jacELToRpRa * (l / omega1) * FEL * nDotDRpDJ + abs(dRaDV) * dRaDV * (dJacELDRa * l / omega1 * FEL * nDotDRpDJ + jacELToRpRa * dLDRa / omega1 * FEL * nDotDRpDJ - jacELToRpRa * l * dOmega1DRa / omega1^2 * FEL * nDotDRpDJ + jacELToRpRa * l / omega1 * (dDFDE * dEDRa + dDFDL * dLDRa) * nDotDRpDJ + jacELToRpRa * l / omega1 * FEL * dnDotDRpDJDRa))
                g0EdgeV = pref * abs(dRpDUloc) * jacELToRpRa * (l / omega1) * FEL * (nDotDRaDJ - nDotDRpDJ)
                dG0EdgeVDU = pref * (sign(dRpDUloc) * d2RpDU2(u) * jacELToRpRa * (l / omega1) * FEL * (nDotDRaDJ - nDotDRpDJ) + abs(dRpDUloc) * dRpDUloc * ((dJacELDRp + dJacELDRa) * l / omega1 * FEL * (nDotDRaDJ - nDotDRpDJ) + jacELToRpRa * (dLDRp + dLDRa) / omega1 * FEL * (nDotDRaDJ - nDotDRpDJ) - jacELToRpRa * l * (dOmega1DRp + dOmega1DRa) / omega1^2 * FEL * (nDotDRaDJ - nDotDRpDJ) + jacELToRpRa * l / omega1 * (dDFDE * (dEDRp + dEDRa) + dDFDL * (dLDRp + dLDRa)) * (nDotDRaDJ - nDotDRpDJ) + jacELToRpRa * l / omega1 * FEL * ((dnDotDRaDJDRp + dnDotDRaDJDRa) - (dnDotDRpDJDRp + dnDotDRpDJDRa))))
                h0 = - nDotOmega
                dH0DU = - dRpDUloc * (dNDotOmegaDRp + dNDotOmegaDRa)
                dH0DV = - dRaDV * dNDotOmegaDRa
                ##########
                tabWMat!(ell, n1,n2,rp,ra,tabWMat[ThrID])
                tabDWMatDRpRa!(ell, n1,n2,rp,ra,tabDWMatDRpRa[ThrID])
                XMatloc = XMat(n1,n2,rp,ra)
                dXMatDRp, dXMatDRa = dXMatDRpRa(n1,n2,rp,ra)
                dXMatDU, dXMatDV = (dXMatDRp + dXMatDRa) * dRpDUloc, dXMatDRa * dRaDV
                ##########
                # Filling in the response matrix
                for np=1:NBasisElements
                    if fillMatrixSingleSquare!(np, ell, u, v, indexInGridU, tabWMat[ThrID],  tabDWMatDRpRa[ThrID], dRpDUloc, dRaDV, XMatloc, dXMatDU, dXMatDV, g0, dG0DU, dG0DV, h0, dH0DU, dH0DV, g0EdgeU, dG0EdgeUDV, g0EdgeV, dG0EdgeVDU, tabMMatLarge[ThrID]) == "error"
                        println("Values of n1, n2, rp, ra: ", n1, " ", n2, " ", rp ," ",ra)
                    end
                end
            end
        end
    end
    tabMMatloc = [zeros(Complex{Float64},NBasisElements,NBasisElements) for j=1:GridTimesSize]
    for i=1:NThrds
        tabMMatloc += tabMMatLarge[i]
    end
    for t=1:GridTimesSize
        tabMMat[t] = tabMMatloc[t]
    end
end
"""
    tabMMatResonances!(ell::Int64, resonanceList, tabMMat::Array{Array{Complex{Float64},2},1})

Function filling the list of matrices tabMMat with all values of the response matrix, for a given value of ell and a given set of resonances resonanceList. The outer dimension of tabMMat corresponds to the different time steps, the inner 2 dimensions to the basis radial orders, np and nq.
This function is parallelised on multiple threads, over the scanning of the action space through the index i.
WARNING: when a resonance (n1,n2) is included, its opposite (-n1,-n2) is included by design.

# Arguments
- `ell`: harmonic number ell.
- `resonanceList`: list of resonance vectors on which the matrix is computed.
- `tabMMat`: list of matrices to fill.

# Output
None
"""
function tabMMatResonances!(ell::Int64, resonanceList, tabMMat::Array{Array{Complex{Float64},2},1})
    tabMMatLarge = [[zeros(Complex{Float64},NBasisElements,NBasisElements) for i=1:GridTimesSize] for j=1:NThrds]
    Threads.@threads for i=1:NUVGrid
        ThrID = Threads.threadid()
        if OnCluster == false
            print("Integration domain: ", i, " / ", NUVGrid, '\r')
        end
        u, v = gridUV[:,i]
        indexInGridU = floor(Int64, (u - UMin) / DeltaUV + 1.0)
        rp, ra = gridRpRa[:,i]
        e, l = eLFromRpRa(rp,ra)
        dRpDUloc = dRpDU(u)
        dRaDV = dRpDU(v)
        dEDRp, dEDRa = dEDRpRa(rp,ra)
        dLDRp, dLDRa = dLDRpRa(rp,ra)
        dRpDE = 1.0 / (dEDRp - dEDRa * dLDRp / dLDRa)
        dRpDL = 1.0 / (dLDRp - dLDRa * dEDRp / dEDRa)
        dRaDE = 1.0 / (dEDRa - dEDRp * dLDRa / dLDRp)
        dRaDL = 1.0 / (dLDRa - dLDRp * dEDRa / dEDRp)
        d2EDRp2, d2EDRpDRaloc, d2EDRa2 = d2EDRpDRa(rp,ra)
        d2LDRp2, d2LDRpDRaloc, d2LDRa2 = d2LDRpDRa(rp,ra)
        dDRpDEDRa = - dRpDE^2 * (d2EDRpDRaloc - d2EDRa2 * dLDRp / dLDRa - dEDRa * d2LDRpDRaloc / dLDRa + dEDRa * dLDRp * d2LDRa2 / dLDRa^2)
        dDRpDLDRa = - dRpDL^2 * (d2LDRpDRaloc - d2LDRa2 * dEDRp / dEDRa - dLDRa * d2EDRpDRaloc / dEDRa + dLDRa * dEDRp * d2EDRa2 / dEDRa^2)
        dDRaDEDRp = - dRaDE^2 * (d2EDRpDRaloc - d2EDRp2 * dLDRa / dLDRp - dEDRp * d2LDRpDRaloc / dLDRp + dEDRp * dLDRa * d2LDRp2 / dLDRp^2)
        dDRaDLDRp = - dRaDL^2 * (d2LDRpDRaloc - d2LDRp2 * dEDRa / dEDRp - dLDRp * d2EDRpDRaloc / dEDRp + dLDRp * dEDRa * d2EDRp2 / dEDRp^2)
        dDRaDEDRa = - dRaDE^2 * (d2EDRa2 - d2EDRpDRaloc * dLDRa / dLDRp - dEDRp * d2LDRa2 / dLDRp + dEDRp * dLDRa * d2LDRpDRaloc / dLDRp^2)
        dDRaDLDRa = - dRaDL^2 * (d2LDRa2 - d2LDRpDRaloc * dEDRa / dEDRp - dLDRp * d2EDRa2 / dEDRp + dLDRp * dEDRa * d2EDRpDRaloc / dEDRp^2)
        dDRpDEDRp = - dRpDE^2 * (d2EDRp2 - d2EDRpDRaloc * dLDRp / dLDRa - dEDRa * d2LDRp2 / dLDRa + dEDRa * dLDRp * d2LDRpDRaloc / dLDRa^2)
        dDRpDLDRp = - dRpDL^2 * (d2LDRp2 - d2LDRpDRaloc * dEDRp / dEDRa - dLDRa * d2EDRp2 / dEDRa + dLDRa * dEDRp * d2EDRpDRaloc / dEDRa^2)
        omega1, omega2 = omega12FromRpRa(rp,ra)
        dOmega1DRp, dOmega1DRa, dOmega2DRp, dOmega2DRa = dOmega12DRpRa(rp,ra)
        jacELToRpRa = jacobianELToRpRa(rp,ra)
        dJacELDRp, dJacELDRa = dJacobianELToRpRaDRpRa(rp,ra)
        jacRpRaToUV = jacobianRpRaToUV(u,v)
        dJacRpRaDU, dJacRpRaDV = dJacobianRpRaToUVDUV(u,v)
        FEL = DF(e,l)
        dDFDE, dDFDL = DDFDEL(e,l)
        d2DFDE2, d2DFDEDL, d2DFDL2 = D2DFDEDL(e,l)
        ##########
        # We benefit from the symmetry W[-n1,-n2]=W[n1,n2] to reduce the number of evaluations of WMat
        # We also note that the resonance (n1,n2)=(0,0) does not contribute to the response matrix, and we make sure not to consider it
        for resonance in resonanceList
            n1,n2 = resonance
            pref = 2.0 * (2.0 * pi)^(3) * CYlm(ell, n2)^(2) / (2.0 * ell + 1.0) # Value of the prefactor
            nDotOmega = n1 * omega1 + n2 * omega2
            dNDotOmegaDRp = n1 * dOmega1DRp + n2 * dOmega2DRp
            dNDotOmegaDRa = n1 * dOmega1DRa + n2 * dOmega2DRa
            nDotDDFDJ = n1 * omega1 * dDFDE + n2 * (dDFDL + omega2 * dDFDE)
            dnDotDDFDJDRp = n1 * (dOmega1DRp * dDFDE + omega1 * (d2DFDE2 * dEDRp + d2DFDEDL * dLDRp)) + n2 * (d2DFDL2 * dLDRp + d2DFDEDL * dEDRp + dOmega2DRp * dDFDE + omega2 * (d2DFDE2 * dEDRp + d2DFDEDL * dLDRp)) 
            dnDotDDFDJDRa = n1 * (dOmega1DRa * dDFDE + omega1 * (d2DFDE2 * dEDRa + d2DFDEDL * dLDRa)) + n2 * (d2DFDL2 * dLDRa + d2DFDEDL * dEDRa + dOmega2DRa * dDFDE + omega2 * (d2DFDE2 * dEDRa + d2DFDEDL * dLDRa)) 
            nDotDRpDJ = n1 * omega1 * dRpDE + n2 * (dRpDL + omega2 * dRpDE)
            dnDotDRpDJDRa = n1 * (dOmega1DRa * dRpDE + omega1 * dDRpDEDRa + n2 * (dDRpDLDRa + dOmega2DRa * dRpDE + omega2 * dDRpDEDRa)) 
            dnDotDRpDJDRp = n1 * (dOmega1DRp * dRpDE + omega1 * dDRpDEDRp + n2 * (dDRpDLDRp + dOmega2DRp * dRpDE + omega2 * dDRpDEDRp))
            nDotDRaDJ = n1 * omega1 * dRaDE + n2 * (dRaDL + omega2 * dRaDE)
            dnDotDRaDJDRp = n1 * (dOmega1DRp * dRaDE + omega1 * dDRaDEDRp + n2 * (dDRaDLDRp + dOmega2DRp * dRaDE + omega2 * dDRaDEDRp)) 
            dnDotDRaDJDRa = n1 * (dOmega1DRa * dRaDE + omega1 * dDRaDEDRa + n2 * (dDRaDLDRa + dOmega2DRa * dRaDE + omega2 * dDRaDEDRa))
            g0 = pref * jacRpRaToUV * jacELToRpRa * (l / omega1) * nDotDDFDJ # Value of the numerator without Wnp and Wnq !! 
            dG0DU = pref * (dJacRpRaDU * jacELToRpRa * (l / omega1) * nDotDDFDJ + jacRpRaToUV * dRpDUloc * ((dJacELDRp + dJacELDRa) * l / omega1 * nDotDDFDJ + jacELToRpRa * (dLDRp + dLDRa) / omega1 * nDotDDFDJ - jacELToRpRa * l * (dOmega1DRp + dOmega1DRa) / omega1^2 * nDotDDFDJ + jacELToRpRa * l / omega1 * (dnDotDDFDJDRp + dnDotDDFDJDRa)))
            dG0DV = pref * (dJacRpRaDV * jacELToRpRa * (l / omega1) * nDotDDFDJ + jacRpRaToUV * dRaDV * (dJacELDRa * l / omega1 * nDotDDFDJ + jacELToRpRa * dLDRa / omega1 * nDotDDFDJ - jacELToRpRa * l * dOmega1DRa / omega1^2 * nDotDDFDJ + jacELToRpRa * l / omega1 * dnDotDDFDJDRa))
            g0EdgeU = pref * abs(dRaDV) * jacELToRpRa * (l / omega1) * FEL * nDotDRpDJ
            dG0EdgeUDV = pref * (sign(dRaDV) * d2RpDU2(v) * jacELToRpRa * (l / omega1) * FEL * nDotDRpDJ + abs(dRaDV) * dRaDV * (dJacELDRa * l / omega1 * FEL * nDotDRpDJ + jacELToRpRa * dLDRa / omega1 * FEL * nDotDRpDJ - jacELToRpRa * l * dOmega1DRa / omega1^2 * FEL * nDotDRpDJ + jacELToRpRa * l / omega1 * (dDFDE * dEDRa + dDFDL * dLDRa) * nDotDRpDJ + jacELToRpRa * l / omega1 * FEL * dnDotDRpDJDRa))
            g0EdgeV = pref * abs(dRpDUloc) * jacELToRpRa * (l / omega1) * FEL * (nDotDRaDJ - nDotDRpDJ)
            dG0EdgeVDU = pref * (sign(dRpDUloc) * d2RpDU2(u) * jacELToRpRa * (l / omega1) * FEL * (nDotDRaDJ - nDotDRpDJ) + abs(dRpDUloc) * dRpDUloc * ((dJacELDRp + dJacELDRa) * l / omega1 * FEL * (nDotDRaDJ - nDotDRpDJ) + jacELToRpRa * (dLDRp + dLDRa) / omega1 * FEL * (nDotDRaDJ - nDotDRpDJ) - jacELToRpRa * l * (dOmega1DRp + dOmega1DRa) / omega1^2 * FEL * (nDotDRaDJ - nDotDRpDJ) + jacELToRpRa * l / omega1 * (dDFDE * (dEDRp + dEDRa) + dDFDL * (dLDRp + dLDRa)) * (nDotDRaDJ - nDotDRpDJ) + jacELToRpRa * l / omega1 * FEL * ((dnDotDRaDJDRp + dnDotDRaDJDRa) - (dnDotDRpDJDRp + dnDotDRpDJDRa))))
            h0 = - nDotOmega
            dH0DU = - dRpDUloc * (dNDotOmegaDRp + dNDotOmegaDRa)
            dH0DV = - dRaDV * dNDotOmegaDRa
            ##########
            tabWMat!(ell, n1,n2,rp,ra,tabWMat[ThrID])
            tabDWMatDRpRa!(ell, n1,n2,rp,ra,tabDWMatDRpRa[ThrID])
            XMatloc = XMat(n1,n2,rp,ra)
            dXMatDRp, dXMatDRa = dXMatDRpRa(n1,n2,rp,ra)
            dXMatDU, dXMatDV = (dXMatDRp + dXMatDRa) * dRpDUloc, dXMatDRa * dRaDV
            ##########
            # Filling in the response matrix
            for np=1:NBasisElements
                fillMatrixSingleSquare!(np, ell, u, v, indexInGridU, tabWMat[ThrID],  tabDWMatDRpRa[ThrID], dRpDUloc, dRaDV, XMatloc, dXMatDU, dXMatDV, g0, dG0DU, dG0DV, h0, dH0DU, dH0DV, g0EdgeU, dG0EdgeUDV, g0EdgeV, dG0EdgeVDU, tabMMatLarge[ThrID])
            end
        end
    end
    tabMMatloc = [zeros(Complex{Float64},NBasisElements,NBasisElements) for j=1:GridTimesSize]
    for i=1:NThrds
        tabMMatloc += tabMMatLarge[i]
    end
    for t=1:GridTimesSize
        tabMMat[t] = tabMMatloc[t]
    end
end


"""
    fillMatrixSingleSquare!(np::Int64, ell::Int64, u::Float64, v::Float64, indexInGridU::Int64, tabWMat::Array{Float64,1},  tabDWMatDRpRa::Array{Float64,2}, dRpDUloc::Float64, dRaDV::Float64, XMatloc::Float64, dXMatDU::Float64, dXMatDV::Float64, g0::Float64, dG0DU::Float64, dG0DV::Float64, h0::Float64, dH0DU::Float64, dH0DV::Float64, g0EdgeU::Float64, dG0EdgeUDV::Float64, g0EdgeV::Float64, dG0EdgeVDU::Float64, matrix)

Function incrementing the line np of the response matrix `matrix`, taken as an argument, with the right square integrals corresponding to a given value of the index i in the u,v plane, and the resonance vector (n1,n2).
We take advantage of the np-nq symmetries in the matrix to restrict the for loop to nq>=np.

# Arguments
- `np`: order of the radial basis element, i.e. index of a line in the response matrix.
- `ell`: harmonic number ell.
- `u,v`: current position in the (u,v) plane. 
- `indexInGridU`: position of u in the full grid of u values.
- `tabWMat`: list of values of WMat for the current (n1,n2,rp,ra) and all possible radial orders.
- `tabDWMatDRpRa`: list of values of (dW/drp,dW/dra) for the current (n1,n2,rp,ra) and all possible radial orders.
- `dRpDUloc`: current value of drp/du.
- `dRaDV`: current value of dra/dv.
- `XMatloc`: value of XMat for the current (n1,n2,rp,ra).
- `dXMatDU,dXMatDV`: derivatives of XMat w.r.t. u and v for the current (n1,n2,rp,ra).
- `g0`: current value of g0, i.e. prefactor * jacobians((Jr,L) -> (u,v)) * n.dDF/dJ.
- `dG0DU,dG0DV`: derivatives of g0 w.r.t. u and v.
- `h0`: current value of h0, i.e. -n.Omega.
- `dH0DU,dH0DV`: derivatives of h0 w.r.t. u and v.
- `g0EdgeU`: value of the version of g0 at the u-edge of the integration domain.
- `dG0EdgeUDV`: derivative of g0EdgeU w.r.t. v.
- `g0EdgeV`: value of the version of g0 at the v-edge of the integration domain.
- `dG0EdgeVDU`: derivative of g0EdgeV w.r.t. u.
- `matrix`: matrix to be filled, i.e. the full response matrix.

# Output
None
"""
function fillMatrixSingleSquare!(np::Int64, ell::Int64, u::Float64, v::Float64, indexInGridU::Int64, tabWMat::Array{Float64,1},  tabDWMatDRpRa::Array{Float64,2}, dRpDUloc::Float64, dRaDV::Float64, XMatloc::Float64, dXMatDU::Float64, dXMatDV::Float64, g0::Float64, dG0DU::Float64, dG0DV::Float64, h0::Float64, dH0DU::Float64, dH0DV::Float64, g0EdgeU::Float64, dG0EdgeUDV::Float64, g0EdgeV::Float64, dG0EdgeVDU::Float64, matrix)
    wnp = tabWMat[np]
    dWnpDRp, dWnpDRa = tabDWMatDRpRa[:,np]
    dWnpDU, dWnpDV = (dWnpDRp + dWnpDRa) * dRpDUloc, dWnpDRa * dRaDV
    wnpDrift = wnp + 4.0 * pi / 3.0 * tabPFactor[np] * XMatloc
    dWnpDUDrift = dWnpDU + 4.0 * pi / 3.0 * tabPFactor[np] * dXMatDU
    dWnpDVDrift = dWnpDV + 4.0 * pi / 3.0 * tabPFactor[np] * dXMatDV
    for nq=np:NBasisElements
        wnq = tabWMat[nq]
        dWnqDRp, dWnqDRa = tabDWMatDRpRa[:,nq]
        dWnqDU, dWnqDV = (dWnqDRp + dWnqDRa) * dRpDUloc, dWnqDRa * dRaDV
        wnqDrift = wnq + 4.0 * pi / 3.0 * tabPFactor[nq] * XMatloc
        dWnqDUDrift = dWnqDU + 4.0 * pi / 3.0 * tabPFactor[nq] * dXMatDU
        dWnqDVDrift = dWnqDV + 4.0 * pi / 3.0 * tabPFactor[nq] * dXMatDV
        # Loop over the considered times 
        for t=2:GridTimesSize
            if ell!=1
                if fillMatrixElements!(g0 * wnp * wnq, dG0DU * wnp * wnq + g0 * dWnpDU * wnq + g0 * wnp * dWnqDU, dG0DV * wnp * wnq + g0 * dWnpDV * wnq + g0 * wnp * dWnqDV, h0 * GridTimes[t], dH0DU * GridTimes[t], dH0DV * GridTimes[t], np, nq, matrix[t]) == "error"
                    println("Values of u, v, g0, h0: ", u, " ", v, " ", g0 ," ",h0)
                    flush(stdout)
                    return "error"
                end
                if (u == DeltaUV * floor(UMax / DeltaUV))
                    fillMatrixEdgeElements!(g0EdgeU * wnp * wnq, dG0EdgeUDV * wnp * wnq + g0EdgeU * dWnpDV * wnq + g0EdgeU * wnp * dWnqDV, h0 * GridTimes[t], dH0DV * GridTimes[t], np, nq, matrix[t])
                end
                if (v == GridVMax[indexInGridU])
                    fillMatrixEdgeElements!(g0EdgeV * wnp * wnq, dG0EdgeVDU * wnp * wnq + g0EdgeV * dWnpDU * wnq + g0EdgeV * wnp * dWnqDU, h0 * GridTimes[t], dH0DU * GridTimes[t], np, nq, matrix[t])
                end
            else
                fillMatrixElementsDrift!(g0 * wnp * wnqDrift, dG0DU * wnp * wnqDrift + g0 * dWnpDU * wnqDrift + g0 * wnp * dWnqDUDrift, dG0DV * wnp * wnqDrift + g0 * dWnpDV * wnqDrift + g0 * wnp * dWnqDVDrift, g0 * wnpDrift * wnq, dG0DU * wnpDrift * wnq + g0 * dWnpDUDrift * wnq + g0 * wnpDrift * dWnqDU, dG0DV * wnpDrift * wnq + g0 * dWnpDVDrift * wnq + g0 * wnpDrift * dWnqDV, h0 * GridTimes[t], dH0DU * GridTimes[t], dH0DV * GridTimes[t], np, nq, matrix[t])
                if (u == DeltaUV * floor(UMax / DeltaUV))
                    fillMatrixEdgeElementsDrift!(g0EdgeU * wnp * wnqDrift, dG0EdgeUDV * wnp * wnqDrift + g0EdgeU * dWnpDV * wnqDrift + g0EdgeU * wnp * dWnqDVDrift, g0EdgeU * wnpDrift * wnq, dG0EdgeUDV * wnpDrift * wnq + g0EdgeU * dWnpDVDrift * wnq + g0EdgeU * wnpDrift * dWnqDV, h0 * GridTimes[t], dH0DV * GridTimes[t], np, nq, matrix[t])
                end
                if (v == GridVMax[indexInGridU])
                    fillMatrixEdgeElementsDrift!(g0EdgeV * wnp * wnqDrift, dG0EdgeVDU * wnp * wnqDrift + g0EdgeV * dWnpDU * wnqDrift + g0EdgeV * wnp * dWnqDUDrift, g0EdgeV * wnpDrift * wnq, dG0EdgeVDU * wnpDrift * wnq + g0EdgeV * dWnpDUDrift * wnq + g0EdgeV * wnpDrift * dWnqDU, h0 * GridTimes[t], dH0DU * GridTimes[t], np, nq, matrix[t])
                end
            end
        end
    end
end


"""
    fillMatrixElements!(a::Float64, b::Float64, c::Float64, d::Float64, e::Float64, f::Float64, np::Int64, nq::Int64, matrix)

Function incrementing the element np,nq of the response matrix `matrix`, taken as an argument, with the right square integral corresponding to a given value of the index i in the u,v plane, of the resonance vector (n1,n2), and of the time step t.
It also increments the symmetric element nq,np.

# Arguments
- `a`: value of the factor before the exponential in the integrand, i.e. g0 * W_np * W_nq
- `b`: derivative of a w.r.t. u.
- `c`: derivative of a w.r.t. v.
- `d`: phase of the exponential factor in the integrand, i.e. h0 * t.
- `e`: derivative of d w.r.t. u.
- `f`: derivative of d w.r.t. v.
- `np,nq`: orders of the radial basis element, indexes of a cell in the response matrix.
- `matrix`: matrix to be filled, i.e. the response matrix at a single time step.

# Output
None
"""
function fillMatrixElements!(a::Float64, b::Float64, c::Float64, d::Float64, e::Float64, f::Float64, np::Int64, nq::Int64, matrix)
    resloc = 2.0 * real(- im * aleph(a,b,c,d,e,f))
    if (np == nq)
        matrix[np,nq] += resloc
    else 
        matrix[np,nq] += resloc
        matrix[nq,np] += resloc
    end
    if isnan(resloc)
        println("resloc is NaN. Values of a, b, c, d, e, f, np, nq: ", a, " ", b, " ", c ," ",d," ",e," ",f," ",np," ",nq)
        flush(stdout)
        return "error"
    else
        return nothing
    end
end
"""
    fillMatrixEdgeElements!(a::Float64, b::Float64, c::Float64, d::Float64, np::Int64, nq::Int64, matrix)

Function incrementing the element np,nq of the response matrix `matrix`, taken as an argument, with the right square integral corresponding to a given value of the index i in the u,v plane located on its edge, of the resonance vector (n1,n2), and of the time step t.
It also increments the symmetric element nq,np.

# Arguments
- `a`: value of the factor before the exponential in the integrand.
- `b`: derivative of a w.r.t. either u or v, depending on the edge.
- `c`: phase of the exponential factor in the integrand, i.e. h0 * t.
- `d`: derivative of c w.r.t either u or v, depending on the edge.
- `np,nq`: orders of the radial basis element, indexes of a cell in the response matrix.
- `matrix`: matrix to be filled, i.e. the response matrix at a single time step.

# Output
None
"""
function fillMatrixEdgeElements!(a::Float64, b::Float64, c::Float64, d::Float64, np::Int64, nq::Int64, matrix)
    resloc = 2.0 * real(- im * alephEdge(a,b,c,d))
    if (np == nq)
        matrix[np,nq] += resloc 
    else 
        matrix[np,nq] += resloc 
        matrix[nq,np] += resloc 
    end
end
"""
    fillMatrixElementsDrift!(a1::Float64, b1::Float64, c1::Float64, a2::Float64, b2::Float64, c2::Float64, d::Float64, e::Float64, f::Float64, np::Int64, nq::Int64, matrix)

Function incrementing the element np,nq of the response matrix `matrix`, taken as an argument, with the right square integral corresponding to a given value of the index i in the u,v plane, of the resonance vector (n1,n2), and of the time step t, in the case where ell=1 and the barycentric drift should be taken into account.
It also increments the symmetric element nq,np. Since the new term introduces an asymmetry in the matrix, separate expressions must be computed for the two symmetric terms.

# Arguments
- `a1`: value of the factor before the exponential in the integrand when the drift concerns the nq dimension, i.e. g0 * W_np * (W_nq + X_nq)
- `b1`: derivative of a1 w.r.t. u.
- `c1`: derivative of a1 w.r.t. v.
- `a2`: value of the factor before the exponential in the integrand when the drift concerns the np dimension, i.e. g0 * (W_np + X_np) * W_nq
- `b2`: derivative of a2 w.r.t. u.
- `c2`: derivative of a2 w.r.t. v.
- `d`: phase of the exponential factor in the integrand, i.e. h0 * t.
- `e`: derivative of d w.r.t. u.
- `f`: derivative of d w.r.t. v.
- `np,nq`: orders of the radial basis element, indexes of a cell in the response matrix.
- `matrix`: matrix to be filled, i.e. the response matrix at a single time step.

# Output
None
"""
function fillMatrixElementsDrift!(a1::Float64, b1::Float64, c1::Float64, a2::Float64, b2::Float64, c2::Float64, d::Float64, e::Float64, f::Float64, np::Int64, nq::Int64, matrix)
    respqloc = 2.0 * real(- im * aleph(a1, b1, c1, d, e, f))
    resqploc = 2.0 * real(- im * aleph(a2, b2, c2, d, e, f))
    if (np == nq)
        matrix[np,nq] += respqloc 
    else 
        matrix[np,nq] += respqloc 
        matrix[nq,np] += resqploc 
    end
end
"""
    fillMatrixEdgeElementsDrift!(a1::Float64, b1::Float64, a2::Float64, b2::Float64, c::Float64, d::Float64, np::Int64, nq::Int64, matrix)

Function incrementing the element np,nq of the response matrix `matrix`, taken as an argument, with the right square integral corresponding to a given value of the index i in the u,v plane located on its edge, of the resonance vector (n1,n2), and of the time step t, in the case where ell=1 and the barycentric drift should be taken into account.
It also increments the symmetric element nq,np. Since the new term introduces an asymmetry in the matrix, separate expressions must be computed for the two symmetric terms.

# Arguments
- `a1`: value of the factor before the exponential in the integrand when the drift concerns the nq dimension.
- `b1`: derivative of a1 w.r.t. either u or v, depending on the edge.
- `a2`: value of the factor before the exponential in the integrand when the drift concerns the np dimension.
- `b2`: derivative of a2 w.r.t. either u or v, depending on the edge.
- `c`: phase of the exponential factor in the integrand, i.e. h0 * t.
- `d`: derivative of c w.r.t either u or v, depending on the edge.
- `np,nq`: orders of the radial basis element, indexes of a cell in the response matrix.
- `matrix`: matrix to be filled, i.e. the response matrix at a single time step.

# Output
None
"""
function fillMatrixEdgeElementsDrift!(a1::Float64, b1::Float64, a2::Float64, b2::Float64, c::Float64, d::Float64, np::Int64, nq::Int64, matrix)
    respqloc = 2.0 * real(- im * alephEdge(a1, b1, c, d))
    resqploc = 2.0 * real(- im * alephEdge(a2, b2, c, d))
    if (np == nq)
        matrix[np,nq] += respqloc 
    else 
        matrix[np,nq] += respqloc 
        matrix[nq,np] += resqploc 
    end
end


"""
    aleph(g::Float64,dGDU::Float64,dGDV::Float64,h::Float64,dHDU::Float64,dHDV::Float64)

Function computing one local orbital contribution.
Corresponds to INT[(g+dgdu*u+dgdv*v)*exp(im*(h+dhdu*u+dhdv*v)),{u,-DeltaUV/2,DeltaUV/2},{v,-DeltaUV/2,DeltaUV/2}]

# Arguments
- `g`: value of the factor before the exponential in the integrand.
- `dGDU, dGDV`: derivative of g w.r.t. u and v.
- `h`: phase of the exponential factor in the integrand, i.e. h0 * t.
- `dHDU,dHDV`: derivative of h w.r.t u and v.

# Output
- value of the integral.
"""
function aleph(g::Float64,dGDU::Float64,dGDV::Float64,h::Float64,dHDU::Float64,dHDV::Float64)
    (g * exp(im * h)) * DeltaUV^2 * alephD(dGDU * DeltaUV / g, dGDV * DeltaUV / g, dHDU * DeltaUV, dHDV * DeltaUV)
end
"""
    alephD(a::Float64,b::Float64,c::Float64,d::Float64)

Analytical function computing one local orbital contribution, with normalised arguments.
Corresponds to INT[(1+ax+by)*exp(im*(cx+dy)),{x,-1/2,1/2},{y,-1/2,1/2}]

# Arguments
- `a`: coefficient before x in the integrand.
- `b`: coefficient before y in the integrand.
- `c`: coefficient before x in the exponential phase of the integrand.
- `d`: coefficient before y in the exponential phase of the integrand.

# Output
- value of the analytical integral.
"""
function alephD(a::Float64,b::Float64,c::Float64,d::Float64)
    scc = sinc(0.5 * c / pi)
    scd = sinc(0.5 * d / pi)
    return scc * scd - im * a / c * scd * (cos(0.5 * c) - scc) - im * b / d * scc * (cos(0.5 * d) - scd)
end
"""
    alephEdge(g::Float64,dGDV::Float64,h::Float64,dHDV::Float64)

Function computing one local orbital contribution on the edge of the integration domain.
Corresponds to INT[(g+dgdv*v)*exp(im*(h+dhdv*v)) ,{v,-DeltaUV/2,DeltaUV/2}]

# Arguments
- `g`: value of the factor before the exponential in the integrand.
- `dGDV`: derivative of g w.r.t. v.
- `h`: phase of the exponential factor in the integrand, i.e. h0 * t.
- `dHDV`: derivative of h w.r.t v.

# Output
- value of the integral.
"""
function alephEdge(g::Float64,dGDV::Float64,h::Float64,dHDV::Float64)
    (g * exp(im * h)) * DeltaUV * alephDEdge(dGDV * DeltaUV / g, dHDV * DeltaUV)
end
"""
    alephDEdge(a::Float64,b::Float64)

Analytical function computing one local orbital contribution, with normalised arguments, on the edge of the integration domain.
Corresponds to INT[(1+ax)*exp(im*bx),{x,-1/2,1/2}]

# Arguments
- `a`: coefficient before x in the integrand.
- `b`: coefficient before x in the exponential phase of the integrand.

# Output
- value of the analytical integral.
"""
function alephDEdge(a::Float64,b::Float64)
    scb = sinc(0.5 * b / pi)
    return scb - im * a / b * (cos(0.5 * b) - scb)
end


