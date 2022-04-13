########################################
# CONSTANT PARAMETERS FOR THE WHOLE COMPUTATION.
########################################
#
# Defining which actions the code will do.
#
const version = "MW_LMC_notAcc" 	# Switch to know which version this is. Possibilities: "MW_LMC_notAcc" "MW_LMC_Acc"
const projectionQ = false 			# Whether the projection should be read or computed. "read", "computed" or false.
const resonanceMatrixQ = false 		# Whether the resonance matrix should be read or computed. "read", "computed" or false.
const responseQ = false 			# Whether the response should be computed. "computed" or false.


#
# Constants setting the dimentionless scales for each possible potential.
#
const G = 1.0						# Newton's constant of gravity.
const MTot = 1.0					# Mass scale.
const BIsochrone = 1.0				# Scale radius of the isochrone sphere.
const BPlummer = 1.0				# Scale radius of the Plummer sphere.
const BHernquist = 1.0				# Scale radius of the Hernquist sphere.
const ANFW = 1.0					# Scale radius of the NFW sphere.
const RMaxNFW = 10.0 * ANFW			# Truncation radius of the NFW sphere.


#
# Code parameters
#
const PotentialType = "Hernquist" 	# Type of background potential from which the orbits are taken. The list of possibilities can be found in "Potentials.jl"
const DFType = "HERNQUISTBaes"		# Type of distribution function sampling the phase space. No requirement that the DF should generate the potential. The DF should be positive. The list of possibilities can be found in "DistributionFunctions.jl"
const BasisType = "Clutton"			# Type of radial basis functions. The list of possibilities can be found in "RadialBasis.jl"
const RMin = 0.1 # 0.1				# Minimum pericentric radius considered.
const RMax = 10.0 # 10.0			# Maximum apocentric radius considered.
const DeltaUV = 1.0					# Step size in the (u,v) grid.
const RBasis = 11.0 # 11.0			# Scale radius of the radial basis.
const NMax = 3						# Maximum order of the radial basis elements.
const Q = 1.0						# Sampling parameter of the (u,v) grid.
const N1Max = 2						# Maximum order of the radial Fourier number n_1.
const Sigma = 0.001					# Width parameter for the construction of the (u,v) grid.
const Ra = 1.0						# Anisotropy radius of the Osipkov-Merrit DFs.
const Beta = 0.01					# Anisotropy parameter of the Baes-Van Hese DFs.
const QPlummer = 0.0				# Anisotropy parameter of the Dejonghe DF of the Plummer sphere.
const EllMax = 1					# Maximum order of the ell harmonic number.
const NStepsWMat = 100				# Sampling size for the Runge-Kutta integration of WMat.
const EpsilonWMat = 1.0e-4			# Sanity edge width for the integration of WMat.
const HWMat = 2.0 * (1.0 - EpsilonWMat) / NStepsWMat # Integration step for the computation of WMat. Avoiding the edges, where vr is ill-defined.
const TMin = 0.0					# Minimum time step. Should remain 0.0 by convention.
const TMax = 19.93994 # Default for the LMC: 19.93994	# Maximum time for computing the evolution.
const DeltaT = 0.996997 # Default: 1/20 of total time = 0.996997	# Time step of the evolution.
