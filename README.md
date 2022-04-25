# LiRGHaM
Linear Response of a Galactic Halo with the Matrix method


## Configuration of the code and first tests

1 - After downloading the project files, the first step is to open `source/MW_LMC_Response/LMC_Main.jl` and replace the string `CodeDirectory` on line 8 with the local path to the folder `source/`.

2 - Launch `julia` and import all required packages (all located in `Packages.jl`) by running 

```
julia> import Pkg
julia> Pkg.add("PackageName")
```

3 - First test: from the terminal, enter the folder `source/MW_LMC_Response/` and launch the `LMC_Main.jl` file, e.g. by running:

```
> julia LMC_Main.jl
```

This will include all files properly, i.e. load packages, generate the action space grid, etc.

Output:

```
Loading packages
  3.097204 seconds (6.87 M allocations: 493.059 MiB, 2.80% gc time, 0.10% compilation time)
Acceleration p factor, timing:
  0.973640 seconds (1.87 M allocations: 115.461 MiB, 2.36% gc time, 99.98% compilation time)
No projection computed.
No matrix computed.
No response computed.
```


4 - Open `source/MW_LMC_Response/MW_LMC_Parameters.jl`. Change the parameters 

```
const projectionQ = "computed"
const resonanceMatrixQ = "computed"
const responseQ = "computed"
```

Then, run the `LMC_Main.jl` file again.

```
> julia LMC_Main.jl
```

This will compute the LMC's projection, the response matrix and the response with a set of small parameters ensuring fast computation. When the code is done running, check that the perturber, matrix and response files were dropped in `source/MW_LMC_Response/Data/`.

Output: 

```
Loading packages
  2.971229 seconds (6.87 M allocations: 493.187 MiB, 3.12% gc time, 0.19% compilation time)
Acceleration p factor, timing:
  0.962040 seconds (1.87 M allocations: 115.461 MiB, 2.17% gc time, 99.98% compilation time)
Computing the trajectory.
  0.055668 seconds (306.67 k allocations: 18.516 MiB, 102.77% compilation time)
Projection of the LMC
Time step = 1; Projection
  0.309196 seconds (4.63 M allocations: 137.284 MiB, 4.97% gc time, 47.68% compilation time)
Time step = 2; Projection
  0.146434 seconds (4.23 M allocations: 116.092 MiB, 3.17% gc time)
```
 
[...]

```
Time step = 20; Projection
  0.119947 seconds (3.40 M allocations: 93.419 MiB, 3.31% gc time)
Time step = 21; Projection
  0.113061 seconds (3.14 M allocations: 86.164 MiB, 3.44% gc time)
  3.102557 seconds (80.79 M allocations: 2.198 GiB, 4.09% gc time, 9.00% compilation time)
Computing the matrix: 
ell = 0; resonance = (1, 0); Response matrix computation
  3.022264 seconds (6.67 M allocations: 365.766 MiB, 3.22% gc time, 99.40% compilation time)
ell = 0; resonance = (2, 0); Response matrix computation
  0.017884 seconds (18.36 k allocations: 500.938 KiB)
```

[...] 

```
ell = 1; resonance = (2, 1); Response matrix computation
  0.020757 seconds (18.36 k allocations: 500.938 KiB)
  3.176705 seconds (6.85 M allocations: 372.725 MiB, 3.06% gc time, 95.40% compilation time)
```

5 - Explore different setups by changing the `MW_LMC_Parameters.jl` file: orbit of the LMC (`version`); anisotropy of the Milky Way halo (`Beta`); set of basis functions (`BasisType`); radial scale of the basis functions (`RBasis`); action space sampling (`RMin`, `RMax`, `DeltaUV`, `Q`); maximum angular harmonic (`EllMax`); maximum radial basis element (`NMax`); maximum angular Fourier number (`N1Max`). 


## Producing figures

The figures are produced using a Wolfram Mathematica script.

1 - Open `plots/ResponsePlots.m` and enter the `figLMCFullResponse` section. Modify the `nbdir` string with the local path to the `plots/` folder. 

2 - Modify the values of the parameters in the `Parameters` subsection. All corresponding data files must have been produced by the Julia code prior to plotting.

3 - From the terminal, run the MathKernel command, e.g. on MacOS through

```
> /Applications/Mathematica.app/Contents/MacOS/MathKernel -script ResponsePlots.m
```

The mathematica script generates plots in the `plots/` folder. 21 of them correspond to the MW halo's response at the 21 integration time steps, one is the colour bar for these plots, and one is a .gif file presenting a movie.


## Example

With the following parameter file `MW_LMC_Parameters.jl`:
```
########################################
# CONSTANT PARAMETERS FOR THE WHOLE COMPUTATION.
########################################
#
# Defining which actions the code will do.
#
const version = "MW_LMC_notAcc" 	# Switch to know which version this is. Possibilities: "MW_LMC_notAcc" "MW_LMC_Acc"
const projectionQ = "computed" 			# Whether the projection should be read or computed. "read", "computed" or false.
const resonanceMatrixQ = "computed" 		# Whether the resonance matrix should be read or computed. "read", "computed" or false.
const responseQ = "computed" 			# Whether the response should be computed. "computed" or false.


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
const NMax = 20						# Maximum order of the radial basis elements.
const Q = 1.0						# Sampling parameter of the (u,v) grid.
const N1Max = 2						# Maximum order of the radial Fourier number n_1.
const Sigma = 0.001					# Width parameter for the construction of the (u,v) grid.
const Ra = 1.0						# Anisotropy radius of the Osipkov-Merrit DFs.
const Beta = 0.01					# Anisotropy parameter of the Baes-Van Hese DFs.
const QPlummer = 0.0				# Anisotropy parameter of the Dejonghe DF of the Plummer sphere.
const EllMax = 4					# Maximum order of the ell harmonic number.
const NStepsWMat = 100				# Sampling size for the Runge-Kutta integration of WMat.
const EpsilonWMat = 1.0e-4			# Sanity edge width for the integration of WMat.
const HWMat = 2.0 * (1.0 - EpsilonWMat) / NStepsWMat # Integration step for the computation of WMat. Avoiding the edges, where vr is ill-defined.
const TMin = 0.0					# Minimum time step. Should remain 0.0 by convention.
const TMax = 19.93994 # Default for the LMC: 19.93994	# Maximum time for computing the evolution.
const DeltaT = 0.996997 # Default: 1/20 of total time = 0.996997	# Time step of the evolution.
```

The `FiducialResponse.gif` movie produced is:

![FiducialResponse](https://user-images.githubusercontent.com/103592382/165072369-4f917577-5476-4356-b3e3-a942beac73a0.gif)
