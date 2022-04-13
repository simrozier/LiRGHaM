########################################
# PHYSICAL QUANTITIES REQUIRED FOR INTEGRATING THE LMC'S TRAJECTORY.
########################################
#
# Physical characteristics of the Milky Way halo, in physical units.
#
const GNewton = 4.30e-6 # (km/s)^2 kpc/Msun. Newton's constant.
const aMW = 40.8 # kpc. Scale radius of the Milky Way halo's Hernquist sphere.
const MMW = 1.57e12 # Msun. Milky Way halo mass.
const RVirMW = 280 # kpc. Virial radius of the Hernquist halo.
#
# Parameters used in the conversion from physical to theoretical units. (G = MMW = aMW = 1)
#
const VUnit = sqrt(GNewton * MMW / aMW) # Convertion from km/s to theoretical units.
const TUnit = aMW / VUnit # Same for the time.
const RhoUnit = MMW / aMW^3 # Same for the density.
#
# Parameters of the LMC.
#
const MLMC = 1.8e11 # Msun
const aLMC = 20.0 # kpc
#
# Parameters for the integration of the LMC's orbit.
# The orbit integration is performed with higher accuracy than what will effectively remain in the matrix computation. So the time step is decreased.
#
const RpLMC = 48.0 # kpc. Pericentric radius of the LMC's orbit around the Milky Way.
const VpLMC = 340.0 # km / s. Velocity of the LMC at pericenter.
const NTStepsIntegration = 621 # Number of time steps for the precise integration of the LMC's orbit.
const NTStepsIntegrationAcc = 1805 # Number of time steps for the precise integration of the LMC's orbit in the case where the Milky Way's motion is taken into account.
const DeltaTIntegration = DeltaT / 10.0 # Integration time step for the precise integration.
const TMaxIntegration = TMin + (NTStepsIntegration - 1) * DeltaTIntegration # Maximum time step of the orbit integration.
const TMaxIntegrationAcc = TMin + (NTStepsIntegrationAcc - 1) * DeltaTIntegration # Maximum time step of the orbit integration, accelerated Milky Way.


