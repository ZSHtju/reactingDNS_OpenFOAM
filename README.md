# reactingDNS_OpenFOAM
Incorporating detailed transport properties with OpenFOAM reactingFoam solver for direct numerical simulation of multi-component reacting flows

In the OpenFOAM's native solver, reactingFoam, the original transport equations for species mass fractions were based on the simple asscumptions of unity Scimidt number. And the dynamic viscosity was calculated using the Sutherland equation. This is acceptable for RANS or LES simulation of turbulent combustion, in which the molecular diffusion being negligible in comparison with the turbulent diffusion. But this is not acceptable for the reaction-diffusion process, e.g. flame propagation. 
Usage:
