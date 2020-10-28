# reactingDNS_OpenFOAM
Incorporating detailed transport properties with OpenFOAM reactingFoam solver for direct numerical simulation of multi-component reacting flows

In the OpenFOAM's native solver, reactingFoam, the original transport equations for species mass fractions were based on the simple asscumptions of unity Scimidt number. And the dynamic viscosity was calculated using the Sutherland equation. This is acceptable for RANS or LES simulation of turbulent combustion, in which the molecular diffusion being negligible in comparison with the turbulent diffusion. But this is not acceptable for the reaction-diffusion process, e.g. laminar flame propagation. 

Current solver use to the third-order logarithm polunomial fitting method to descirbe the detailed transport properties (viscosity, thermal conducitivity and binary diffusion coefficient) for each species [1]. And the mixture-average transport model is used to calculate the mixture properties.

Usage:

Step 1:

Using ANSYS CHEMKIN to fit the transport data. The tran.out can be outputted when you selected to process Transport properties with Fit with Verbose Output. 

Step 2:

Using the formatConvert python code to convert tran.out to four files: binaryDiff, speciesLambda, speciesMu and thermoDiff.

Only H and H2 are considered to have the thermoDiff.

Step 3:

Put the binaryDiff, speciesLambda, speciesMu, thermoDiff and trandat into the constant directory.

Step 4:

Trun on/off the switch of differetialDiffusion and thermalDiffusion. If differetialDiffusion is off, the unity Lewis number will be adopted for all species.
