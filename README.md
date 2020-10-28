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

Attention:

When the species name starts with number in trandat, you need to change it to word by add double quotation marks, e.g. 1H2 to "1H2".

Current method is simple and can be easily implemented in other solvers in OpenFOAM. When you are using this solver to publish paper, please kindly consider to cite following papers, in which n-heptane/air laminar premixed flames in a low-temperature ignition regime [2], spherical laminar premixed flame with intrinsic flame instability [3] and turbulent flame propagation and ignition [4, 5] were stuided.

[1] Kee, R. J.; Dixon-Lewis, G.; Warnatz, J.; Coltrin, M. E.; Miller, J. A. A Fortran com-puter code package for the evaluation of gas-phase multicomponent transport proper-ties.Sandia National Laboratories Report SAND86-82461986,13, 80401–1887.
[2] Zhong, S.; Zhang, F.; Jangi, M.; Bai, X.-S.; Yao, M.; Peng, Z. Structure and propagationof n-heptane/air premixed flame in low temperature ignition regime.Applied  Energy2020,275, 115320.
[3] Zhang, N.; Zhang, F.; Zhong, S.; Peng, Z.; Yu, J.; Liu, H.; Xu, C. Numerical and theo-retical investigation of ethanol/air flame instability.Combustion Theory and Modelling2020, 1–22.
[4] Zhong, S.;  Zhang, F.;  Peng, Z.;  Bai, F.;  Du, Q. Roles of CO2 and H2O in premixedturbulent oxy-fuel combustion.Fuel2018,234, 1044–1054.
[5] Zhong, S.; Zhang, F.; Du, Q.; Peng, Z. Characteristics of reactivity controlled combus-tion with n-heptane low temperature reforming products.Fuel2020,275, 117980.
