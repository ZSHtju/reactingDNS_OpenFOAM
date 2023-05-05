# reactingDNS_OpenFOAM
Incorporating detailed transport properties with OpenFOAM reactingFoam solver for direct numerical simulation of multi-component reacting flows

In the OpenFOAM's native solver, reactingFoam, the original transport equations for species mass fractions were based on the simple asscumptions of unity Schmidt number. And the dynamic viscosity was calculated using the Sutherland equation. This is acceptable for RANS or LES simulation of turbulent combustion, in which the molecular diffusion being negligible in comparison with the turbulent diffusion. But this is not acceptable for the reaction-diffusion process, e.g. laminar flame propagation. 

Current solver reactingDNS use the third-order logarithm polynomial fitting method to descirbe the detailed transport properties (viscosity, thermal conducitivity and binary diffusion coefficient) of each species [1]. And the mixture-average transport model is used to calculate the mixture properties.

Usage:

Step 1:

Using ANSYS CHEMKIN premixed freely propagation laminar flame models to fit the transport data. The tran.out can be outputted when you selecte to process Transport properties with Fit with Verbose Output. 

Step 2:

Using the formatConvert python code to convert tran.out to four files: binaryDiff, speciesLambda, speciesMu and thermoDiff.

Only H and H2 are considered to have the thermoDiff.

Step 3:

Put the binaryDiff, speciesLambda, speciesMu, thermoDiff and trandat into the constant directory.

Step 4:

Trun on/off the switch of differentialDiffusion and thermalDiffusion (Soret diffusion). If differetialDiffusion is off, the unity Lewis number will be adopted for all species.

Attention:

When the species name starts with number in trandat, you need to change it to word by add double quotation marks, e.g. 1H2 to "1H2".

Please ensure there is no extra useless lines at the end of the tran.out when you use the writetofoam.

Please ensure you have a blank line at the beginning of the trandat when you put it into the case.

Current method is simple and can be easily implemented in other solvers in OpenFOAM. When you are using this solver to publish paper, please kindly consider to cite following papers, in which n-heptane/air laminar premixed flames in a low-temperature ignition regime [2], spherical laminar premixed flame with intrinsic flame instability [3], turbulent flame propagation and ignition [4, 5] and ammonia/methane/air premixed laminar burner flame [6] were stuided.

[1] Kee, R. J.; Dixon-Lewis, G.; Warnatz, J.; Coltrin, M. E.; Miller, J. A. A Fortran com-puter code package for the evaluation of gas-phase multicomponent transport properties.Sandia National Laboratories Report SAND86-8246, 1986, 13, 80401–1887.

[2] Zhong, S.; Zhang, F.; Jangi, M.; Bai, X.-S.; Yao, M.; Peng, Z. Structure and propagation of n-heptane/air premixed flame in low temperature ignition regime. Applied  Energy, 2020, 275, 115320.

[3] Zhang, N.; Zhang, F.; Zhong, S.; Peng, Z.; Yu, J.; Liu, H.; Xu, C. Numerical and theoretical investigation of ethanol/air flame instability. Combustion Theory and Modelling, 2020, 1–22.

[4] Zhong, S.;  Zhang, F.;  Peng, Z.;  Bai, F.;  Du, Q. Roles of CO2 and H2O in premixed turbulent oxy-fuel combustion. Fuel, 2018, 234, 1044–1054.

[5] Zhong, S.; Zhang, F.; Du, Q.; Peng, Z. Characteristics of reactivity controlled combustion with n-heptane low temperature reforming products. Fuel, 2020, 275, 117980.

[6] Rodolfo C. Rocha, Shenghui Zhong, Leilei Xu, Xue-Song Bai, Mário Costa, Xiao Cai, Haisol Kim, Christian Brackmann, Zhongshan Li, and Marcus Aldén.
Energy & Fuels 2021, 35 (9), 7179-7192. DOI: 10.1021/acs.energyfuels.0c03520.


Other published papers using reactingDNS solver:

[1] Tahsin Berk Kıymaz, Emre Böncü, Dilay Güleryüz, Mehmet Karaca, Barış Yılmaz, Christophe Allouis, İskender Gökalp. Numerical investigations on flashback dynamics of premixed methane-hydrogen-air laminar flames. International Journal of Hydrogen Energy, 2022, 47(59), 25022-25033. https://doi.org/10.1016/j.ijhydene.2022.05.230.

[2] Shervin Karimkashi, Ossi Kaario, Ville Vuorinen. Effects of hydrogen enrichment and turbulence intensity on the combustion mode in locally stratified dual-fuel mixtures of n-dodecane/methane. Applications in Energy and Combustion Science, 2022,  10, 100072. https://doi.org/10.1016/j.jaecs.2022.100072.

[3] Shervin Karimkashi, Parsa Tamadonfar, Ossi Kaario, Ville Vuorinen. A Numerical Investigation on Effects of Hydrogen Enrichment and Turbulence on NO Formation Pathways in Premixed Ammonia/Air Flames. Combustion Science and Technology, 2023,  https://doi.org/10.1080/00102202.2023.2180634.

[4] Parsa Tamadonfar, Shervin Karimkashi, Ossi Kaario, Ville Vuorinen. A study of flame dynamics and structure in premixed turbulent planar NH3/H2/air flames. International Journal of Engine Research, 2023, https://doi.org/10.1177/14680874231170657.

