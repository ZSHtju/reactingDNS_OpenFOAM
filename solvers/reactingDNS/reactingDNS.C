/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactingDNS

Description
	Solver for combustion with chemical reactions and direct numerical simulation (Laminar turbulence model).
	The mixture-averaged transport properties model is hard-coded into reactingFoam solver.
	The detailed transport properties for each species are evaluated by a third-order logarithm polynomial fittting method.
	(Kee, R. J; et al. Sandia National Laboratories Report SAND86-8246 1986, 13, 80401-1887.)
	
	Author: Shenghui Zhong, State Key Laboratory of Engines, Tianjin Univerisity, Tianjin, China.
	zhongshenghui@tju.edu.cn
	
	This solver had been successfully applied in following papers.
	Zhong, S; et al. Fuel 2018, 234, 1044-1054.
	Zhong, S; et al. Fuel 2020, 275, 117980.
	Zhong, S; et al. Applied Energy 2020, 275, 115320.
	Zhang, N; et al. Combustion Theory and Modelling, 2020 1-22.
	
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

/* #include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "perfectFluid.H"
#include "adiabaticPerfectFluid.H"
#include "rhoConst.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include "sutherlandTransport.H"
#include "constTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "IFstream.H"
#include "OFstream.H"
#include "Switch.H"
#include "Random.H"
#include "mathematicalConstants.H" */
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    #include "readChemistryProperties.H"

	#include "createDiffFields.H"
	#include "createTDiffFields.H"
	#include "createLambdaFields.H"
	#include "createMuFields.H"

    turbulence->validate();


    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	Info<< "\nStarting time loop\n" << endl;
	label nStep = 0;

    #include "readTranData.H"
    #include "readBinaryDiff.H"
	#include "readThermalDiff.H"
    #include "readSpeciesLambda.H"
    #include "readSpeciesMu.H"
    //#include "getElementData.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"

		#include "compressibleCourantNo.H"
		#include "setDeltaT.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        nStep +=1;
        runTime++;

        Info<< "Time = " << runTime.timeName() << " nStep = "<< nStep << nl << endl;

        #include "rhoEqn.H"


        while (pimple.loop())
        {
            #include "updateTransProperties.H"
            #include "UEqn.H"

			//chemistry.solve(runTime.deltaT().value());
			reaction->correct();
			if(differentialDiffusion)
			{
				#include "Y-hEqn_DD.H"
			}
			else
			{
				#include "Y-hEqn.H"
			}

             // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

		rho = thermo.rho();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
       Info<< "End\n" << endl;

       return 0;
}


// ************************************************************************* //;
