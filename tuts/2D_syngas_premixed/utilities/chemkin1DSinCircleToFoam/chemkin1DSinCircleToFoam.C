/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "BasicChemistryModel.H"
#include "chemistrySolver.H"
#include "basicSpecieMixture.H"
#include "multiComponentMixture.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

//using namespace std;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	timeSelector::addOptions();
	#include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

	Info<< nl << "Reading thermophysicalProperties" << endl;

	autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
	psiReactionThermo& thermo = pThermo();
	thermo.validate(args.executable(), "h", "e");

	autoPtr< BasicChemistryModel<psiReactionThermo>> pChemistry
	(
	  BasicChemistryModel<psiReactionThermo>::New(thermo)
	);

	basicSpecieMixture& composition = thermo.composition();
	PtrList<volScalarField>& Y = composition.Y();

	const word inertSpecie(thermo.lookup("inertSpecie"));
	if (!composition.species().found(inertSpecie))
	{
		FatalIOErrorIn(args.executable().c_str(), thermo)
			<< "Inert specie " << inertSpecie << " not found in available species "
			<< composition.species() << exit(FatalIOError);
	}

	volScalarField& T = const_cast<volScalarField&>(thermo.T());

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI);

		Info<< "Time = " << runTime.timeName() << endl;

		Info<< "Reading CHEMKIN 1D data file" << endl;

		Info<< "Please make sure the scalar in SI unit" << endl;
		ifstream fileRow("CHEMKIN1D.txt");
		string tmp;
		label count = 0;

		while(getline(fileRow,tmp,'\n'))
		{
			if(tmp.size() > 0)
			count++;
		}
		label NN = count;
		Info << count << endl;
		scalarList d(NN,0.);
		scalarList TT(NN,0.);
		PtrList<scalarList> YY(Y.size()+2);
		for(label i = 0; i < (Y.size()+2) ; i++)
		{
			YY.set
			(
				i,
				new scalarList(NN,0.)
			);
		}
		
		ifstream myfile("CHEMKIN1D.txt");
		if (myfile.fail())
		{
			Info << " file not found, the name of the file shoule be CHEMKIN1D.txt" << endl;
		}
		else
		{
			char sline[300];
			for(label i = 0; i < NN ; i++)
			{
				for(label j = 0; j < (Y.size()+2); j++)
				{
					myfile  >>  YY[j][i];
				}
				myfile.getline(sline,300);
			}

			d = YY[0];
			TT = YY[1];

			Info << " distance "<< d << endl;
			Info << " Temperature " << TT << endl;

			forAll(Y,i)
			{
				Info << Y[i].name() << YY[i+2] << endl;
			}
		}

		IOdictionary flameBallProperties
		(
			IOobject
			(
			"flameBallProperties",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE,
			false
			)
		);

		scalar CCx =  readScalar(flameBallProperties.lookup("CCx"));
		scalar CCy =  readScalar(flameBallProperties.lookup("CCy"));
		scalar CCz =  readScalar(flameBallProperties.lookup("CCz"));

		dimensionedScalar innerBunrtRadius("innerBunrtRadius", dimLength , 0.0);
		flameBallProperties.lookup("innerBunrtRadius") >> innerBunrtRadius;
		dimensionedScalar outerBunrtRadius("outerBunrtRadius", dimLength , 0.0);
		flameBallProperties.lookup("outerBunrtRadius") >> outerBunrtRadius;

		scalar amplitudeFactor = readScalar(flameBallProperties.lookup("amplitudeFactor"));
		scalar frequency= readScalar(flameBallProperties.lookup("frequency"));

		scalar r = 0;
		scalar theta = 0;

		scalar delta_r = 0;
		volScalarField y
		(
			IOobject
			(
				"y",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("y", dimless, 0.0)
		);

		forAll(mesh.C(),celli)
		{
			scalar r1 = innerBunrtRadius.value();
			scalar r2 = outerBunrtRadius.value();
			r = Foam::sqrt(Foam::sqr(mesh.C()[celli].x()-CCx) + Foam::sqr(mesh.C()[celli].y()-CCy) +Foam::sqr(mesh.C()[celli].z()-CCz)); 
			theta = Foam::asin(Foam::mag(mesh.C()[celli].y()-CCy)/r);	
			if(mesh.C()[celli].x() > CCx && mesh.C()[celli].y() < CCy )  
			{
				y[celli] = theta*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() < CCx && mesh.C()[celli].y() < CCy ) 
			{
				y[celli]  = (Foam::constant::mathematical::pi - theta)*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() < CCx && mesh.C()[celli].y() > CCy ) 
			{
				y[celli]  = (Foam::constant::mathematical::pi + theta)*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() > CCx && mesh.C()[celli].y() > CCy ) 
			{
				y[celli]  = (2*Foam::constant::mathematical::pi - theta)*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() ==CCx && mesh.C()[celli].y() < CCy ) 
			{
				y[celli]  = Foam::constant::mathematical::pi*(r1+r2)/4;
			}
			else  if (mesh.C()[celli].x() == CCx && mesh.C()[celli].y() > CCy ) 
			{
				y[celli]  = Foam::constant::mathematical::pi*(r1+r2)*3/4;
			}
			else if (mesh.C()[celli].x() > CCx && mesh.C()[celli].y() == CCy ) 
			{
				y[celli]  = 0;
			}
			else if(mesh.C()[celli].x() < CCx && mesh.C()[celli].y() == CCy ) 
			{
				y[celli]  = Foam::constant::mathematical::pi*(r1+r2)*2/4;
			}

			delta_r = ((r2-r1)/amplitudeFactor)*Foam::sin( y[celli] *frequency*2/(r1+r2) );	
			r1 = r1 + delta_r;
			r2 = r2 + delta_r;

			if (r < r1)
			{
				T[celli] = TT[NN-1];
				label inertIndex1 = -1;
				scalar Yt1 = 0;
				forAll(Y,j)
				{
					if(Y[j].name() != inertSpecie)
					{
						Y[j][celli] =  YY[j+2][NN-1];
						Yt1 += Y[j][celli];
					}
					else
					{
						inertIndex1 = j;
					}
			}

			Y[inertIndex1][celli] = 1 - Yt1;
			}

			if(r>=r1 && r<=r2)
			{
				label MM = 0;
				for(label i = 0; i < NN; i++)
				{
					if( (r2-r) >= (d[i]-d[0]) )
					{
						if ((r2-r) <= (d[i+1])-d[0] )
						{
							MM = i;
							break;
						}
					}
				}

				T[celli] = TT[MM]+(TT[MM+1]-TT[MM])*(r2-r-d[MM]+d[0])/(d[MM+1]-d[MM]);
				scalar Yt2 = 0;
				label inertIndex2 = -1;
				forAll(Y,j)
				{
				if(Y[j].name() != inertSpecie)
				{
				Y[j][celli] = YY[j+2][MM]+(YY[j+2][MM+1]-YY[j+2][MM])*(r2-r-d[MM]+d[0])/(d[MM+1]-d[MM]);
				Yt2 += Y[j][celli];
				}
				else
				{
				inertIndex2 = j;
				}
				}

				Y[inertIndex2][celli] = 1 - Yt2;
			}
		}

		forAll(Y,i)
		{
			Y[i].write();
		}
			T.write();
		}

		Info<< "End\n" << endl;
		return 0;
}


// ************************************************************************* //
