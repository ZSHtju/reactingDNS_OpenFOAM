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

		Info<< "Reading  1D Flame data file" << endl;

		Info<< "Please make sure the scalar in SI unit" << endl;
		ifstream fileRow("1DFlame.dat");
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
		
		ifstream myfile("1DFlame.dat");
		if (myfile.fail())
		{
			Info << " file not found, the name of the file shoule be 1DFlame.dat" << endl;
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

		IOdictionary initialFlameProperties
		(
			IOobject
			(
				"initialFlameProperties",
				runTime.constant(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		word type(initialFlameProperties.lookup("type"));
		
		Switch reverse(initialFlameProperties.lookupOrDefault("reverse",false));
		dimensionedVector startPoint(initialFlameProperties.lookup("startPoint"));			
		scalar amplitudeFactor = readScalar(initialFlameProperties.lookup("amplitudeFactor"));
		scalar frequency= readScalar(initialFlameProperties.lookup("frequency"));

	if(type == "plane")
	{
		Switch twin(initialFlameProperties.lookupOrDefault("twin",false));		
		dimensionedVector flameNormalDirection(initialFlameProperties.lookup("flameNormalDirection")); 	
		dimensionedVector flameFluctuationDirection(initialFlameProperties.lookup("flameFluctuationDirection"));
		dimensionedScalar distance(initialFlameProperties.lookup("distance"));	
		

		forAll(mesh.C(),celli)
		{
			
			scalar r1 = distance.value();
			scalar r2 = distance.value() + d[NN-1]-d[1];
			scalar r = (mesh.C()[celli] - startPoint.value())&flameNormalDirection.value();
			
			if(twin)
			{
				r = std::abs((mesh.C()[celli] - startPoint.value())&flameNormalDirection.value());
			}

			scalar delta_r = ((d[NN-1]-d[1])*amplitudeFactor)*Foam::sin((mesh.C()[celli]&flameNormalDirection.value()) *frequency/(d[NN-1]-d[0]));
			
			r1 = r1 + delta_r;
			r2 = r2 + delta_r;
			
			if(!reverse)
			{
				if (r < r1)
				{
					T[celli] = TT[0];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][0];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;
				}
				else if(r>=r1 && r<=r2)
				{
					label MM = 0;
					for(label i = 0; i < NN-1; i++)
					{
						if((r-r1) >= (d[i]-d[1]) && (r-r1) <= (d[i+1]-d[1]))
						{
								MM = i;
								break;
						}
					}

					T[celli] = TT[MM]+(TT[MM+1]-TT[MM])*(r-r1-d[MM]+d[1])/(d[MM+1]-d[MM]);
					scalar Yt = 0;
					label inertIndex2 = -1;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] = YY[j+2][MM]+(YY[j+2][MM+1]-YY[j+2][MM])*(r-r1-d[MM]+d[1])/(d[MM+1]-d[MM]);
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex2 = j;
						}
					}
					Y[inertIndex2][celli] = 1 - Yt;
				}
				else
				{
					T[celli] = TT[NN-1];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][NN-1];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;				
				}
			}
			else
			{
				if (r < r1)
				{
					T[celli] = TT[NN-1];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][NN-1];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;
				}
				else if(r>=r1 && r<=r2)
				{
					label MM = 0;
					for(label i = 0; i < NN-1; i++)
					{
						if( (r2-r) >= (d[i]-d[1]) )
						{
							if ((r2-r) <= (d[i+1])-d[1] )
							{
								MM = i;
								break;
							}
						}
					}
					T[celli] = TT[MM]+(TT[MM+1]-TT[MM])*(r2-r-d[MM]+d[1])/(d[MM+1]-d[MM]);
					scalar Yt = 0;
					label inertIndex2 = -1;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] = YY[j+2][MM]+(YY[j+2][MM+1]-YY[j+2][MM])*(r2-r-d[MM]+d[1])/(d[MM+1]-d[MM]);
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex2 = j;
						}
					}
					Y[inertIndex2][celli] = 1 - Yt;
				}
				else
				{
					T[celli] = TT[0];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][0];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;				
				}				
			}
		}

	}
	
	
	if(type == "circle")
	{
		dimensionedScalar radius(initialFlameProperties.lookup("radius"));
		
		forAll(mesh.C(),celli)
		{
			scalar r1 = radius.value();
			scalar r2 = radius.value()+d[NN-1]-d[1];
			scalar r = Foam::sqrt(Foam::sqr(mesh.C()[celli].x()-startPoint.value().x()) + Foam::sqr(mesh.C()[celli].y()-startPoint.value().y()) +Foam::sqr(mesh.C()[celli].z()-startPoint.value().z())); 
			scalar theta = Foam::asin(Foam::mag(mesh.C()[celli].y()-startPoint.value().y())/r);	
			scalar y = 0;
			if(mesh.C()[celli].x() > startPoint.value().x() && mesh.C()[celli].y() < startPoint.value().y() )  
			{
				y = theta*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() < startPoint.value().x() && mesh.C()[celli].y() < startPoint.value().y() ) 
			{
				y  = (Foam::constant::mathematical::pi - theta)*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() < startPoint.value().x() && mesh.C()[celli].y() > startPoint.value().y() ) 
			{
				y  = (Foam::constant::mathematical::pi + theta)*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() > startPoint.value().x() && mesh.C()[celli].y() > startPoint.value().y() ) 
			{
				y  = (2*Foam::constant::mathematical::pi - theta)*(r1+r2)/2;
			}
			else if (mesh.C()[celli].x() ==startPoint.value().x() && mesh.C()[celli].y() < startPoint.value().y() ) 
			{
				y  = Foam::constant::mathematical::pi*(r1+r2)/4;
			}
			else  if (mesh.C()[celli].x() == startPoint.value().x() && mesh.C()[celli].y() > startPoint.value().y() ) 
			{
				y  = Foam::constant::mathematical::pi*(r1+r2)*3/4;
			}
			else if (mesh.C()[celli].x() > startPoint.value().x() && mesh.C()[celli].y() == startPoint.value().y() ) 
			{
				y  = 0;
			}
			else if(mesh.C()[celli].x() < startPoint.value().x() && mesh.C()[celli].y() == startPoint.value().y() ) 
			{
				y  = Foam::constant::mathematical::pi*(r1+r2)*2/4;
			}

			scalar delta_r = ((d[NN-2]-d[1])*amplitudeFactor)*Foam::sin(y*frequency*2/(r1+r2));
			r1 = r1 + delta_r;
			r2 = r2 + delta_r;
			if(!reverse)
			{
				if (r < r1)
				{
					T[celli] = TT[NN-1];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][NN-1];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;
				}
				else if(r>=r1 && r<=r2)
				{
					label MM = 0;
					for(label i = 0; i < NN-1; i++)
					{
						if( (r2-r) >= (d[i]-d[1]) )
						{
							if ((r2-r) <= (d[i+1])-d[1] )
							{
								MM = i;
								break;
							}
						}
					}
					T[celli] = TT[MM]+(TT[MM+1]-TT[MM])*(r2-r-d[MM]+d[1])/(d[MM+1]-d[MM]);
					scalar Yt = 0;
					label inertIndex2 = -1;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] = YY[j+2][MM]+(YY[j+2][MM+1]-YY[j+2][MM])*(r2-r-d[MM]+d[1])/(d[MM+1]-d[MM]);
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex2 = j;
						}
					}
					Y[inertIndex2][celli] = 1 - Yt;
				}
				else
				{
					T[celli] = TT[0];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][0];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;				
				}
			}
			else
			{
				if (r < r1)
				{
					T[celli] = TT[0];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][0];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;
				}
				else if(r>=r1 && r<=r2)
				{
					label MM = 0;
					for(label i = 0; i < NN-1; i++)
					{
						if((r-r1) >= (d[i]-d[1]) && (r-r1) <= (d[i+1]-d[1]))
						{
								MM = i;
								break;
						}
					}

					T[celli] = TT[MM]+(TT[MM+1]-TT[MM])*(r-r1-d[MM]+d[1])/(d[MM+1]-d[MM]);
					scalar Yt = 0;
					label inertIndex2 = -1;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] = YY[j+2][MM]+(YY[j+2][MM+1]-YY[j+2][MM])*(r-r1-d[MM]+d[1])/(d[MM+1]-d[MM]);
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex2 = j;
						}
					}
					Y[inertIndex2][celli] = 1 - Yt;
				}
				else
				{
					T[celli] = TT[NN-1];
					label inertIndex1 = -1;
					scalar Yt = 0;
					forAll(Y,j)
					{
						if(Y[j].name() != inertSpecie)
						{
							Y[j][celli] =  YY[j+2][NN-1];
							Yt += Y[j][celli];
						}
						else
						{
							inertIndex1 = j;
						}
					}
					Y[inertIndex1][celli] = 1 - Yt;				
				}								
			}
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
