/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::CCMChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    CCMChemistryModelI.H
    CCMChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef CCMChemistryModel_H
#define CCMChemistryModel_H

#include "StandardChemistryModel.H"
#include "Reaction.H"
#include "volFields.H"
#include "simpleMatrix.H"

#include "OFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class CCMChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class CCMChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
    // Private Member Functions

        //- Disallow copy constructor
        CCMChemistryModel(const CCMChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const CCMChemistryModel&);

        // override for CCM
       // template<class DeltaTType>
       // scalar solve(const DeltaTType& deltaT);

protected:
	// Protected members for CCM
	label nActiveCell_;
	scalar J_H_Ox_;
	scalar J_H_fu_;
	scalar J_C_Ox_;
	scalar J_C_fu_;
	scalar J_O_Ox_;
	scalar J_O_fu_;
	PtrList<scalarField> spInfo_;
	label	nMCCM_;
	scalarField ZoneSpan_;
	scalarField MZmin_;
	scalarField MZmax_;
	bool expChi_;
	bool chemistryLoadBalancing_;	
	bool pressureMapping_;

public:

    //- Runtime type information
    TypeName("CCM");


    // Constructors

        //- Construct from thermo
        CCMChemistryModel(ReactionThermo& thermo);


    //- Destructor
    virtual ~CCMChemistryModel();




    // CCM memebr function
          scalar solveCCM
         (
			const scalar deltaT,
            const label nBeg,
            const label nLast,
            const PtrList<scalarField>& FieldCCM,
            PtrList<scalarField>& RRCCM,
			const labelList cellIndexTmp
         );

         void readCCMProperties();  
		 
		virtual scalar CCMIndexMapSolve
		(
			const scalar deltaT,
			const PtrList<scalarField>& Yfield,
			const scalarField& Tfield,
			const scalarField& pfield,
			const scalarField& Zfield, // only for index
			const scalarField& Chifield, // only for index
			const scalarField& hsfield, // density is recalculated by enthalpy
			scalarField& CAtomNumber, // for baba correct
			scalarField& HAtomNumber, // for baba correct
			scalarField& ZoneNumber, // transfer RRCCM to RR
			label& gnActiveCell,
			PtrList<scalarField>& RRCCM
		);
		
		scalar  nActiveCellCCM() const;   
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "CCMChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CCMChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
