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

\*---------------------------------------------------------------------------*/

#include "CCMChemistryModel.H"
//#include "SpecieMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"

#include "atomicWeights.H"
#include "thermoPhysicsTypes.H"
#include "chemistryReader.H"
#include "clockTime.H"
#include "HashTableSH.H"
#include "Random.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::CCMChemistryModel<ReactionThermo, ThermoType>::CCMChemistryModel
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo),	
	// CCM protected members
	nActiveCell_(0),
	J_H_Ox_(0),
	J_H_fu_(0),
	J_C_Ox_(0),
	J_C_fu_(0),
	J_O_Ox_(0),
	J_O_fu_(0),
	spInfo_(this->nSpecie_),
	nMCCM_(3),
	ZoneSpan_(this->nSpecie_+4, -1),
	MZmin_(this->nSpecie_+4, -great),
	MZmax_(this->nSpecie_+4, great),
	expChi_(true),
	chemistryLoadBalancing_(true),
    pressureMapping_(true)	
	{
		const label nSizeSpInfo = 8;
		forAll(spInfo_,i)
		{
			spInfo_.set
			(
				i,
				new scalarField(nSizeSpInfo, 0.0)
			);
		}
	}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template<class ReactionThermo, class ThermoType>
Foam::CCMChemistryModel<ReactionThermo, ThermoType>::
~CCMChemistryModel()
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::CCMChemistryModel<ReactionThermo, ThermoType>::solveCCM
(
	const scalar deltaT,
	const label nBeg,
	const label nLast,
	const PtrList<scalarField>& FieldCCM,
	PtrList<scalarField>& RRCCM,
	const labelList cellIndexTmp
)
{

  //  BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);

    for(label icell=nBeg; icell < nLast; icell ++)
    {
		label celli = cellIndexTmp[icell];
		
        scalar Ti = FieldCCM[0][celli];
		if (Ti  > this->Treact_)
        {
        scalar pi = FieldCCM[1][celli];
        const scalar rhoi = FieldCCM[2][celli];

            for (label i=0; i<this->nSpecie_; i++)
            {
                c[i] = rhoi*FieldCCM[i+5][celli]/this->specieThermo_[i].W();
                c0[i] = c[i];
            }
		
            // Initialise time progress
            scalar timeLeft = deltaT;
            scalar deltaTChemIni = this->deltaTChemIni_;

            // Calculate the chemical source terms
            while (timeLeft > SMALL)
            {
                scalar dt = timeLeft;
                this->solve(c, Ti, pi, dt, deltaTChemIni);
                timeLeft -= dt;
            }

            deltaTMin = min(deltaTChemIni, deltaTMin);
           // Info << "deltaTMin = " << deltaTMin << endl;
		   // Info << "deltaTChemIni_ = " << this->deltaTChemIni_ << endl;
			
            for (label i=0; i<this->nSpecie_; i++)
            {
               RRCCM[i][celli] = (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT/rhoi;
            }
		}
		else
        {
            for (label i=0; i<this->nSpecie_; i++)
            {
                RRCCM[i][celli] = 0;
            }
        }
    }

    return deltaTMin;

}


// ************************************************************************* //
template<class ReactionThermo, class ThermoType>
void Foam::CCMChemistryModel<ReactionThermo, ThermoType>::readCCMProperties()
{
	IOdictionary CCMProperties
	(
		IOobject
		(
			"CCMProperties",
			this->time().constant(),
			this->mesh(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE,
			false
		)
	);

        List<List<specieElement> > specieComp(this->nSpecie_);
	    speciesTable speciesTab = this->thermo().composition().species();
        autoPtr<chemistryReader<ThermoType>> tchemRead
        (
        chemistryReader<ThermoType>::New(this->thermo(),speciesTab)
        );


        const HashTable<List<specieElement> >& specComp(tchemRead().specieComposition());

	forAll(specieComp,i)
	{
		specieComp[i] = specComp[this->Y_[i].name()];
	}

	wordList elementNames(4);
	{
		elementNames[0] = "H";
		elementNames[1] = "C";
		elementNames[2] = "O";
		elementNames[3] = "N";

	}

	forAll(this->Y_, i)
	{
		Info << "specie = " << this->Y_[i].name() << "\t";
		forAll(elementNames, en)
		{
				forAll(specieComp[i], kk)
				{
					
					if(specieComp[i][kk].name() == elementNames[en])
					{
						Info << specieComp[i][kk].name() << ", nAtoms = " << specieComp[i][kk].nAtoms() << "\t";
						spInfo_[i][en] = specieComp[i][kk].nAtoms();
					}
				}

		}

		Info << this->specieThermo_[i].W() << endl;
		spInfo_[i][4] = this->specieThermo_[i].W();
	}


	CCMProperties.lookup("J_H_Ox:J_H_fu")>> J_H_Ox_ >> J_H_fu_;
	CCMProperties.lookup("J_C_Ox:J_C_fu")>> J_C_Ox_ >> J_C_fu_;
	CCMProperties.lookup("J_O_Ox:J_O_fu")>> J_O_Ox_ >> J_O_fu_;

	Info <<" J_H_Ox= "<< J_H_Ox_ << " J_H_fu= " << J_H_fu_<< endl;
	Info <<" J_C_Ox= "<< J_C_Ox_ << " J_C_fu= " << J_C_fu_<< endl;
	Info <<" J_O_Ox= "<< J_O_Ox_ << " J_O_fu= " << J_O_fu_<< endl;


	CCMProperties.lookup("min:max:SpanZoneJe"  ) >> MZmin_[0] >> MZmax_[0] >> ZoneSpan_[0];
	CCMProperties.lookup("min:max:SpanZoneT"   ) >> MZmin_[1] >> MZmax_[1] >> ZoneSpan_[1];
	CCMProperties.lookup("min:max:SpanZoneXi"  ) >> MZmin_[2] >> MZmax_[2] >> ZoneSpan_[2];

	CCMProperties.lookup("expChi" ) >> expChi_;
    CCMProperties.lookup("pressureMapping" ) >> pressureMapping_;


	IFstream speciesInfo(this->time().constant()/"CCMProperties");
	Foam::string   headerStr;
	Foam::word     thisSpecies;
	scalar curSpSpan=-1.;
	scalar maxSpMass=-1.;
	scalar minSpMass=-1.;


	Foam::string CCMkeyword=' ';
	Foam::string CCMKEYWORD="CCMspecies";
	while(CCMkeyword != CCMKEYWORD && !speciesInfo.eof() )
	{
		speciesInfo.Foam::ISstream::getLine(headerStr);
		CCMkeyword = headerStr.operator()(CCMKEYWORD.size() );
	}

	if(CCMkeyword != CCMKEYWORD)
	{
		FatalErrorIn("CCMreader cannot find keyword \"CCMspecies\" in your CCMproperties ")
		<< abort(FatalError);
	}

	speciesInfo.Foam::ISstream::getLine(headerStr);

	while( !speciesInfo.eof() )
	{
		speciesInfo >> thisSpecies;
		if(!speciesInfo.bad() )
		{
			speciesInfo >> minSpMass >> maxSpMass >>curSpSpan;
			speciesInfo.Foam::ISstream::getLine(headerStr);
			Info <<" this Species "<< thisSpecies << " is added to the CCM phase space with resolution of: " << curSpSpan << endl;
			for(label ii=0; ii< this->nSpecie_ ; ii++)
			{
				if(thisSpecies == this->Y_[ii].name())
				{ 
					spInfo_[ii][5]=curSpSpan;
					spInfo_[ii][6]=minSpMass;
					spInfo_[ii][7]=maxSpMass;
					if(spInfo_[ii][5]>0)
					{
						Info <<  this->Y_[ii].name() << " dimention is added to CCM " << endl;
						nMCCM_ ++;
					}
					curSpSpan=-1;
				}

			}
		}
	}

	for(label ii=0; ii<this->nSpecie_; ii++)
	{
		 Info << this->Y_[ii].name() <<tab
		 <<"H " << spInfo_[ii][0] <<tab  <<"C " <<  spInfo_[ii][1] <<tab
		 <<"O " << spInfo_[ii][2] << tab <<"N  "<<  spInfo_[ii][3] <<tab
		 <<"MW "<< spInfo_[ii][4]<< tab << "speciesSpan "<< spInfo_[ii][5]<< tab
		 <<" minSp "<< spInfo_[ii][6] << tab<<" maxSp "<< spInfo_[ii][7] <<"\n"<< endl;
	}

	label mccm_YM =3;
	for(label i=0; i < this->nSpecie_; i++)
	{
	  if(spInfo_[i][5] > 0)
	  {
		ZoneSpan_[mccm_YM] = spInfo_[i][5];
		MZmin_[mccm_YM]    = spInfo_[i][6];
		MZmax_[mccm_YM]    = spInfo_[i][7];
		mccm_YM ++;
	  }
	}

	for(label mccm=0; mccm < nMCCM_; mccm ++)
	{
		scalar N_Zone =(MZmax_[mccm]-MZmin_[mccm])/ZoneSpan_[mccm]+1;
		Info << mccm <<" min(MZ)= " << MZmin_[mccm] << " max(MZ)= " << MZmax_[mccm] << " ZoneSpan= "<< ZoneSpan_[mccm] << " Nu. Zone in CCMMesh "  << N_Zone << endl;
	}

	if(pressureMapping_)
	{
		nMCCM_++; //add new dimension for pressure
		CCMProperties.lookup("min:max:SpanZoneP"  ) >> MZmin_[nMCCM_-1] >> MZmax_[nMCCM_-1] >> ZoneSpan_[nMCCM_-1];
	}
}

template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::CCMChemistryModel<ReactionThermo, ThermoType>::CCMIndexMapSolve
(
const scalar deltaT,
const PtrList<scalarField>& Yfield,
const scalarField& Tfield,
const scalarField& pfield,
const scalarField& Zfield,    // only for index
const scalarField& Chifield, // only for index
const scalarField& hsfield, // density is recalculated by enthalpy
scalarField& CAtomNumber,  // for baba correct
scalarField& HAtomNumber,   // for baba correct
scalarField& ZoneNumber, // transfer RRCCM to RR
label& gnActiveCell,
PtrList<scalarField>& RRCCM
)
{
	    const clockTime clockTime_= clockTime();
        clockTime_.timeIncrement();
		scalar timeForCCMindex = 0.0;
		scalar timeForCCMmap = 0.0;
		scalar timeForCCMreaction = 0.0;
		nActiveCell_	= 0;
		label fieldSize = Tfield.size();
		scalarField ZoneIndex(fieldSize,0);
		//scalarField ZoneNumber(fieldSize,0); // declare in solver
		
		List<scalarField> MZ(nMCCM_);
		forAll(MZ, nz)  MZ[nz].setSize(fieldSize ,0.);

		forAll(MZ, mccm) MZ[mccm] *= 0;

		MZ[0] = Zfield;
			
		MZ[1] = Tfield;

		MZ[2] = Chifield;
		

		for(label icell= 0; icell < fieldSize; icell++)
		{
			if(expChi_) MZ[2][icell] = Foam::exp(-MZ[2][icell]);
		}

		label mccm_YM =3;
		for(label i=0; i < this->nSpecie_; i++)
		{
		 if(spInfo_[i][5] > 0)
		 {

			for(label icell= 0; icell < fieldSize; icell++)
			{
			  MZ[mccm_YM][icell] = Yfield[i][icell];
			}
			  mccm_YM ++;
		 }
		}

		scalarField allMinMZ(nMCCM_,0);
		scalarField allMaxMZ(nMCCM_,0);
		
		for(label mccm=0; mccm < MZ.size(); mccm ++)
		{
			allMinMZ[mccm] =  min(MZ[mccm]);
			allMaxMZ[mccm] =  max(MZ[mccm]);
			reduce(allMaxMZ[mccm], maxOp<scalar>());
			reduce(allMinMZ[mccm], minOp<scalar>());
		}
		
		for(label mccm=0; mccm < MZ.size(); mccm ++)
		{
			Info <<" mccm "<< mccm << " max/min MZ "<< allMaxMZ[mccm] <<" "<<  allMinMZ[mccm] << endl;
		}
		const labelList HashSizeList = {1, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
		label  curHashSize = fieldSize*0.5;
		label  whichSize = 10;
		for(int i = 0; i < (HashSizeList.size() - 1); i++)
		{
			if( curHashSize < HashSizeList[i+1] && curHashSize >= HashSizeList[i])
			{
				whichSize = i+1;
			}
		}
		HashTableSH HashMapP(HashSizeList[whichSize]); //HashMapP(1280)
		int nActiveCell = 0;
		scalarField ZoneNumberIndex;
		for(label icell= 0; icell < fieldSize; icell++)
		{
		  ZoneIndex[icell]   = 0;
		  scalar Zindex=0;
		  scalarField nZone(MZ.size(),0);							
		 for(label mccm=0; mccm < MZ.size(); mccm ++)
		 {
			//Map in MZ direction
			scalar   left = 0;
			nZone[mccm] = (MZmax_[mccm]-MZmin_[mccm])/ZoneSpan_[mccm]+1;			//scalar   nZone = (MZmax_[mccm]-MZmin_[mccm])/ZoneSpan_[mccm]+1;
			scalar   right = nZone[mccm];											//scalar   right = nZone;
			while(right>left+1)
			{
				label  mid = (left+right)/2;
				if(MZ[mccm][icell] >= mid*ZoneSpan_[mccm]+MZmin_[mccm])
				{
					left = mid;
				}
				else
				{
					right = mid;
				}
			}

			Zindex = left;
			scalar termIndex =1;
			if(mccm > 0 )
			{
				for(label _mccm = mccm; _mccm > 0; _mccm --) termIndex *= nZone[_mccm-1]; //for(label _mccm = mccm; _mccm > 0; _mccm --)termIndex *=nZone;
			}
		    ZoneIndex[icell] += Zindex * termIndex;
		}
			// processor index
			long int index = ZoneIndex[icell];
            int ZN = -1;
			if(HashMapP.set(index, nActiveCell,true, ZN))
			{
				ZoneNumber[icell] = nActiveCell;
				ZoneNumberIndex.resize(nActiveCell+1);
				ZoneNumberIndex[nActiveCell] = index;
				nActiveCell ++;
			}
			else
			{
				ZoneNumber[icell] = ZN;
			}			 
		}
/////////////////////////////Shijie Xu, 20190508; Shenghui Zhong 20190805 
// global index
	label ZoneSize = ZoneNumberIndex.size();
	label gZoneSize = ZoneSize;
	reduce(gZoneSize, sumOp<scalar>());
	scalarField gZoneNumberIndex(gZoneSize, 0);
	scalarField gCCMZoneNumber(gZoneSize, 0); 
	scalarField CCMZoneNumber(ZoneSize,0);


	label nProc=Pstream::nProcs();
	label curProc = Pstream::myProcNo();


	labelList allZoneSize(nProc,0);
	allZoneSize[curProc] = ZoneSize;
	forAll(allZoneSize,iproc)
	{
		reduce(allZoneSize[iproc], sumOp<scalar>());
	}

	label nZoneStart    = 0;
	for(label iproc = 0; iproc < curProc; iproc++) nZoneStart += allZoneSize[iproc];


	for(label icell = 0; icell < ZoneSize; icell++) 
	{
	   label gicell = icell + nZoneStart;
	   gZoneNumberIndex[gicell] = ZoneNumberIndex[icell]; 
	}
	reduce(gZoneNumberIndex, sumOp<scalarField>());
			
	//Shenghui & Shijie  20200331 redistribute the CCM gZoneIndex based on the remainder
	scalarField pZoneNumberIndex;
	label nCount = 0;
	forAll(gZoneNumberIndex,i)
	{
		long int index = gZoneNumberIndex[i];
		label remainder  =  index%nProc;
		if(remainder == curProc)
		{
			pZoneNumberIndex.resize(nCount+1);
			pZoneNumberIndex[nCount] = gZoneNumberIndex[i];			
			nCount++;			
		}
	}
		curHashSize = nCount*0.8;
		whichSize = 10;
		for(int i = 0; i < (HashSizeList.size() - 1); i++)
		{
			if( curHashSize < HashSizeList[i+1] && curHashSize >= HashSizeList[i])
			{
				whichSize = i+1;
			}
		}
		HashTableSH HashMapGP(HashSizeList[whichSize]); //HashMapP(1280)	
		int GnActiveCell = 0;
		int GZN = -1;
		scalarField gpCCMZoneNumber(nCount,0);
		for(label gicell= 0; gicell < nCount; gicell++)
		{
				long int index = pZoneNumberIndex[gicell];
				if(HashMapGP.set(index, GnActiveCell,true, GZN))
				{
					//gpCCMZoneNumber.resize(GnActiveCell+1);
					gpCCMZoneNumber[gicell]  = GnActiveCell;
					GnActiveCell ++;
				}
				else
				{
					gpCCMZoneNumber[gicell] = GZN;
				}
		}
				
		labelList gnActiveCellProcessor(nProc,0);
		gnActiveCellProcessor[curProc] = GnActiveCell;
		forAll(gnActiveCellProcessor,iproc)
		{
			reduce(gnActiveCellProcessor[iproc], sumOp<scalar>());
		}
				
		label 	deltaGZoneNumber = 0;
		for(label iproc = 0; iproc < curProc; iproc++)
		{
			deltaGZoneNumber +=  gnActiveCellProcessor[iproc];			
		}
		for(label gicell= 0; gicell < nCount; gicell++)
		{
			gpCCMZoneNumber[gicell] += deltaGZoneNumber;
		}
				
	    label NCount = 0;
		forAll(gZoneNumberIndex,i)
		{
			long int index = gZoneNumberIndex[i];
			label remainder  =  index%nProc;
			if(remainder == curProc)
			{
				gCCMZoneNumber[i] = gpCCMZoneNumber[NCount];
				NCount++;
			}
		}
	  		reduce(gCCMZoneNumber, sumOp<scalarField>());
				
	for(label icell = 0; icell < ZoneSize; icell++) 
	{
	   label gicell = icell + nZoneStart;
	   CCMZoneNumber[icell] = gCCMZoneNumber[gicell]; 
	}
	gnActiveCell = GnActiveCell;
	reduce(gnActiveCell, sumOp<scalar>());
	for(label icell = 0; icell < fieldSize; icell++) 
	{
	   label gicell = ZoneNumber[icell];
	   ZoneNumber[icell] = CCMZoneNumber[gicell]; 
	}
		timeForCCMindex = clockTime_.timeIncrement();
		Info << "try hard, do better"  <<  endl;
		Info << "CCM index done "  << "("  << timeForCCMindex << " s)" <<  "."  << endl;
	//


/* 		HashTableSH HashMapG(1280);
		int GnActiveCell = 0;
		int GZN = -1;
		for(label gicell= 0; gicell < gZoneSize; gicell++)
		{
				long int index = gZoneNumberIndex[gicell];
				if(HashMapG.set(index, GnActiveCell,true, GZN))
				{
					gCCMZoneNumber[gicell]  = GnActiveCell;
					//ZoneNumberIndex.resize(nActiveCell+1);
					//ZoneNumberIndex[nActiveCell] = index;
					GnActiveCell ++;
				}
				else
				{
					gCCMZoneNumber[gicell] = GZN;
				}
		}
        gnActiveCell = GnActiveCell;

		
	for(label icell = 0; icell < ZoneSize; icell++) 
	{
	   label gicell = icell + nZoneStart;
	   CCMZoneNumber[icell] = gCCMZoneNumber[gicell]; 
	}	
	
	for(label icell = 0; icell < fieldSize; icell++) 
	{
	   label gicell = ZoneNumber[icell];
	   ZoneNumber[icell] = CCMZoneNumber[gicell]; 
	}

		timeForCCMindex = clockTime_.timeIncrement();
		Info << "try hard, do better"  <<  endl;
		Info << "CCM index done "  << "("  << timeForCCMindex << " s)" <<  "."  << endl; */

	PtrList<scalarField> gFieldCCM(this->nSpecie_+5);
	forAll(gFieldCCM, i)
	{
		gFieldCCM.set
		(
			i,
			new scalarField(gnActiveCell,0.)
		);	
	}

	scalarField gnCell_CCM(gnActiveCell, 0.0);
   // scalarField gTimeOfOde_CCM(gnActiveCell, 0.0);
	for(label icell= 0; icell < fieldSize; icell++)
	{
		label izone            = ZoneNumber[icell];
		gFieldCCM[0][izone]    += Tfield[icell];
		gFieldCCM[1][izone]    += pfield[icell];
		gFieldCCM[3][izone]    += hsfield[icell];
     //   gTimeOfOde_CCM[izone] += timeOfOde[icell];
		for (label i=0; i<this->nSpecie_; i++)
		{
			gFieldCCM[5+i][izone]   += Yfield[i][icell];
		}
		gnCell_CCM[izone] ++;
	}

	forAll(gFieldCCM, i)
	{
		reduce(gFieldCCM[i], sumOp<scalarField>());
	}
	
	reduce(gnCell_CCM, sumOp<scalarField>());
//	reduce(gTimeOfOde_CCM, sumOp<scalarField>());
	


// FOR CHEMISTRY
		 
/* 	List <scalarField> RRCCM(this->nSpecie_);
	forAll(RRCCM,i)
	{
		RRCCM[i].resize(gnActiveCell);
	}
	PtrList<scalarField> gRRCCM(this->nSpecie_ );*/
	forAll(RRCCM, i)
	{
				RRCCM.set
				(
					i,
					new scalarField(gnActiveCell,0.)
				);
	} 
	
	label nBeg  = curProc*(gnActiveCell/nProc); //250717
	label nLast = nBeg + gnActiveCell/nProc;

	if(curProc == nProc - 1 ) nLast = gnActiveCell;

	 //Shenghui Zhong , 20200301, random load balance
	 //Random access to mesh cell
	 labelList cellIndexTmp = identity(gnActiveCell);
	 Random randGenerator(unsigned(time(NULL)));
	 //Pout << unsigned(time(NULL)) << endl;
	 //Random randGenerator(gnActiveCell);
	 label j;
	 for(label i = 0; i < gnActiveCell; i++)
	 {
		 j = randGenerator.sampleAB<int>(i,gnActiveCell);
		 label tmp = cellIndexTmp[i];
		 cellIndexTmp[i] = cellIndexTmp[j];
		 cellIndexTmp[j] = tmp;			 
	 } 
	 //Pout << "cellIndexTmp = " << cellIndexTmp << endl;
	 //
//
    for(label zonei =nBeg; zonei < nLast; zonei ++)//for(label izone =nBeg; izone < nLast; izone ++)
    {
		label izone = cellIndexTmp[zonei];
        gFieldCCM[0][izone] /= gnCell_CCM[izone];
        gFieldCCM[1][izone] /= gnCell_CCM[izone];
        gFieldCCM[3][izone] /= gnCell_CCM[izone];
      //  gTimeOfOde_CCM[izone] /= gnCell_CCM[izone];
        scalar Ytot=0;
        for (label i=0; i<this->nSpecie_; i++)
        {
			if(this->Y_[i].name() == "FN2") // faker species (reset to 0) cannot affect the real mixture composition, only for reform gas case, Shenghui Zhong
			{
			 gFieldCCM[5+i][izone] =  0.0;
			 Ytot +=gFieldCCM[5+i][izone];		
			}
			else
			{
			 gFieldCCM[5+i][izone] /= gnCell_CCM[izone];
			 Ytot +=gFieldCCM[5+i][izone];				
			}
		}
		for (label i=0; i<this->nSpecie_; i++)	gFieldCCM[5+i][izone] /= Ytot; 
        
    }


    for(label zonei =nBeg; zonei < nLast; zonei ++)//for(label izone =nBeg; izone < nLast; izone ++)
    {
		label izone = cellIndexTmp[zonei];
        ThermoType mixture(0.0*this->specieThermo_[0]);

        gFieldCCM[4][izone] = 0;
        for(label kk=0; kk<this->nSpecie_; kk++)
        {
            const scalar spHci = this->specieThermo_[kk].Hc();
            gFieldCCM[4][izone] += spHci*gFieldCCM[5+kk][izone] ;
        }
        for(label kk=0; kk < this->nSpecie_; kk++)
        {
            mixture += gFieldCCM[5+kk][izone] *this->specieThermo_[kk];    
        }
		scalar totH = gFieldCCM[3][izone] + gFieldCCM[4][izone];
	    gFieldCCM[0][izone]   = mixture.THa(totH, gFieldCCM[1][izone], gFieldCCM[0][izone]);
        gFieldCCM[2][izone]   = mixture.rho(gFieldCCM[1][izone], gFieldCCM[0][izone]);
    }

// for correctCH
    CAtomNumber.setSize(gnActiveCell,0);
	HAtomNumber.setSize(gnActiveCell,0);
    for(label zonei =nBeg; zonei < nLast; zonei ++)//for(label izone =nBeg; izone < nLast; izone ++)
    {
		 label izone = cellIndexTmp[zonei];
		 for(label i=0; i<this->nSpecie_; i++)
		 {
			if(this->Y_[i].name() != "H2O" && this->Y_[i].name() != "h2o")
			HAtomNumber[izone] += atomicWeights["H"]*gFieldCCM[i+5][izone]*spInfo_[i][0]/spInfo_[i][4];

			if(this->Y_[i].name() != "CO2" && this->Y_[i].name() != "co2")
			CAtomNumber[izone] += atomicWeights["C"]*gFieldCCM[i+5][izone]*spInfo_[i][1]/spInfo_[i][4];
		}
	}	
	reduce(CAtomNumber, sumOp<scalarField>());
	reduce(HAtomNumber, sumOp<scalarField>());

/////////////////////////////Shijie Xu, 20190508 And Shenghui Zhong , 20190723
	     timeForCCMmap = clockTime_.timeIncrement();
		 Info << "CCM map done "  <<  "("  << timeForCCMmap << " s)" <<  "."  << endl;
		 //

			this->solveCCM
			(
				deltaT,
				nBeg,
				nLast,
				gFieldCCM,
				RRCCM,
				cellIndexTmp  
			);
// Shenghui Zhong 2019-7-12      
			forAll(RRCCM,i)
			{
				reduce(RRCCM[i], sumOp<scalarField>());
			}
			
	timeForCCMreaction = clockTime_.timeIncrement();
    Info << "CCM reaction done "  <<  "("  << timeForCCMreaction << " s)" <<  "."  << endl;

	return gnActiveCell;
}
// ************************************************************************* //
