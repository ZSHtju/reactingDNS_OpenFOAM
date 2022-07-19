mixMu *= 0;
mixLambdaA *= 0 ; // reset the number to zero;
mixLambdaB *= 0 ;
DiffError *= 0 ; 
// get mixture mole mass
MW *= 0;
forAll(Y,i)
{
    MW += Y[i]/composition.Wi(i);
}
MW = 1/MW;
// get species mole fraction
forAll(Y,i)
{
    
    X[i] = Y[i]*MW/composition.Wi(i);
}

if(differentialDiffusion)
{
	forAll(groupSpecies,i)
	{
		if(groupSpecies[i][0] != -1)
		{
			for(label j = i; j < groupSpecies.size();j ++)
			{
				if(groupSpecies[j][0] != -1)
				{
					label m = max(groupSpecies[i][0],groupSpecies[j][0]);
					label n = min(groupSpecies[i][0],groupSpecies[j][0]);	
					label k = Y.size()*n+m-0.5*n*(n+1);
					volScalarField tmpd = 1e-4 * exp(Diff1[k]*pow(log(T/TUnit),0) + Diff2[k]*pow(log(T/TUnit),1) + Diff3[k]*pow(log(T/TUnit),2) + Diff4[k]*pow(log(T/TUnit),3)) / (p / 101325 /pUnit);
					forAll(groupSpecies[i],h)
					{
						forAll(groupSpecies[j],z)
						{
								label m = max(groupSpecies[i][h],groupSpecies[j][z]);
								label n = min(groupSpecies[i][h],groupSpecies[j][z]);	
								label k = Y.size()*n+m-0.5*n*(n+1);
								d[k] = tmpd;
						}
					}
				}
			}
		}
	}
/* 	for (label k = 0; k < (0.5*Y.size()*(Y.size()+1)); k++)
    {
        d[k] = 1e-4 * exp(Diff1[k]*pow(log(T/TUnit),0) + Diff2[k]*pow(log(T/TUnit),1) + Diff3[k]*pow(log(T/TUnit),2) + Diff4[k]*pow(log(T/TUnit),3)) / (p / 101325 /pUnit);// 1e-4 cm2/s to m2/s   p 1 dyne/cm2 = 0.1 Pa //
    } */
}

forAll(groupSpecies,i)
{
	if(groupSpecies[i][0] != -1)
	{
		Lambda[groupSpecies[i][0]] = lambdaUnit*1e-5*exp(Lambda1[groupSpecies[i][0]]*pow(log(T/TUnit),0) + Lambda2[groupSpecies[i][0]]*pow(log(T/TUnit),1) + Lambda3[groupSpecies[i][0]]*pow(log(T/TUnit),2) + Lambda4[groupSpecies[i][0]]*pow(log(T/TUnit),3));	
		Mu[groupSpecies[i][0]] = 0.1*MuUnit*exp(Mu1[groupSpecies[i][0]]*pow(log(T/TUnit),0) + Mu2[groupSpecies[i][0]]*pow(log(T/TUnit),1) + Mu3[groupSpecies[i][0]]*pow(log(T/TUnit),2) + Mu4[groupSpecies[i][0]]*pow(log(T/TUnit),3));		
		forAll(groupSpecies[i],j)
		{
			if(j > 0)
			{
				Lambda[groupSpecies[i][j]] = Lambda[groupSpecies[i][0]];
				Mu[groupSpecies[i][j]] = Mu[groupSpecies[i][0]];
			}
		}		
	}
}

forAll(Y,i)
{
//volScalarField& 	Xi = X[i];
	if(differentialDiffusion)
	{
		volScalarField& 	Yi = Y[i];
		DiNumerator *= 0 ;//reset to 0
		DiDenominator *= 0;
		for (label j = 0; j< Y.size(); j++)
		{
			label m = max(i,j);
			label n = min(i,j);
			volScalarField& Xj = X[j];

			if(j != i)
			{
			 label k = Y.size()*n+m-0.5*n*(n+1);
			 //d[k] = 1e-4 * exp(Diff1[k]*pow(log(T/TUnit),0) + Diff2[k]*pow(log(T/TUnit),1) + Diff3[k]*pow(log(T/TUnit),2) + Diff4[k]*pow(log(T/TUnit),3)) / (p / 101325 /pUnit);// 1e-4 cm2/s to m2/s   p 1 dyne/cm2 = 0.1 Pa //
			 DiNumerator += (Xj+1e-12) * composition.Wi(j);
			 DiDenominator += MW*(Xj+1e-12)/d[k];
			} 
		}
		Diff[i] = DiNumerator*DUnit/ DiDenominator; //*********** Di m2/s*******************//
		//DiffError += (rho*Diff[i]*fvc::grad(Yi) + Yi*rho*Diff[i]*fvc::grad(MW)/MW);
	}
	//Lambdai = lambdaUnit*1e-5*exp(Lambda1[i]*pow(log(T/TUnit),0) + Lambda2[i]*pow(log(T/TUnit),1) + Lambda3[i]*pow(log(T/TUnit),2) + Lambda4[i]*pow(log(T/TUnit),3));
	// mixLambdaA += Xi*Lambda[i];
	// mixLambdaB += Xi/Lambda[i];
	//		
	//Mu[i] = 0.1*MuUnit*exp(Mu1[i]*pow(log(T/TUnit),0) + Mu2[i]*pow(log(T/TUnit),1) + Mu3[i]*pow(log(T/TUnit),2) + Mu4[i]*pow(log(T/TUnit),3));
}
	
forAll(Y,i)
{
     volScalarField& Xi = X[i];
     volScalarField& Mui = Mu[i];
     yhi *= 0;
     forAll (Y,z)
	 {
		volScalarField& Xz = X[z];
		volScalarField& Muz = Mu[z];
		yhi += Xz*0.353553*pow(1.0+Foam::sqrt(Mui/Muz)*pow025(composition.Wi(z)/composition.Wi(i)),2)/Foam::sqrt(1.0+composition.Wi(i)/composition.Wi(z));
	 }
     mixMu += (Xi)*Mui/yhi;
	 mixLambdaA += Xi*Lambda[i];
     mixLambdaB += Xi/Lambda[i];
}

mixAlpha = 0.5*(mixLambdaA+(1/mixLambdaB))/thermo.Cp();
// for soret effect, thermal diffusion for light species
if(thermalDiffusion)
{
	TDRatio_H *= 0;
	TDRatio_H2 *= 0;
	forAll(Y,i)
	{
		if(Y[i].name() == "H")
		{
			forAll(Y,j)
			{
				if(j != i)
				{
					TDRatio_H += X[i]*X[j]*(TDiff1[j]*pow(T/TUnit,0) + TDiff2[j]*pow(T/TUnit,1) + TDiff3[j]*pow(T/TUnit,2) + TDiff4[j]*pow(T/TUnit,3));
				}
			}		
			//DiffError += 	rho*Diff[i]*TDRatio_H*Y[i]*fvc::grad(T)/(T*(X[i]+SMALL));	
		}
		
		if(Y[i].name() == "H2")
		{
			forAll(Y,j)
			{
				if(j != i)
				{
					TDRatio_H2 += X[i]*X[j]*(TDiff1[j+Y.size()]*pow(T/TUnit,0) + TDiff2[j+Y.size()]*pow(T/TUnit,1) + TDiff3[j+Y.size()]*pow(T/TUnit,2) + TDiff4[j+Y.size()]*pow(T/TUnit,3));
				}
			}
			//DiffError += 	rho*Diff[i]*TDRatio_H2*Y[i]*fvc::grad(T)/(T*(X[i]+SMALL));				
		}		
	}
}


if(differentialDiffusion)
{
	if(thermalDiffusion)
	{
		Info << " WITH differential diffusion and  thermalDiffusion"  << endl;
	}
	else
	{
		Info << " WITH differential diffusion"  << endl;
	}
}
else
{
    Info << " WITHOUT differential diffusion "  << endl;
}

Info << "update transport properties" << endl;
Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;
//Info << "code by zsh" << endl;
