tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);
label tmpDebug = solverPerformance::debug; 
solverPerformance::debug = 0; 
vk *= 0;
diffsh *= 0;
// step : advection-diffusoin
	//label inertIndex = -1;
	volScalarField Yt(0.0*Y[0]);
    surfaceScalarField  phiUc = linearInterpolate(DiffError) & mesh.Sf();
    for(label i=0; i < Y.size(); i++)
    {
		volScalarField& Yi = Y[i];
		volScalarField& Xi = X[i];
		volScalarField& Diffi = Diff[i];
		volScalarField  ThermalDiff = TDRatio_H*0;
		if(Y[i].name() == "H")
		{
		   ThermalDiff = TDRatio_H;
		}
		else if(Y[i].name() == "H2")
		{
		   ThermalDiff = TDRatio_H2;
		}		   
		else
		{
		   ThermalDiff *= 0;
		}			   
        
	   for(label icell = 0; icell < rho.size(); icell++)
	   {
		 Hsi[icell] = composition.Hs(i, p[icell], T[icell]);
	   }

		volScalarField::Boundary& HsiBf = Hsi.boundaryFieldRef();

		forAll(HsiBf, patchi)
		{
			scalarField& HsiPatch = HsiBf[patchi];
			const scalarField& pp = p.boundaryField()[patchi];
			const scalarField& Tp = T.boundaryField()[patchi];

			forAll(HsiPatch, facei)
			{
				HsiPatch[facei] = composition.Hs(i, pp[facei], Tp[facei]);
			}
		}
		vk += Hsi*(rho*Diffi*fvc::grad(Yi)+Yi*rho*Diffi*fvc::grad(MW)/MW-Yi*DiffError); //vk += -rho*Hai*(Diffi+Lei*mixAlpha/rho)*(fvc::grad(Yi));
		diffsh += fvc::laplacian(mixAlpha*Hsi,Yi);// override temperature to enthalpy
   
		fvScalarMatrix YEqn
		(
			fvm::ddt(rho, Yi)
		  + mvConvection->fvmDiv(phi, Yi)
		  - fvm::laplacian(rho*Diffi+turbulence->mut()/Sct, Yi)   //fvm::laplacian(rho*(Diffi+Lei*mixAlpha/rho), Yi)
		  + mvConvection->fvmDiv(phiUc,Yi)              
		  - fvc::laplacian(rho*Diffi*Yi/MW,MW)
		  - fvc::laplacian(rho*Diffi*ThermalDiff*Yi/(T*(Xi+SMALL)),T) // soret effect
		  ==
			reaction->R(Yi)
		  + fvOptions(rho, Yi)
		);
		YEqn.relax();

		fvOptions.constrain(YEqn);

		YEqn.solve("Yi");

		fvOptions.correct(Yi);

		Yi.max(0.0);
		Yt += Yi;
    }
// for N2			
    Y[inertIndex] = Y[inertIndex] + scalar(1) - Yt;
    Y[inertIndex].max(0.0);        
    chemistryhsSource = reaction->Qdot();

// Energy Equation	

    volScalarField& he = thermo.he();

	fvScalarMatrix EEqn
	(
		fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
	  + fvc::ddt(rho, K) + fvc::div(phi, K)
	  + (
			he.name() == "e"
		  ? fvc::div
			(
				fvc::absolute(phi/fvc::interpolate(rho), U),
				p,
				"div(phiv,p)"
			)
		  : -dpdt
		)
	  - fvm::laplacian(mixAlpha+turbulence->mut()/Prt, he)
	  + diffsh
	 ==
        reaction->Qdot()
	  + fvc::div(vk)
	  + fvOptions(rho, he)
	);

	EEqn.relax();

	fvOptions.constrain(EEqn);

	EEqn.solve();

	fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;


Info << " chemistryhsSource max/min: "<<  max(chemistryhsSource).value() <<" " << min(chemistryhsSource).value() << endl;
solverPerformance::debug = tmpDebug; 
