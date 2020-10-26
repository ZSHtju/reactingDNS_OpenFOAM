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
	//surfaceScalarField  phiUc = linearInterpolate(DiffError) & mesh.Sf();
    for(label i=0; i < Y.size(); i++)
    {
		volScalarField& Yi = Y[i];
       
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
   	//Hsi.correctBoundaryConditions();
           vk += Hsi*mixAlpha*fvc::grad(Yi); //vk += -rho*Hai*(Diffi+Lei*mixAlpha/rho)*(fvc::grad(Yi));
           diffsh += fvc::laplacian(mixAlpha*Hsi,Yi);
	   //volScalarField& Ymean = Y[i];
	    fvScalarMatrix YEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(mixAlpha, Yi)   //fvm::laplacian(rho*(Diffi+Lei*mixAlpha/rho), Yi)             
              ==
                chemistry.RR(i)
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
    chemistryhsSource = chemistry.Qdot()();

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
      - fvm::laplacian(mixAlpha, he)
      + diffsh
     ==
        chemistry.Qdot()()
      + fvc::div(vk)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

if(boundTemperature)
{
	forAll(rho, icell)
	{ 
		gasHThermoPhysics mixture(0.0*specieData[0]);
		for(label kk=0; kk < Y.size(); kk++)
		{
			mixture += Y[kk][icell]*specieData[kk];
		}
		he[icell]   = max(mixture.Hs(p[icell],Tmin.value()),he[icell]);
		he[icell]   = min(mixture.Hs(p[icell],Tmax.value()),he[icell]);
	}


	forAll(he.boundaryField(), patchi)
	{
		forAll(he.boundaryField()[patchi],facei)
		{
			gasHThermoPhysics mixture(0.0*specieData[0]);
			for(label kk=0; kk < Y.size(); kk++)
			{
				mixture += Y[kk].boundaryField()[patchi][facei]*specieData[kk];
			}
			scalar pi = p.boundaryField()[patchi][facei];
			scalar& hei = he.boundaryFieldRef()[patchi][facei]; //250717
			he.boundaryFieldRef()[patchi][facei]   = max(mixture.Hs(pi,Tmin.value()),hei);
			he.boundaryFieldRef()[patchi][facei]   = min(mixture.Hs(pi,Tmax.value()),hei);

		}
	}

}

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;


Info << " chemistryhsSource max/min: "<<  max(chemistryhsSource).value() <<" " << min(chemistryhsSource).value() << endl;
solverPerformance::debug = tmpDebug; 
