#include "createRDeltaT.H"

Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<fluidReactionThermo> pThermo(fluidReactionThermo::New(mesh));
fluidReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

/* autoPtr< BasicChemistryModel<psiReactionThermo>> pChemistry
(
  BasicChemistryModel<psiReactionThermo>::New(thermo)
);
 */
basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const dictionary thermoDict(thermo.properties());
const word inertSpecie(thermoDict.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermoDict)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureReference pressureReference(p,rho , pimple.dict(), false);

mesh.schemes().setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidReactionThermophysicalTransportModel>thermophysicalTransport
(
    fluidReactionThermophysicalTransportModel::New(turbulence(), thermo)
);


Info<< "Creating reaction model\n" << endl;
autoPtr<combustionModel> reaction
(
    combustionModel::New(thermo, turbulence())
);

dimensioned<scalar> unburntRho(dimDensity,1.0);



Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);


Info << "Creating displacement speed\n" << endl;
volScalarField displacementSpeed
(
    IOobject
    (
        "displacementSpeed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVelocity, 0)
);

Info << "Creating field reactionSourceTerm\n" << endl;
volScalarField reactionSourceTerm
(
    IOobject
    (
        "reactionSourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimMass/dimVolume/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

// for print reaction rates for each species
bool writeRR = runTime.controlDict().lookupOrDefault<bool>("writeRR",false);
Info << "writeRR: " << writeRR << endl;

PtrList<volScalarField> RYi(Y.size());

// 
if (writeRR)
{
    forAll(RYi,i)
    {
        RYi.set
        (
            i,
            new volScalarField 
            (
                IOobject
                (
                    "R_" + Y[i].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless, 0)
            )
        );
    }
}
    // RYii.field() = -reaction->R(Y[i])->source();


//thermo.correct();

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"