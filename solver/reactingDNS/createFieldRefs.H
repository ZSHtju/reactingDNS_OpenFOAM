const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
const label inertIndex(composition.species()[inertSpecie]);
//BasicChemistryModel<psiReactionThermo>& chemistry = pChemistry();
//volScalarField& h = thermo.he();
//scalar dtChem = refCast<const psiChemistryModel>(chemistry).deltaTChem()[0];

volScalarField& psiToWrite(const_cast<volScalarField&>(thermo.psi()));
volScalarField& kappaToWrite(const_cast<volScalarField&>(thermo.kappa()));
volScalarField& muToWrite(const_cast<volScalarField&>(thermo.mu()()));

psiToWrite.writeOpt() = IOobject::AUTO_WRITE;
kappaToWrite.writeOpt() = IOobject::AUTO_WRITE;
muToWrite.writeOpt() = IOobject::AUTO_WRITE;