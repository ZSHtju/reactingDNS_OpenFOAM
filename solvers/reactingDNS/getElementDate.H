label  nSpecies = composition.Y().size();
List<List<specieElement> > specieComp(Y.size());
speciesTable speciesTab = thermo.composition().species();
autoPtr<chemistryReader<gasHThermoPhysics>> tchemRead
(
    chemistryReader<gasHThermoPhysics>::New(thermo,speciesTab)
);

const HashTable<List<specieElement> >& specComp(tchemRead().specieComposition());
forAll(specieComp,i)
{
    specieComp[i] = specComp[Y[i].name()];

}

wordList elementNames(4);

{
    elementNames[0] = "H";
    elementNames[1] = "C";
    elementNames[2] = "O";
    elementNames[3] = "N";
}

List<scalarField> spInfo(nSpecies);
forAll(spInfo,i)  spInfo[i].setSize(8,0);

forAll(Y, i)
{
    Info << "specie = " << Y[i].name() << "\t";
    forAll(elementNames, en)
    {
            forAll(specieComp[i], kk)
            {
                if(specieComp[i][kk].name() == elementNames[en])
                {
                    Info << specieComp[i][kk].name() << ", nAtoms = " << specieComp[i][kk].nAtoms() << "\t";
                    spInfo[i][en] = specieComp[i][kk].nAtoms();
                }
            }

    }

    Info << specieData[i].W() << endl;
    spInfo[i][4] = specieData[i].W();
}
