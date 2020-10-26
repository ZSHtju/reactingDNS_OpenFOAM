Info<< "Reading chemistry properties\n" << endl;

IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

//Switch prescibed_Le(chemistryProperties.lookup("prescibed_Le"));
//Switch writeAvgPressure(chemistryProperties.lookup("writeAvgPressure"));
Switch differentialDiffusion(chemistryProperties.lookup("differentialDiffusion"));
Switch thermalDiffusion(chemistryProperties.lookup("thermalDiffusion"));
