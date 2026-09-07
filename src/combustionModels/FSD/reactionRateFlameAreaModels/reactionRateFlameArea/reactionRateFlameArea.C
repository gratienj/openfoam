/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "reactionRateFlameArea.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(reactionRateFlameArea, 0);
    defineRunTimeSelectionTable(reactionRateFlameArea, dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::reactionRateFlameArea::reactionRateFlameArea
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
:
    coeffDict_
    (
        dict.optionalSubDict(modelType + "Coeffs")
    ),

    mesh_(mesh),

    combModel_(combModel),

    fuel_
    (
        dict.lookupOrDefault<word>("fuel", "H2")
    ),

    omega_
    (
        IOobject
        (
            "FSDomega",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "FSDomega",
            dimensionSet(1, -2, -1, 0, 0, 0, 0),
            Zero
        )
    )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::reactionRateFlameArea>
Foam::reactionRateFlameArea::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
{
    const word modelType
    (
        dict.lookup("reactionRateFlameArea")
    );

    Info<< "Selecting reaction rate flame area correlation "
        << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "reactionRateFlameArea",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    const word className =
        modelType.substr(0, modelType.find('<'));

    return autoPtr<reactionRateFlameArea>
    (
        ctorPtr(className, dict, mesh, combModel)
    );
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::reactionRateFlameArea::~reactionRateFlameArea()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::reactionRateFlameArea::read(const dictionary& dict)
{
    fuel_ = dict.lookupOrDefault<word>("fuel", "H2");

    return true;
}


// ************************************************************************* //
