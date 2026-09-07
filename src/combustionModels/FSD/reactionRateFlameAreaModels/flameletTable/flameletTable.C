/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Direct lookup/interpolation of the Cantera flamelet omega(K) table.
\*---------------------------------------------------------------------------*/

#include "flameletTable.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "Tuple2.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace reactionRateFlameAreaModels
{
    defineTypeNameAndDebug(flameletTable, 0);
    addToRunTimeSelectionTable
    (
        reactionRateFlameArea,
        flameletTable,
        dictionary
    );
}
}


Foam::reactionRateFlameAreaModels::flameletTable::flameletTable
(
    const word modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
:
    reactionRateFlameArea(modelType, dict, mesh, combModel),
    K_(),
    omegaTable_(),
    omegaMin_(coeffDict_.getOrDefault<scalar>("omegaMin", 0.0)),
    belowMin_(coeffDict_.getOrDefault<word>("belowMin", "clamp")),
    aboveMax_(coeffDict_.getOrDefault<word>("aboveMax", "clamp")),
    tableFile_(coeffDict_.getOrDefault<word>("tableFile", "FSD_H2_flameletTable"))
{
    readTable();
}


Foam::reactionRateFlameAreaModels::flameletTable::~flameletTable()
{}


void Foam::reactionRateFlameAreaModels::flameletTable::readTable()
{
    IOdictionary tableDict
    (
        IOobject
        (
            tableFile_,
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const List<Tuple2<scalar, scalar>> tableData
    (
        tableDict.lookup("table")
    );

    if (tableData.size() < 2)
    {
        FatalIOErrorInFunction(tableDict)
            << "The flamelet table '" << tableFile_
            << "' must contain at least two (K omega) entries."
            << exit(FatalIOError);
    }

    K_.setSize(tableData.size());
    omegaTable_.setSize(tableData.size());

    forAll(tableData, i)
    {
        K_[i] = tableData[i].first();
        omegaTable_[i] = max(tableData[i].second(), omegaMin_);

        if (K_[i] <= 0.0)
        {
            FatalIOErrorInFunction(tableDict)
                << "Invalid flamelet-table K at index " << i
                << ": K = " << K_[i] << " s^-1."
                << exit(FatalIOError);
        }

        if (i > 0 && K_[i] <= K_[i - 1])
        {
            FatalIOErrorInFunction(tableDict)
                << "Flamelet-table K values must be strictly increasing. "
                << "Found K[" << i - 1 << "] = " << K_[i - 1]
                << " and K[" << i << "] = " << K_[i] << "."
                << exit(FatalIOError);
        }
    }

    if
    (
        belowMin_ != "clamp"
     && belowMin_ != "zero"
    )
    {
        FatalIOErrorInFunction(tableDict)
            << "belowMin must be 'clamp' or 'zero', but is '"
            << belowMin_ << "'."
            << exit(FatalIOError);
    }

    if
    (
        aboveMax_ != "clamp"
     && aboveMax_ != "zero"
    )
    {
        FatalIOErrorInFunction(tableDict)
            << "aboveMax must be 'clamp' or 'zero', but is '"
            << aboveMax_ << "'."
            << exit(FatalIOError);
    }

    Info<< "flameletTable: read " << K_.size()
        << " Cantera omega(K) points from " << tableFile_ << nl
        << "    K range       = " << K_.first() << " .. " << K_.last()
        << " s^-1" << nl
        << "    omega range  = " << min(omegaTable_)
        << " .. " << max(omegaTable_)
        << " kg/m2/s" << nl
        << "    belowMin     = " << belowMin_ << nl
        << "    aboveMax     = " << aboveMax_ << nl;
}


Foam::scalar
Foam::reactionRateFlameAreaModels::flameletTable::interpolate
(
    const scalar K
) const
{
    if (K <= K_.first())
    {
        return
        (
            belowMin_ == "zero"
          ? 0.0
          : omegaTable_.first()
        );
    }

    if (K >= K_.last())
    {
        return
        (
            aboveMax_ == "zero"
          ? 0.0
          : omegaTable_.last()
        );
    }

    label lo = 0;
    label hi = K_.size() - 1;

    while (hi - lo > 1)
    {
        const label mid = (lo + hi)/2;

        if (K_[mid] <= K)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    const scalar dK = K_[hi] - K_[lo];
    const scalar f = (K - K_[lo])/dK;

    return max
    (
        omegaTable_[lo]
      + f*(omegaTable_[hi] - omegaTable_[lo]),
        omegaMin_
    );
}


void Foam::reactionRateFlameAreaModels::flameletTable::correct
(
    const volScalarField& sigma
)
{
    volScalarField::Internal& iOmega = omega_;

    forAll(iOmega, celli)
    {
        iOmega[celli] = interpolate(max(sigma[celli], scalar(0)));
    }

    volScalarField::Boundary& bOmega = omega_.boundaryFieldRef();

    forAll(bOmega, patchi)
    {
        forAll(bOmega[patchi], facei)
        {
            bOmega[patchi][facei] = interpolate
            (
                max
                (
                    sigma.boundaryField()[patchi][facei],
                    scalar(0)
                )
            );
        }
    }
}


bool Foam::reactionRateFlameAreaModels::flameletTable::read
(
    const dictionary& dict
)
{
    if (!reactionRateFlameArea::read(dict))
    {
        return false;
    }

    coeffDict_ = dict.optionalSubDict
    (
        typeName + "Coeffs"
    );

    coeffDict_.readIfPresent("omegaMin", omegaMin_);
    coeffDict_.readIfPresent("belowMin", belowMin_);
    coeffDict_.readIfPresent("aboveMax", aboveMax_);
    coeffDict_.readIfPresent("tableFile", tableFile_);

    readTable();

    return true;
}


// ************************************************************************* //

