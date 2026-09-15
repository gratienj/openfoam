/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Direct lookup/interpolation of a Cantera CSV flamelet omega(K) table.
\*---------------------------------------------------------------------------*/

#include "flameletTableCSV.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "IOobject.H"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cstdlib>

namespace Foam
{
namespace reactionRateFlameAreaModels
{

    defineTypeNameAndDebug(flameletTableCSV, 0);

    addToRunTimeSelectionTable
    (
        reactionRateFlameArea,
        flameletTableCSV,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

std::string
Foam::reactionRateFlameAreaModels::flameletTableCSV::trim
(
    const std::string& s
)
{
    const std::string whitespace = " \t\r\n";

    const std::size_t first = s.find_first_not_of(whitespace);

    if (first == std::string::npos)
    {
        return "";
    }

    const std::size_t last = s.find_last_not_of(whitespace);

    return s.substr(first, last - first + 1);
}


// ************************************************************************* //

std::vector<std::string>
Foam::reactionRateFlameAreaModels::flameletTableCSV::splitCSV
(
    const std::string& line
)
{
    std::vector<std::string> fields;

    std::stringstream ss(line);
    std::string field;

    while (std::getline(ss, field, ','))
    {
        fields.push_back(trim(field));
    }

    // Handle a trailing comma
    if (!line.empty() && line.back() == ',')
    {
        fields.push_back("");
    }

    return fields;
}


// ************************************************************************* //

Foam::label
Foam::reactionRateFlameAreaModels::flameletTableCSV::findColumn
(
    const std::vector<std::string>& header,
    const word& columnName
)
{
    const std::string target = columnName.c_str();

    for
    (
        std::size_t i = 0;
        i < header.size();
        ++i
    )
    {
        if (header[i] == target)
        {
            return static_cast<label>(i);
        }
    }

    return -1;
}


// ************************************************************************* //

Foam::reactionRateFlameAreaModels::flameletTableCSV::flameletTableCSV
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

    omegaMin_
    (
        coeffDict_.getOrDefault<scalar>
        (
            "omegaMin",
            0.0
        )
    ),

    belowMin_
    (
        coeffDict_.getOrDefault<word>
        (
            "belowMin",
            "clamp"
        )
    ),

    aboveMax_
    (
        coeffDict_.getOrDefault<word>
        (
            "aboveMax",
            "clamp"
        )
    ),

    tableFile_
    (
        coeffDict_.getOrDefault<fileName>
        (
            "tableFile",
            "H2_FSD_V11_flamelet.csv"
        )
    ),

    KColumn_
    (
        coeffDict_.getOrDefault<word>
        (
            "KColumn",
            "K"
        )
    ),

    omegaColumn_
    (
        coeffDict_.getOrDefault<word>
        (
            "omegaColumn",
            "omega0_H2"
        )
    )
{
    readTable();
}


// ************************************************************************* //

Foam::reactionRateFlameAreaModels::flameletTableCSV::~flameletTableCSV()
{}


// * * * * * * * * * * * * * Table Reading  * * * * * * * * * * * * * * * //

void
Foam::reactionRateFlameAreaModels::flameletTableCSV::readTable()
{
    K_.clear();
    omegaTable_.clear();


    // ---------------------------------------------------------------------
    // Locate CSV in constant/
    // ---------------------------------------------------------------------

    IOobject tableIO
    (
        tableFile_,
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    const fileName csvPath = tableIO.objectPath();

    Info<< "flameletTableCSV: reading CSV file "
        << csvPath << nl;


    // ---------------------------------------------------------------------
    // Open CSV
    // ---------------------------------------------------------------------

    std::ifstream csv(csvPath.c_str());

    if (!csv.good())
    {
        FatalIOErrorInFunction(csvPath)
            << "Cannot open flamelet CSV file: "
            << csvPath
            << exit(FatalIOError);
    }


    // ---------------------------------------------------------------------
    // Read header
    // ---------------------------------------------------------------------

    std::string line;

    bool headerFound = false;

    std::vector<std::string> header;

    while (std::getline(csv, line))
    {
        line = trim(line);

        if (line.empty())
        {
            continue;
        }

        // Skip comments
        if (line[0] == '#')
        {
            continue;
        }

        header = splitCSV(line);
        headerFound = true;

        break;
    }

    if (!headerFound)
    {
        FatalIOErrorInFunction(csvPath)
            << "The flamelet CSV file '" << csvPath
            << "' contains no header."
            << exit(FatalIOError);
    }


    // ---------------------------------------------------------------------
    // Find requested columns
    // ---------------------------------------------------------------------

    const label KIndex =
        findColumn(header, KColumn_);

    const label omegaIndex =
        findColumn(header, omegaColumn_);


    if (KIndex < 0)
    {
        FatalIOErrorInFunction(csvPath)
            << "Column '" << KColumn_
            << "' was not found in flamelet CSV file '"
            << csvPath << "'."
            << nl
            << "Available columns:" << nl;

        forAll(header, i)
        {
            FatalIOErrorInFunction(csvPath)
                << "    [" << i << "] "
                << header[i] << nl;
        }

        exit(FatalIOError);
    }


    if (omegaIndex < 0)
    {
        FatalIOErrorInFunction(csvPath)
            << "Column '" << omegaColumn_
            << "' was not found in flamelet CSV file '"
            << csvPath << "'."
            << nl
            << "Available columns:" << nl;

        forAll(header, i)
        {
            FatalIOErrorInFunction(csvPath)
                << "    [" << i << "] "
                << header[i] << nl;
        }

        exit(FatalIOError);
    }


    Info<< "    K column       = " << KColumn_
        << " [" << KIndex << "]" << nl
        << "    omega column   = " << omegaColumn_
        << " [" << omegaIndex << "]" << nl;


    // ---------------------------------------------------------------------
    // Read data
    // ---------------------------------------------------------------------

    label lineNumber = 1;

    while (std::getline(csv, line))
    {
        ++lineNumber;

        line = trim(line);

        if (line.empty())
        {
            continue;
        }

        // Skip comments
        if (line[0] == '#')
        {
            continue;
        }

        const std::vector<std::string> fields = splitCSV(line);

        const label nFields =
            static_cast<label>(fields.size());

        if
        (
            KIndex >= nFields
         || omegaIndex >= nFields
        )
        {
            WarningInFunction
                << "Skipping CSV line " << lineNumber
                << ": not enough fields."
                << nl;

            continue;
        }


        const std::string& KString =
            fields[KIndex];

        const std::string& omegaString =
            fields[omegaIndex];


        // -------------------------------------------------------------
        // Ignore non-numeric / NaN rows
        // -------------------------------------------------------------

        if
        (
            KString.empty()
         || omegaString.empty()
         || KString == "nan"
         || KString == "NaN"
         || KString == "NAN"
         || omegaString == "nan"
         || omegaString == "NaN"
         || omegaString == "NAN"
        )
        {
            continue;
        }


        char* KEnd = nullptr;
        char* omegaEnd = nullptr;

        const scalar K =
            std::strtod
            (
                KString.c_str(),
                &KEnd
            );

        const scalar omega =
            std::strtod
            (
                omegaString.c_str(),
                &omegaEnd
            );


        // Check conversion
        if
        (
            KEnd == KString.c_str()
         || omegaEnd == omegaString.c_str()
        )
        {
            WarningInFunction
                << "Skipping non-numeric CSV line "
                << lineNumber << "."
                << nl;

            continue;
        }


        // -------------------------------------------------------------
        // Validate K
        // -------------------------------------------------------------

        if (K <= 0.0)
        {
            FatalIOErrorInFunction(csvPath)
                << "Invalid K at CSV line "
                << lineNumber
                << ": K = " << K
                << " s^-1."
                << exit(FatalIOError);
        }


        // -------------------------------------------------------------
        // Check monotonicity
        // -------------------------------------------------------------

        if
        (
            !K_.empty()
         && K <= K_.last()
        )
        {
            FatalIOErrorInFunction(csvPath)
                << "K values must be strictly increasing."
                << nl
                << "At CSV line " << lineNumber
                << ": previous K = " << K_.last()
                << ", current K = " << K
                << exit(FatalIOError);
        }


        K_.append(K);

        omegaTable_.append
        (
            max(omega, omegaMin_)
        );
    }


    csv.close();


    // ---------------------------------------------------------------------
    // Check table
    // ---------------------------------------------------------------------

    if (K_.size() < 2)
    {
        FatalIOErrorInFunction(csvPath)
            << "The flamelet CSV table '" << csvPath
            << "' must contain at least two valid "
            << "(K, omega0_H2) entries."
            << exit(FatalIOError);
    }


    // ---------------------------------------------------------------------
    // Check extrapolation behaviour
    // ---------------------------------------------------------------------

    if
    (
        belowMin_ != "clamp"
     && belowMin_ != "zero"
    )
    {
        FatalIOErrorInFunction(csvPath)
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
        FatalIOErrorInFunction(csvPath)
            << "aboveMax must be 'clamp' or 'zero', but is '"
            << aboveMax_ << "'."
            << exit(FatalIOError);
    }


    // ---------------------------------------------------------------------
    // Information
    // ---------------------------------------------------------------------

    Info<< "flameletTableCSV: successfully read "
        << K_.size()
        << " Cantera flamelet points from "
        << csvPath << nl
        << "    K range       = "
        << K_.first() << " .. "
        << K_.last()
        << " s^-1" << nl
        << "    omega range   = "
        << min(omegaTable_) << " .. "
        << max(omegaTable_)
        << " kg/m2/s" << nl
        << "    belowMin      = "
        << belowMin_ << nl
        << "    aboveMax      = "
        << aboveMax_ << nl;
}


// * * * * * * * * * * * * * Interpolation  * * * * * * * * * * * * * * //

Foam::scalar
Foam::reactionRateFlameAreaModels::flameletTableCSV::interpolate
(
    const scalar K
) const
{
    // Below table
    if (K <= K_.first())
    {
        return
        (
            belowMin_ == "zero"
          ? 0.0
          : omegaTable_.first()
        );
    }


    // Above table
    if (K >= K_.last())
    {
        return
        (
            aboveMax_ == "zero"
          ? 0.0
          : omegaTable_.last()
        );
    }


    // Binary search
    label lo = 0;
    label hi = K_.size() - 1;

    while (hi - lo > 1)
    {
        const label mid =
            (lo + hi)/2;

        if (K_[mid] <= K)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }


    const scalar dK =
        K_[hi] - K_[lo];

    const scalar f =
        (K - K_[lo])/dK;


    return max
    (
        omegaTable_[lo]
      + f*(omegaTable_[hi] - omegaTable_[lo]),
        omegaMin_
    );
}


// * * * * * * * * * * * * * Correct  * * * * * * * * * * * * * * * * * //

void
Foam::reactionRateFlameAreaModels::flameletTableCSV::correct
(
    const volScalarField& sigma
)
{
    volScalarField::Internal& iOmega = omega_;

    forAll(iOmega, celli)
    {
        iOmega[celli] =
            interpolate
            (
                max
                (
                    sigma[celli],
                    scalar(0)
                )
            );
    }


    volScalarField::Boundary& bOmega =
        omega_.boundaryFieldRef();


    forAll(bOmega, patchi)
    {
        forAll(bOmega[patchi], facei)
        {
            bOmega[patchi][facei] =
                interpolate
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


// * * * * * * * * * * * * * Read  * * * * * * * * * * * * * * * * * * * //

bool
Foam::reactionRateFlameAreaModels::flameletTableCSV::read
(
    const dictionary& dict
)
{
    if (!reactionRateFlameArea::read(dict))
    {
        return false;
    }


    coeffDict_ =
        dict.optionalSubDict
        (
            typeName + "Coeffs"
        );


    coeffDict_.readIfPresent
    (
        "omegaMin",
        omegaMin_
    );

    coeffDict_.readIfPresent
    (
        "belowMin",
        belowMin_
    );

    coeffDict_.readIfPresent
    (
        "aboveMax",
        aboveMax_
    );

    coeffDict_.readIfPresent
    (
        "tableFile",
        tableFile_
    );

    coeffDict_.readIfPresent
    (
        "KColumn",
        KColumn_
    );

    coeffDict_.readIfPresent
    (
        "omegaColumn",
        omegaColumn_
    );


    readTable();

    return true;
}


// ************************************************************************* //