/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Direct lookup/interpolation of a CSV flamelet omega(K) table.

    CSV file:

        <case>/constant/H2_FSD_V11_flamelet.csv

    Expected columns by default:

        K
        omega0_H2

    Other columns are ignored.
\*---------------------------------------------------------------------------*/

#include "flameletTableCSV.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cctype>

#include "Pstream.H"
#include "mathematicalConstants.H"

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cerrno>
#include <cmath>


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


// * * * * * * * * * * * * * * Private functions * * * * * * * * * * * * * //


std::string
Foam::reactionRateFlameAreaModels::flameletTableCSV::trim
(
    const std::string& value
)
{
    const std::string whitespace(" \t\r\n");

    const std::string::size_type first =
        value.find_first_not_of(whitespace);

    if (first == std::string::npos)
    {
        return std::string();
    }

    const std::string::size_type last =
        value.find_last_not_of(whitespace);

    return value.substr(first, last - first + 1);
}


std::vector<std::string>
Foam::reactionRateFlameAreaModels::flameletTableCSV::splitCSV
(
    const std::string& line
)
{
    std::vector<std::string> fields;

    std::string field;
    bool insideQuotes = false;

    for (std::string::size_type i = 0; i < line.size(); ++i)
    {
        const char c = line[i];

        if (c == '"')
        {
            insideQuotes = !insideQuotes;
        }
        else if (c == ',' && !insideQuotes)
        {
            fields.push_back(trim(field));
            field.clear();
        }
        else
        {
            field += c;
        }
    }

    fields.push_back(trim(field));

    return fields;
}


Foam::label
Foam::reactionRateFlameAreaModels::flameletTableCSV::findColumn
(
    const std::vector<std::string>& header,
    const word& columnName
)
{
    const std::string target(columnName.c_str());

    for
    (
        label i = 0;
        i < static_cast<label>(header.size());
        ++i
    )
    {
        if (trim(header[i]) == target)
        {
            return i;
        }
    }

    return -1;
}


Foam::scalar
Foam::reactionRateFlameAreaModels::flameletTableCSV::readScalar
(
    const std::string& value,
    const label lineNumber,
    const word& columnName
)
{
    const std::string cleaned = trim(value);

    if (cleaned.empty())
    {
        FatalErrorInFunction
            << "Empty value in column '" << columnName
            << "' at CSV line " << lineNumber
            << exit(FatalError);
    }

    char* endPtr = nullptr;

    errno = 0;

    const double result =
        std::strtod(cleaned.c_str(), &endPtr);

    if
    (
        endPtr == cleaned.c_str()
     || *endPtr != '\0'
     || errno == ERANGE
     || !std::isfinite(result)
    )
    {
        FatalErrorInFunction
            << "Invalid numerical value '" << cleaned
            << "' in column '" << columnName
            << "' at CSV line " << lineNumber
            << exit(FatalError);
    }

    return scalar(result);
}


// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //


Foam::reactionRateFlameAreaModels::flameletTableCSV::flameletTableCSV
(
    const word modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const combustionModel& combModel
)
:
    reactionRateFlameArea
    (
        modelType,
        dict,
        mesh,
        combModel
    ),

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


Foam::reactionRateFlameAreaModels::flameletTableCSV::~flameletTableCSV()
{}


// * * * * * * * * * * * * * Read CSV table * * * * * * * * * * * * * * * //


void Foam::reactionRateFlameAreaModels::flameletTableCSV::readTable()
{
    const fileName csvPath =
        mesh_.time().globalPath() / "constant" / tableFile_;

    Info<< "flameletTableCSV: reading table from \""
        << csvPath << "\"" << nl
        << "    K column       = " << KColumn_ << nl
        << "    omega column   = " << omegaColumn_ << nl;

    std::ifstream csv(csvPath.c_str());

    if (!csv.good())
    {
        FatalErrorInFunction
            << "Cannot open flamelet CSV file: "
            << csvPath << nl
            << exit(FatalError);
    }

    std::string line;

    // ---------------------------------------------------------------------
    // Read header
    // ---------------------------------------------------------------------

    if (!std::getline(csv, line))
    {
        FatalErrorInFunction
            << "Cannot read header from CSV file: "
            << csvPath << nl
            << exit(FatalError);
    }

    // Remove possible UTF-8 BOM
    if (line.size() >= 3
     && static_cast<unsigned char>(line[0]) == 0xEF
     && static_cast<unsigned char>(line[1]) == 0xBB
     && static_cast<unsigned char>(line[2]) == 0xBF)
    {
        line.erase(0, 3);
    }

    // ---------------------------------------------------------------------
    // Split header
    // ---------------------------------------------------------------------

    std::vector<std::string> header;

    {
        std::stringstream ss(line);
        std::string field;

        while (std::getline(ss, field, ','))
        {
            // trim
            while (!field.empty()
                && std::isspace(static_cast<unsigned char>(field.front())))
            {
                field.erase(field.begin());
            }

            while (!field.empty()
                && std::isspace(static_cast<unsigned char>(field.back())))
            {
                field.pop_back();
            }

            header.push_back(field);
        }
    }

    label KIndex = -1;
    label omegaIndex = -1;

    for (label i = 0; i < static_cast<label>(header.size()); ++i)
    {
        if (header[i] == KColumn_)
        {
            KIndex = i;
        }

        if (header[i] == omegaColumn_)
        {
            omegaIndex = i;
        }
    }

    if (KIndex < 0)
    {
        FatalErrorInFunction
            << "Column \"" << KColumn_
            << "\" not found in CSV header." << nl
            << exit(FatalError);
    }

    if (omegaIndex < 0)
    {
        FatalErrorInFunction
            << "Column \"" << omegaColumn_
            << "\" not found in CSV header." << nl
            << exit(FatalError);
    }

    Info<< "    K index         = " << KIndex << nl
        << "    omega index     = " << omegaIndex << nl;

    // ---------------------------------------------------------------------
    // Read data
    // ---------------------------------------------------------------------

    std::vector<scalar> Kvalues;
    std::vector<scalar> omegaValues;

    label lineNumber = 1;

    while (std::getline(csv, line))
    {
        ++lineNumber;

        // Skip empty lines
        if (line.empty())
        {
            continue;
        }

        std::stringstream ss(line);
        std::string field;

        std::vector<std::string> fields;

        while (std::getline(ss, field, ','))
        {
            fields.push_back(field);
        }

        const label nFields = fields.size();

        if (KIndex >= nFields || omegaIndex >= nFields)
        {
            WarningInFunction
                << "Skipping malformed CSV line "
                << lineNumber << ": " << line << nl;

            continue;
        }

        try
        {
            const scalar K =
                std::stod(fields[KIndex]);

            const scalar omega =
                std::stod(fields[omegaIndex]);

            if (!std::isfinite(K) || !std::isfinite(omega))
            {
                WarningInFunction
                    << "Skipping non-finite CSV line "
                    << lineNumber << ": " << line << nl;

                continue;
            }

            Kvalues.push_back(K);
            omegaValues.push_back(omega);
        }
        catch (...)
        {
            WarningInFunction
                << "Skipping non-numeric CSV line "
                << lineNumber << ": " << line << nl;

            continue;
        }
    }

    csv.close();

    // ---------------------------------------------------------------------
    // Check table
    // ---------------------------------------------------------------------

    if (Kvalues.size() < 2)
    {
        FatalErrorInFunction
            << "CSV table contains fewer than 2 valid data points." << nl
            << "File: " << csvPath << nl
            << exit(FatalError);
    }

    for (label i = 1; i < static_cast<label>(Kvalues.size()); ++i)
    {
        if (Kvalues[i] <= Kvalues[i-1])
        {
            FatalErrorInFunction
                << "K column is not strictly increasing." << nl
                << "At CSV line " << i + 2 << ":" << nl
                << "    K[" << i-1 << "] = " << Kvalues[i-1] << nl
                << "    K[" << i   << "] = " << Kvalues[i] << nl
                << exit(FatalError);
        }
    }

    // ---------------------------------------------------------------------
    // Transfer to OpenFOAM containers
    // ---------------------------------------------------------------------

    K_.setSize(Kvalues.size());
    omegaTable_.setSize(omegaValues.size());

    for (label i = 0; i < static_cast<label>(Kvalues.size()); ++i)
    {
        K_[i] = Kvalues[i];
        omegaTable_[i] = omegaValues[i];
    }

    omegaMin_ = gMin(omegaTable_);
    omegaMax_ = gMax(omegaTable_);

    Info<< "    Number of points = " << K_.size() << nl
        << "    K min            = " << K_.first() << nl
        << "    K max            = " << K_.last() << nl
        << "    omega min        = " << omegaMin_ << nl
        << "    omega max        = " << omegaMax_ << nl
        << endl;
}


// * * * * * * * * * * * * * Interpolation * * * * * * * * * * * * * * * //


Foam::scalar
Foam::reactionRateFlameAreaModels::flameletTableCSV::interpolate
(
    const scalar K
) const
{
    // ---------------------------------------------------------------------
    // Below table
    // ---------------------------------------------------------------------

    if (K <= K_.first())
    {
        if (belowMin_ == "zero")
        {
            return 0.0;
        }

        return max(omegaTable_.first(), omegaMin_);
    }


    // ---------------------------------------------------------------------
    // Above table
    // ---------------------------------------------------------------------

    if (K >= K_.last())
    {
        if (aboveMax_ == "zero")
        {
            return 0.0;
        }

        return max(omegaTable_.last(), omegaMin_);
    }


    // ---------------------------------------------------------------------
    // Binary search
    // ---------------------------------------------------------------------

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


// * * * * * * * * * * * * * Correct * * * * * * * * * * * * * * * * * * //


void
Foam::reactionRateFlameAreaModels::flameletTableCSV::correct
(
    const volScalarField& sigma
)
{
    volScalarField::Internal& iOmega = omega_;


    // ---------------------------------------------------------------------
    // Internal cells
    // ---------------------------------------------------------------------

    forAll(iOmega, celli)
    {
        const scalar K =
            max
            (
                sigma[celli],
                scalar(0)
            );


        iOmega[celli] =
            interpolate(K);
    }


    // ---------------------------------------------------------------------
    // Boundary faces
    // ---------------------------------------------------------------------

    volScalarField::Boundary& bOmega =
        omega_.boundaryFieldRef();


    forAll(bOmega, patchi)
    {
        forAll(bOmega[patchi], facei)
        {
            const scalar K =
                max
                (
                    sigma.boundaryField()[patchi][facei],
                    scalar(0)
                );


            bOmega[patchi][facei] =
                interpolate(K);
        }
    }
}


// * * * * * * * * * * * * * Read dictionary * * * * * * * * * * * * * * //


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
