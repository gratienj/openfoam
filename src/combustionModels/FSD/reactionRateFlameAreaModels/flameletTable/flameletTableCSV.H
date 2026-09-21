/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Direct lookup/interpolation of a CSV flamelet omega(K) table.

    The CSV file is read from:

        <case>/constant/<tableFile>

    This is compatible with parallel runs since the path is resolved
    from the global case directory rather than processorN/constant/.
\*---------------------------------------------------------------------------*/

#ifndef flameletTableCSV_H
#define flameletTableCSV_H

#include "reactionRateFlameArea.H"
#include "scalarList.H"
#include "word.H"
#include "fileName.H"

#include <string>
#include <vector>


namespace Foam
{
namespace reactionRateFlameAreaModels
{

class flameletTableCSV
:
    public reactionRateFlameArea
{
    // Private data

        //- Strain-rate / flamelet coordinate
        scalarList K_;

        //- Flame consumption rate corresponding to K_
        scalarList omegaTable_;

        //- Minimum allowed omega
        scalar omegaMin_;

        //- Maximum allowed omega
        scalar omegaMax_;

        //- Behaviour below minimum K
        word belowMin_;

        //- Behaviour above maximum K
        word aboveMax_;

        //- CSV file name, relative to <case>/constant/
        fileName tableFile_;

        //- Name of the K column in CSV
        word KColumn_;

        //- Name of the omega column in CSV
        word omegaColumn_;


    // Private Member Functions

        //- Read the CSV flamelet table
        void readTable();

        //- Interpolate omega(K)
        scalar interpolate(const scalar K) const;

        //- Split one CSV line into fields
        static std::vector<std::string> splitCSV
        (
            const std::string& line
        );

        //- Remove leading/trailing spaces
        static std::string trim
        (
            const std::string& value
        );

        //- Convert string to scalar
        static scalar readScalar
        (
            const std::string& value,
            const label lineNumber,
            const word& columnName
        );

        //- Find a column in the CSV header
        static label findColumn
        (
            const std::vector<std::string>& header,
            const word& columnName
        );


public:

    //- Runtime type information
    TypeName("flameletTableCSV");


    // Constructors

        //- Construct from dictionary
        flameletTableCSV
        (
            const word modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
        );


    //- Destructor
    virtual ~flameletTableCSV();


    // Member Functions

        //- Correct reaction-rate flame area correlation
        virtual void correct
        (
            const volScalarField& sigma
        );

        //- Read dictionary
        virtual bool read
        (
            const dictionary& dict
        );
};


} // End namespace reactionRateFlameAreaModels
} // End namespace Foam

#endif
