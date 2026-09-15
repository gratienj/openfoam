/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Direct lookup/interpolation of a Cantera CSV flamelet table.

    Expected CSV columns:
        K
        omega0_H2

    Other CSV columns are allowed and ignored.

    Example:
        K,Uin,mdot,omega0_H2,Sc_H2,Sc_heat,Tmax,...
        177.285279,2.363..., ...,7.03596e-03,...
        ...
\*---------------------------------------------------------------------------*/

#ifndef flameletTableCSV_H
#define flameletTableCSV_H

#include "reactionRateFlameArea.H"
#include "scalarList.H"
#include "word.H"
#include "fileName.H"

namespace Foam
{
namespace reactionRateFlameAreaModels
{

class flameletTableCSV
:
    public reactionRateFlameArea
{
    // Private data

        //- Flamelet abscissa: K [s^-1]
        scalarList K_;

        //- H2 consumption rate: omega0_H2 [kg/m2/s]
        scalarList omegaTable_;

        //- Minimum allowed omega
        scalar omegaMin_;

        //- Behaviour below first table value
        word belowMin_;

        //- Behaviour above last table value
        word aboveMax_;

        //- CSV file name
        fileName tableFile_;

        //- CSV column names
        word KColumn_;
        word omegaColumn_;


    // Private Member Functions

        //- Read CSV flamelet table
        void readTable();

        //- Interpolate omega(K)
        scalar interpolate(const scalar K) const;

        //- Trim whitespace from a string
        static std::string trim(const std::string& s);

        //- Split one CSV line
        static std::vector<std::string> splitCSV(const std::string& line);

        //- Find a column in CSV header
        static label findColumn
        (
            const std::vector<std::string>& header,
            const word& columnName
        );


public:

    //- Runtime type information
    TypeName("flameletTableCSV");


    // Constructors

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

        //- Correct reaction rate flame area
        virtual void correct(const volScalarField& sigma);

        //- Read model dictionary
        virtual bool read(const dictionary& dict);
};

} // End namespace reactionRateFlameAreaModels
} // End namespace Foam

#endif

// ************************************************************************* //