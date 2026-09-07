/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Direct flamelet-table reaction-rate model for the FSD combustion model.

    The table contains pairs:

        K [s^-1]     omega0_H2 [kg m^-2 s^-1]

    The local FSD strain-rate sigma is used directly as the table abscissa.
    Linear interpolation is used between the actual Cantera continuation
    points produced by analyze_flamelet_v9_openfoam_fixed.py.
\*---------------------------------------------------------------------------*/

#ifndef flameletTable_H
#define flameletTable_H

#include "reactionRateFlameArea.H"
#include "List.H"
#include "word.H"

namespace Foam
{
namespace reactionRateFlameAreaModels
{

/*---------------------------------------------------------------------------*\\
                         Class flameletTable Declaration
\*---------------------------------------------------------------------------*/

class flameletTable
:
    public reactionRateFlameArea
{
    // Private data

        //- Flamelet strain-rate abscissa [s^-1]
        List<scalar> K_;

        //- Flamelet H2 consumption rate per unit flame area [kg m^-2 s^-1]
        List<scalar> omegaTable_;

        //- Minimum allowed omega
        scalar omegaMin_;

        //- Behaviour below the first tabulated K: clamp or zero
        word belowMin_;

        //- Behaviour above the last tabulated K: clamp or zero
        word aboveMax_;

        //- File read from the constant directory
        word tableFile_;


    // Private member functions

        //- Read the external table
        void readTable();

        //- Linear interpolation of omega(K)
        scalar interpolate(const scalar K) const;

        //- No copy construct
        flameletTable(const flameletTable&) = delete;

        //- No copy assignment
        void operator=(const flameletTable&) = delete;


public:

    //- Runtime type information
    TypeName("flameletTable");


    // Constructors

        //- Construct from components
        flameletTable
        (
            const word modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
        );


    // Destructor

        virtual ~flameletTable();


    // Member functions

        //- Update omega from the local FSD strain-rate field
        virtual void correct(const volScalarField& sigma) override;

        //- Update model properties from dictionary
        virtual bool read(const dictionary& dict) override;
};

} // End namespace reactionRateFlameAreaModels
} // End namespace Foam

#endif

// ************************************************************************* //
