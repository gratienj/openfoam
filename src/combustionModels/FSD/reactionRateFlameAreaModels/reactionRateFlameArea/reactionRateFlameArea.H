/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Abstract reaction-rate-per-flame-area model for the FSD combustion model.

    This version initializes FSDomega if it is not present at time 0 and
    contains the run-time selector implementation used by the derived models.
\*---------------------------------------------------------------------------*/

#ifndef reactionRateFlameArea_H
#define reactionRateFlameArea_H

#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "volFields.H"
#include "combustionModel.H"

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\\
                            Class reactionRateFlameArea Declaration
\*---------------------------------------------------------------------------*/

class reactionRateFlameArea
{
protected:

    // Protected data

        //- Dictionary containing the model coefficients
        dictionary coeffDict_;

        //- Mesh reference
        const fvMesh& mesh_;

        //- Combustion model owner
        const combustionModel& combModel_;

        //- Fuel name
        word fuel_;

        //- Fuel consumption rate per unit flame area [kg m^-2 s^-1]
        volScalarField omega_;


private:

    // Private member functions

        //- No copy construct
        reactionRateFlameArea(const reactionRateFlameArea&) = delete;

        //- No copy assignment
        void operator=(const reactionRateFlameArea&) = delete;


public:

    //- Runtime type information
    TypeName("reactionRateFlameArea");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            reactionRateFlameArea,
            dictionary,
            (
                const word modelType,
                const dictionary& dict,
                const fvMesh& mesh,
                const combustionModel& combModel
            ),
            (modelType, dict, mesh, combModel)
        );


    // Constructors

        //- Construct from components
        reactionRateFlameArea
        (
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
        );


    // Selector

        static autoPtr<reactionRateFlameArea> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const combustionModel& combModel
        );


    // Destructor

        virtual ~reactionRateFlameArea();


    // Member functions

        //- Return fuel consumption rate per unit flame area
        const volScalarField& omega() const
        {
            return omega_;
        }

        //- Correct omega from the local strain-rate field
        virtual void correct(const volScalarField& sigma) = 0;

        //- Update from dictionary
        virtual bool read(const dictionary& dictProperties);
};

} // End namespace Foam

#endif

// ************************************************************************* //
