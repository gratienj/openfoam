/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::yPlus

Group
    grpFieldFunctionObjects

Description
    This function object evaluates and outputs turbulence y+ for turbulence
    models. The field is stored on the mesh database so that it can be
    retrieved and used for other applications.

    Example of function object specification to calculate the y+ (LES):
    \verbatim
    yPlus1
    {
        type        yPlus;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: yPlus        | yes         |
        phiName      | Name of flux field      | no          | phi
        resultName   | Name of y+ field        | no          | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


SeeAlso
    Foam::functionObject
    Foam::functionObjects::writeFiles
    Foam::functionObjects::timeControl

SourceFiles
    yPlus.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_yPlus_H
#define functionObjects_yPlus_H

#include "writeFiles.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class turbulenceModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class yPlus Declaration
\*---------------------------------------------------------------------------*/

class yPlus
:
    public writeFiles
{
    // Private data

        //- Name of mass/volume flux field (optional, default = phi)
        word phiName_;

        //- Result name
        word resultName_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate y+
        void calcYPlus
        (
            const turbulenceModel& turbModel,
            const fvMesh& mesh,
            volScalarField& yPlus
        );

        //- Disallow default bitwise copy construct
        yPlus(const yPlus&);

        //- Disallow default bitwise assignment
        void operator=(const yPlus&);


public:

    //- Runtime type information
    TypeName("yPlus");


    // Constructors

        //- Construct from Time and dictionary
        yPlus
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~yPlus();


    // Member Functions

        //- Read the yPlus data
        virtual bool read(const dictionary&);

        //- Calculate the yPlus field
        virtual bool execute(const bool postProcess = false);

        //- Write the yPlus field
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
