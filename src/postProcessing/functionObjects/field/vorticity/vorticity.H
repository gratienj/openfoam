/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::functionObjects::vorticity

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates and outputs the vorticity, the curl of
    the velocity as a volvectorField.  The field is stored on the mesh
    database so that it can be retrieved and used for other applications.

    Example of function object specification to calculate the vorticity:
    \verbatim
    vorticity1
    {
        type        vorticity;
        functionObjectLibs ("libutilityFunctionObjects.so");
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: vorticity    | yes         |
        UName        | Name of velocity field  | no          | U
        resultName   | Name of Courant number field | no    | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable


SeeAlso
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    vorticity.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_vorticity_H
#define functionObjects_vorticity_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class vorticity Declaration
\*---------------------------------------------------------------------------*/

class vorticity
:
    public fieldExpression
{

public:

    //- Runtime type information
    TypeName("vorticity");


    // Constructors

        //- Construct from Time and dictionary
        vorticity
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~vorticity();


    // Member Functions

        //- Execute, currently does nothing
        virtual bool execute(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
