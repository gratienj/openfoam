/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2019-2019 DLR
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

\*---------------------------------------------------------------------------*/

#include "impSphere.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunction
    {
        defineTypeNameAndDebug(impSphere, 0);
        addToRunTimeSelectionTable(implicitFunctions, impSphere, dict);
    }

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunction::impSphere::impSphere
(
    const point& origin,
    const scalar radius,
    const scalar scale
)
:
    origin_(origin),
    radius_(radius),
    scale_(scale)
{

}


Foam::implicitFunction::impSphere::impSphere
(
    const dictionary& dict
)
:
    origin_(dict.lookup("origin")),
    radius_(readScalar(dict.lookup("radius"))),
    scale_(dict.lookupOrDefault<scalar>("scale",1))
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::implicitFunction::impSphere::~impSphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
