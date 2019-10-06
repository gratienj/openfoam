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

#include "impSin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace implicitFunction
    {
        defineTypeNameAndDebug(impSin, 0);
        addToRunTimeSelectionTable(implicitFunctions, impSin, dict);
    }

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunction::impSin::impSin
(
    const scalar period,
    const scalar phase,
    const scalar amplitude,
    const vector direction,
    const vector up,
    const vector origin
)
:
    period_(period),
    phase_(phase),
    amplitude_(amplitude),
    up_(up),
    direction_(direction),
    origin_(origin)
{

}


Foam::implicitFunction::impSin::impSin
(
    const dictionary& dict
)
:
    period_(readScalar(dict.lookup("period"))),
    phase_(dict.lookupOrDefault<scalar>("phase",0.0)),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    up_(dict.lookup("up")),
    direction_(dict.lookup("direction")),
    origin_(dict.lookup("origin"))
{
    direction_ /= mag(direction_);
    up_ /= mag(up_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::implicitFunction::impSin::~impSin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
