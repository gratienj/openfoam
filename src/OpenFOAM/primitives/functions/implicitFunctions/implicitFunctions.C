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

#include "implicitFunctions.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(implicitFunctions, 0);
    defineRunTimeSelectionTable(implicitFunctions, dict);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::implicitFunctions> Foam::implicitFunctions::New
(
    const word& implicitFunctionsType,
    const dictionary& dict
)
{
    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_->find(implicitFunctionsType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown implicitFunctions type " << implicitFunctionsType
            << endl << endl
            << "Valid implicitFunctions types : " << endl
            << dictConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<implicitFunctions>(cstrIter()( dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::implicitFunctions::implicitFunctions()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::implicitFunctions::~implicitFunctions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
