/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                isoAdvector | Copyright (C) 2016-2017 DHI
              Modified work | Copyright (C) 2018-2019 Johan Roenby
              Modified work | Copyright (C) 2019 DLR
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

#include "advectionSchemes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::advectionSchemes>
Foam::advectionSchemes::New
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
{
    const dictionary& solversDict = alpha1.mesh().solution::subDict("solvers");
    const word schemeType
    (
        solversDict.subDict(alpha1.name()).lookup("advectionScheme")
    );

    Info<< "Selecting advectionSchemes: " << schemeType << endl;

    auto cstrIter = componentsConstructorTablePtr_->find(schemeType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown advectionSchemes type "
            << schemeType << endl << endl
            << "Valid  advectionSchemess are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<advectionSchemes>(cstrIter()( alpha1, phi,U));
}


// ************************************************************************* //
