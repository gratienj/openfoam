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

Class
    Foam::implicitFunction::impSin

Description

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    impSin.C

\*---------------------------------------------------------------------------*/

#ifndef impSin_H
#define impSin_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impSin Declaration
\*---------------------------------------------------------------------------*/

class impSin
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- origin point
        const scalar period_;

        //- Radius
        const scalar phase_;

        const scalar amplitude_;

        vector up_;

        vector direction_;

        vector origin_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        impSin(const impSin&);

        //- Disallow default bitwise assignment
        void operator=(const impSin&);


public:

    //- Runtime type information
    TypeName("sin");


    // Constructors

        //- Construct from components
        impSin
        (
            const scalar period,
            const scalar phase,
            const scalar amplitude,
            vector direction,
            vector up,
            vector origin
        );

        //- Construct from dictionary (used by implicitFunctions)
        impSin
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impSin();

    virtual scalar value(const vector p) 
    {

        scalar x = (p - origin_) & direction_;
        scalar z = (p - origin_) & -up_;

        return z - amplitude_*Foam::sin(2*constant::mathematical::pi*x/period_ + phase_);
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
