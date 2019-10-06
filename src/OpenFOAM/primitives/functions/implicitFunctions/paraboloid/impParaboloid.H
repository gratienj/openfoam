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
    Foam::implicitFunction::impParaboloid

Description
    creates a paraboloid

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    impParaboloid.C

\*---------------------------------------------------------------------------*/

#ifndef impParaboloid_H
#define impParaboloid_H

#include "implicitFunctions.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impParaboloid Declaration
\*---------------------------------------------------------------------------*/

class impParaboloid
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- ax^2 + bx*y + cy^2
        vector coeffs_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        impParaboloid(const impParaboloid&);

        //- Disallow default bitwise assignment
        void operator=(const impParaboloid&);


public:

    //- Runtime type information
    TypeName("paraboloid");


    // Constructors

        //- Construct from components
        impParaboloid(const vector coeffs);

        //- Construct from dictionary (used by implicitFunctions)
        impParaboloid
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impParaboloid();

    virtual scalar value(const vector p) 
    {
          return coeffs_.x()*sqr(p.x()) + coeffs_.y()*p.x()*p.y()  + coeffs_.z()*sqr(p.y()) - p.z();
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
