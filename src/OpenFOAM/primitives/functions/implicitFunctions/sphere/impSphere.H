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
    Foam::implicitFunction::impSphere

Description
    creates a sin wave

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    impSphere.C

\*---------------------------------------------------------------------------*/

#ifndef impSphere_H
#define impSphere_H

#include "implicitFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impSphere Declaration
\*---------------------------------------------------------------------------*/

class impSphere
:
    public implicitFunctions
{
private:

    // Private Member Data

        //- origin point
        const point origin_;

        //- Radius
        const scalar radius_;

        const scalar scale_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        impSphere(const impSphere&);

        //- Disallow default bitwise assignment
        void operator=(const impSphere&);


public:

    //- Runtime type information
    TypeName("sphere");


    // Constructors

        //- Construct from components
        impSphere( const point&, const scalar radius,const scalar scale);

        //- Construct from dictionary (used by implicitFunctions)
        impSphere
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impSphere();

    // Member Functions

    virtual scalar value(const vector p) 
    {
          return (-mag(p-origin_)+radius_)*scale_;
    }

    virtual vector grad(const vector p) 
    {
        return (origin_-p)*scale_;
    }

    virtual scalar distanceToSurfaces(const vector p) 
    {
        return mag(mag(p-origin_)-radius_)*scale_;
    }





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
