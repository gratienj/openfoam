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
    Foam::implicitFunction::impCylinder

Description
    creates a infintite long cylinder

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    impCylinder.C

\*---------------------------------------------------------------------------*/

#ifndef impCylinder_H
#define impCylinder_H

#include "implicitFunctions.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace implicitFunction
{
// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class impCylinder Declaration
\*---------------------------------------------------------------------------*/

class impCylinder
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

        vector direction_;

        tensor project_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        impCylinder(const impCylinder&);

        //- Disallow default bitwise assignment
        void operator=(const impCylinder&);


public:

    //- Runtime type information
    TypeName("cylinder");


    // Constructors

        //- Construct from components
        impCylinder( const point&, const scalar radius,const scalar scale,const vector direction);

        //- Construct from dictionary (used by implicitFunctions)
        impCylinder
        (
            const dictionary& dict
        );

    //- Destructor
    virtual ~impCylinder();

    virtual scalar value(const vector p) 
    {
        return (-mag(project_ & (p-origin_))+radius_)*scale_;
    }

    virtual vector grad(const vector p) 
    {
        return -(project_ & (p-origin_))*scale_;
    }

    virtual scalar distanceToSurfaces(const vector p) 
    {
        return mag(mag(project_ & (p-origin_))-radius_)*scale_;
    }


    // Member Functions


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace implicitFunction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
