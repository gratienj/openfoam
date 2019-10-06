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
    Foam::leastSquareFitParabolid

Description
    Fits paraboloid in the a point cloud; fit is in the local
    coordinate system. 

    Original code supplied by Henning Scheufler, DLR (2019)

SourceFiles
    leastSquareFitParabolid.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquareFitParabolid_H
#define leastSquareFitParabolid_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"



#include "multiDimPolyFitter.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class leastSquareFitParabolid Declaration
\*---------------------------------------------------------------------------*/

class leastSquareFitParabolid
{



private:
    // Private data
    multiDimPolyFitter<scalar> polyFitter_;

    Vector<label> geomDir_;

    label nDims_;

    cartesianCS makeLocalCoordSystem
    (
        const point& basePoint,
        const vector& normal
    );


public:

    // Constructors

        //- Construct from fvMesh and globalIndex
        leastSquareFitParabolid
        (
            const Vector<label> geomDir,
            const Vector<label> explicitDim
        );


    // Member functions

        //- return the fit parameters
        scalarField fitParaboloid
        (
            const point& centre,
            const vector& normal,
            const vectorField& positions
        );

        //- return the fit parameters
        scalarField fitParaboloid
        (
            const point& centre,
            const vector& normal,
            const vectorField& positions,
            const scalarField& weight
        );

        label nCoeffs()
        {
            return polyFitter_.nCoeffs();
        } 




};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
