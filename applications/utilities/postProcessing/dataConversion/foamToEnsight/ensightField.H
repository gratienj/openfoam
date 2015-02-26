/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

InApplication
    foamToEnsight

Description

SourceFiles
    ensightField.C

\*---------------------------------------------------------------------------*/

#ifndef ensightField_H
#define ensightField_H

#include "ensightMesh.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Wrapper to get hold of the field or the subsetted field
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
volField
(
    const Foam::fvMeshSubset&,
    const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& vf
);


//- Ensight component ordering (different for symmTensor)
template<class Type>
class ensightTraits
{

public:

    static Foam::labelList initOrder()
    {
        return Foam::identity(Foam::pTraits<Type>::nComponents);
    }

    static const Foam::labelList order;
};


template<>
Foam::labelList ensightTraits<Foam::symmTensor>::initOrder()
{
    Foam::labelList order(6);
    // Ensight order for symmTensor is xx yy zz xy xz yz
    order[0] = 0;
    order[1] = 3;
    order[2] = 5;
    order[3] = 1;
    order[4] = 2;
    order[5] = 4;
    return order;
}


template<class Type>
const Foam::labelList ensightTraits<Type>::order
(
    ensightTraits<Type>::initOrder()
);


template<class Type>
void ensightField
(
    const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>& vf,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    const bool binary,
    const bool nodeValues,
    Foam::Ostream& ensightCaseFile
);


template<class Type>
void writePatchField
(
    const Foam::word& fieldName,
    const Foam::Field<Type>& pf,
    const Foam::word& patchName,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    Foam::Ostream& ensightCaseFile
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ensightField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
