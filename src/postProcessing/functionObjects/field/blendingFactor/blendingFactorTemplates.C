/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "gaussConvectionScheme.H"
#include "boundedConvectionScheme.H"
#include "blendedSchemeBase.H"
#include "fvcCellReduce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::blendingFactor::calc()
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    if (!foundField<FieldType>(fieldName_))
    {
        return false;
    }

    const FieldType& field = lookupField<FieldType>(fieldName_);

    const word divScheme("div(" + phiName_ + ',' + fieldName_ + ')');
    ITstream& its = mesh_.divScheme(divScheme);

    const surfaceScalarField& phi = lookupField<surfaceScalarField>(phiName_);

    tmp<fv::convectionScheme<Type>> cs =
        fv::convectionScheme<Type>::New(mesh_, phi, its);

    const fv::convectionScheme<Type>& cs = tcs();

    if (isA<fv::boundedConvectionScheme<Type>>(cs))
    {
        const fv::boundedConvectionScheme<Type>& bcs =
            refCast<const fv::boundedConvectionScheme<Type>>(cs);

    // Retrieve the face-based blending factor
    const blendedSchemeBase<Type>& blendedScheme =
        refCast<const blendedSchemeBase<Type>>(interpScheme);
    tmp<surfaceScalarField> factorf(blendedScheme.blendingFactor(field));

    // Convert into vol field whose values represent the local face maxima
    return store
    (
        resultName_,
        fvc::cellReduce(factorf, maxEqOp<scalar>())
    );
}


// ************************************************************************* //
