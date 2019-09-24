/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "weightedFlux.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class Type>
void Foam::weightedFlux<Type>::clearOut()
{
    deleteDemandDrivenData(oDelta_);
    deleteDemandDrivenData(nDelta_);
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::weightedFlux<Type>::~weightedFlux()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::weightedFlux<Type>::makeDeltas() const
{

    const fvMesh& mesh = this->mesh();

    oDelta_ = new surfaceScalarField
    (
        IOobject
        (
            "oDelta",
            mesh.pointsInstance(),
            mesh
        ),
        mesh,
        dimLength
    );
    surfaceScalarField& oDelta = *oDelta_;

    nDelta_ = new surfaceScalarField
    (
        IOobject
        (
            "nDelta",
            mesh.pointsInstance(),
            mesh
        ),
        mesh,
        dimLength
    );
    surfaceScalarField& nDelta = *nDelta_;

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const surfaceVectorField n = mesh.Sf()/mesh.magSf();

    const vectorField& C = mesh.cellCentres();
    const vectorField& Cf = mesh.faceCentres();

    // all distances are NORMAL to the face,
    // as in the weighting factors in surfaceInterpolation.C
    forAll(owner, facei)
    {
        oDelta[facei] =
            mag(n[facei] & (C[owner[facei]] - Cf[facei]));
        nDelta[facei] =
            mag(n[facei] & (C[neighbour[facei]] - Cf[facei]));
    }

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        // patch normal vector
        const vectorField nPatch = currPatch.Sf()/currPatch.magSf();

        // processor patch
        if (currPatch.coupled())
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
            const vectorField& pCf = mesh.Cf().boundaryField()[patchi];

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                // all distances are NORMAL to the face
                oDelta.boundaryFieldRef()[patchi][facei] =
                    mag(nPatch[facei] & (pCf[facei] - C[own]));
            }

            // weight = delta_neighbour / delta in ORTHOGONAL direction,
            nDelta.boundaryFieldRef()[patchi] =
                currPatch.weights()*oDelta.boundaryFieldRef()[patchi]
               /(1.0 - currPatch.weights());
        }
        else
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
            const vectorField& pCf = mesh.Cf().boundaryField()[patchi];

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];
                // all distances are NORMAL to the face!
                oDelta.boundaryFieldRef()[patchi][facei] =
                    mag(nPatch[facei] & (pCf[facei] - C[own]));
                nDelta.boundaryFieldRef()[patchi][facei] =
                    mag(nPatch[facei] & (pCf[facei] - C[own]));
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::weightedFlux<Type>::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = vf.mesh();
    const surfaceScalarField& oDelta = weightedFlux<Type>::oDelta();
    const surfaceScalarField& nDelta = weightedFlux<Type>::nDelta();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tvff
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "weightedFlux::interpolate(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& vff = tvff.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(vff, facei)
    {
        scalar sigmaDeltaO = sigma_[owner[facei]]/oDelta[facei];
        scalar sigmaDeltaN = sigma_[neighbour[facei]]/nDelta[facei];

        vff[facei] =
            (vf[owner[facei]]*sigmaDeltaO + vf[neighbour[facei]]*sigmaDeltaN)
           /(sigmaDeltaO + sigmaDeltaN);
    }

    // interpolate patches
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bvff =
        vff.boundaryFieldRef();

    forAll(bvff, patchi)
    {
        fvsPatchField<Type>& pvff = bvff[patchi];
        // if not coupled - simply copy the boundary values of the field
        if (!pvff.coupled())
        {
            pvff = vf.boundaryField()[patchi];
        }
        // e.g. processor patches have to calculated separately
        else
        {
            const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
            scalarField sigmaN =
                sigma_.boundaryField()[patchi].patchNeighbourField();

            Field<Type> vfO = vf.boundaryField()[patchi].patchInternalField();
            Field<Type> vfN = vf.boundaryField()[patchi].patchNeighbourField();

            forAll(pOwner, facei)
            {
                label own = pOwner[facei];

                scalar sigmaDeltaO =
                    sigma_[own]/oDelta.boundaryField()[patchi][facei];
                scalar sigmaDeltaN =
                    sigmaN[facei]/nDelta.boundaryField()[patchi][facei];

                pvff[facei] =
                    (vfO[facei]*sigmaDeltaO + vfN[facei]*sigmaDeltaN)
                   /(sigmaDeltaO + sigmaDeltaN);
            }
        }
    }
    return tvff;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(weightedFlux)
}

// ************************************************************************* //
