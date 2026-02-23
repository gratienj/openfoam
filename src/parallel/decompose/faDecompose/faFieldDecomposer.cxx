/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2026 OpenCFD Ltd.
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

#include "faFieldDecomposer.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faFieldDecomposer::patchFieldDecomposer::patchFieldDecomposer
(
    const label sizeBeforeMapping,
    const labelUList& addressingSlice,
    const label addressingOffset
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    directAddressing_(addressingSlice)
{
    forAll(directAddressing_, i)
    {
// if constexpr (withTurningIndex_)
        // Subtract one to align addressing.
        // directAddressing_[i] -= addressingOffset + 1;
        // ZT, 12/Nov/2010
        directAddressing_[i] -= addressingOffset;
    }
}


Foam::faFieldDecomposer::processorAreaPatchFieldDecomposer::
processorAreaPatchFieldDecomposer
(
    const label nTotalFaces,
    const labelUList& owner,  // == mesh.edgeOwner()
    const labelUList& neigh,  // == mesh.edgeNeighbour()
    const labelUList& addressingSlice,
    const bitSet& flip
)
:
    sizeBeforeMapping_(nTotalFaces),
    directAddressing_(addressingSlice.size())
{
    forAll(directAddressing_, i)
    {
        // Subtract one to align addressing.
        label ai = addressingSlice[i];
// if constexpr (withTurningIndex_)
//         label ai = mag(addressingSlice[i]) - 1;

        if (ai < neigh.size())
        {
            // This is a regular edge. it has been an internal edge
            // of the original mesh and now it has become a edge
            // on the parallel boundary

// if constexpr (withTurningIndex_)
            // if (addressingSlice[i] >= 0)
            if (!flip.test(i))
            {
                // We are the owner side so use the neighbour value
                directAddressing_[i] = neigh[ai];
            }
            else
            {
                // We are the neighbour side so use the owner value
                directAddressing_[i] = owner[ai];
            }
        }
        else
        {
            // This is a edge that used to be on a cyclic boundary
            // but has now become a parallel patch edge. I cannot
            // do the interpolation properly (I would need to look
            // up the different (edge) list of data), so I will
            // just grab the value from the owner face

            directAddressing_[i] = owner[ai];
        }
    }
}


Foam::faFieldDecomposer::processorEdgePatchFieldDecomposer::
processorEdgePatchFieldDecomposer
(
    label sizeBeforeMapping,
    const labelUList& addressingSlice
)
:
    sizeBeforeMapping_(sizeBeforeMapping),
    addressing_(addressingSlice.size()),
    weights_(addressingSlice.size())
{
    forAll(addressing_, i)
    {
        addressing_[i].resize(1);
        weights_[i].resize(1);

// if constexpr (withTurningIndex_)
//      addressing_[i][0] = mag(addressingSlice[i]) - 1;
//      weights_[i][0] = sign(addressingSlice[i]);

        addressing_[i][0] = mag(addressingSlice[i]);
        weights_[i][0] = 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faFieldDecomposer::faFieldDecomposer
(
    Foam::zero,
    const faMesh& procMesh,
    const labelUList& edgeAddressing,
    const labelUList& faceAddressing,
    const labelUList& boundaryAddressing
)
:
    procMesh_(procMesh),
    edgeAddressing_(edgeAddressing),
    faceAddressing_(faceAddressing),
    boundaryAddressing_(boundaryAddressing)
{}


Foam::faFieldDecomposer::faFieldDecomposer
(
    const faMesh& completeMesh,
    const faMesh& procMesh,
    const labelUList& edgeAddressing,
    const labelUList& faceAddressing,
    const labelUList& boundaryAddressing
)
:
    faFieldDecomposer
    (
        Foam::zero{},
        procMesh,
        edgeAddressing,
        faceAddressing,
        boundaryAddressing
    )
{
    reset(completeMesh);
}


Foam::faFieldDecomposer::faFieldDecomposer
(
    const label nTotalFaces,
    const UList<labelRange>& boundaryRanges,
    const labelUList& edgeOwner,
    const labelUList& edgeNeigbour,

    const faMesh& procMesh,
    const labelUList& edgeAddressing,
    const labelUList& faceAddressing,
    const labelUList& boundaryAddressing
)
:
    faFieldDecomposer
    (
        Foam::zero{},
        procMesh,
        edgeAddressing,
        faceAddressing,
        boundaryAddressing
    )
{
    reset(nTotalFaces, boundaryRanges, edgeOwner, edgeNeigbour);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faFieldDecomposer::empty() const noexcept
{
    return patchFieldDecomposers_.empty();
}


void Foam::faFieldDecomposer::clear()
{
    patchFieldDecomposers_.clear();
    processorAreaPatchFieldDecomposers_.clear();
    processorEdgePatchFieldDecomposers_.clear();
    edgeSigns_.clear();
}


void Foam::faFieldDecomposer::reset
(
    const label nTotalFaces,
    const UList<labelRange>& boundaryRanges,
    const labelUList& edgeOwner,
    const labelUList& edgeNeigbour
)
{
    const label nMappers = procMesh_.boundary().size();

    patchFieldDecomposers_.resize_null(nMappers);
    processorAreaPatchFieldDecomposers_.resize_null(nMappers);
    processorEdgePatchFieldDecomposers_.resize_null(nMappers);
    edgeSigns_.resize_null(nMappers);

    bitSet flipMap;

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];
        const faPatch& fap = procMesh_.boundary()[patchi];
        const labelSubList localPatchSlice(fap.patchSlice(edgeAddressing_));

        if (oldPatchi >= 0)
        {
            patchFieldDecomposers_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    boundaryRanges[oldPatchi].size(),
                    localPatchSlice,
                    boundaryRanges[oldPatchi].start()
                )
            );
        }
        else
        {
            // No oldPatch - is processor patch. edgeAddressing_ does
            // not have 'flip' sign so use the face map to see which side
            // we've got.
            flipMap.clear();
            flipMap.resize(localPatchSlice.size());

            forAll(localPatchSlice, i)
            {
                const label ownFacei = faceAddressing_[fap.edgeFaces()[i]];
                flipMap.set(i, (edgeOwner[localPatchSlice[i]] != ownFacei));
            }

            processorAreaPatchFieldDecomposers_.set
            (
                patchi,
                new processorAreaPatchFieldDecomposer
                (
                    nTotalFaces,
                    edgeOwner,
                    edgeNeigbour,
                    localPatchSlice,
                    flipMap
                )
            );

            processorEdgePatchFieldDecomposers_.set
            (
                patchi,
                new processorEdgePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    localPatchSlice
                )
            );

            auto& s = edgeSigns_.emplace_set(patchi, localPatchSlice.size());
            forAll(localPatchSlice, i)
            {
                s[i] = (flipMap.test(i) ? -1 : 1);
            }
        }
    }
}


void Foam::faFieldDecomposer::reset(const faMesh& completeMesh)
{
    const label nMappers = procMesh_.boundary().size();

    patchFieldDecomposers_.resize_null(nMappers);
    processorAreaPatchFieldDecomposers_.resize_null(nMappers);
    processorEdgePatchFieldDecomposers_.resize_null(nMappers);
    edgeSigns_.resize_null(nMappers);

    // Create weightings now - needed for proper parallel synchronization
    //// (void)completeMesh.weights();
    // Disabled the above (2022-04-04)
    // Use weights if they already exist, otherwise simply ignore

    // faPatches don't have their own start() - so these are invariant
    const labelList completePatchStarts
    (
        completeMesh.boundary().patchStarts()
    );

    const auto& edgeOwner = completeMesh.edgeOwner();
    const auto& edgeNeighbour = completeMesh.edgeNeighbour();

    bitSet flipMap;

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];
        const faPatch& fap = procMesh_.boundary()[patchi];
        const labelSubList localPatchSlice(fap.patchSlice(edgeAddressing_));

        if (oldPatchi >= 0)
        {
            patchFieldDecomposers_.set
            (
                patchi,
                new patchFieldDecomposer
                (
                    completeMesh.boundary()[oldPatchi].size(),
                    localPatchSlice,
                    completePatchStarts[oldPatchi]
                )
            );
        }
        else
        {
            // No oldPatch - is processor patch. edgeAddressing_ does
            // not have 'flip' sign so use the face map to see which side
            // we've got.
            flipMap.clear();
            flipMap.resize(localPatchSlice.size());

            forAll(localPatchSlice, i)
            {
                const label ownFacei = faceAddressing_[fap.edgeFaces()[i]];
                flipMap.set(i, (edgeOwner[localPatchSlice[i]] != ownFacei));
            }

            processorAreaPatchFieldDecomposers_.set
            (
                patchi,
                new processorAreaPatchFieldDecomposer
                (
                    completeMesh.nFaces(),
                    edgeOwner,
                    edgeNeighbour,
                    localPatchSlice,
                    flipMap
                )
            );

            processorEdgePatchFieldDecomposers_.set
            (
                patchi,
                new processorEdgePatchFieldDecomposer
                (
                    procMesh_.boundary()[patchi].size(),
                    localPatchSlice
                )
            );

            auto& s = edgeSigns_.emplace_set(patchi, localPatchSlice.size());
            forAll(localPatchSlice, i)
            {
                s[i] = (flipMap.test(i) ? -1 : 1);
            }
        }
    }
}


// ************************************************************************* //
