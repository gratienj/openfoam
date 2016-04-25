/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "regionSplit.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "globalIndex.H"
#include "syncTools.H"
#include "FaceCellWave.H"
#include "minData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(regionSplit, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionSplit::calcNonCompactRegionSplit
(
    const label facei,
    const label otherFaceI,

    labelList& cellRegion
) const
{
    if (faceRegion[facei] >= 0)
    {
        if (faceRegion[otherFaceI] == -1)
        {
            faceRegion[otherFaceI] = faceRegion[facei];
            newChangedFaces.append(otherFaceI);
        }
        else if (faceRegion[otherFaceI] == -2)
        {
            // otherFaceI blocked but facei is not. Is illegal for coupled
            // faces, not for explicit connections.
        }
        else if (faceRegion[otherFaceI] != faceRegion[facei])
        {
            FatalErrorInFunction
                  << "Problem : coupled face " << facei
                  << " on patch " << mesh().boundaryMesh().whichPatch(facei)
                  << " has region " << faceRegion[facei]
                  << " but coupled face " << otherFaceI
                  << " has region " << faceRegion[otherFaceI]
                  << endl
                  << "Is your blocked faces specification"
                  << " synchronized across coupled boundaries?"
                  << abort(FatalError);
        }
    }
    else if (faceRegion[facei] == -1)
    {
        if (blockedFace.size() && blockedFace[faceI])
        {
            faceRegion[facei] = faceRegion[otherFaceI];
            newChangedFaces.append(facei);
        }
        else
        {
            // otherFaceI blocked but facei is not. Is illegal for coupled
            // faces, not for explicit connections.
        }
    }

    // Seed unblocked faces
    labelList seedFaces(nUnblocked);
    List<minData> seedData(nUnblocked);
    nUnblocked = 0;


    forAll(faceData, faceI)
    {
        label facei = cFaces[i];

        if (faceRegion[facei] == -1)
        {
            faceRegion[facei] = markValue;
            changedFaces[nFaces++] = facei;
        }
    }


    // Propagate information inwards
    FaceCellWave<minData> deltaCalc
    (
        mesh(),
        explicitConnections,
        false,  // disable walking through cyclicAMI for backwards compatibility
        seedFaces,
        seedData,
        faceData,
        cellData,
        mesh().globalData().nTotalCells()+1
    );


    // And extract
    cellRegion.setSize(mesh().nCells());
    forAll(cellRegion, cellI)
    {
        if (cellData[cellI].valid(deltaCalc.data()))
        {
            label facei = changedFaces[i];

            label own = mesh().faceOwner()[facei];

            if (cellRegion[own] == -1)
            {
                cellRegion[own] = markValue;
                changedCells.append(own);
            }

            if (mesh().isInternalFace(facei))
            {
                label nei = mesh().faceNeighbour()[facei];

                if (cellRegion[nei] == -1)
                {
                    cellRegion[nei] = markValue;
                    changedCells.append(nei);
                }
            }
        }
        else
        {
            label celli = changedCells[i];

            const cell& cFaces = mesh().cells()[celli];

            if (blockedFace.size() && !blockedFace[faceI])
            {
                label facei = cFaces[cFaceI];

                if (faceRegion[facei] == -1)
                {
                    faceRegion[facei] = markValue;
                    newChangedFaces.append(facei);
                }
            }
        }


        //if (debug)
        //{
        //    Pout<< "regionSplit::fillSeedMask : changedFaces before sync:"
        //        << changedFaces.size() << endl;
        //}


        // Check for changes to any locally coupled face.
        // Global connections are done later.

        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if
            (
                isA<cyclicPolyPatch>(pp)
             && refCast<const cyclicPolyPatch>(pp).owner()
            )
            {
                // Transfer from neighbourPatch to here or vice versa.

                const cyclicPolyPatch& cycPatch =
                    refCast<const cyclicPolyPatch>(pp);

                label facei = cycPatch.start();

                forAll(cycPatch, i)
                {
                    label otherFaceI = cycPatch.transformGlobalFace(facei);

                    transferCoupledFaceRegion
                    (
                        facei,
                        otherFaceI,
                        faceRegion,
                        newChangedFaces
                    );

                    facei++;
                }
            }
            cellRegion[cellI] = globalFaces.toGlobal(faceI);
        }
    }
}


Foam::autoPtr<Foam::globalIndex> Foam::regionSplit::calcRegionSplit
(
    const bool doGlobalRegions,
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,

    labelList& cellRegion
) const
{
    if (debug)
    {
        if (blockedFace.size())
        {
            // Check that blockedFace is synced.
            boolList syncBlockedFace(blockedFace);
            syncTools::swapFaceList(mesh(), syncBlockedFace);

            forAll(syncBlockedFace, facei)
            {
                if (syncBlockedFace[facei] != blockedFace[facei])
                {
                    FatalErrorInFunction
                        << "Face " << facei << " not synchronised. My value:"
                        << blockedFace[facei] << "  coupled value:"
                        << syncBlockedFace[facei]
                        << abort(FatalError);
                }
            }
        }
    }

    // Region per face.
    // -1 unassigned
    // -2 blocked
    labelList faceRegion(mesh().nFaces(), -1);

    if (blockedFace.size())
    {
        forAll(blockedFace, facei)
        {
            if (blockedFace[facei])
            {
                faceRegion[facei] = -2;
            }
        }
    }


    // Assign local regions
    // ~~~~~~~~~~~~~~~~~~~~

    // Start with region 0
    label nLocalRegions = 0;


    if (!doGlobalRegions)
    {
        // Block all parallel faces to avoid comms across
        boolList coupledOrBlockedFace(blockedFace);
        const polyBoundaryMesh& pbm = mesh().boundaryMesh();

        if (coupledOrBlockedFace.size())
        {
            forAll(pbm, patchI)
            {
                const polyPatch& pp = pbm[patchI];
                if (isA<processorPolyPatch>(pp))
                {
                    label faceI = pp.start();
                    forAll(pp, i)
                    {
                        coupledOrBlockedFace[faceI++] = true;
                    }
                }
            }
        }

        // Create dummy (local only) globalIndex
        labelList offsets(Pstream::nProcs()+1, 0);
        for (label i = Pstream::myProcNo()+1; i < offsets.size(); i++)
        {
            offsets[i] = mesh().nFaces();
        }
        const globalIndex globalRegions(offsets.xfer());

        // Minimise regions across connected cells
        // Note: still uses global decisions so all processors are running
        //       in lock-step, i.e. slowest determines overall time.
        //       To avoid this we could switch off Pstream::parRun.
        calcNonCompactRegionSplit
        (
            globalRegions,
            coupledOrBlockedFace,
            explicitConnections,
            cellRegion
        );

        // Current unsetCell has now been handled. Go to next region.
        nLocalRegions++;
        unsetCellI++;
    }
    while (true);


    if (debug)
    {
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] < 0)
            {
                FatalErrorInFunction
                    << "cell:" << celli << " region:" << cellRegion[celli]
                    << abort(FatalError);
            }
            else
            {
                globalRegion = fnd();
            }
            cellRegion[cellI] = globalRegion;
        }

        forAll(faceRegion, facei)
        {
            if (faceRegion[facei] == -1)
            {
                FatalErrorInFunction
                    << "face:" << facei << " region:" << faceRegion[facei]
                    << abort(FatalError);
            }
        }

        return autoPtr<globalIndex>(new globalIndex(compactOffsets.xfer()));
    }



    // Initial global region numbers
    const globalIndex globalRegions(mesh().nFaces());

    // Minimise regions across connected cells (including parallel)
    calcNonCompactRegionSplit
    (
        globalRegions,
        blockedFace,
        explicitConnections,
        cellRegion
    );

    if (!doGlobalRegions)
    {
        return autoPtr<globalIndex>(new globalIndex(nLocalRegions));
    }


    // 2. Assign global regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Offset local regions to create unique global regions.

    globalIndex globalRegions(nLocalRegions);


    // Convert regions to global ones
    forAll(cellRegion, celli)
    {
        cellRegion[celli] = globalRegions.toGlobal(cellRegion[celli]);
    }

    // Now our cellRegion will have
    // - non-local regions (i.e. originating from other processors)
    // - non-compact locally originating regions
    // so we'll need to compact

    // 4a: count per originating processor the number of regions
    labelList nOriginating(Pstream::nProcs(), 0);
    {
        labelHashSet haveRegion(mesh().nCells()/8);

        labelList nbrRegion(mesh().nFaces()-mesh().nInternalFaces(), -1);
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            // Count originating processor. Use isLocal as efficiency since
            // most cells are locally originating.
            if (globalRegions.isLocal(region))
            {
                if (haveRegion.insert(region))
                {
                    label facei = pp.start()+i;
                    if (!blockedFace.size() || !blockedFace[facei])
                    {
                        patchNbrRegion[i] = cellRegion[patchCells[i]];
                    }
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh(), nbrRegion);

        Map<label> globalToMerged(mesh().nFaces()-mesh().nInternalFaces());

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                label procI = globalRegions.whichProcID(region);
                if (haveRegion.insert(region))
                {
                    label facei = pp.start()+i;

                    if (!blockedFace.size() || !blockedFace[facei])
                    {
                        if (patchNbrRegion[i] < cellRegion[patchCells[i]])
                        {
                            //Pout<< "on patch:" << pp.name()
                            //    << " cell:" << patchCells[i]
                            //    << " at:"
                            // << mesh().cellCentres()[patchCells[i]]
                            //    << " was:" << cellRegion[patchCells[i]]
                            //    << " nbr:" << patchNbrRegion[i]
                            //    << endl;

                            globalToMerged.insert
                            (
                                cellRegion[patchCells[i]],
                                patchNbrRegion[i]
                            );
                        }
                    }
                }
            }
        }


        label nMerged = returnReduce(globalToMerged.size(), sumOp<label>());

        if (debug)
        {
            Pout<< "nMerged:" << nMerged << endl;
        }

        if (nMerged == 0)
        {
            break;
        }

        // Renumber the regions according to the globalToMerged
        forAll(cellRegion, celli)
        {
            label regionI = cellRegion[celli];
            Map<label>::const_iterator iter = globalToMerged.find(regionI);
            if (iter != globalToMerged.end())
            {
                 cellRegion[celli] = iter();
            }
        }
    }


    // Now our cellRegion will have non-local elements in it. So compact
    // it.

    // 4a: count. Use a labelHashSet to count regions only once.
    label nCompact = 0;
    {
        labelHashSet localRegion(mesh().nFaces()-mesh().nInternalFaces());
        forAll(cellRegion, celli)
        {
            if
            (
                globalRegions.isLocal(cellRegion[celli])
             && localRegion.insert(cellRegion[celli])
            )
            {
                nCompact++;
            }
        }
    }

    if (debug)
    {
        Pout<< "Counted " << nOriginating[Pstream::myProcNo()]
            << " local regions." << endl;
    }


    // Global numbering for compacted local regions
    autoPtr<globalIndex> globalCompactPtr
    (
        new globalIndex(nOriginating[Pstream::myProcNo()])
    );
    const globalIndex& globalCompact = globalCompactPtr();


    // 4b: renumber
    // Renumber into compact indices. Note that since we've already made
    // all regions global we now need a Map to store the compacting information
    // instead of a labelList - otherwise we could have used a straight
    // labelList.

    // Local compaction map
    Map<label> globalToCompact(2*nOriginating[Pstream::myProcNo()]);
    // Remote regions we want the compact number for
    List<labelHashSet> nonLocal(Pstream::nProcs());
    forAll(nonLocal, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            nonLocal[procI].resize(2*nOriginating[procI]);
        }
    }

    forAll(cellRegion, celli)
    {
        label region = cellRegion[celli];
        if (globalRegions.isLocal(region))
        {
            // Insert new compact region (if not yet present)
            globalToCompact.insert
            (
                region,
                globalCompact.toGlobal(globalToCompact.size())
            );
        }
        else
        {
            nonLocal[globalRegions.whichProcID(region)].insert(region);
        }
    }


    // Now we have all the local regions compacted. Now we need to get the
    // non-local ones from the processors to whom they are local.
    // Convert the nonLocal (labelHashSets) to labelLists.

    labelListList sendNonLocal(Pstream::nProcs());
    forAll(sendNonLocal, procI)
    {
        sendNonLocal[procI] = nonLocal[procI].toc();
    }

    if (debug)
    {
        forAll(sendNonLocal, procI)
        {
            Pout<< "    from processor " << procI
                << " want " << sendNonLocal[procI].size()
                << " region numbers."
                << endl;
        }
        Pout<< endl;
    }


    // Get the wanted region labels into recvNonLocal
    labelListList recvNonLocal;
    Pstream::exchange<labelList, label>(sendNonLocal, recvNonLocal);

    // Now we have the wanted compact region labels that procI wants in
    // recvNonLocal[procI]. Construct corresponding list of compact
    // region labels to send back.

    labelListList sendWantedLocal(Pstream::nProcs());
    forAll(recvNonLocal, procI)
    {
        const labelList& nonLocal = recvNonLocal[procI];
        sendWantedLocal[procI].setSize(nonLocal.size());

        forAll(nonLocal, i)
        {
            sendWantedLocal[procI][i] = globalToCompact[nonLocal[i]];
        }
    }


    // Send back (into recvNonLocal)
    recvNonLocal.clear();
    Pstream::exchange<labelList, label>(sendWantedLocal, recvNonLocal);
    sendWantedLocal.clear();

    // Now recvNonLocal contains for every element in setNonLocal the
    // corresponding compact number. Insert these into the local compaction
    // map.

    forAll(recvNonLocal, procI)
    {
        const labelList& wantedRegions = sendNonLocal[procI];
        const labelList& compactRegions = recvNonLocal[procI];

        forAll(wantedRegions, i)
        {
            globalToCompact.insert(wantedRegions[i], compactRegions[i]);
        }
    }

    // Finally renumber the regions
    forAll(cellRegion, celli)
    {
        cellRegion[celli] = globalToCompact[cellRegion[celli]];
    }

    return globalCompactPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSplit::regionSplit(const polyMesh& mesh, const bool doGlobalRegions)
:
    MeshObject<polyMesh, Foam::TopologicalMeshObject, regionSplit>(mesh),
    labelList(mesh.nCells(), -1)
{
    globalNumberingPtr_ = calcRegionSplit
    (
        doGlobalRegions,    //do global regions
        boolList(0, false), //blockedFaces
        List<labelPair>(0), //explicitConnections,
        *this
    );
}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const bool doGlobalRegions
)
:
    MeshObject<polyMesh, Foam::TopologicalMeshObject, regionSplit>(mesh),
    labelList(mesh.nCells(), -1)
{
    globalNumberingPtr_ = calcRegionSplit
    (
        doGlobalRegions,
        blockedFace,        //blockedFaces
        List<labelPair>(0), //explicitConnections,
        *this
    );
}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,
    const bool doGlobalRegions
)
:
    MeshObject<polyMesh, Foam::TopologicalMeshObject, regionSplit>(mesh),
    labelList(mesh.nCells(), -1)
{
    globalNumberingPtr_ = calcRegionSplit
    (
        doGlobalRegions,
        blockedFace,            //blockedFaces
        explicitConnections,    //explicitConnections,
        *this
    );
}


// ************************************************************************* //
