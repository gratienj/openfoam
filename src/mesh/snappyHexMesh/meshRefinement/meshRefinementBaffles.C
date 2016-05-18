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

#include "meshRefinement.H"
#include "fvMesh.H"
#include "syncTools.H"
#include "Time.H"
#include "refinementSurfaces.H"
#include "faceSet.H"
#include "indirectPrimitivePatch.H"
#include "polyTopoChange.H"
#include "meshTools.H"
#include "polyModifyFace.H"
#include "polyModifyCell.H"
#include "polyAddFace.H"
#include "polyRemoveFace.H"
#include "localPointRegion.H"
#include "duplicatePoints.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "unitConversion.H"
#include "OBJstream.H"
#include "patchFaceOrientation.H"
#include "PatchEdgeFaceWave.H"
#include "patchEdgeFaceRegion.H"
#include "polyMeshAdder.H"
#include "IOmanip.H"
#include "refinementParameters.H"

#include "zeroGradientFvPatchFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshRefinement::createBaffle
(
    const label facei,
    const label ownPatch,
    const label neiPatch,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[facei];
    label zoneID = mesh_.faceZones().whichZone(facei);
    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            f,                          // modified face
            facei,                      // label of face
            mesh_.faceOwner()[facei],   // owner
            -1,                         // neighbour
            false,                      // face flip
            ownPatch,                   // patch for face
            false,                      // remove from zone
            zoneID,                     // zone for face
            zoneFlip                    // face flip in zone
        )
    );


    label dupFaceI = -1;

    if (mesh_.isInternalFace(facei))
    {
        if (neiPatch == -1)
        {
            FatalErrorInFunction
                << "No neighbour patch for internal face " << facei
                << " fc:" << mesh_.faceCentres()[facei]
                << " ownPatch:" << ownPatch << abort(FatalError);
        }

        bool reverseFlip = false;
        if (zoneID >= 0)
        {
            reverseFlip = !zoneFlip;
        }

        dupFaceI = meshMod.setAction
        (
            polyAddFace
            (
                f.reverseFace(),            // modified face
                mesh_.faceNeighbour()[facei],// owner
                -1,                         // neighbour
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                facei,                      // masterFaceID,
                true,                       // face flip
                neiPatch,                   // patch for face
                zoneID,                     // zone for face
                reverseFlip                 // face flip in zone
            )
        );
    }
    return dupFaceI;
}


void Foam::meshRefinement::getBafflePatches
(
    const labelList& surfacesToTest,
    const pointField& neiCc,
    labelList& ownPatch,
    labelList& neiPatch
) const
{
    autoPtr<OBJstream> str;
    if (debug&OBJINTERSECTIONS)
    {
        mkDir(mesh_.time().path()/timeName());
        str.reset
        (
            new OBJstream
            (
                mesh_.time().path()/timeName()/"intersections.obj"
            )
        );

        Pout<< "getIntersections : Writing surface intersections to file "
            << str().name() << nl << endl;
    }


    globalRegion1.setSize(mesh_.nFaces());
    globalRegion1 = -1;
    globalRegion2.setSize(mesh_.nFaces());
    globalRegion2 = -1;


    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        label own = mesh_.faceOwner()[facei];

        if (mesh_.isInternalFace(facei))
        {
            start[i] = cellCentres[own];
            end[i] = cellCentres[mesh_.faceNeighbour()[facei]];
        }
        else
        {
            start[i] = cellCentres[own];
            end[i] = neiCc[facei-mesh_.nInternalFaces()];
        }
    }

    // Extend segments a bit
    {
        labelList minLevel;
        calcCellCellRays
        (
            neiCc,
            labelList(neiCc.size(), -1),
            testFaces,
            start,
            end,
            minLevel
        );
    }


    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    surfaces_.findNearestIntersection
    (
        surfacesToTest,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
    );


    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (hit1[i].hit() && hit2[i].hit())
        {
            if (str.valid())
            {
                str().write(linePointRef(start[i], hit1[i].rawPoint()));
                str().write
                (
                    linePointRef(hit1[i].rawPoint(), hit2[i].rawPoint())
                );
                str().write(linePointRef(hit2[i].rawPoint(), end[i]));
            }

            // Pick up the patches
            ownPatch[facei] = globalToMasterPatch
            [
                surfaces_.globalRegion(surface1[i], region1[i])
            ];
            neiPatch[facei] = globalToMasterPatch
            [
                surfaces_.globalRegion(surface2[i], region2[i])
            ];

            if (ownPatch[facei] == -1 || neiPatch[facei] == -1)
            {
                FatalErrorInFunction
                    << "problem." << abort(FatalError);
            }
        }
    }
}


void Foam::meshRefinement::getBafflePatches
(
    const labelList& globalToMasterPatch,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,

    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& ownPatch,
    labelList& neiPatch
) const
{
    // This determines the patches for the intersected faces to
    // - remove the outside of the mesh
    // - introduce baffles for (non-faceZone) intersections
    // Any baffles for faceZones (faceType 'baffle'/'boundary') get introduced
    // later


    // 1. Determine intersections with unnamed surfaces and cell zones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList cellToZone;
    labelList unnamedRegion1;
    labelList unnamedRegion2;
    labelList namedSurfaceIndex;
    {
        PackedBoolList posOrientation;
        zonify
        (
            true,               // allowFreeStandingZoneFaces
            -2,                 // zone to put unreached cells into
            locationsInMesh,
            zonesInMesh,

            cellToZone,
            unnamedRegion1,
            unnamedRegion2,
            namedSurfaceIndex,
            posOrientation
        );
    }


    // 2. Baffle all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // The logic is that all intersections with unnamed surfaces become
    // baffles except where we're inbetween a cellZone and background
    // or inbetween two different cellZones.

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    const labelList testFaces(intersectedFaces());

    ownPatch.setSize(mesh_.nFaces());
    ownPatch = -1;
    neiPatch.setSize(mesh_.nFaces());
    neiPatch = -1;
    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        if (unnamedRegion1[faceI] != -1 || unnamedRegion2[faceI] != -1)
        {
            label ownMasterPatch = -1;
            if (unnamedRegion1[faceI] != -1)
            {
                ownMasterPatch = globalToMasterPatch[unnamedRegion1[faceI]];
            }
            label neiMasterPatch = -1;
            if (unnamedRegion2[faceI] != -1)
            {
                neiMasterPatch = globalToMasterPatch[unnamedRegion2[faceI]];
            }


            // Now we always want to produce a baffle except if
            // one side is a valid cellZone

            label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
            label neiZone = -1;

            if (mesh_.isInternalFace(faceI))
            {
                neiZone = cellToZone[mesh_.faceNeighbour()[faceI]];
            }
            else
            {
                neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];
            }


            if
            (
                (ownZone != neiZone)
             && (
                    (ownZone >= 0 && neiZone != -2)
                 || (neiZone >= 0 && ownZone != -2)
                )
             && (
                    namedSurfaceIndex.size() == 0
                 || namedSurfaceIndex[faceI] == -1
                )
            )
            {
                // One side is a valid cellZone and the neighbour is a different
                // one (or -1 but not -2). Do not baffle unless the user has
                // put both an unnamed and a named surface there. In that
                // case assume that the user wants both: baffle and faceZone.
            }
            else
            {
                ownPatch[faceI] = ownMasterPatch;
                neiPatch[faceI] = neiMasterPatch;
            }
        }
    }

    // No need to parallel sync since intersection data (surfaceIndex_ etc.)
    // already guaranteed to be synced...
    // However:
    // - owncc and neicc are reversed on different procs so might pick
    //   up different regions reversed? No problem. Neighbour on one processor
    //   might not be owner on the other processor but the neighbour is
    //   not used when creating baffles from proc faces.
    // - tolerances issues occasionally crop up.
    syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    syncTools::syncFaceList(mesh_, neiPatch, maxEqOp<label>());
}


Foam::Map<Foam::labelPair> Foam::meshRefinement::getZoneBafflePatches
(
    const bool allowBoundary,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
) const
{
    Map<labelPair> bafflePatch(mesh_.nFaces()/1000);

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    const faceZoneMesh& fZones = mesh_.faceZones();

    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            // Get zone
            label zoneI = fZones.findZoneID(faceZoneName);

            const faceZone& fZone = fZones[zoneI];

            // Get patch allocated for zone
            label globalRegionI = surfaces_.globalRegion(surfI, 0);
            labelPair zPatches
            (
                globalToMasterPatch[globalRegionI],
                globalToSlavePatch[globalRegionI]
            );

            Info<< "For zone " << fZone.name() << " found patches "
                << mesh_.boundaryMesh()[zPatches[0]].name() << " and "
                << mesh_.boundaryMesh()[zPatches[1]].name()
                << endl;

            forAll(fZone, i)
            {
                label facei = fZone[i];

                if (allowBoundary || mesh_.isInternalFace(facei))
                {
                    labelPair patches = zPatches;
                    if (fZone.flipMap()[i])
                    {
                       patches = reverse(patches);
                    }

                    if (!bafflePatch.insert(facei, patches))
                    {
                        FatalErrorInFunction
                            << "Face " << facei
                            << " fc:" << mesh_.faceCentres()[facei]
                            << " in zone " << fZone.name()
                            << " is in multiple zones!"
                            << abort(FatalError);
                    }
                }
            }
        }
    }
    return bafflePatch;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createBaffles
(
    const labelList& ownPatch,
    const labelList& neiPatch
)
{
    if
    (
        ownPatch.size() != mesh_.nFaces()
     || neiPatch.size() != mesh_.nFaces()
    )
    {
        FatalErrorInFunction
            << "Illegal size :"
            << " ownPatch:" << ownPatch.size()
            << " neiPatch:" << neiPatch.size()
            << ". Should be number of faces:" << mesh_.nFaces()
            << abort(FatalError);
    }

    if (debug)
    {
        labelList syncedOwnPatch(ownPatch);
        syncTools::syncFaceList(mesh_, syncedOwnPatch, maxEqOp<label>());
        labelList syncedNeiPatch(neiPatch);
        syncTools::syncFaceList(mesh_, syncedNeiPatch, maxEqOp<label>());

        forAll(syncedOwnPatch, facei)
        {
            if
            (
                (ownPatch[facei] == -1 && syncedOwnPatch[facei] != -1)
             || (neiPatch[facei] == -1 && syncedNeiPatch[facei] != -1)
            )
            {
                FatalErrorInFunction
                    << "Non synchronised at face:" << facei
                    << " on patch:" << mesh_.boundaryMesh().whichPatch(facei)
                    << " fc:" << mesh_.faceCentres()[facei] << endl
                    << "ownPatch:" << ownPatch[facei]
                    << " syncedOwnPatch:" << syncedOwnPatch[facei]
                    << " neiPatch:" << neiPatch[facei]
                    << " syncedNeiPatch:" << syncedNeiPatch[facei]
                    << abort(FatalError);
            }
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nBaffles = 0;

    forAll(ownPatch, facei)
    {
        if (ownPatch[facei] != -1)
        {
            // Create baffle or repatch face. Return label of inserted baffle
            // face.
            createBaffle
            (
                facei,
                ownPatch[facei],   // owner side patch
                neiPatch[facei],   // neighbour side patch
                meshMod
            );
            nBaffles++;
        }
    }
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        label coupledPatchI = -1;
        if
        (
            pp.coupled()
        && !refCast<const coupledPolyPatch>(pp).separated()
        )
        {
            coupledPatchI = patchI;
        }

        forAll(pp, i)
        {
            label faceI = pp.start()+i;

            if (ownPatch[faceI] != -1)
            {
                createBaffle
                (
                    faceI,
                    ownPatch[faceI],    // owner side patch
                    neiPatch[faceI],    // neighbour side patch
                    meshMod
                );

                if (coupledPatchI != -1)
                {
                    faceToCoupledPatch_.insert(faceI, coupledPatchI);
                }

                nBaffles++;
            }
        }
    }


    autoPtr<mapPolyMesh> map;
    if (returnReduce(nBaffles, sumOp<label>()))
    {
        // Change the mesh (no inflation, parallel sync)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }


        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

    // Pick up owner side of baffle
    forAll(ownPatch, oldFaceI)
    {
        label facei = reverseFaceMap[oldFaceI];

        if (ownPatch[oldFaceI] != -1 && facei >= 0)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[facei]];

            if (ownPatch[oldFaceI] != -1 && faceI >= 0)
            {
                const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

                forAll(ownFaces, i)
                {
                    baffledFacesSet.insert(ownFaces[i]);
                }
            }
        }
    }
    // Pick up neighbour side of baffle (added faces)
    forAll(faceMap, facei)
    {
        label oldFaceI = faceMap[facei];

        if (oldFaceI >= 0 && reverseFaceMap[oldFaceI] != facei)
        {
            const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[facei]];

            if (oldFaceI >= 0 && reverseFaceMap[oldFaceI] != faceI)
            {
                const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

                forAll(ownFaces, i)
                {
                    baffledFacesSet.insert(ownFaces[i]);
                }
            }
        }
        baffledFacesSet.sync(mesh_);

        updateMesh(map, baffledFacesSet.toc());
    }

    return map;
}


Foam::labelList Foam::meshRefinement::getZones
(
    const List<surfaceZonesInfo::faceZoneType>& fzTypes
) const
{
    const faceZoneMesh& fZones = mesh_.faceZones();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (hasInfo && findIndex(fzTypes, fzType) != -1)
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label zoneI = fZones.whichZone(facei);

                if (zoneI != -1)
                {
                    FatalErrorInFunction
                        << "face:" << facei << " on patch " << pp.name()
                        << " is in zone " << fZones[zoneI].name()
                        << exit(FatalError);
                }
            }
        }
    }

    return newBaffles;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createZoneBaffles
(
    const labelList& zoneIDs,
    List<labelPair>& baffles,
    labelList& originatingFaceZone
)
{
    autoPtr<mapPolyMesh> map;

    if (zoneIDs.size() > 0)
    {
        const faceZoneMesh& faceZones = mesh_.faceZones();

        // Split internal faces on interface surfaces
        Info<< "Converting zoned faces into baffles ..." << endl;

        // Per (internal) face the patch it should go into
        labelList ownPatch(mesh_.nFaces(), -1);
        labelList neiPatch(mesh_.nFaces(), -1);
        labelList faceZoneID(mesh_.nFaces(), -1);

        labelList nBaffles(zoneIDs.size(), 0);

        forAll(zoneIDs, j)
        {
            label zoneI = zoneIDs[j];

            const faceZone& fz = faceZones[zoneI];

            const word& masterName = faceZoneToMasterPatch_[fz.name()];
            label masterPatchI = mesh_.boundaryMesh().findPatchID(masterName);
            const word& slaveName = faceZoneToSlavePatch_[fz.name()];
            label slavePatchI = mesh_.boundaryMesh().findPatchID(slaveName);

            if (masterPatchI == -1 || slavePatchI == -1)
            {
                FatalErrorInFunction
                    << "Problem: masterPatchI:" << masterPatchI
                    << " slavePatchI:" << slavePatchI << exit(FatalError);
            }

            forAll(fz, i)
            {
                label faceI = fz[i];
                if (mesh_.isInternalFace(faceI))
                {
                    if (fz.flipMap()[i])
                    {
                        ownPatch[faceI] = slavePatchI;
                        neiPatch[faceI] = masterPatchI;
                    }
                    else
                    {
                        ownPatch[faceI] = masterPatchI;
                        neiPatch[faceI] = slavePatchI;
                    }
                    faceZoneID[faceI] = zoneI;

                    nBaffles[j]++;
                }
            }
        }

        label nLocalBaffles = sum(nBaffles);


        label nTotalBaffles = returnReduce(nLocalBaffles, sumOp<label>());

        if (nTotalBaffles > 0)
        {
            Pstream::listCombineGather(nBaffles, plusEqOp<label>());
            Pstream::listCombineScatter(nBaffles);

            Info<< nl
                << setf(ios_base::left)
                << setw(30) << "FaceZone"
                << setw(10) << "FaceType"
                << setw(10) << "nBaffles"
                << nl
                << setw(30) << "--------"
                << setw(10) << "--------"
                << setw(10) << "--------"
                << endl;

            forAll(zoneIDs, j)
            {
                label zoneI = zoneIDs[j];
                const faceZone& fz = faceZones[zoneI];

                label mpI, spI;
                surfaceZonesInfo::faceZoneType fzType;
                bool hasInfo = getFaceZoneInfo(fz.name(), mpI, spI, fzType);

                if (hasInfo)
                {
                    Info<< setf(ios_base::left)
                        << setw(30) << fz.name()
                        << setw(10)
                        << surfaceZonesInfo::faceZoneTypeNames[fzType]
                        << setw(10) << nBaffles[j]
                        << nl;
                }
            }
            Info<< endl;

            // Create baffles.
            map = createBaffles(ownPatch, neiPatch);

            // Get pairs of faces created.
            // Just loop over faceMap and store baffle if we encounter a slave
            // face.

            baffles.setSize(nLocalBaffles);
            originatingFaceZone.setSize(nLocalBaffles);
            label baffleI = 0;

            const labelList& faceMap = map().faceMap();
            const labelList& reverseFaceMap = map().reverseFaceMap();

            forAll(faceMap, facei)
            {
                label oldFaceI = faceMap[facei];

                // Does face originate from face-to-patch
                Map<labelPair>::const_iterator iter = faceToPatch.find
                (
                    oldFaceI
                );

                if (iter != faceToPatch.end())
                {
                    label masterFaceI = reverseFaceMap[oldFaceI];
                    if (facei != masterFaceI)
                    {
                        baffles[baffleI++] = labelPair(masterFaceI, facei);
                    }
                }
            }

            if (baffleI != baffles.size())
            {
                FatalErrorInFunction
                    << "Had " << faceToPatch.size() << " patches to create "
                    << " but encountered " << baffleI
                    << " slave faces originating from patcheable faces."
                    << abort(FatalError);
            }

            if (debug&MESH)
            {
                const_cast<Time&>(mesh_.time())++;
                Pout<< "Writing zone-baffled mesh to time " << timeName()
                    << endl;
                write
                (
                    debugType(debug),
                    writeType(writeLevel() | WRITEMESH),
                    mesh_.time().path()/"baffles"
                );
            }
        }
        Info<< "Created " << nTotalBaffles << " baffles in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
    return map;
}


Foam::List<Foam::labelPair> Foam::meshRefinement::freeStandingBaffles
(
    const List<labelPair>& couples,
    const scalar planarAngle
) const
{
    // Done by counting the number of baffles faces per mesh edge. If edge
    // has 2 boundary faces and both are baffle faces it is the edge of a baffle
    // region.

    // All duplicate faces on edge of the patch are to be merged.
    // So we count for all edges of duplicate faces how many duplicate
    // faces use them.
    labelList nBafflesPerEdge(mesh_.nEdges(), 0);


    // This algorithm is quite tricky. We don't want to use edgeFaces and
    // also want it to run in parallel so it is now an algorithm over
    // all (boundary) faces instead.
    // We want to pick up any edges that are only used by the baffle
    // or internal faces but not by any other boundary faces. So
    // - increment count on an edge by 1 if it is used by any (uncoupled)
    //   boundary face.
    // - increment count on an edge by 1000000 if it is used by a baffle face
    // - sum in parallel
    //
    // So now any edge that is used by baffle faces only will have the
    // value 2*1000000+2*1.


    const label baffleValue = 1000000;


    // Count number of boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        // Count number of boundary faces. Discard coupled boundary faces.
        if (!pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                const labelList& fEdges = mesh_.faceEdges(facei);

                forAll(fEdges, fEdgeI)
                {
                    nBafflesPerEdge[fEdges[fEdgeI]]++;
                }
                facei++;
            }
        }
    }


    DynamicList<label> fe0;
    DynamicList<label> fe1;


    // Count number of duplicate boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(couples, i)
    {
        {
            label f0 = couples[i].first();
            const labelList& fEdges0 = mesh_.faceEdges(f0, fe0);
            forAll(fEdges0, fEdgeI)
            {
                nBafflesPerEdge[fEdges0[fEdgeI]] += baffleValue;
            }
        }

        {
            label f1 = couples[i].second();
            const labelList& fEdges1 = mesh_.faceEdges(f1, fe1);
            forAll(fEdges1, fEdgeI)
            {
                nBafflesPerEdge[fEdges1[fEdgeI]] += baffleValue;
            }
        }
    }

    // Add nBaffles on shared edges
    syncTools::syncEdgeList
    (
        mesh_,
        nBafflesPerEdge,
        plusEqOp<label>(),  // in-place add
        label(0)            // initial value
    );


    // Baffles which are not next to other boundaries and baffles will have
    // nBafflesPerEdge value 2*baffleValue+2*1 (from 2 boundary faces which
    // are both baffle faces)

    List<labelPair> filteredCouples(couples.size());
    label filterI = 0;

    forAll(couples, i)
    {
        const labelPair& couple = couples[i];

        if
        (
            patches.whichPatch(couple.first())
         == patches.whichPatch(couple.second())
        )
        {
            const labelList& fEdges = mesh_.faceEdges(couple.first());

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (nBafflesPerEdge[edgeI] == 2*baffleValue+2*1)
                {
                    filteredCouples[filterI++] = couple;
                    break;
                }
            }
        }
    }
    filteredCouples.setSize(filterI);


    label nFiltered = returnReduce(filteredCouples.size(), sumOp<label>());

    Info<< "freeStandingBaffles : detected "
        << nFiltered
        << " free-standing baffles out of "
        << returnReduce(couples.size(), sumOp<label>())
        << nl << endl;


    if (nFiltered > 0)
    {
        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(filteredCouples.size());
        pointField end(filteredCouples.size());

        const pointField& cellCentres = mesh_.cellCentres();

        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];
            start[i] = cellCentres[mesh_.faceOwner()[couple.first()]];
            end[i] = cellCentres[mesh_.faceOwner()[couple.second()]];
        }

        // Extend segments a bit
        {
            const vectorField smallVec(ROOTSMALL*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do test for intersections
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;

        surfaces_.findNearestIntersection
        (
            identity(surfaces_.surfaces().size()),
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );

        //mkDir(mesh_.time().path()/timeName());
        //OBJstream str
        //(
        //    mesh_.time().path()/timeName()/"flatBaffles.obj"
        //);

        const scalar planarAngleCos = Foam::cos(degToRad(planarAngle));

        label filterI = 0;
        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];

            // Note: for a baffle-surface we do not want to merge the baffle.
            // We could either check for hitting the same triangle (but you
            // might hit same point on neighbouring triangles due to tolerance
            // issues) or better just to compare the hit point.
            // This might still go wrong for a ray in the plane of the triangle
            // which might hit two different points on the same triangle due
            // to tolerances...

            if
            (
                hit1[i].hit()
             && hit2[i].hit()
             && mag(hit1[i].hitPoint()-hit2[i].hitPoint()) > mergeDistance_
            )
            {
                // Two different hits. Check angle.
                //str.write
                //(
                //    linePointRef(hit1[i].hitPoint(), hit2[i].hitPoint()),
                //    normal1[i],
                //    normal2[i]
                //);

                if ((normal1[i]&normal2[i]) > planarAngleCos)
                {
                    // Both normals aligned
                    vector n = end[i]-start[i];
                    scalar magN = mag(n);
                    if (magN > VSMALL)
                    {
                        filteredCouples[filterI++] = couple;
                    }
                }
            }
            else if (hit1[i].hit() || hit2[i].hit())
            {
                // Single hit. Do not include in freestanding baffles.
            }
        }

        filteredCouples.setSize(filterI);

        Info<< "freeStandingBaffles : detected "
            << returnReduce(filterI, sumOp<label>())
            << " planar (within " << planarAngle
            << " degrees) free-standing baffles out of "
            << nFiltered
            << nl << endl;
    }

    return filteredCouples;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeBaffles
(
    const List<labelPair>& couples,
    const Map<label>& faceToPatch
)
{
    autoPtr<mapPolyMesh> map;

    if (returnReduce(couples.size()+faceToPatch.size(), sumOp<label>()))
    {
        // Mesh change engine
        polyTopoChange meshMod(mesh_);

        const faceList& faces = mesh_.faces();
        const labelList& faceOwner = mesh_.faceOwner();
        const faceZoneMesh& faceZones = mesh_.faceZones();

        forAll(couples, i)
        {
            label face0 = couples[i].first();
            label face1 = couples[i].second();

            // face1 < 0 signals a coupled face that has been converted to
            // baffle

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (face1 < 0 || own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = faceZones.whichZone(face0);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                }

                label nei = (face1 < 0 ? -1 : own1);

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        nei,                    // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
            else
            {
                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(face1);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                }

                meshMod.setAction(polyRemoveFace(face0));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face1],           // modified face
                        face1,                  // label of face being modified
                        own1,                   // owner
                        own0,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }

        forAllConstIter(Map<label>, faceToPatch, iter)
        {
            label faceI = iter.key();
            label patchI = iter();

            if (!mesh_.isInternalFace(faceI))
            {
                FatalErrorInFunction
                    << "problem: face:" << faceI
                    << " at:" << mesh_.faceCentres()[faceI]
                    << "(wanted patch:" << patchI
                    << ") is an internal face" << exit(FatalError);
            }

            label zoneID = faceZones.whichZone(faceI);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    faces[faceI],           // modified face
                    faceI,                  // label of face being modified
                    faceOwner[faceI],       // owner
                    -1,                     // neighbour
                    false,                  // face flip
                    patchI,                 // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }


        // Change the mesh (no inflation)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Recalculate intersections on merged faces since
        // this seems to give problems? Note: should not be necessary since
        // baffles preserve intersections from when they were created.
        labelList newExposedFaces(2*couples.size());
        label newI = 0;

        forAll(couples, i)
        {
            label newFace0 = map().reverseFaceMap()[couples[i].first()];
            if (newFace0 != -1)
            {
                newExposedFaces[newI++] = newFace0;
            }

            label newFace1 = map().reverseFaceMap()[couples[i].second()];
            if (newFace1 != -1)
            {
                newExposedFaces[newI++] = newFace1;
            }
        }
        newExposedFaces.setSize(newI);
        updateMesh(map, newExposedFaces);
    }

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeZoneBaffles
(
    const bool doInternalZones,
    const bool doBaffleZones
)
{
    labelList zoneIDs;
    {
        DynamicList<surfaceZonesInfo::faceZoneType> fzTypes;
        if (doInternalZones)
        {
            fzTypes.append(surfaceZonesInfo::INTERNAL);
        }
        if (doBaffleZones)
        {
            fzTypes.append(surfaceZonesInfo::BAFFLE);
        }
        zoneIDs = getZones(fzTypes);
    }

    List<labelPair> zoneBaffles
    (
        subsetBaffles
        (
            mesh_,
            zoneIDs,
            localPointRegion::findDuplicateFacePairs(mesh_)
        )
    );

    autoPtr<mapPolyMesh> mapPtr;
    if (returnReduce(zoneBaffles.size(), sumOp<label>()))
    {
        mapPtr = mergeBaffles(zoneBaffles, Map<label>(0));
    }
    return mapPtr;
}


// Finds region per cell for cells inside closed named surfaces
void Foam::meshRefinement::findCellZoneGeometric
(
    const pointField& neiCc,
    const labelList& closedNamedSurfaces,   // indices of closed surfaces
    labelList& namedSurfaceIndex,           // per face index of named surface
    const labelList& surfaceToCellZone,     // cell zone index per surface

    labelList& cellToZone
) const
{
    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Check if cell centre is inside
    labelList insideSurfaces;
    surfaces_.findInside
    (
        closedNamedSurfaces,
        cellCentres,
        insideSurfaces
    );

    forAll(insideSurfaces, celli)
    {
        if (cellToZone[celli] == -2)
        {
            label surfI = insideSurfaces[celli];

        if (surfI != -1)
        {
            if (cellToZone[cellI] == -2)
            {
                cellToZone[celli] = surfaceToCellZone[surfI];
            }
            else if (cellToZone[cellI] == -1)
            {
                cellToZone[celli] = surfaceToCellZone[surfI];
            }
        }
    }


    // Some cells with cell centres close to surface might have
    // had been put into wrong surface. Recheck with perturbed cell centre.


    // 1. Collect points

    // Count points to test.
    label nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            if (mesh_.isInternalFace(facei))
            {
                nCandidates += 2;
            }
            else
            {
                nCandidates += 1;
            }
        }
    }

    // Collect points.
    pointField candidatePoints(nCandidates);
    nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            label own = faceOwner[facei];
            const point& ownCc = cellCentres[own];

            if (mesh_.isInternalFace(facei))
            {
                label nei = faceNeighbour[facei];
                const point& neiCc = cellCentres[nei];
                // Perturbed cc
                const vector d = 1e-4*(neiCc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
                candidatePoints[nCandidates++] = neiCc+d;
            }
            else
            {
                //const point& neiFc = mesh_.faceCentres()[facei];
                const point& neiFc = neiCc[facei-mesh_.nInternalFaces()];

                // Perturbed cc
                const vector d = 1e-4*(neiFc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
            }
        }
    }


    // 2. Test points for inside

    surfaces_.findInside
    (
        closedNamedSurfaces,
        candidatePoints,
        insideSurfaces
    );


    // 3. Update zone information

    nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            label own = faceOwner[facei];

            if (mesh_.isInternalFace(facei))
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }

                label neiSurfI = insideSurfaces[nCandidates++];
                if (neiSurfI != -1)
                {
                    label nei = faceNeighbour[facei];

                    cellToZone[nei] = surfaceToCellZone[neiSurfI];
                }
            }
            else
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }
            }
        }
    }


    // Adapt the namedSurfaceIndex
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // for if any cells were not completely covered.

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label ownZone = cellToZone[mesh_.faceOwner()[facei]];
        label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];

        if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
        {
            // Give face the zone of max cell zone
            namedSurfaceIndex[facei] = findIndex
            (
                surfaceToCellZone,
                max(ownZone, neiZone)
            );
        }
    }

            label minZone;
            if (ownZone == -1)
            {
                minZone = neiZone;
            }
            else if (neiZone == -1)
            {
                minZone = ownZone;
            }
            else
            {
                minZone = min(ownZone, neiZone);
            }

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

            if (geomSurfI != -1)
            {
                label facei = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                neiCellZone[facei-mesh_.nInternalFaces()] = ownZone;
            }
        }
    }

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
                {
                    // Give face the max cell zone
                    namedSurfaceIndex[facei] = findIndex
                    (
                        surfaceToCellZone,
                        max(ownZone, neiZone)
                    );
                }
            }
        }
    }

    // Sync
    syncTools::syncFaceList(mesh_, namedSurfaceIndex, maxEqOp<label>());
}


void Foam::meshRefinement::findCellZoneInsideWalk
(
    const pointField& locationsInMesh,
    const labelList& zonesInMesh,
    const labelList& faceToZone, // per face -1 or some index >= 0

    labelList& cellToZone
) const
{
    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());
    //selectSeparatedCoupledFaces(blockedFace);

    forAll(namedSurfaceIndex, facei)
    {
        if (namedSurfaceIndex[facei] == -1)
        {
            blockedFace[facei] = false;
        }
        else
        {
            blockedFace[facei] = true;
        }
    }
    // No need to sync since faceToZone already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];
        label zoneID = zonesInMesh[i];

        // Find the region containing the insidePoint
        label keepRegionI = findRegion
        (
            mesh_,
            cellRegion,
            mergeDistance_*vector(1,1,1),
            insidePoint
        );

        Info<< "For cellZone "
            <<  (
                    zoneID == -1
                  ? "none"
                  : mesh_.cellZones()[zoneID].name()
                )
            << " found point " << insidePoint
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        // Set all cells with this region
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == keepRegionI)
            {
                if (cellToZone[celli] == -2)
                {
                    cellToZone[celli] = surfaceToCellZone[surfI];
                }
                else if (cellToZone[celli] != surfaceToCellZone[surfI])
                {
                    WarningInFunction
                        << "Cell " << celli
                        << " at " << mesh_.cellCentres()[celli]
                        << " is inside surface " << surfaces_.names()[surfI]
                        << " but already marked as being in zone "
                        << cellToZone[celli] << endl
                        << "This can happen if your surfaces are not"
                        << " (sufficiently) closed."
                        << endl;
                }
            }
        }
    }
}


void Foam::meshRefinement::findCellZoneInsideWalk
(
    const pointField& locationsInMesh,
    const wordList& zoneNamesInMesh,
    const labelList& faceToZone,    // per face -1 or some index >= 0

    labelList& cellToZone
) const
{
    const cellZoneMesh& czs = mesh_.cellZones();

    labelList zoneIDs(zoneNamesInMesh.size());
    forAll(zoneNamesInMesh, i)
    {
        zoneIDs[i] = czs.findZoneID(zoneNamesInMesh[i]);
    }
    findCellZoneInsideWalk
    (
        locationsInMesh,
        zoneIDs,
        faceToZone,
        cellToZone
    );
}


bool Foam::meshRefinement::calcRegionToZone
(
    const label backgroundZoneID,
    const label surfZoneI,
    const label ownRegion,
    const label neiRegion,

    labelList& regionToCellZone
) const
{
    bool changed = false;

    // Check whether inbetween different regions
    if (ownRegion != neiRegion)
    {
        // Jump. Change one of the sides to my type.

        // 1. Interface between my type and unset region.
        // Set region to keepRegion

        if (regionToCellZone[ownRegion] == -2)
        {
            if (surfZoneI == -1)
            {
                // Special: face is -on faceZone  -not real boundary
                //          -not on cellZone
                // so make regions same on either side
                if (regionToCellZone[neiRegion] != -2)
                {
                    regionToCellZone[ownRegion] = regionToCellZone[neiRegion];
                    changed = true;
                }
            }
            else if (regionToCellZone[neiRegion] == surfZoneI)
            {
                // Face between unset and my region. Put unset
                // region into background region
                //MEJ: see comment in findCellZoneTopo
                if (backgroundZoneID != -2)
                {
                    regionToCellZone[ownRegion] = backgroundZoneID;
                    changed = true;
                }
            }
            else if (regionToCellZone[neiRegion] != -2)
            {
                // Face between unset and other region.
                // Put unset region into my region
                regionToCellZone[ownRegion] = surfZoneI;
                changed = true;
            }
        }
        else if (regionToCellZone[neiRegion] == -2)
        {
            if (surfZoneI == -1)
            {
                // Special: face is -on faceZone  -not real boundary
                //          -not on cellZone
                // so make regions same on either side
                regionToCellZone[neiRegion] = regionToCellZone[ownRegion];
                changed = true;
            }
            else if (regionToCellZone[ownRegion] == surfZoneI)
            {
                // Face between unset and my region. Put unset
                // region into background region
                if (backgroundZoneID != -2)
                {
                    regionToCellZone[neiRegion] = backgroundZoneID;
                    changed = true;
                }
            }
            else if (regionToCellZone[ownRegion] != -2)
            {
                // Face between unset and other region.
                // Put unset region into my region
                regionToCellZone[neiRegion] = surfZoneI;
                changed = true;
            }
        }
    }
    return changed;
}


void Foam::meshRefinement::findCellZoneTopo
(
    const label backgroundZoneID,
    const pointField& locationsInMesh,
    const labelList& unnamedSurfaceRegion,
    const labelList& namedSurfaceIndex,
    const labelList& surfaceToCellZone,
    labelList& cellToZone
) const
{
    // This routine fixes small problems with left over unassigned regions
    // (after all off the unreachable bits of the mesh have been removed).
    // This routine splits the mesh into regions, based on the intersection
    // with a surface. The problem is that we know the surface which
    // (intersected) face belongs to (in namedSurfaceIndex) but we don't
    // know which side of the face it relates to. So all we are doing here
    // is get the correspondence between surface/cellZone and regionSplit
    // region. See the logic in calcRegionToZone.
    // Basically it looks at the neighbours of a face on a named surface.
    // If one side is in the cellZone corresponding to the surface and
    // the other side is unassigned (-2) it sets this to the background zone.
    // So the zone to set these unmarked cells to is provided as argument:
    // - backgroundZoneID = -2 : do not change so remove cells
    // - backgroundZoneID = -1 : put into background

    // Assumes:
    // - region containing keepPoint does not go into a cellZone
    // - all other regions can be found by crossing faces marked in
    //   namedSurfaceIndex.

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());

    forAll(namedSurfaceIndex, facei)
    {
        if (namedSurfaceIndex[facei] == -1)
        {
            blockedFace[facei] = false;
        }
        else
        {
            blockedFace[facei] = true;
        }
    }
    // No need to sync since namedSurfaceIndex already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Per mesh region the zone the cell should be put in.
    // -2   : not analysed yet
    // -1   : keepPoint region. Not put into any cellzone.
    // >= 0 : index of cellZone
    labelList regionToCellZone(cellRegion.nRegions(), -2);

    // See which cells already are set in the cellToZone (from geometric
    // searching) and use these to take over their zones.
    // Note: could be improved to count number of cells per region.
    forAll(cellToZone, celli)
    {
        if (cellToZone[celli] != -2)
        {
            regionToCellZone[cellRegion[celli]] = cellToZone[celli];
        }
    }

    // Synchronise regionToCellZone.
    // Note:
    // - region numbers are identical on all processors
    // - keepRegion is identical ,,
    // - cellZones are identical ,,
    Pstream::listCombineGather(regionToCellZone, maxEqOp<label>());
    Pstream::listCombineScatter(regionToCellZone);


    // Find the region containing the keepPoint
    forAll(locationsInMesh, i)
    {
        FatalErrorInFunction
            << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.bounds()
            << exit(FatalError);
    }

        Info<< "Found point " << keepPoint
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << keepPoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        // Mark default region with zone -1. Note that all regions should
        // already be matched to a cellZone through the loop over cellToZone.
        // This is just to mop up any remaining regions. Not sure whether
        // this is needed?
        if (regionToCellZone[keepRegionI] == -2)
        {
            regionToCellZone[keepRegionI] = -1;
        }
    }


    // Find correspondence between cell zone and surface
    // by changing cell zone every time we cross a surface.
    while (true)
    {
        // Synchronise regionToCellZone.
        // Note:
        // - region numbers are identical on all processors
        // - keepRegion is identical ,,
        // - cellZones are identical ,,
        // This done at top of loop to account for geometric matching
        // not being synchronised.
        Pstream::listCombineGather(regionToCellZone, maxEqOp<label>());
        Pstream::listCombineScatter(regionToCellZone);


        bool changed = false;

        // Internal faces

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label surfI = namedSurfaceIndex[facei];

            // Connected even if no cellZone defined for surface
            if (unnamedSurfaceRegion[faceI] == -1 && surfI != -1)
            {
                // Calculate region to zone from cellRegions on either side
                // of internal face.
                bool changedCell = calcRegionToZone
                (
                    backgroundZoneID,
                    surfaceToCellZone[surfI],
                    cellRegion[mesh_.faceOwner()[facei]],
                    cellRegion[mesh_.faceNeighbour()[facei]],
                    regionToCellZone
                );

                changed = changed | changedCell;
            }
        }

        // Boundary faces

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Get coupled neighbour cellRegion
        labelList neiCellRegion(mesh_.nFaces()-mesh_.nInternalFaces());
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    neiCellRegion[facei-mesh_.nInternalFaces()] =
                        cellRegion[mesh_.faceOwner()[facei]];
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiCellRegion);

        // Calculate region to zone from cellRegions on either side of coupled
        // face.
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;

                    label surfI = namedSurfaceIndex[facei];

                    // Connected even if no cellZone defined for surface
                    if (unnamedSurfaceRegion[faceI] == -1 && surfI != -1)
                    {
                        bool changedCell = calcRegionToZone
                        (
                            backgroundZoneID,
                            surfaceToCellZone[surfI],
                            cellRegion[mesh_.faceOwner()[facei]],
                            neiCellRegion[facei-mesh_.nInternalFaces()],
                            regionToCellZone
                        );

                        changed = changed | changedCell;
                    }
                }
            }
        }

        if (!returnReduce(changed, orOp<bool>()))
        {
            break;
        }
    }


    if (debug)
    {
        forAll(regionToCellZone, regionI)
        {
            Pout<< "Region " << regionI
                << " becomes cellZone:" << regionToCellZone[regionI]
                << endl;
        }
    }

    // Rework into cellToZone
    forAll(cellToZone, celli)
    {
        cellToZone[celli] = regionToCellZone[cellRegion[celli]];
    }
}


void Foam::meshRefinement::makeConsistentFaceIndex
(
    const labelList& surfaceMap,
    const labelList& cellToZone,
    labelList& namedSurfaceIndex
) const
{
    // Make namedSurfaceIndex consistent with cellToZone - clear out any
    // blocked faces inbetween same cell zone (or background (=-1))
    // Do not do this for surfaces relating to 'pure' faceZones i.e.
    // faceZones without a cellZone. Note that we cannot check here
    // for different cellZones on either side but no namedSurfaceIndex
    // since cellZones can now originate from locationsInMesh as well
    // (instead of only through named surfaces)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label ownZone = cellToZone[faceOwner[facei]];
        label neiZone = cellToZone[faceNeighbour[facei]];

        if (ownZone == neiZone && namedSurfaceIndex[facei] != -1)
        {
            namedSurfaceIndex[facei] = -1;
        }
        else if (ownZone != neiZone && namedSurfaceIndex[facei] == -1)
        {
            FatalErrorInFunction
                << "Different cell zones on either side of face " << facei
                << " at " << mesh_.faceCentres()[facei]
                << " but face not marked with a surface."
                << abort(FatalError);
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get coupled neighbour cellZone
    labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces());
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                neiCellZone[facei-mesh_.nInternalFaces()] =
                    cellToZone[mesh_.faceOwner()[facei]];
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

    // Use coupled cellZone to do check
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                if (ownZone == neiZone && namedSurfaceIndex[facei] != -1)
                {
                    namedSurfaceIndex[facei] = -1;
                }
                else if (ownZone != neiZone && namedSurfaceIndex[facei] == -1)
                {
                    FatalErrorInFunction
                        << "Different cell zones on either side of face "
                        << facei << " at " << mesh_.faceCentres()[facei]
                        << " but face not marked with a surface."
                        << abort(FatalError);
                }
            }
        }
        else
        {
            // Unzonify boundary faces
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                namedSurfaceIndex[facei] = -1;
            }

            makeConsistentFaceIndex
            (
                zoneToNamedSurface,
                cellToZone,
                namedSurfaceIndex
            );
        }
    }
}


void Foam::meshRefinement::getIntersections
(
    const labelList& surfacesToTest,
    const pointField& neiCc,
    const labelList& testFaces,

    labelList& namedSurfaceIndex,
    PackedBoolList& posOrientation
) const
{
    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    namedSurfaceIndex.setSize(mesh_.nFaces());
    namedSurfaceIndex = -1;

    posOrientation.setSize(mesh_.nFaces());
    posOrientation = false;

    // Statistics: number of faces per faceZone
    labelList nSurfFaces(surfaces_.surfZones().size(), 0);

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        forAll(facePatch, facei)
        {
            if (facePatch[facei] != -1)
            {
                problemFaces.insert(facei);
            }
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(ROOTSMALL*(end-start));
        start -= smallVec;
        end += smallVec;
    }


    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note that we intersect all intersected faces again. Could reuse
    // the information already in surfaceIndex_.

    labelList surface1;
    List<pointIndexHit> hit1;
    vectorField normal1;
    labelList surface2;
    List<pointIndexHit> hit2;
    vectorField normal2;
    {
        labelList region1;
        labelList region2;
        surfaces_.findNearestIntersection
        (
            surfacesToTest,
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );
    }

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];
        const vector& area = mesh_.faceAreas()[faceI];

        if (surface1[i] != -1)
        {
            // If both hit should probably choose 'nearest'
            if
            (
                surface2[i] != -1
             && (
                    magSqr(hit2[i].hitPoint())
                  < magSqr(hit1[i].hitPoint())
                )
            )
            {
                namedSurfaceIndex[faceI] = surface2[i];
                posOrientation[faceI] = ((area&normal2[i]) > 0);
                nSurfFaces[surface2[i]]++;
            }
            else
            {
                namedSurfaceIndex[faceI] = surface1[i];
                posOrientation[faceI] = ((area&normal1[i]) > 0);
                nSurfFaces[surface1[i]]++;
            }
        }
        else if (surface2[i] != -1)
        {
            namedSurfaceIndex[faceI] = surface2[i];
            posOrientation[faceI] = ((area&normal2[i]) > 0);
            nSurfFaces[surface2[i]]++;
        }
    }


    // surfaceIndex might have different surfaces on both sides if
    // there happen to be a (obviously thin) surface with different
    // regions between the cell centres. If one is on a named surface
    // and the other is not this might give problems so sync.
    syncTools::syncFaceList
    (
        mesh_,
        namedSurfaceIndex,
        maxEqOp<label>()
    );

    // Print a bit
    if (debug)
    {
        forAll(nSurfFaces, surfI)
        {
            Pout<< "Surface:"
                << surfaces_.names()[surfI]
                << "  nZoneFaces:" << nSurfFaces[surfI] << nl;
        }
        Pout<< endl;
    }
}


void Foam::meshRefinement::zonify
(
    const bool allowFreeStandingZoneFaces,
    const label backgroundZoneID,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,

    labelList& cellToZone,
    labelList& unnamedRegion1,
    labelList& unnamedRegion2,
    labelList& namedSurfaceIndex,
    PackedBoolList& posOrientation
) const
{
    // Determine zones for cells and faces
    // cellToZone:
    // -2  : unset
    // -1  : not in any zone (zone 'none' or background zone)
    // >=0 : zoneID
    // namedSurfaceIndex, posOrientation:
    // -1  : face not intersected by named surface
    // >=0 : index of named surface
    //       (and posOrientation: surface normal v.s. face normal)

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);


    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));
    labelList unnamedSurfaces(surfaceZonesInfo::getUnnamedSurfaces(surfZones));

    // Get map from surface to cellZone (or -1)
    labelList surfaceToCellZone;
    if (namedSurfaces.size())
    {
        // Get/add cellZones corresponding to surface names
        surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
    }


    cellToZone.setSize(mesh_.nCells());
    cellToZone = -2;

    namedSurfaceIndex.clear();
    posOrientation.clear();



    // 1. Test all (unnamed & named) surfaces

    // Unnamed surfaces. Global surface region of intersection (or -1)
    getIntersections
    (
        unnamedSurfaces,
        neiCc,
        intersectedFaces(),
        unnamedRegion1,
        unnamedRegion2
    );


    if (namedSurfaces.size())
    {
        getIntersections
        (
            namedSurfaces,
            neiCc,
            intersectedFaces(),
            namedSurfaceIndex,
            posOrientation
        );
    }


    // 2. Walk from locationsInMesh. Hard set cellZones.
    //    Note: walk through faceZones! (these might get overridden later)

    if (locationsInMesh.size())
    {
        Info<< "Setting cellZones according to locationsInMesh:" << endl;

        labelList locationsZoneIDs(zonesInMesh.size(), -1);
        forAll(locationsInMesh, i)
        {
            const word& name = zonesInMesh[i];

            Info<< "Location : " << locationsInMesh[i] << nl
                << "    cellZone : " << name << endl;

            if (name != "none")
            {
                label zoneID = mesh_.cellZones().findZoneID(name);
                if (zoneID == -1)
                {
                    FatalErrorInFunction << "problem" << abort(FatalError);
                }
                locationsZoneIDs[i] = zoneID;
            }
        }
        Info<< endl;


        // Assign cellZone according to seed points. Walk through named surfaces
        findCellZoneInsideWalk
        (
            locationsInMesh,    // locations
            locationsZoneIDs,   // index of cellZone (or -1)
            unnamedRegion1,     // per face -1 (unblocked) or >= 0 (blocked)
            cellToZone
        );
    }


    // 3. Mark named-surfaces-with-insidePoint. Hard set cellZones.

    labelList locationSurfaces
    (
        surfaceZonesInfo::getInsidePointNamedSurfaces(surfZones)
    );

    if (locationSurfaces.size())
    {
        Info<< "Found " << locationSurfaces.size()
            << " named surfaces with a provided inside point."
            << " Assigning cells inside these surfaces"
            << " to the corresponding cellZone."
            << nl << endl;

        // Collect per surface the -insidePoint -cellZone
        pointField insidePoints(locationSurfaces.size());
        labelList insidePointCellZoneIDs(locationSurfaces.size(), -1);
        forAll(locationSurfaces, i)
        {
            label surfI = locationSurfaces[i];
            insidePoints[i] = surfZones[surfI].zoneInsidePoint();

            const word& name = surfZones[surfI].cellZoneName();
            if (name != "none")
            {
                label zoneID = mesh_.cellZones().findZoneID(name);
                if (zoneID == -1)
                {
                    FatalErrorInFunction
                        << "problem"
                        << abort(FatalError);
                }
                insidePointCellZoneIDs[i] = zoneID;
            }
        }


        // Stop at unnamed or named surface
        labelList allRegion1(mesh_.nFaces(), -1);
        forAll(namedSurfaceIndex, faceI)
        {
            allRegion1[faceI] = max
            (
                unnamedRegion1[faceI],
                namedSurfaceIndex[faceI]
            );
        }

        findCellZoneInsideWalk
        (
            insidePoints,           // locations
            insidePointCellZoneIDs, // index of cellZone
            allRegion1,             // per face -1 (unblocked) or >= 0 (blocked)
            cellToZone
        );
    }


    // 4. Mark named-surfaces-with-geometric faces. Do geometric test. Soft set
    // cellZones. Correct through making consistent.

    // Closed surfaces with cellZone specified.
    labelList closedNamedSurfaces
    (
        surfaceZonesInfo::getClosedNamedSurfaces
        (
            surfZones,
            surfaces_.geometry(),
            surfaces_.surfaces()
        )
    );

    if (closedNamedSurfaces.size())
    {
        Info<< "Found " << closedNamedSurfaces.size()
            << " closed, named surfaces. Assigning cells in/outside"
            << " these surfaces to the corresponding cellZone."
            << nl << endl;

        findCellZoneGeometric
        (
            neiCc,
            closedNamedSurfaces,    // indices of closed surfaces
            namedSurfaceIndex,      // per face index of named surface
            surfaceToCellZone,      // cell zone index per surface

            cellToZone
        );
    }


    // 5. Find any unassigned regions (from regionSplit)

    if (namedSurfaces.size())
    {
        Info<< "Walking from known cellZones; crossing a faceZone "
            << "face changes cellZone" << nl << endl;

        findCellZoneTopo
        (
            backgroundZoneID,
            pointField(0),
            unnamedRegion1,      // Intersections with unnamed surfaces
            namedSurfaceIndex,   // Intersections with named surfaces
            surfaceToCellZone,
            cellToZone
        );

        // Make sure namedSurfaceIndex is unset inbetween same cell zones.
        if (!allowFreeStandingZoneFaces)
        {
            Info<< "Only keeping zone faces inbetween different cellZones."
                << nl << endl;

            // Surfaces with faceZone but no cellZone
            labelList standaloneNamedSurfaces
            (
                surfaceZonesInfo::getStandaloneNamedSurfaces
                (
                    surfZones
                )
            );

            // Construct map from surface index to index in
            // standaloneNamedSurfaces (or -1)
            labelList surfaceMap(surfZones.size(), -1);
            forAll(standaloneNamedSurfaces, i)
            {
                surfaceMap[standaloneNamedSurfaces[i]] = i;
            }

            makeConsistentFaceIndex
            (
                surfaceMap,
                cellToZone,
                namedSurfaceIndex
            );
        }
    }


    // Some stats
    if (debug)
    {
        label nZones = gMax(cellToZone)+1;

        label nUnvisited = 0;
        label nBackgroundCells = 0;
        labelList nZoneCells(nZones, 0);
        forAll(cellToZone, cellI)
        {
            label zoneI = cellToZone[cellI];
            if (zoneI >= 0)
            {
                nZoneCells[zoneI]++;
            }
            else if (zoneI == -1)
            {
                nBackgroundCells++;
            }
            else if (zoneI == -2)
            {
                nUnvisited++;
            }
            else
            {
                FatalErrorInFunction
                    << "problem" << exit(FatalError);
            }
        }
        reduce(nUnvisited, sumOp<label>());
        reduce(nBackgroundCells, sumOp<label>());
        forAll(nZoneCells, zoneI)
        {
            reduce(nZoneCells[zoneI], sumOp<label>());
        }
        Info<< "nUnvisited      :" << nUnvisited << endl;
        Info<< "nBackgroundCells:" << nBackgroundCells << endl;
        Info<< "nZoneCells      :" << nZoneCells << endl;
    }
    if (debug&MESH)
    {
        const_cast<Time&>(mesh_.time())++;
        Pout<< "Writing cell zone allocation on mesh to time "
            << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            mesh_.time().path()/"cell2Zone"
        );
        volScalarField volCellToZone
        (
            IOobject
            (
                "cellToZone",
                timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(cellToZone, cellI)
        {
            volCellToZone[cellI] = cellToZone[cellI];
        }
        volCellToZone.write();
    }
}


void Foam::meshRefinement::handleSnapProblems
(
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
{
    Info<< nl
        << "Introducing baffles to block off problem cells" << nl
        << "----------------------------------------------" << nl
        << endl;

    labelList facePatch;
    if (useTopologicalSnapDetection)
    {
        facePatch = markFacesOnProblemCells
        (
            motionDict,
            removeEdgeConnectedCells,
            perpendicularAngle,
            globalToMasterPatch
        );
    }
    else
    {
        facePatch = markFacesOnProblemCellsGeometric
        (
            snapParams,
            motionDict,
            globalToMasterPatch,
            globalToSlavePatch
        );
    }
    Info<< "Analyzed problem cells in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug&MESH)
    {
        faceSet problemFaces(mesh_, "problemFaces", mesh_.nFaces()/100);

        forAll(facePatch, facei)
        {
            if (facePatch[facei] != -1)
            {
                problemFaces.insert(facei);
            }
        }
        problemFaces.instance() = timeName();
        Pout<< "Dumping " << problemFaces.size()
            << " problem faces to " << problemFaces.objectPath() << endl;
        problemFaces.write();
    }

    Info<< "Introducing baffles to delete problem cells." << nl << endl;

    if (debug)
    {
        runTime++;
    }

    // Create baffles with same owner and neighbour for now.
    createBaffles(facePatch, facePatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles");

    if (debug&MESH)
    {
        Pout<< "Writing extra baffled mesh to time "
            << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"extraBaffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


Foam::labelList Foam::meshRefinement::freeStandingBaffleFaces
(
    const labelList& faceToZone,
    const labelList& cellToZone,
    const labelList& neiCellZone
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    // We want to pick up the faces to orient. These faces come in
    // two variants:
    // - faces originating from stand-alone faceZones
    //   (these will most likely have no cellZone on either side so
    //    ownZone and neiZone both -1)
    // - sticky-up faces originating from a 'bulge' in a outside of
    //   a cellZone. These will have the same cellZone on either side.
    //   How to orient these is not really clearly defined so do them
    //   same as stand-alone faceZone faces for now. (Normally these will
    //   already have been removed by the 'allowFreeStandingZoneFaces=false'
    //   default setting)

    // Note that argument neiCellZone will have -1 on uncoupled boundaries.

    DynamicList<label> faceLabels(mesh_.nFaces()/100);

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceToZone[facei] != -1)
        {
            // Free standing baffle?
            label ownZone = cellToZone[faceOwner[facei]];
            label neiZone = cellToZone[faceNeighbour[facei]];
            if (ownZone == neiZone)
            {
                faceLabels.append(facei);
            }
        }
    }
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, i)
        {
            label facei = pp.start()+i;
            if (faceToZone[facei] != -1)
            {
                // Free standing baffle?
                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                if (ownZone == neiZone)
                {
                    faceLabels.append(facei);
                }
            }
        }
    }
    return faceLabels.shrink();
}


void Foam::meshRefinement::calcPatchNumMasterFaces
(
    const PackedBoolList& isMasterFace,
    const indirectPrimitivePatch& patch,
    labelList& nMasterFacesPerEdge
) const
{
    // Number of (master)faces per edge
    nMasterFacesPerEdge.setSize(patch.nEdges());
    nMasterFacesPerEdge = 0;

    forAll(patch.addressing(), facei)
    {
        const label meshFaceI = patch.addressing()[facei];

        if (isMasterFace[meshFaceI])
        {
            const labelList& fEdges = patch.faceEdges()[facei];
            forAll(fEdges, fEdgeI)
            {
                nMasterFacesPerEdge[fEdges[fEdgeI]]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        patch.meshEdges(mesh_.edges(), mesh_.pointEdges()),
        nMasterFacesPerEdge,
        plusEqOp<label>(),
        label(0)
    );
}


Foam::label Foam::meshRefinement::markPatchZones
(
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    labelList& faceToZone
) const
{
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());


    // Protect all non-manifold edges
    {
        label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = -2;
                nProtected++;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " non-manifold edges" << nl << endl;
    }


    // Hand out zones

    DynamicList<label> changedEdges;
    DynamicList<patchEdgeFaceRegion> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchEdgeFaceRegion
    >::propagationTol();

    int dummyTrackData;

    const globalIndex globalFaces(patch.size());

    label facei = 0;

    label currentZoneI = 0;

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        for (; facei < allFaceInfo.size(); facei++)
        {
            if (!allFaceInfo[facei].valid(dummyTrackData))
            {
                globalSeed = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label procI = globalFaces.whichProcID(globalSeed);
        label seedFaceI = globalFaces.toLocal(procI, globalSeed);

        //Info<< "Seeding zone " << currentZoneI
        //    << " from processor " << procI << " face " << seedFaceI
        //    << endl;

        if (procI == Pstream::myProcNo())
        {
            patchEdgeFaceRegion& faceInfo = allFaceInfo[seedFaceI];


            // Set face
            faceInfo = currentZoneI;

            // .. and seed its edges
            const labelList& fEdges = patch.faceEdges()[seedFaceI];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchEdgeFaceRegion& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFaceI,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduce(changedEdges.size(), sumOp<label>()) == 0)
        {
            break;
        }


        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchEdgeFaceRegion
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );

        currentZoneI++;
    }


    faceToZone.setSize(patch.size());
    forAll(allFaceInfo, facei)
    {
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            FatalErrorInFunction
                << "Problem: unvisited face " << facei
                << " at " << patch.faceCentres()[facei]
                << exit(FatalError);
        }
        faceToZone[facei] = allFaceInfo[facei].region();
    }

    return currentZoneI;
}


void Foam::meshRefinement::consistentOrientation
(
    const PackedBoolList& isMasterFace,
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    const labelList& faceToZone,
    const Map<label>& zoneToOrientation,
    PackedBoolList& meshFlipMap
) const
{
    const polyBoundaryMesh& bm = mesh_.boundaryMesh();

    // Data on all edges and faces
    List<patchFaceOrientation> allEdgeInfo(patch.nEdges());
    List<patchFaceOrientation> allFaceInfo(patch.size());

    // Make sure we don't walk through
    // - slaves of coupled faces
    // - non-manifold edges
    {
        label nProtected = 0;

        forAll(patch.addressing(), facei)
        {
            const label meshFaceI = patch.addressing()[facei];
            const label patchi = bm.whichPatch(meshFaceI);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFaceI]
            )
            {
                // Slave side. Mark so doesn't get visited.
                allFaceInfo[facei] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " slaves of coupled faces" << nl << endl;
    }
    {
        label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " non-manifold edges" << nl << endl;
    }



    DynamicList<label> changedEdges;
    DynamicList<patchFaceOrientation> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchFaceOrientation
    >::propagationTol();

    int dummyTrackData;

    globalIndex globalFaces(patch.size());

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        forAll(allFaceInfo, facei)
        {
            if (allFaceInfo[facei] == orientedSurface::UNVISITED)
            {
                globalSeed = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label procI = globalFaces.whichProcID(globalSeed);
        label seedFaceI = globalFaces.toLocal(procI, globalSeed);

        //Info<< "Seeding from processor " << procI << " face " << seedFaceI
        //    << endl;

        if (procI == Pstream::myProcNo())
        {
            // Determine orientation of seedFace

            patchFaceOrientation& faceInfo = allFaceInfo[seedFaceI];

            // Start off with correct orientation
            faceInfo = orientedSurface::NOFLIP;

            if (zoneToOrientation[faceToZone[seedFaceI]] < 0)
            {
                faceInfo.flip();
            }


            const labelList& fEdges = patch.faceEdges()[seedFaceI];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchFaceOrientation& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFaceI,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduce(changedEdges.size(), sumOp<label>()) == 0)
        {
            break;
        }



        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchFaceOrientation
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );
    }


    // Push master zone info over to slave (since slave faces never visited)
    {
        labelList neiStatus
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            orientedSurface::UNVISITED
        );

        forAll(patch.addressing(), i)
        {
            const label meshFaceI = patch.addressing()[i];
            if (!mesh_.isInternalFace(meshFaceI))
            {
                neiStatus[meshFaceI-mesh_.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiStatus);

        forAll(patch.addressing(), i)
        {
            const label meshFaceI = patch.addressing()[i];
            const label patchi = bm.whichPatch(meshFaceI);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFaceI]
            )
            {
                // Slave side. Take flipped from neighbour
                label bFaceI = meshFaceI-mesh_.nInternalFaces();

                if (neiStatus[bFaceI] == orientedSurface::NOFLIP)
                {
                    allFaceInfo[i] = orientedSurface::FLIP;
                }
                else if (neiStatus[bFaceI] == orientedSurface::FLIP)
                {
                    allFaceInfo[i] = orientedSurface::NOFLIP;
                }
                else
                {
                    FatalErrorInFunction
                        << "Incorrect status for face " << meshFaceI
                        << abort(FatalError);
                }
            }
        }
    }


    // Convert to meshFlipMap and adapt faceZones

    meshFlipMap.setSize(mesh_.nFaces());
    meshFlipMap = false;

    forAll(allFaceInfo, faceI)
    {
        label meshFaceI = patch.addressing()[faceI];

        if (allFaceInfo[faceI] == orientedSurface::NOFLIP)
        {
            meshFlipMap[meshFaceI] = false;
        }
        else if (allFaceInfo[faceI] == orientedSurface::FLIP)
        {
            meshFlipMap[meshFaceI] = true;
        }
        else
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << faceI
                << " centre:" << mesh_.faceCentres()[meshFaceI]
                << abort(FatalError);
        }
    }
}


void Foam::meshRefinement::zonify
(
    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList& isMasterFace,
    const labelList& cellToZone,
    const labelList& neiCellZone,
    const labelList& faceToZone,
    const PackedBoolList& meshFlipMap,
    polyTopoChange& meshMod
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label faceZoneI = faceToZone[faceI];

        if (faceZoneI != -1)
        {
            // Orient face zone to have slave cells in min cell zone.
            // Note: logic to use flipMap should be consistent with logic
            //       to pick up the freeStandingBaffleFaces!

            label ownZone = cellToZone[faceOwner[faceI]];
            label neiZone = cellToZone[faceNeighbour[faceI]];

            bool flip;

            if (ownZone == neiZone)
            {
                // free-standing face. Use geometrically derived orientation
                flip = meshFlipMap[faceI];
            }
            else
            {
                flip =
                (
                    ownZone == -1
                 || (neiZone != -1 && ownZone > neiZone)
                );
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[faceI],           // modified face
                    faceI,                          // label of face
                    faceOwner[faceI],               // owner
                    faceNeighbour[faceI],           // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    faceZoneI,                      // zone for face
                    flip                            // face flip in zone
                )
            );
        }
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Set owner as no-flip
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        label faceI = pp.start();

        forAll(pp, i)
        {
            const label meshFaceI = patch.addressing()[i];
            const label patchi = bm.whichPatch(meshFaceI);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFaceI]
            )
            {
                label ownZone = cellToZone[faceOwner[faceI]];
                label neiZone = neiCellZone[faceI-mesh_.nInternalFaces()];

                bool flip;

                if (ownZone == neiZone)
                {
                    // free-standing face. Use geometrically derived orientation
                    flip = meshFlipMap[faceI];
                }
                else
                {
                    flip =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[faceI],           // modified face
                        faceI,                          // label of face
                        faceOwner[faceI],               // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchI,                         // patch for face
                        false,                          // remove from zone
                        faceZoneI,                      // zone for face
                        flip                            // face flip in zone
                    )
                );
            }
            faceI++;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(allFaceInfo, facei)
    {
        label meshFaceI = patch.addressing()[facei];

        if (allFaceInfo[facei] == orientedSurface::NOFLIP)
        {
            meshMod.setAction
            (
                polyModifyCell
                (
                    cellI,
                    false,          // removeFromZone
                    zoneI
                )
            );
        }
        else if (allFaceInfo[facei] == orientedSurface::FLIP)
        {
            Swap(key.first(), key.second());
        }

        HashTable<word, labelPair, labelPair::Hash<> >::
        const_iterator zoneFnd = zoneIDsToFaceZone.find
        (
            key
        );

        if (zoneFnd == zoneIDsToFaceZone.end())
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << facei
                << " centre:" << mesh_.faceCentres()[meshFaceI]
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshRefinement::baffleAndSplitMesh
(
    const bool doHandleSnapProblems,
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,

    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh
)
{
    // Introduce baffles
    // ~~~~~~~~~~~~~~~~~

    // Split the mesh along internal faces wherever there is a pierce between
    // two cell centres.

    Info<< "Introducing baffles for "
        << returnReduce(countHits(), sumOp<label>())
        << " faces that are intersected by the surface." << nl << endl;

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    labelList ownPatch, neiPatch;
    getBafflePatches
    (
        globalToMasterPatch,

        locationsInMesh,
        zonesInMesh,

        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    createBaffles(ownPatch, neiPatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }

    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles");

    if (debug&MESH)
    {
        Pout<< "Writing baffled mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"baffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }


    // Introduce baffles to delete problem cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Create some additional baffles where we want surface cells removed.

    if (doHandleSnapProblems)
    {
        handleSnapProblems
        (
            snapParams,
            useTopologicalSnapDetection,
            removeEdgeConnectedCells,
            perpendicularAngle,
            motionDict,
            runTime,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Removing additional cells might have created disconnected bits
        // so re-do the surface intersections
        {
            // Swap neighbouring cell centres and cell level
            neiLevel.setSize(mesh_.nFaces()-mesh_.nInternalFaces());
            neiCc.setSize(mesh_.nFaces()-mesh_.nInternalFaces());
            calcNeighbourData(neiLevel, neiCc);

            labelList ownPatch, neiPatch;
            getBafflePatches
            (
                globalToMasterPatch,

                locationsInMesh,
                zonesInMesh,

                neiLevel,
                neiCc,

                ownPatch,
                neiPatch
            );

            createBaffles(ownPatch, neiPatch);
        }

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            checkData();
        }
    }


    // Select part of mesh
    // ~~~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Remove unreachable sections of mesh" << nl
        << "-----------------------------------" << nl
        << endl;

    if (debug)
    {
        runTime++;
    }

    splitMeshRegions
    (
        globalToMasterPatch,
        globalToSlavePatch,
        locationsInMesh,
        locationsOutsideMesh
    );

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Split mesh in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After subsetting");

    if (debug&MESH)
    {
        Pout<< "Writing subsetted mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


void Foam::meshRefinement::mergeFreeStandingBaffles
(
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const scalar planarAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh
)
{
    // Merge baffles
    // ~~~~~~~~~~~~~

    Info<< nl
        << "Merge free-standing baffles" << nl
        << "---------------------------" << nl
        << endl;


    // List of pairs of freestanding baffle faces.
    List<labelPair> couples
    (
        freeStandingBaffles    // filter out freestanding baffles
        (
            localPointRegion::findDuplicateFacePairs(mesh_),
            planarAngle
        )
    );

    label nCouples = couples.size();
    reduce(nCouples, sumOp<label>());

    Info<< "Detected free-standing baffles : " << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        mergeBaffles(couples, Map<label>(0));

        // Detect any problem cells resulting from merging of baffles
        // and delete them
        handleSnapProblems
        (
            snapParams,
            useTopologicalSnapDetection,
            removeEdgeConnectedCells,
            perpendicularAngle,
            motionDict,
            runTime,
            globalToMasterPatch,
            globalToSlavePatch
        );

        // Very occasionally removing a problem cell might create a disconnected
        // region so re-check

        Info<< nl
            << "Remove unreachable sections of mesh" << nl
            << "-----------------------------------" << nl
            << endl;

        if (debug)
        {
            runTime++;
        }

        splitMeshRegions
        (
            globalToMasterPatch,
            globalToSlavePatch,
            locationsInMesh,
            locationsOutsideMesh
        );


        if (debug)
        {
            // Debug:test all is still synced across proc patches
            checkData();
        }
    }
    Info<< "Merged free-standing baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMesh
(
    const label nBufferLayers,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,

    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    const pointField& locationsOutsideMesh
)
{
    // Determine patches to put intersections into
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);

    // Find intersections with all unnamed(!) surfaces
    labelList ownPatch, neiPatch;
    getBafflePatches
    (
        globalToMasterPatch,

        locationsInMesh,
        zonesInMesh,

        neiLevel,
        neiCc,

        ownPatch,
        neiPatch
    );

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces(), false);

    forAll(ownPatch, facei)
    {
        if (ownPatch[facei] != -1 || neiPatch[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());


    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Set unreachable cells to -1
    findRegions
    (
        mesh_,
        mergeDistance_*vector(1,1,1),   // perturbVec
        locationsInMesh,
        locationsOutsideMesh,
        cellRegion.nRegions(),
        cellRegion
    );

    Info<< "Found point " << keepPoint
        << " in global region " << keepRegionI
        << " out of " << cellRegion.nRegions() << " regions." << endl;

    if (keepRegionI == -1)
    {
        FatalErrorInFunction
            << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh_.bounds()
            << exit(FatalError);
    }


    // Walk out nBufferlayers from region boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (modifies cellRegion, ownPatch)
    // Takes over face patch onto points and then back to faces and cells
    // (so cell-face-point walk)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Patch for exposed faces for lack of anything sensible.
    label defaultPatch = 0;
    if (globalToMasterPatch.size())
    {
        defaultPatch = globalToMasterPatch[0];
    }

    for (label i = 0; i < nBufferLayers; i++)
    {
        // 1. From cells (via faces) to points

        labelList pointBaffle(mesh_.nPoints(), -1);

        forAll(faceNeighbour, facei)
        {
            const face& f = mesh_.faces()[facei];

            label ownRegion = cellRegion[faceOwner[facei]];
            label neiRegion = cellRegion[faceNeighbour[facei]];

            if (ownRegion == -1 && neiRegion != -1)
            {
                // Note max(..) since possibly regionSplit might have split
                // off extra unreachable parts of mesh. Note: or can this only
                // happen for boundary faces?
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[facei]);
                }
            }
            else if (ownRegion != -1 && neiRegion == -1)
            {
                label newPatchI = neiPatch[facei];
                if (newPatchI == -1)
                {
                    newPatchI = max(defaultPatch, ownPatch[facei]);
                }
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = newPatchI;
                }
            }
        }
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            const face& f = mesh_.faces()[facei];

            label ownRegion = cellRegion[faceOwner[facei]];

            if (ownRegion == -1)
            {
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[facei]);
                }
            }
        }

        // Sync
        syncTools::syncPointList
        (
            mesh_,
            pointBaffle,
            maxEqOp<label>(),
            label(-1)           // null value
        );


        // 2. From points back to faces

        const labelListList& pointFaces = mesh_.pointFaces();

        forAll(pointFaces, pointi)
        {
            if (pointBaffle[pointi] != -1)
            {
                const labelList& pFaces = pointFaces[pointi];

                forAll(pFaces, pFaceI)
                {
                    label facei = pFaces[pFaceI];

                    if (ownPatch[facei] == -1)
                    {
                        ownPatch[facei] = pointBaffle[pointi];
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());


        // 3. From faces to cells (cellRegion) and back to faces (ownPatch)

        labelList newOwnPatch(ownPatch);

        forAll(ownPatch, facei)
        {
            if (ownPatch[facei] != -1)
            {
                label own = faceOwner[facei];

                if (cellRegion[own] == -1)
                {
                    cellRegion[own] = labelMax;

                    const cell& ownFaces = mesh_.cells()[own];
                    forAll(ownFaces, j)
                    {
                        if (ownPatch[ownFaces[j]] == -1)
                        {
                            newOwnPatch[ownFaces[j]] = ownPatch[facei];
                        }
                    }
                }
                if (mesh_.isInternalFace(facei))
                {
                    label nei = faceNeighbour[facei];

                    if (cellRegion[nei] == -1)
                    {
                        cellRegion[nei] = labelMax;

                        const cell& neiFaces = mesh_.cells()[nei];
                        forAll(neiFaces, j)
                        {
                            if (ownPatch[neiFaces[j]] == -1)
                            {
                                newOwnPatch[neiFaces[j]] = ownPatch[facei];
                            }
                        }
                    }
                }
            }
        }

        ownPatch.transfer(newOwnPatch);

        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    }



    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] != keepRegionI)
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label nCellsToKeep = mesh_.nCells() - cellsToRemove.size();
    reduce(nCellsToKeep, sumOp<label>());

    Info<< "Keeping all cells containing points " << locationsInMesh << endl
        << "Selected for keeping : " << nCellsToKeep << " cells." << endl;


    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatches(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        if (ownPatch[facei] != -1)
        {
            exposedPatches[i] = ownPatch[facei];
        }
        else
        {
            WarningInFunction
                << "For exposed face " << facei
                << " fc:" << mesh_.faceCentres()[facei]
                << " found no patch." << endl
                << "    Taking patch " << defaultPatch
                << " instead." << endl;
            exposedPatches[i] = defaultPatch;
        }
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints
(
    const localPointRegion& regionSide
)
{
    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nNonManifPoints = returnReduce
    (
        regionSide.meshPointMap().size(),
        sumOp<label>()
    );

    Info<< "dupNonManifoldPoints : Found : " << nNonManifPoints
        << " non-manifold points (out of "
        << mesh_.globalData().nTotalPoints()
        << ')' << endl;


    autoPtr<mapPolyMesh> map;

    if (nNonManifPoints)
    {
        // Topo change engine
        duplicatePoints pointDuplicator(mesh_);

        // Insert changes into meshMod
        pointDuplicator.setRefinement(regionSide, meshMod);

        // Change the mesh (no inflation, parallel sync)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Is mapping only (no faces created, positions
        // stay same) so no need to recalculate intersections.
        updateMesh(map, labelList(0));
    }

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints()
{
    // Analyse which points need to be duplicated
    localPointRegion regionSide(mesh_);

    return dupNonManifoldPoints(regionSide);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::zonify
(
    const labelList& pointToDuplicate
)
{
    label nPointPairs = 0;
    forAll(pointToDuplicate, pointI)
    {
        label otherPointI = pointToDuplicate[pointI];
        if (otherPointI != -1)
        {
            nPointPairs++;
        }
    }

    autoPtr<mapPolyMesh> map;

    if (returnReduce(nPointPairs, sumOp<label>()))
    {
        Map<label> pointToMaster(2*nPointPairs);
        forAll(pointToDuplicate, pointI)
        {
            label otherPointI = pointToDuplicate[pointI];
            if (otherPointI != -1)
            {
                // Slave point
                pointToMaster.insert(pointI, otherPointI);
            }
        }

        // Topochange container
        polyTopoChange meshMod(mesh_);

        // Insert changes
        polyMeshAdder::mergePoints(mesh_, pointToMaster, meshMod);

        // Change the mesh (no inflation, parallel sync)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Is mapping only (no faces created, positions
        // stay same) so no need to recalculate intersections.
        updateMesh(map, labelList(0));
    }

    return map;
}


// Duplicate points on 'boundary' zones. Do not duplicate points on
// 'internal' or 'baffle' zone. Whether points are on normal patches does
// not matter
Foam::autoPtr<Foam::mapPolyMesh>
Foam::meshRefinement::dupNonManifoldBoundaryPoints()
{
    const labelList boundaryFaceZones
    (
        getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::BOUNDARY
            )
        )
    );
    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = getZones(fzTypes);
    }



    // 0 : point used by normal, unzoned boundary faces
    // 1 : point used by 'boundary' zone
    // 2 : point used by internal/baffle zone
    PackedList<2> pointStatus(mesh_.nPoints(), 0u);

    forAll(boundaryFaceZones, j)
    {
        const faceZone& fZone = mesh_.faceZones()[boundaryFaceZones[j]];
        forAll(fZone, i)
        {
            const face& f = mesh_.faces()[fZone[i]];
            forAll(f, fp)
            {
                pointStatus[f[fp]] = max(pointStatus[f[fp]], 1u);
            }
        }
    }
    forAll(internalOrBaffleFaceZones, j)
    {
        const faceZone& fZone = mesh_.faceZones()[internalOrBaffleFaceZones[j]];
        forAll(fZone, i)
        {
            const face& f = mesh_.faces()[fZone[i]];
            forAll(f, fp)
            {
                pointStatus[f[fp]] = max(pointStatus[f[fp]], 2u);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pointStatus,
        maxEqOp<unsigned int>(),    // combine op
        0u                          // null value
    );

    // Pick up points on boundary zones that are not on internal/baffle zones
    label n = 0;
    forAll(pointStatus, pointI)
    {
        if (pointStatus[pointI] == 1u)
        {
            n++;
        }
    }

    label globalNPoints = returnReduce(n, sumOp<label>());
    Info<< "Duplicating " << globalNPoints << " points on"
        << " faceZones of type "
        << surfaceZonesInfo::faceZoneTypeNames[surfaceZonesInfo::BOUNDARY]
        << endl;

    autoPtr<mapPolyMesh> map;

    if (globalNPoints)
    {
        labelList candidatePoints(n);
        n = 0;
        forAll(pointStatus, pointI)
        {
            if (pointStatus[pointI] == 1u)
            {
                candidatePoints[n++] = pointI;
            }
        }
        localPointRegion regionSide(mesh_, candidatePoints);
        map = dupNonManifoldPoints(regionSide);
    }
    return map;
}


// Zoning
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::zonify
(
    const bool allowFreeStandingZoneFaces,
    const pointField& locationsInMesh,
    const wordList& zonesInMesh,
    wordPairHashTable& zonesToFaceZone
)
{
    if (locationsInMesh.size() != zonesInMesh.size())
    {
        FatalErrorInFunction << "problem" << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();


    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(neiLevel, neiCc);


    // Add any faceZones, cellZones originating from surface to the mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList surfaceToCellZone;
    labelList surfaceToFaceZone;

    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));
    if (namedSurfaces.size())
    {
        Info<< "Setting cellZones according to named surfaces:" << endl;
        forAll(namedSurfaces, i)
        {
<<<<<<< HEAD
            const vectorField smallVec(ROOTSMALL*(end-start));
            start -= smallVec;
            end += smallVec;
=======
            label surfI = namedSurfaces[i];

            Info<< "Surface : " << surfaces_.names()[surfI] << nl
                << "    faceZone : " << surfZones[surfI].faceZoneName() << nl
                << "    cellZone : " << surfZones[surfI].cellZoneName() << endl;
>>>>>>> ENH: snappyHexMesh: various improvements. See below or the default snappyHexMeshDict.
        }
        Info<< endl;

        // Add zones to mesh
        surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
        surfaceToFaceZone = surfaceZonesInfo::addFaceZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
    }



    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Zone per cell:
    // -2 : unset
    // -1 : not in any zone (zone 'none')
    // >=0: zoneID
    // namedSurfaceIndex:
    // -1  : face not intersecting a named surface
    // >=0 : index of named surface
    labelList cellToZone;
    labelList namedSurfaceIndex;
    PackedBoolList posOrientation;
    {
        labelList unnamedRegion1;
        labelList unnamedRegion2;

        zonify
        (
            allowFreeStandingZoneFaces,
            -1,             // Set all cells with cellToZone -2 to -1
            locationsInMesh,
            zonesInMesh,

            cellToZone,
            unnamedRegion1,
            unnamedRegion2,
            namedSurfaceIndex,
            posOrientation
        );
    }


    // Convert namedSurfaceIndex (index of named surfaces) to
    // actual faceZone index

    //- Per face index of faceZone or -1
    labelList faceToZone(mesh_.nFaces(), -1);

    forAll(namedSurfaceIndex, faceI)
    {
        label surfI = namedSurfaceIndex[faceI];
        if (surfI != -1)
        {
            label facei = testFaces[i];

            if (mesh_.isInternalFace(facei))
            {
                start[i] = cellCentres[faceOwner[facei]];
                end[i] = cellCentres[faceNeighbour[facei]];
            }
            else
            {
                start[i] = cellCentres[faceOwner[facei]];
                end[i] = neiCc[facei-mesh_.nInternalFaces()];
            }
        }
    }



    // Allocate and assign faceZones from cellZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // 1. Detect inter-region face and allocate names

        forAll(testFaces, i)
        {
            label facei = testFaces[i];
            const vector& area = mesh_.faceAreas()[facei];

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            if (faceToZone[faceI] == -1)    // Not named surface
            {
                // Face not yet in a faceZone. (it might already have been
                // done so by a 'named' surface). Check if inbetween different
                // cellZones
                allocateInterRegionFaceZone
                (
                    surface2[i] != -1
                 && (
                        magSqr(hit2[i].hitPoint())
                      < magSqr(hit1[i].hitPoint())
                    )
                )
                {
                    namedSurfaceIndex[facei] = surface2[i];
                    posOrientation[facei] = ((area&normal2[i]) > 0);
                    nSurfFaces[surface2[i]]++;
                }
                else
                {
                    namedSurfaceIndex[facei] = surface1[i];
                    posOrientation[facei] = ((area&normal1[i]) > 0);
                    nSurfFaces[surface1[i]]++;
                }
            }
            else if (surface2[i] != -1)
            {
                namedSurfaceIndex[facei] = surface2[i];
                posOrientation[facei] = ((area&normal2[i]) > 0);
                nSurfFaces[surface2[i]]++;
            }
        }

        labelList neiCellZone;
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

        forAll(neiCellZone, bFaceI)
        {
            label faceI = bFaceI + mesh_.nInternalFaces();
            if (faceToZone[faceI] == -1)
            {
                allocateInterRegionFaceZone
                (
                    cellToZone[mesh_.faceOwner()[faceI]],
                    neiCellZone[bFaceI],
                    zonesToFaceZone,
                    zoneIDsToFaceZone
                );
            }
        }


        // 2.Combine faceZoneNames allocated on different processors

        Pstream::mapCombineGather(zonesToFaceZone, eqOp<word>());
        Pstream::mapCombineScatter(zonesToFaceZone);


        // 3. Allocate faceZones from (now synchronised) faceZoneNames
        //    Note: the faceZoneNames contain the same data but in different
        //          order. We could sort the contents but instead just loop
        //          in sortedToc order.

        Info<< "Setting faceZones according to neighbouring cellZones:"
            << endl;

        // From cellZone indices to faceZone index
        HashTable<label, labelPair, labelPair::Hash<> > fZoneLookup
        (
            zonesToFaceZone.size()
        );

        const cellZoneMesh& cellZones = mesh_.cellZones();

        {
            List<Pair<word> > czs(zonesToFaceZone.sortedToc());

            forAll(czs, i)
            {
                const Pair<word>& cz = czs[i];
                const word& fzName = zonesToFaceZone[cz];

                Info<< indent<< "cellZones : "
                    << cz[0] << ' ' << cz[1] << nl
                    << "    faceZone : " << fzName << endl;

                label faceZoneI = surfaceZonesInfo::addFaceZone
                (
                    fzName,                 // name
                    labelList(0),           // addressing
                    boolList(0),            // flipMap
                    mesh_
                );

                label cz0 = cellZones.findZoneID(cz[0]);
                label cz1 = cellZones.findZoneID(cz[1]);

                fZoneLookup.insert(labelPair(cz0, cz1), faceZoneI);
            }
        }


        // 4. Set faceToZone with new faceZones


        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            if (faceToZone[faceI] == -1)
            {
                // Face not yet in a faceZone. (it might already have been
                // done so by a 'named' surface). Check if inbetween different
                // cellZones

                label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
                label neiZone = cellToZone[mesh_.faceNeighbour()[faceI]];
                if (ownZone != neiZone)
                {
                    bool swap =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                    labelPair key(ownZone, neiZone);
                    if (swap)
                    {
                        Swap(key.first(), key.second());
                    }
                    faceToZone[faceI] = fZoneLookup[key];
                }
            }
        }
        forAll(neiCellZone, bFaceI)
        {
            label faceI = bFaceI + mesh_.nInternalFaces();
            if (faceToZone[faceI] == -1)
            {
                label ownZone = cellToZone[mesh_.faceOwner()[faceI]];
                label neiZone = neiCellZone[bFaceI];
                if (ownZone != neiZone)
                {
                    bool swap =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                    labelPair key(ownZone, neiZone);
                    if (swap)
                    {
                        Swap(key.first(), key.second());
                    }
                    faceToZone[faceI] = fZoneLookup[key];
                }
            }
        }
        Info<< endl;
    }


    //- Per face index of faceZone or -1
    labelList faceToZone(mesh_.nFaces(), -1);

    // Convert namedSurfaceIndex (index of named surfaces) to
    // actual faceZone index

    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];
        if (surfI != -1)
        {
            faceToZone[facei] = surfaceToFaceZone[surfI];
        }
    }


    // Topochange container
    polyTopoChange meshMod(mesh_);



    // Get coupled neighbour cellZone. Set to -1 on non-coupled patches.
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!pp.coupled())
        {
            label bFaceI = pp.start()-mesh_.nInternalFaces();
            forAll(pp, i)
            {
                neiCellZone[bFaceI++] = -1;
            }
        }
    }



    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));


    // faceZones
    // ~~~~~~~~~
    // Faces on faceZones come in two variants:
    // - faces on the outside of a cellZone. They will be oriented to
    //   point out of the maximum cellZone.
    // - free-standing faces. These will be oriented according to the
    //   local surface normal. We do this in a two step algorithm:
    //      - do a consistent orientation
    //      - check number of faces with consistent orientation
    //      - if <0 flip the whole patch
    PackedBoolList meshFlipMap(mesh_.nFaces(), false);
    {
        // Collect all data on zone faces without cellZones on either side.
        const indirectPrimitivePatch patch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                freeStandingBaffleFaces
                (
                    faceToZone,
                    cellToZone,
                    neiCellZone
                )
            ),
            mesh_.points()
        );

        label nFreeStanding = returnReduce(patch.size(), sumOp<label>());
        if (nFreeStanding > 0)
        {
            Info<< "Detected " << nFreeStanding << " free-standing zone faces"
                << endl;

            if (debug)
            {
                OBJstream str(mesh_.time().path()/"freeStanding.obj");
                str.write(patch.localFaces(), patch.localPoints(), false);
            }


            // Detect non-manifold edges
            labelList nMasterFacesPerEdge;
            calcPatchNumMasterFaces(isMasterFace, patch, nMasterFacesPerEdge);

            // Mark zones. Even a single original surface might create multiple
            // disconnected/non-manifold-connected zones
            labelList faceToConnectedZone;
            const label nZones = markPatchZones
            (
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone
            );

            Map<label> nPosOrientation(2*nZones);
            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                nPosOrientation.insert(zoneI, 0);
            }

            // Make orientations consistent in a topological way. This just
            // checks  the first face per zone for whether nPosOrientation
            // is negative (which it never is at this point)
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );

            // Count per region the number of orientations (taking the new
            // flipMap into account)
            forAll(patch.addressing(), facei)
            {
                label meshFaceI = patch.addressing()[facei];

                if (isMasterFace[meshFaceI])
                {
                    label n = 1;
                    if
                    (
                        bool(posOrientation[meshFaceI])
                     == meshFlipMap[meshFaceI]
                    )
                    {
                        n = -1;
                    }

                    nPosOrientation.find(faceToConnectedZone[facei])() += n;
                }
            }
            Pstream::mapCombineGather(nPosOrientation, plusEqOp<label>());
            Pstream::mapCombineScatter(nPosOrientation);


            Info<< "Split " << nFreeStanding << " free-standing zone faces"
                << " into " << nZones << " disconnected regions with size"
                << " (negative denotes wrong orientation) :"
                << endl;

            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                Info<< "    " << zoneI << "\t" << nPosOrientation[zoneI]
                    << endl;
            }
            Info<< endl;


            // Re-apply with new counts (in nPosOrientation). This will cause
            // zones with a negative count to be flipped.
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );
        }
    }


    // Put the faces into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label faceZoneI = faceToZone[facei];

        if (faceZoneI != -1)
        {
            // Orient face zone to have slave cells in max cell zone.
            // Note: logic to use flipMap should be consistent with logic
            //       to pick up the freeStandingBaffleFaces!

            label ownZone = cellToZone[faceOwner[facei]];
            label neiZone = cellToZone[faceNeighbour[facei]];

            bool flip;

            if (ownZone == neiZone)
            {
                // free-standing face. Use geometrically derived orientation
                flip = meshFlipMap[facei];
            }
            else
            {
                flip =
                (
                    ownZone == -1
                 || (neiZone != -1 && ownZone > neiZone)
                );
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[facei],           // modified face
                    facei,                          // label of face
                    faceOwner[facei],               // owner
                    faceNeighbour[facei],           // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    faceZoneI,                      // zone for face
                    flip                            // face flip in zone
                )
            );
        }
    }


    // Set owner as no-flip
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        label facei = pp.start();

        forAll(pp, i)
        {
            label faceZoneI = faceToZone[facei];

            if (faceZoneI != -1)
            {
                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                bool flip;

                if (ownZone == neiZone)
                {
                    // free-standing face. Use geometrically derived orientation
                    flip = meshFlipMap[facei];
                }
                else
                {
                    flip =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[facei],           // modified face
                        facei,                          // label of face
                        faceOwner[facei],               // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchi,                         // patch for face
                        false,                          // remove from zone
                        faceZoneI,                      // zone for face
                        flip                            // face flip in zone
                    )
                );
            }
            facei++;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToZone, celli)
    {
        label zoneI = cellToZone[celli];

        if (zoneI >= 0)
        {
            meshMod.setAction
            (
                polyModifyCell
                (
                    celli,
                    false,          // removeFromZone
                    zoneI
                )
            );
        }
    }

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Print some stats (note: zones are synchronised)
    if (mesh_.cellZones().size() > 0)
    {
        Info<< "CellZones:" << endl;
        forAll(mesh_.cellZones(), zoneI)
        {
            const cellZone& cz = mesh_.cellZones()[zoneI];
            Info<< "    " << cz.name()
                << "\tsize:" << returnReduce(cz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }
    if (mesh_.faceZones().size() > 0)
    {
        Info<< "FaceZones:" << endl;
        forAll(mesh_.faceZones(), zoneI)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            Info<< "    " << fz.name()
                << "\tsize:" << returnReduce(fz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }

    // None of the faces has changed, only the zones. Still...
    updateMesh(map, labelList());

    return map;
}


// ************************************************************************* //
