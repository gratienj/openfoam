/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2026 OpenCFD Ltd.
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

#include "faMeshDecomposition.H"
#include "Time.H"
#include "dictionary.H"
#include "labelIOList.H"
#include "Map.H"
#include "ListOps.H"
#include "globalMeshData.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMeshDecomposition::distributeFaces()
{
    const word& polyMeshRegionName = faMesh::mesh().name();

    Info<< "\nCalculating distribution of finite-area faces ["
        << polyMesh::regionName(areaName_) << "]" << endl;

    cpuTime decompositionTime;

    for (label proci = 0; proci < nProcs(); ++proci)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(proci)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        IOobject ioFvAddr
        (
            "procAddressing",
            "constant",
            polyMesh::meshSubDir,
            procFvMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );


        // faceProcAddressing (polyMesh)
        ioFvAddr.resetHeader("faceProcAddressing");
        const labelList fvFaceProcAddressing
        (
            labelIOList::readContents(ioFvAddr)
        );

        labelHashSet faceProcAddressingHash;
        faceProcAddressingHash.reserve(fvFaceProcAddressing.size());

        // If faMesh's fvPatch is a part of the global face zones, faces of that
        // patch will be present on all processors. Because of that, looping
        // through faceProcAddressing will decompose global faMesh faces to the
        // very last processor regardless of where fvPatch is really decomposed.
        // Since global faces which do not belong to specific processor are
        // located at the end of the faceProcAddressing, cutting it at
        // i = owner.size() will correctly decompose faMesh faces.
        // Vanja Skuric, 2016-04-21
        if (hasGlobalFaceZones_)
        {
            // owner (polyMesh)
            ioFvAddr.resetHeader("owner");
            const label ownerSize = labelIOList::readContentsSize(ioFvAddr);

            for (label i = 0; i < ownerSize; ++i)
            {
                faceProcAddressingHash.insert(fvFaceProcAddressing[i]);
            }
        }
        else
        {
            faceProcAddressingHash.insert(fvFaceProcAddressing);
        }

        forAll(faceLabels(), facei)
        {
            // With +1 for lookup in faceMap with flip encoding
            const label index = (faceLabels()[facei] + 1);

            if (faceProcAddressingHash.contains(index))
            {
                faceToProc_[facei] = proci;
            }
        }
    }

    Info<< "\nFinished decomposition in "
        << decompositionTime.elapsedCpuTime()
        << " s" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshDecomposition::faMeshDecomposition
(
    const word& areaName,
    const polyMesh& mesh,
    const label nProcessors,
    const dictionary& params
)
:
    faMesh(areaName, mesh),
    areaName_(areaName.empty() ? polyMesh::defaultRegion : areaName),
    nProcs_(nProcessors),
    distributed_(false),
    hasGlobalFaceZones_(false),
    cyclicParallel_(false),
    faceToProc_(faMesh::nFaces()),
    procFaceLabels_(nProcs_),
    procMeshEdgesMap_(nProcs_),
    procNInternalEdges_(nProcs_, Zero),
    procPatchEdgeLabels_(nProcs_),
    procPatchPointAddressing_(nProcs_),
    procPatchEdgeAddressing_(nProcs_),
    procEdgeAddressing_(nProcs_),
    procFaceAddressing_(nProcs_),
    procBoundaryAddressing_(nProcs_),
    procPatchSize_(nProcs_),
    procPatchStartIndex_(nProcs_),
    procNeighbourProcessors_(nProcs_),
    procProcessorPatchSize_(nProcs_),
    procProcessorPatchStartIndex_(nProcs_),
    globallySharedPoints_()
{
    updateParameters(params);
}


Foam::faMeshDecomposition::faMeshDecomposition
(
    const polyMesh& mesh,
    const label nProcessors,
    const dictionary& params
)
:
    faMeshDecomposition
    (
        polyMesh::defaultRegion,
        mesh,
        nProcessors,
        params
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMeshDecomposition::updateParameters
(
    const dictionary& params
)
{
    params.readIfPresent("distributed", distributed_);
    if (params.found("globalFaceZones"))
    {
        hasGlobalFaceZones_ = true;
    }
}


void Foam::faMeshDecomposition::decomposeMesh()
{
    // Decide which cell goes to which processor
    distributeFaces();

    const word& polyMeshRegionName = faMesh::mesh().name();

    Info<< "\nDistributing faces to processors ["
        << polyMesh::regionName(areaName_) << "]" << endl;

    labelList nLocalFaces(nProcs_, Zero);

    // Pass 1: determine local sizes, sanity check

    forAll(faceToProc_, facei)
    {
        const label proci = faceToProc_[facei];

        if (proci < 0 || proci >= nProcs_)
        {
            FatalErrorInFunction
                << "Invalid processor label " << proci
                << " for face " << facei << nl
                << abort(FatalError);
        }
        else
        {
            ++nLocalFaces[proci];
        }
    }

    // Adjust lengths
    forAll(nLocalFaces, proci)
    {
        procFaceAddressing_[proci].resize(nLocalFaces[proci]);
        nLocalFaces[proci] = 0;  // restart list
    }

    // Pass 2: fill in local lists
    forAll(faceToProc_, facei)
    {
        const label proci = faceToProc_[facei];
        const label localFacei = nLocalFaces[proci];
        ++nLocalFaces[proci];

        procFaceAddressing_[proci][localFacei] = facei;
    }


    // Find processor mesh faceLabels and ...

    for (label procI = 0; procI < nProcs(); procI++)
    {
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        IOobject ioAddr
        (
            "procAddressing",
            "constant",
            polyMesh::meshSubDir,
            procFvMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );


        // pointProcAddressing (polyMesh)
        ioAddr.resetHeader("pointProcAddressing");
        const labelList fvPointProcAddressing
        (
            labelIOList::readContents(ioAddr)
        );

        Map<label> fvFaceProcAddressingHash;

        {
            // faceProcAddressing (polyMesh)
            ioAddr.resetHeader("faceProcAddressing");
            const labelList fvFaceProcAddressing
            (
                labelIOList::readContents(ioAddr)
            );
            fvFaceProcAddressingHash = invertToMap(fvFaceProcAddressing);
        }


        const labelList& curProcFaceAddressing = procFaceAddressing_[procI];

        labelList& curFaceLabels = procFaceLabels_[procI];

        curFaceLabels.resize_fill(curProcFaceAddressing.size(), -1);

        forAll(curProcFaceAddressing, faceI)
        {
            curFaceLabels[faceI] =
                fvFaceProcAddressingHash.find
                (
                    faceLabels()[curProcFaceAddressing[faceI]] + 1
                ).val();
        }

        // Create processor finite-area mesh
        faMesh procMesh
        (
            areaName_,
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );

        const uindirectPrimitivePatch& patch = this->patch();
        const Map<label>& map = patch.meshPointMap();

        EdgeMap<label> edgesHash;
        edgesHash.reserve(patch.nEdges());

        const label nIntEdges = patch.nInternalEdges();

        for (label edgei = 0; edgei < nIntEdges; ++edgei)
        {
            edgesHash.insert(patch.edges()[edgei], edgesHash.size());
        }

        for (const auto& fap : faMesh::boundary())
        {
            // Also include emptyFaPatch etc
            for (const label edgei : fap.edgeLabels())
            {
                edgesHash.insert(patch.edges()[edgei], edgesHash.size());
            }
        }


        const uindirectPrimitivePatch& procPatch = procMesh.patch();
        const labelUList& procMeshPoints = procPatch.meshPoints();
        const edgeList& procEdges = procPatch.edges();

        labelList& curPatchPointAddressing = procPatchPointAddressing_[procI];
        curPatchPointAddressing.resize(procMeshPoints.size(), -1);

        forAll(procMeshPoints, pointi)
        {
            curPatchPointAddressing[pointi] =
                map[fvPointProcAddressing[procMeshPoints[pointi]]];
        }

        procNInternalEdges_[procI] = procPatch.nInternalEdges();

        auto& curPatchEdgeAddressing = procPatchEdgeAddressing_[procI];
        curPatchEdgeAddressing.resize(procEdges.size(), -1);

        Map<label>& curMap = procMeshEdgesMap_[procI];
        curMap.clear();
        curMap.reserve(procEdges.size());

        forAll(procEdges, edgeI)
        {
            edge curGlobalEdge(curPatchPointAddressing, procEdges[edgeI]);

            if (auto iter = edgesHash.cfind(curGlobalEdge); iter.good())
            {
                // The edgeID (not edgeLabel) in serial
                auto globalEdgeId = iter.val();

                // For each proc edgeLabel, the serial edgeID.
                curPatchEdgeAddressing[edgeI] = globalEdgeId;

                // (key = serial edgeID; val = proc edgeLabel)
                curMap.insert(globalEdgeId, edgeI);
            }
            else
            {
                // Only fails if the polyMesh pointProcAddressing is corrupt
                FatalErrorInFunction
                    << "Failed edge lookup of " << curGlobalEdge << endl
                    << exit(FatalError);
            }
        }
    }


    Info << "\nDistributing edges to processors" << endl;

    // Loop through all internal edges and decide which processor they
    // belong to. First visit all internal edges.

    // set references to the original mesh
    const faBoundaryMesh& patches = boundary();
    const edgeList& edges = this->edges();
    const labelUList& owner = edgeOwner();
    const labelUList& neighbour = edgeNeighbour();

    // Memory management
    {
        List<DynamicList<label>> procEdgeList(nProcs());

        forAll(procEdgeList, procI)
        {
            const auto& procEdgeAddr = procPatchEdgeAddressing_[procI];
            const auto nProcInternalEdges = procNInternalEdges_[procI];

            // Copy the internal edges ids
            procEdgeList[procI].reserve(procEdgeAddr.size());
            procEdgeList[procI].push_back
            (
                procEdgeAddr.slice(0, nProcInternalEdges)
            );
        }


        // Detect inter-processor boundaries
        // Track processor boundaries as (neighbour rank, edgeLabels)
        // lists for each subdomain
        List
        <
            DynamicList<std::pair<label, DynamicList<label>>>
        > interProcBoundaries(nProcs());

        forAll(neighbour, edgeI)
        {
            const label ownProc = faceToProc_[owner[edgeI]];
            const label neiProc = faceToProc_[neighbour[edgeI]];

            if (ownProc != neiProc)
            {
                // inter - processor patch edge found. Go through the list of
                // inside boundaries for the owner processor and try to find
                // this inter-processor patch.

                bool interProcBouFound = false;

                for
                (
                    auto& [ownNbrProc, ownProcEdges]
                  : interProcBoundaries[ownProc]
                )
                {
                    if (ownNbrProc == neiProc)
                    {
                        // the inter - processor boundary exists
                        interProcBouFound = true;

                        ownProcEdges.push_back(edgeI);

                        bool neighbourFound = false;

                        for
                        (
                            auto& [neiNbrProc, neiProcEdges]
                          : interProcBoundaries[neiProc]
                        )
                        {
                            if (neiNbrProc == ownProc)
                            {
                                // boundary found. Add the face
                                neighbourFound = true;

                                neiProcEdges.push_back(edgeI);
                            }

                            if (neighbourFound) break;
                        }

                        if (interProcBouFound && !neighbourFound)
                        {
                            FatalErrorInFunction
                                << "Inconsistency in inter-processor"
                                << " boundary lists for processors "
                                << ownProc << " and " << neiProc
                                << abort(FatalError);
                        }
                    }

                    if (interProcBouFound) break;
                }

                if (!interProcBouFound)
                {
                    // inter - processor boundaries do not exist
                    // and need to be created

                    // owner -> neighbour
                    auto& [ownNbrProc, ownProcEdges] =
                        interProcBoundaries[ownProc].emplace_back();

                    ownNbrProc = neiProc;
                    ownProcEdges.push_back(edgeI);

                    // neighbour -> owner
                    auto& [neiNbrProc, neiProcEdges] =
                        interProcBoundaries[neiProc].emplace_back();

                    neiNbrProc = ownProc;
                    neiProcEdges.push_back(edgeI);
                }
            }
        }


        // Loop through patches. For cyclic boundaries detect inter-processor
        // edges; for all other, add edges to the edge list and remember start
        // and size of all patches.

        // for all processors, set the size of start index and patch size
        // lists to the number of patches in the mesh
        forAll(procPatchSize_, procI)
        {
            procPatchSize_[procI].setSize(patches.size());
            procPatchStartIndex_[procI].setSize(patches.size());
        }

        forAll(patches, patchI)
        {
            const faPatch& fap = patches[patchI];

            // Reset size and start index for all processors
            forAll(procPatchSize_, procI)
            {
                procPatchSize_[procI][patchI] = 0;
                procPatchStartIndex_[procI][patchI] =
                    procEdgeList[procI].size();
            }

            const label patchStart = fap.start();

//             if (!isA<cyclicFaPatch>(patches[patchI]))
            if (true)
            {
                // Normal patch. Add edges to processor where the face
                // next to the edge lives

                const labelUList& patchEdgeLabels = fap.edgeLabels();

                forAll(patchEdgeLabels, patchEdgei)
                {
                    const label edgeLabel = patchEdgeLabels[patchEdgei];
                    const label facei = patch().edgeOwner(edgeLabel);
                    const label curProc = faceToProc_[facei];

                    // Add to the list of edges
                    procEdgeList[curProc].push_back(patchStart + patchEdgei);

                    // increment the number of edges for this patch
                    procPatchSize_[curProc][patchI]++;
                }
            }
            else
            {
                // Cyclic patch special treatment

                const faPatch& cPatch = patches[patchI];

                const label cycOffset = cPatch.size()/2;

                // Set reference to faceCells for both patches
                const auto firstEdgeFaces
                (
                    cPatch.edgeFaces().slice(0, cycOffset)
                );

                const auto secondEdgeFaces
                (
                    cPatch.edgeFaces().slice(cycOffset)
                );

                const auto startPart1 = (patchStart);
                const auto startPart2 = (patchStart + cycOffset);

                forAll(firstEdgeFaces, edgeI)
                {
                    const label ownProc = faceToProc_[firstEdgeFaces[edgeI]];
                    const label neiProc = faceToProc_[secondEdgeFaces[edgeI]];

                    const label firstEdgei  = (startPart1 + edgeI);
                    const label secondEdgei = (startPart2 + edgeI);

                    if (ownProc != neiProc)
                    {
                        // This edge becomes an inter-processor boundary edge
                        // inter - processor patch edge found. Go through
                        // the list of inside boundaries for the owner
                        // processor and try to find this inter-processor
                        // patch.

                        cyclicParallel_ = true;

                        bool interProcBouFound = false;

                        for
                        (
                            auto& [ownNbrProc, ownProcEdges]
                          : interProcBoundaries[ownProc]
                        )
                        {
                            if (ownNbrProc == neiProc)
                            {
                                // the inter - processor boundary exists.
                                interProcBouFound = true;

                                ownProcEdges.push_back(firstEdgei);

                                bool neighbourFound = false;

                                for
                                (
                                    auto& [neiNbrProc, neiProcEdges]
                                  : interProcBoundaries[neiProc]
                                )
                                {
                                    if (neiNbrProc == ownProc)
                                    {
                                        // boundary found. Add the face
                                        neighbourFound = true;

                                        neiProcEdges.push_back(secondEdgei);
                                    }

                                    if (neighbourFound) break;
                                }

                                if (interProcBouFound && !neighbourFound)
                                {
                                    FatalErrorInFunction
                                        << "Inconsistency in inter-processor"
                                        << " boundary lists for processors "
                                        << ownProc << " and " << neiProc
                                        << " in cyclic boundary matching"
                                        << abort(FatalError);
                                }
                            }

                            if (interProcBouFound) break;
                        }

                        if (!interProcBouFound)
                        {
                            // inter - processor boundaries do not exist
                            // and need to be created

                            // owner -> neighbour
                            auto& [ownNbrProc, ownProcEdges] =
                                interProcBoundaries[ownProc].emplace_back();

                            ownNbrProc = neiProc;
                            ownProcEdges.push_back(firstEdgei);

                            // neighbour -> owner
                            auto& [neiNbrProc, neiProcEdges] =
                                interProcBoundaries[neiProc].emplace_back();

                            neiNbrProc = ownProc;
                            neiProcEdges.push_back(secondEdgei);
                        }
                    }
                    else
                    {
                        // This cyclic edge remains on the processor

                        // add the first edge
                        procEdgeList[ownProc].push_back(firstEdgei);

                        // increment the number of edges for this patch
                        procPatchSize_[ownProc][patchI]++;

                        // Note: I cannot add the other side of the cyclic
                        // boundary here because this would violate the order.
                        // They will be added in a separate loop below
                    }
                }

                // Ordering in cyclic boundaries is important.
                // Add the other half of cyclic edges for cyclic boundaries
                // that remain on the processor
                forAll(secondEdgeFaces, edgeI)
                {
                    const label ownProc = faceToProc_[firstEdgeFaces[edgeI]];
                    const label neiProc = faceToProc_[secondEdgeFaces[edgeI]];

                    //const label firstEdgei  = (startPart1 + edgeI);
                    const label secondEdgei = (startPart2 + edgeI);

                    if (ownProc == neiProc)
                    {
                        // This cyclic edge remains on the processor

                        // add the second edge
                        procEdgeList[ownProc].push_back(secondEdgei);

                        // increment the number of edges for this patch
                        procPatchSize_[ownProc][patchI]++;
                    }
                }
            }
        }


        // Sort processor connections according to neighbProcNo.
        // Not needed for functionality, but gives consistently
        // ordered boundaries.

        for (auto& procBoundaries : interProcBoundaries)
        {
            Foam::sort(procBoundaries);
        }

        // Add inter-processor boundaries and remember start indices
        forAll(procEdgeList, procI)
        {
            // Get internal and regular boundary processor faces
            const auto& curProcEdges = procEdgeList[procI];

            // Get reference to processor edge addressing
            labelList& curProcEdgeAddressing = procEdgeAddressing_[procI];

            labelList& curProcNeighbourProcessors =
                procNeighbourProcessors_[procI];

            labelList& curProcProcessorPatchStartIndex =
                procProcessorPatchStartIndex_[procI];

            labelList& curProcProcessorPatchSize =
                procProcessorPatchSize_[procI];

            // Number of internal and non-processor edges
            const label nNonProcessorEdges = curProcEdges.size();

            const auto& curInterProcBoundaries = interProcBoundaries[procI];

            const auto nProcPatches = interProcBoundaries[procI].size();

            // Flattened values for processor patches
            curProcNeighbourProcessors.resize_nocopy(nProcPatches);
            curProcProcessorPatchStartIndex.resize_nocopy(nProcPatches);
            curProcProcessorPatchSize.resize_nocopy(nProcPatches);

            // Number of processor edges
            label nProcessorEdges = 0;

            for (label procPatchi = 0; procPatchi < nProcPatches; ++procPatchi)
            {
                const auto& [bndNbrProc, bndProcEdges] =
                    curInterProcBoundaries[procPatchi];

                curProcNeighbourProcessors[procPatchi] = bndNbrProc;
                curProcProcessorPatchSize[procPatchi] = bndProcEdges.size();

                // Starts after all previous
                curProcProcessorPatchStartIndex[procPatchi] =
                (
                    nNonProcessorEdges + nProcessorEdges
                );

                nProcessorEdges += bndProcEdges.size();
            }

            // Resize addressing
            curProcEdgeAddressing.resize(nNonProcessorEdges + nProcessorEdges);

            // Fill in the list. Calculate turning index.
            // Turning index will be -1 only for some edges on processor
            // boundaries, i.e. the ones where the current processor ID
            // is in the face which is a edge neighbour.
            // Turning index is stored as the sign of the edge addressing list

            label nEdges = 0;

            // Add internal and boundary edges
            // Remember to increment the index by one such that the
            // turning index works properly.
            for (const label procEdgei : curProcEdges)
            {
                curProcEdgeAddressing[nEdges] = procEdgei;
//                 curProcEdgeAddressing[nEdges] = procEdgei + 1;
                ++nEdges;
            }

            // Processor boundaries
            for (label procPatchi = 0; procPatchi < nProcPatches; ++procPatchi)
            {
                const auto& [bndNbrProc, bndProcEdges] =
                    curInterProcBoundaries[procPatchi];

                for (label edgei : bndProcEdges)
                {
                    // Remember to increment the index by one such that the
                    // turning index works properly.
                    if (faceToProc_[owner[edgei]] == procI)
                    {
                        curProcEdgeAddressing[nEdges] = edgei;
//                      curProcEdgeAddressing[nEdges] = edgei + 1;
                    }
                    else
                    {
                        // turning edge
                        curProcEdgeAddressing[nEdges] = edgei;
//                      curProcEdgeAddressing[nEdges] = -(edgei + 1);
                    }
                    ++nEdges;
                }

                // Debug and suppress [[maybe_unused]] warning
                if (debug & 2)
                {
                    Info<< "proc" << procI << "to" << bndNbrProc
                        << " edgeLabels "; bndProcEdges.writeList(Info) << endl;
                }
            }
        }
    }

    Info << "\nCalculating processor boundary addressing" << endl;
    // For every patch: the original patch index.
    // - identity for non-processor patches (ie, globally identical)
    // - '-1' for processor patches
    forAll(procPatchSize_, proci)
    {
        const auto nNonProcessorPatches = procPatchSize_[proci].size();
        const auto nProcPatches = procProcessorPatchSize_[proci].size();

        auto& curBoundaryAddressing = procBoundaryAddressing_[proci];

        curBoundaryAddressing =
            Foam::identity(nNonProcessorPatches + nProcPatches);

        // Mark processor patches
        curBoundaryAddressing.slice(nNonProcessorPatches) = -1;
    }


    // Gather data about globally shared points

    labelList globallySharedPoints_;

    // Memory management
    {
        labelList pointsUsage(nPoints(), Zero);

        // Globally shared points are the ones used by more than 2 processors
        // Size the list approximately and gather the points
        labelHashSet gSharedPoints;
        gSharedPoints.reserve(Foam::min(128, nPoints()/1000));

        // Loop through all the processors and mark up points used by
        // processor boundaries.  When a point is used twice, it is a
        // globally shared point

        for (label procI = 0; procI < nProcs(); procI++)
        {
            // Get list of edge labels
            const labelList& curEdgeLabels = procEdgeAddressing_[procI];

            // Get start of processor faces
            const labelList& curProcessorPatchStarts =
                procProcessorPatchStartIndex_[procI];

            const labelList& curProcessorPatchSizes =
                procProcessorPatchSize_[procI];

            // Reset the lookup list
            pointsUsage = 0;

            forAll(curProcessorPatchStarts, patchI)
            {
                const auto patchEdgeLabels = curEdgeLabels.slice
                (
                    curProcessorPatchStarts[patchI],
                    curProcessorPatchSizes[patchI]
                );

                for (const label edgei : patchEdgeLabels)
                {
                    // Mark the original edge as used
                    // Remember to decrement the index by one (turning index)

                    const edge& e = edges[edgei];

                    forAll(e, pointI)
                    {
                        if (pointsUsage[e[pointI]] == 0)
                        {
                            // Point not previously used
                            pointsUsage[e[pointI]] = patchI + 1;
                        }
                        else if (pointsUsage[e[pointI]] != patchI + 1)
                        {
                            // Point used by some other patch = global point!
                            gSharedPoints.insert(e[pointI]);
                        }
                    }
                }
            }
        }

        // Grab the result from the hash list
        globallySharedPoints_ = gSharedPoints.sortedToc();
    }


    // Edge label for faPatches

    for (label procI = 0; procI < nProcs(); procI++)
    {
        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );


        // Read volume mesh
        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        // Create processor finite-area mesh
        faMesh procMesh
        (
            areaName_,
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );

        // The serial edgeId values for the local mesh edges
        const labelList& curEdgeAddressing = procEdgeAddressing_[procI];

        const labelList& curPatchStartIndex = procPatchStartIndex_[procI];
        const labelList& curPatchSize = procPatchSize_[procI];

        const labelList& curProcessorPatchStartIndex =
            procProcessorPatchStartIndex_[procI];

        const labelList& curProcessorPatchSize =
            procProcessorPatchSize_[procI];

        const label nNonProcessorPatches = curPatchSize.size();
        const label nProcPatches = curProcessorPatchSize.size();

        labelListList& curPatchEdgeLabels = procPatchEdgeLabels_[procI];
        curPatchEdgeLabels.resize_nocopy(nNonProcessorPatches + nProcPatches);

        for (label patchi = 0; patchi < nNonProcessorPatches; ++patchi)
        {
            // [output] : edgeLabels for the patch
            auto& curEdgeLabels = curPatchEdgeLabels[patchi];

            // Copy the local mesh edge ids for the patch
            curEdgeLabels = curEdgeAddressing.slice
            (
                curPatchStartIndex[patchi],
                curPatchSize[patchi]
            );

            // Inplace renumber with Map
            // (key = serial edgeID; val = proc edgeLabel)

            inplaceRenumber(procMeshEdgesMap_[procI], curEdgeLabels);
        }

        for (label procPatchi = 0; procPatchi < nProcPatches; ++procPatchi)
        {
            // [output] : edgeLabels for the patch
            auto& curEdgeLabels =
                curPatchEdgeLabels[nNonProcessorPatches + procPatchi];

            // Copy the local mesh edge ids for the patch
            curEdgeLabels = curEdgeAddressing.slice
            (
                curProcessorPatchStartIndex[procPatchi],
                curProcessorPatchSize[procPatchi]
            );

            // Inplace renumber with Map
            // (key = serial edgeID; val = proc edgeLabel)
            inplaceRenumber(procMeshEdgesMap_[procI], curEdgeLabels);
        }
    }
}


bool Foam::faMeshDecomposition::writeDecomposition() const
{
    const word& polyMeshRegionName = faMesh::mesh().name();

    Info<< "\nConstructing processor FA meshes" << endl;

    // Make a lookup map for globally shared points
    Map<label> sharedPointLookup(invertToMap(globallySharedPoints_));


    label maxProcFaces = 0, totProcFaces = 0;
    label maxProcEdges = 0, totProcEdges = 0;
    label maxProcPatches = 0, totProcPatches = 0;

    // Write out the meshes
    for (label procI = 0; procI < nProcs(); procI++)
    {
        // Create processor mesh without a boundary

        // create a database
        Time processorDb
        (
            Time::controlDictName,
            time().rootPath(),
            time().caseName()/("processor" + Foam::name(procI)),
            false,  // No function objects
            false   // No extra controlDict libs
        );

        // Read volume mesh
        polyMesh procFvMesh
        (
            IOobject
            (
                polyMeshRegionName,
                processorDb.timeName(),
                processorDb
            )
        );

        IOobject ioFvAddr
        (
            "procAddressing",
            "constant",
            polyMesh::meshSubDir,
            procFvMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // boundaryProcAddressing (polyMesh)
        ioFvAddr.resetHeader("boundaryProcAddressing");
        const labelList fvBoundaryProcAddressing
        (
            labelIOList::readContents(ioFvAddr)
        );


        // Create processor finite-area mesh
        faMesh procMesh
        (
            areaName_,
            procFvMesh,
            labelList(procFaceLabels_[procI])
        );

        // Create processor boundary patches
        const labelList& curBoundaryAddressing =
            procBoundaryAddressing_[procI];

        const labelList& curPatchSizes = procPatchSize_[procI];

        const labelList& curNeighbourProcessors =
            procNeighbourProcessors_[procI];

        const labelList& curProcessorPatchSizes =
            procProcessorPatchSize_[procI];

        const labelListList& curPatchEdgeLabels =
            procPatchEdgeLabels_[procI];

        const faPatchList& meshPatches = boundary();

        faPatchList procPatches
        (
            curPatchSizes.size() + curProcessorPatchSizes.size()
        );

        label nPatches = 0;

        forAll(curPatchSizes, patchi)
        {
            const labelList& curEdgeLabels = curPatchEdgeLabels[nPatches];

            const label neiPolyPatchId =
                fvBoundaryProcAddressing.find
                (
                    meshPatches[curBoundaryAddressing[patchi]]
                    .ngbPolyPatchIndex()
                );

            procPatches.set
            (
                nPatches,
                meshPatches[curBoundaryAddressing[patchi]].clone
                (
                    procMesh.boundary(),
                    curEdgeLabels,
                    nPatches,
                    neiPolyPatchId
                )
            );
            ++nPatches;
        }

        forAll(curProcessorPatchSizes, procPatchI)
        {
            const labelList& curEdgeLabels = curPatchEdgeLabels[nPatches];

            procPatches.set
            (
                nPatches,
                new processorFaPatch
                (
                    curEdgeLabels,
                    nPatches,
                    procMesh.boundary(),
                    -1,
                    procI,
                    curNeighbourProcessors[procPatchI]
                )
            );

            ++nPatches;
        }

        // Add boundary patches
        procMesh.addFaPatches(procPatches);

        // More precision (for points data)
        IOstream::minPrecision(10);

        procMesh.write();

        // Statistics
        Info<< nl << "Processor " << procI;

        if (procMesh.nFaces())
        {
            Info<< nl << "    ";
        }
        else
        {
            Info<< ": ";
        }

        Info<< "Number of faces = " << procMesh.nFaces() << nl;

        if (procMesh.nFaces())
        {
            Info<< "    Number of points = " << procMesh.nPoints() << nl;
        }

        totProcFaces += procMesh.nFaces();
        maxProcFaces = max(maxProcFaces, procMesh.nFaces());

        label nBoundaryEdges = 0;
        label nProcPatches = 0;
        label nProcEdges = 0;

        for (const faPatch& fap : procMesh.boundary())
        {
            if (const auto* ppp = isA<processorFaPatch>(fap); ppp)
            {
                const auto& procPatch = *ppp;

                Info<< "    Number of edges shared with processor "
                    << procPatch.neighbProcNo() << " = "
                    << procPatch.size() << nl;

                nProcEdges += procPatch.size();
                ++nProcPatches;
            }
            else
            {
                nBoundaryEdges += fap.size();
            }
        }

        if (procMesh.nFaces() && (nBoundaryEdges || nProcEdges))
        {
            Info<< "    Number of processor patches = " << nProcPatches << nl
                << "    Number of processor edges = " << nProcEdges << nl
                << "    Number of boundary edges = " << nBoundaryEdges << nl;
        }

        totProcEdges += nProcEdges;
        totProcPatches += nProcPatches;
        maxProcEdges = Foam::max(maxProcEdges, nProcEdges);
        maxProcPatches = Foam::max(maxProcPatches, nProcPatches);

        // Write the addressing information
        IOobject ioAddr
        (
            "procAddressing",
            "constant",
            faMesh::meshSubDir,
            procMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // pointProcAddressing
        ioAddr.rename("pointProcAddressing");
        labelIOList::writeContents(ioAddr, procPatchPointAddressing_[procI]);

        // edgeProcAddressing
        ioAddr.rename("edgeProcAddressing");
        labelIOList::writeContents(ioAddr, procEdgeAddressing_[procI]);

        // faceProcAddressing
        ioAddr.rename("faceProcAddressing");
        labelIOList::writeContents(ioAddr, procFaceAddressing_[procI]);

        // boundaryProcAddressing
        ioAddr.rename("boundaryProcAddressing");
        labelIOList::writeContents(ioAddr, procBoundaryAddressing_[procI]);
    }


    // Summary stats
    Info<< nl
        << "Number of processor edges = " << (totProcEdges/2) << nl
        << "Max number of faces = " << maxProcFaces;

    if (maxProcFaces != totProcFaces)
    {
        scalar avgValue = scalar(totProcFaces)/nProcs_;

        Info<< " (" << 100.0*(maxProcFaces-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of processor patches = " << maxProcPatches;
    if (totProcPatches)
    {
        scalar avgValue = scalar(totProcPatches)/nProcs_;

        Info<< " (" << 100.0*(maxProcPatches-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of edges between processors = " << maxProcEdges;
    if (totProcEdges)
    {
        scalar avgValue = scalar(totProcEdges)/nProcs_;

        Info<< " (" << 100.0*(maxProcEdges-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl << endl;

    return true;
}


// ************************************************************************* //
