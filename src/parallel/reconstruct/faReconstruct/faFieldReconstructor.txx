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

#include "faFieldReconstructor.H"
#include "Time.H"
#include "emptyFaPatch.H"
#include "faPatchFields.H"
#include "faePatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::faFieldReconstructor::reconstructField
(
    const IOobject& fieldObject,
    const UPtrList<GeometricField<Type, faPatchField, areaMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nFaces());

    // Create the patch fields
    PtrList<faPatchField<Type>> patchFields(mesh_.boundary().size());

    // The patch starts (global mesh)
    const labelList gStarts(mesh_.boundary().patchStarts());

    forAll(procMeshes_, proci)
    {
        const auto& procField = procFields[proci];
        const auto& procMesh = procMeshes_[proci];

        // The (edge,face,boundary)ProcAddressing for the current procMesh
        const auto& edgeProcAddr = edgeProcAddressing_[proci];
        const auto& faceProcAddr = faceProcAddressing_[proci];
        const auto& boundaryProcAddr = boundaryProcAddressing_[proci];

        // Set the face values in the reconstructed field
        internalField.rmap(procField.internalField(), faceProcAddr);


        // The patch starts (local mesh)
        const labelList starts(procMesh.boundary().patchStarts());

        // Set the boundary patch values in the reconstructed field

        forAll(boundaryProcAddr, patchI)
        {
            // Get addressing slice for this patch
//             const auto cp =
//                 procMesh.boundary()[patchI].patchSlice(edgeProcAddr);
            const auto cp =
                edgeProcAddr.slice
                (
                    starts[patchI],
                    procMesh.boundary()[patchI].size()
                );


            // Get patch index of the original patch,
            // check if the boundary patch is not a processor patch
            if
            (
                const auto tgtPatchi = boundaryProcAddr[patchI];
                (tgtPatchi >= 0)
            )
            {
                // Regular patch. Fast looping

                if (!patchFields.test(tgtPatchi))
                {
                    patchFields.set
                    (
                        tgtPatchi,
                        faPatchField<Type>::New
                        (
                            procField.boundaryField()[patchI],
                            mesh_.boundary()[tgtPatchi],
                            faPatchField<Type>::Internal::null(),
                            faPatchFieldReconstructor
                            (
                                mesh_.boundary()[tgtPatchi].size(),
                                procField.boundaryField()[patchI].size()
                            )
                        )
                    );
                }

                const label tgtPatchStart = gStarts[tgtPatchi];
//                     mesh_.boundary()[tgtPatchi].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, edgeI)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
//                     reverseAddressing[edgeI] = cp[edgeI] - 1 - tgtPatchStart;
                    reverseAddressing[edgeI] = cp[edgeI] - tgtPatchStart;
                }

                patchFields[tgtPatchi].rmap
                (
                    procField.boundaryField()[patchI],
                    reverseAddressing
                );
            }
            else
            {
                // Processor patch

                const Field<Type>& curPatchField =
                    procField.boundaryField()[patchI];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, edgeI)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
//                     label tgtEdgei = cp[edgeI] - 1;
                    label tgtEdgei = cp[edgeI];

                    // The target edge
                    if (tgtEdgei < 0)
                    {
                        // Edge is incorrectly flipped - should not happen
                    }
                    else if (tgtEdgei < mesh_.nInternalEdges())
                    {
                        // Target edge is an internal edge - ignore
                    }
                    else
                    {
                        // Target edge is a boundary, find which one.

//                     label tgtPatchi =
//                         mesh_.boundary().whichPatch(tgtEdgei);

                        // Binary search in patch starts (with +1 to
                        // include the start in the comparison)

                        const label tgtPatchi =
                            Foam::findLower(gStarts, (tgtEdgei+1));

                        if (tgtPatchi < 0)
                        {
                            FatalErrorInFunction
                                << "Edge " << tgtEdgei
                                << " not found in any of the patches" << nl
                                << "The patches appear to be inconsistent"
                                   " with the mesh :" << endl
                                << abort(FatalError);
                        }

                        if (!patchFields.test(tgtPatchi))
                        {
                            patchFields.set
                            (
                                tgtPatchi,
                                faPatchField<Type>::New
                                (
                                    mesh_.boundary()[tgtPatchi].type(),
                                    mesh_.boundary()[tgtPatchi],
                                    faPatchField<Type>::Internal::null()
                                )
                            );
                        }

                        // add the edge
//                         label tgtPatchEdgei =
//                             mesh_.boundary()[tgtPatchi].whichEdge(tgtEdgei);

                        label tgtPatchEdgei = (tgtEdgei - gStarts[tgtPatchi]);

                        patchFields[tgtPatchi][tgtPatchEdgei] =
                            curPatchField[edgeI];
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchI)
    {
        // add empty patches
        if
        (
            isA<emptyFaPatch>(mesh_.boundary()[patchI])
         && !patchFields.test(patchI)
        )
        {
            patchFields.set
            (
                patchI,
                faPatchField<Type>::New
                (
                    faPatchFieldBase::emptyType(),
                    mesh_.boundary()[patchI],
                    faPatchField<Type>::Internal::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    auto tfield = tmp<GeometricField<Type, faPatchField, areaMesh>>::New
    (
        fieldObject,
        mesh_,
        procFields[0].dimensions(),
        internalField,
        patchFields
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;

}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::faFieldReconstructor::reconstructField
(
    const IOobject& fieldObject,
    const UPtrList<GeometricField<Type, faePatchField, edgeMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nInternalEdges());

    // Create the patch fields
    PtrList<faePatchField<Type>> patchFields(mesh_.boundary().size());


    // The patch starts (global mesh)
    const labelList gStarts(mesh_.boundary().patchStarts());

    // HACK: until we add flip information there is no other way
    // to track when the processor boundary value has flipped

    const auto& edgeOwner = mesh_.edgeOwner();

    Field<scalar> boundaryEdgeSigns;
    {
        label maxBndEdges = 0;
        for (const auto& m : procMeshes_)
        {
            maxBndEdges = Foam::max(maxBndEdges, m.nBoundaryEdges());
        }

        boundaryEdgeSigns.resize(maxBndEdges, scalar(1));
    }

    forAll(procMeshes_, proci)
    {
        const auto& procField = procFields[proci];
        const auto& procMesh = procMeshes_[proci];

        // The (edge,face,boundary)ProcAddressing for the current procMesh
        const auto& edgeProcAddr = edgeProcAddressing_[proci];
        const auto& faceProcAddr = faceProcAddressing_[proci];
        const auto& boundaryProcAddr = boundaryProcAddressing_[proci];

        // Set the edge values in the reconstructed field

        // It is necessary to create a copy of the addressing array to
        // take care of the face direction offset trick.
        //
        {
            labelList curAddr(edgeProcAddr);

//             forAll(curAddr, addrI)
//             {
//                 curAddr[addrI] -= 1;
//             }

            // Set the edge values in the reconstructed field
            internalField.rmap(procField.internalField(), curAddr);

            // HACK: track sign flips without any edge flip information!!

            label bndEdgei = 0;
            for
            (
                label edgei = procMesh.nInternalEdges();
                edgei < procMesh.nEdges();
                ++edgei
            )
            {
                // The corresponding owner face in the serial mesh:
                auto serialEdgei = edgeProcAddr[edgei];
                auto ownFacei = faceProcAddr[procMesh.edgeOwner()[edgei]];

                boundaryEdgeSigns[bndEdgei] =
                (
                    (edgeOwner[serialEdgei] == ownFacei) ? 1 : -1
                );

                ++bndEdgei;
            }
        }

        // The patch starts (local mesh)
        const labelList starts(procMesh.boundary().patchStarts());

        // Set the boundary patch values in the reconstructed field

        forAll(boundaryProcAddr, patchI)
        {
            // Get addressing slice for this patch
//             const auto cp =
//                 procMesh.boundary()[patchI].patchSlice(edgeProcAddr);

            const auto cp =
                edgeProcAddr.slice
                (
                    starts[patchI],
                    procMesh.boundary()[patchI].size()
                );

            // Get patch index of the original patch,
            // check if the boundary patch is not a processor patch
            if
            (
                const auto tgtPatchi = boundaryProcAddr[patchI];
                (tgtPatchi >= 0)
            )
            {
                // Regular patch. Fast looping

                if (!patchFields.test(tgtPatchi))
                {
                    patchFields.set
                    (
                        tgtPatchi,
                        faePatchField<Type>::New
                        (
                            procField.boundaryField()[patchI],
                            mesh_.boundary()[tgtPatchi],
                            faePatchField<Type>::Internal::null(),
                            faPatchFieldReconstructor
                            (
                                mesh_.boundary()[tgtPatchi].size(),
                                procField.boundaryField()[patchI].size()
                            )
                        )
                    );
                }

                const label tgtPatchStart = gStarts[tgtPatchi];
//                     mesh_.boundary()[tgtPatchi].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, edgeI)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
//                     reverseAddressing[faceI] = cp[faceI] - 1 - tgtPatchStart;
                    reverseAddressing[edgeI] = cp[edgeI] - tgtPatchStart;
                }

                patchFields[tgtPatchi].rmap
                (
                    procField.boundaryField()[patchI],
                    reverseAddressing
                );
            }
            else
            {
                // Processor patch

                const Field<Type>& curPatchField =
                    procField.boundaryField()[patchI];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop

                // For the loop:
                // - track boundaryEdgei on the proc-local mesh as addressing
                //   into boundaryEdgeSigns.
                // - track patchEdgei on the proc-local mesh patch.
                for
                (
                    label boundaryEdgei
                  = (starts[patchI] - procMesh.nInternalEdges()),
                    patchEdgei = 0;
                    (patchEdgei < cp.size());
                  ++patchEdgei, ++boundaryEdgei
                )
                {
                    // The target edge
//                     label tgtEdgei = cp[patchEdgei] - 1;
                    label tgtEdgei = cp[patchEdgei];

                    if (tgtEdgei < 0)
                    {
                        // Edge is incorrectly flipped - should not happen
                    }
                    else if (tgtEdgei < mesh_.nInternalEdges())
                    {
                        // Target edge is an internal edge

                        // Processor patch -> internal face
                        // TBD: avoid copying with sign change, which
                        // would preserve the owner side and ignore
                        // (de-duplicate) the neigbour side
                        internalField[tgtEdgei] =
                        (
                            curPatchField[patchEdgei]
                            // HACK: sign flip without edge flip info!
                          * boundaryEdgeSigns[boundaryEdgei]
                        );
                    }
                    else
                    {
                        // Target edge is a boundary, find which one.

//                         label tgtPatchi =
//                             mesh_.boundary().whichPatch(tgtEdgei);

                        // Binary search in patch starts (with +1 to
                        // include the start in the comparison)

                        const label tgtPatchi =
                            Foam::findLower(gStarts, (tgtEdgei+1));

                        if (tgtPatchi < 0)
                        {
                            FatalErrorInFunction
                                << "Edge " << tgtEdgei
                                << " not found in any of the patches" << nl
                                << "The patches appear to be inconsistent"
                                   " with the mesh :" << endl
                                << abort(FatalError);
                        }

                        if (!patchFields.test(tgtPatchi))
                        {
                            patchFields.set
                            (
                                tgtPatchi,
                                faePatchField<Type>::New
                                (
                                    mesh_.boundary()[tgtPatchi].type(),
                                    mesh_.boundary()[tgtPatchi],
                                    faePatchField<Type>::Internal::null()
                                )
                            );
                        }

                        // add the value
//                         label tgtPatchEdgei =
//                             mesh_.boundary()[tgtPatchi].whichEdge(tgtEdgei);

                        label tgtPatchEdgei(tgtEdgei - gStarts[tgtPatchi]);

                        patchFields[tgtPatchi][tgtPatchEdgei] =
                            curPatchField[patchEdgei];
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchI)
    {
        // add empty patches
        if
        (
            isA<emptyFaPatch>(mesh_.boundary()[patchI])
         && !patchFields.test(patchI)
        )
        {
            patchFields.set
            (
                patchI,
                faePatchField<Type>::New
                (
                    faePatchFieldBase::emptyType(),
                    mesh_.boundary()[patchI],
                    faePatchField<Type>::Internal::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    auto tfield = tmp<GeometricField<Type, faePatchField, edgeMesh>>::New
    (
        fieldObject,
        mesh_,
        procFields[0].dimensions(),
        internalField,
        patchFields
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::faFieldReconstructor::reconstructAreaField
(
    const IOobject& fieldObject
)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, faPatchField, areaMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.emplace_set
        (
            proci,
            IOobject
            (
                fieldObject.name(),
                procMeshes_[proci].thisDb().time().timeName(),
                procMeshes_[proci].thisDb(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            procMeshes_[proci]
        );
    }

    return reconstructField
    (
        IOobject
        (
            fieldObject.name(),
            mesh_.thisDb().time().timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procFields
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::faFieldReconstructor::reconstructEdgeField
(
    const IOobject& fieldObject
)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, faePatchField, edgeMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.emplace_set
        (
            proci,
            IOobject
            (
                fieldObject.name(),
                procMeshes_[proci].thisDb().time().timeName(),
                procMeshes_[proci].thisDb(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            procMeshes_[proci]
        );
    }

    return reconstructField
    (
        IOobject
        (
            fieldObject.name(),
            mesh_.thisDb().time().timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procFields
    );
}


template<class Type>
Foam::label Foam::faFieldReconstructor::reconstructAreaFields
(
    const UPtrList<const IOobject>& fieldObjects
)
{
    typedef GeometricField<Type, faPatchField, areaMesh> fieldType;

    label nFields = 0;

    for (const IOobject& io : fieldObjects)
    {
        if (io.isHeaderClass<fieldType>())
        {
            if (verbose_)
            {
                if (!nFields)
                {
                    Info<< "    Reconstructing "
                        << fieldType::typeName << "s\n" << nl;
                }
                Info<< "        " << io.name() << endl;
            }
            ++nFields;

            reconstructAreaField<Type>(io)().write();
            ++nReconstructed_;
        }
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::faFieldReconstructor::reconstructEdgeFields
(
    const UPtrList<const IOobject>& fieldObjects
)
{
    typedef GeometricField<Type, faePatchField, edgeMesh> fieldType;

    label nFields = 0;

    for (const IOobject& io : fieldObjects)
    {
        if (io.isHeaderClass<fieldType>())
        {
            if (verbose_)
            {
                if (!nFields)
                {
                    Info<< "    Reconstructing "
                        << fieldType::typeName << "s\n" << nl;
                }
                Info<< "        " << io.name() << endl;
            }
            ++nFields;

            reconstructEdgeField<Type>(io)().write();
            ++nReconstructed_;
        }
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::faFieldReconstructor::reconstructAreaFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef GeometricField<Type, faPatchField, areaMesh> fieldType;

    return reconstructAreaFields<Type>
    (
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    );
}


template<class Type>
Foam::label Foam::faFieldReconstructor::reconstructEdgeFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef GeometricField<Type, faePatchField, edgeMesh> fieldType;

    return reconstructEdgeFields<Type>
    (
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    );
}


// ************************************************************************* //
