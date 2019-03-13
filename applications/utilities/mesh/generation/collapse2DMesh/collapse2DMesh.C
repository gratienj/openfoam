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

Application
    collapse2DMesh

Group
    grpMeshGenerationUtilities

Description
    Takes 3D mesh with empty front and back and generate 2D version

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "topoSet.H"
#include "processorMeshes.H"
#include "emptyPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Converts 2D polyMesh into polyMesh2D"
    );

    //argList::addArgument("patch", "The (empty) patch giving the geometry");
    argList::validArgs.append("patch");

    #include "addOverwriteOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const word patchName = args[1];
    //const bool overwrite = args.found("overwrite");
    const bool overwrite = args.optionFound("overwrite");

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const label keepPatchi = pbm.findPatchID(patchName);
    if (keepPatchi == -1 || !isA<emptyPolyPatch>(pbm[keepPatchi]))
    {
        FatalErrorInFunction
            << "Patch " << patchName << " does not exist or is not empty"
            << exit(FatalError);
    }

    // Constructing a 2D mesh
    polyTopoChange meshMod(mesh, true, 2);

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        face f(2);
        f[0] = mesh.faces()[facei][0];
        f[1] = mesh.faces()[facei][1];
        meshMod.modifyFace
        (
            f,
            facei,
            own[facei],
            nei[facei],
            false,
            -1,
            -1,
            false
        );
    }
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        if (!isA<emptyPolyPatch>(pp))   // || patchi == keepPatchi)
        {
            forAll(pp, patchFacei)
            {
                face f(2);
                f[0] = pp[patchFacei][0];
                f[1] = pp[patchFacei][1];
                label meshFacei = pp.start()+patchFacei;
                meshMod.modifyFace
                (
                    f,
                    meshFacei,
                    own[meshFacei],
                    -1,
                    false,
                    patchi,
                    -1,
                    false
                );
            }
        }
        else
        {
            forAll(pp, patchFacei)
            {
                meshMod.removeFace(pp.start()+patchFacei, -1);
            }
        }
    }

    // Create a mesh from topo changes.
    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    mesh.updateMesh(morphMap);

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
        mesh.setInstance(oldInstance);
    }

    Info<< "\nWriting collapsed mesh to time = " << runTime.timeName()
        << nl << endl;

    mesh.write();
    //topoSet::removeFiles(mesh);   // only (boundary)faceSets are wrong
    processorMeshes::removeFiles(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
