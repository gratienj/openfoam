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

#include "patchProbes.H"
#include "volFields.H"
#include "IOmanip.H"
#include "mappedPatchBase.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchProbes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchProbes,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchProbes::findElements(const fvMesh& mesh)
{
    (void)mesh.tetBasePtIs();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    label patchi = bm.findPatchID(patchName_);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << " Unknown patch name "
            << patchName_ << endl
            << exit(FatalError);
    }

    const labelList patchIDs(bm.patchSet(patchNames_).sortedToc());

    const polyPatch& pp = bm[patchi];

    if (nFaces > 0)
    {
        // Collect mesh faces and bounding box
        labelList bndFaces(nFaces);
        treeBoundBox overallBb(treeBoundBox::invertedBox);

        nFaces = 0;
        forAll(patchIDs, i)
        {
            const polyPatch& pp = bm[patchIDs[i]];
            forAll(pp, i)
            {
                bndFaces[nFaces++] = pp.start()+i;
                const face& f = pp[i];
                forAll(f, fp)
                {
                    const point& pt = pp.points()[f[fp]];
                    overallBb.min() = min(overallBb.min(), pt);
                    overallBb.max() = max(overallBb.max(), pt);
                }
            }
        }

        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1e-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        const indexedOctree<treeDataFace> boundaryTree
        (
            treeDataFace    // all information needed to search faces
            (
                false,                      // do not cache bb
                mesh,
                bndFaces                    // patch faces only
            ),
            overallBb,                      // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );


        forAll(probeLocations(), probei)
        {
            const point sample = probeLocations()[probei];

            scalar span = boundaryTree.bb().mag();

            pointIndexHit info = boundaryTree.findNearest
            (
                sample,
                Foam::sqr(span)
            );

            if (!info.hit())
            {
                info = boundaryTree.findNearest(sample, Foam::sqr(GREAT));
            }

            label facei = boundaryTree.shapes().faceLabels()[info.index()];

            const label patchi = bm.whichPatch(facei);

            if (isA<emptyPolyPatch>(bm[patchi]))
            {
                WarningInFunction
                << " The sample point: " << sample
                << " belongs to " << patchi
                << " which is an empty patch. This is not permitted. "
                << " This sample will not be included "
                << endl;
            }
            else if (info.hit())
            {
                const point& fc = mesh.faceCentres()[facei];

                mappedPatchBase::nearInfo sampleInfo;

                sampleInfo.first() = pointIndexHit
                (
                    true,
                    fc,
                    facei
                );

                sampleInfo.second().first() = magSqr(facePt-sample);
                sampleInfo.second().second() = Pstream::myProcNo();

                nearest[probei]= sampleInfo;
            }
        }
    }


    // Find nearest.
    Pstream::listCombineGather(nearest, mappedPatchBase::nearestEqOp());
    Pstream::listCombineScatter(nearest);


    // Update actual probe locations
    forAll(nearest, sampleI)
    {
        operator[](sampleI) = nearest[sampleI].first().rawPoint();
    }


    if (debug)
    {
        InfoInFunction << endl;
        forAll(nearest, sampleI)
        {
            label procI = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< operator[](sampleI)
                << " found on processor:" << procI
                << " in local face:" << localI
                << " with location:" << nearest[sampleI].first().rawPoint()
                << endl;
        }
    }


    // Extract any local faces to sample
    elementList_.setSize(nearest.size());
    elementList_ = -1;
    faceList_.setSize(nearest.size());
    faceList_ = -1;

    forAll(nearest, sampleI)
    {
        if (nearest[sampleI].second().second() == Pstream::myProcNo())
        {
            // Store the face to sample
            faceList_[sampleI] = nearest[sampleI].first().index();
        }
    }
}


void Foam::patchProbes::readDict(const dictionary& dict)
{
    if (!dict.readIfPresent("patches", patchNames_))
    {
        word patchName(dict.lookup("patchName"));
        patchNames_ = wordReList(1, wordRe(patchName));
    }
    probes::readDict(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchProbes::patchProbes
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    probes(name, t, dict)
{
    // When constructing probes above it will have called the
    // probes::findElements (since the virtual mechanism not yet operating).
    // Not easy to workaround (apart from feeding through flag into constructor)
    // so clear out any cells found for now.
    elementList_.clear();
    faceList_.clear();

    read(dict);
}


Foam::patchProbes::patchProbes
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool doFindElements
)
:
    probes(name, obr, dict)
{
    readDict(dict);

    if (doFindElements)
    {
        // Find the elements
        findElements(mesh_);

        // Open the probe streams
        prepare();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchProbes::~patchProbes()
{}


bool Foam::patchProbes::write()
{
    if (this->size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


bool Foam::patchProbes::read(const dictionary& dict)
{
    dict.lookup("patchName") >> patchName_;
    return probes::read(dict);
}


// ************************************************************************* //
