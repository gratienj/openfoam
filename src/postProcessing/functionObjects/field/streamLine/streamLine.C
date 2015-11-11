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

\*---------------------------------------------------------------------------*/

#include "streamLine.H"
#include "fvMesh.H"
#include "streamLineParticleCloud.H"
#include "sampledSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(streamLine, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::streamLine::track()
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(obr_);

    IDLList<streamLineParticle> initialParticles;
    streamLineParticleCloud particles
    (
        mesh,
        cloudName_,
        initialParticles
    );

    const sampledSet& seedPoints = sampledSetPtr_();

    forAll(seedPoints, i)
    {
        particles.addParticle
        (
            new streamLineParticle
            (
                mesh,
                seedPoints[i],
                seedPoints.cells()[i],
                lifeTime_               // lifetime
            )
        );
    }

    label nSeeds = returnReduce(particles.size(), sumOp<label>());

    if (log_) Info<< "    seeded " << nSeeds << " particles" << endl;

    // Read or lookup fields
    PtrList<volScalarField> vsFlds;
    PtrList<interpolation<scalar> > vsInterp;
    PtrList<volVectorField> vvFlds;
    PtrList<interpolation<vector> > vvInterp;
    label UIndex = -1;

    if (loadFromFiles_)
    {
        IOobjectList allObjects(mesh, runTime.timeName());

        IOobjectList objects(2*fields_.size());
        forAll(fields_, i)
        {
            objects.add(*allObjects[fields_[i]]);
        }

        ReadFields(mesh, objects, vsFlds);
        vsInterp.setSize(vsFlds.size());
        forAll(vsFlds, i)
        {
            vsInterp.set
            (
                i,
                interpolation<scalar>::New
                (
                    interpolationScheme_,
                    vsFlds[i]
                )
            );
        }
        ReadFields(mesh, objects, vvFlds);
        vvInterp.setSize(vvFlds.size());
        forAll(vvFlds, i)
        {
            vvInterp.set
            (
                i,
                interpolation<vector>::New
                (
                    interpolationScheme_,
                    vvFlds[i]
                )
            );
        }
    }
    else
    {
        label nScalar = 0;
        label nVector = 0;

        forAll(fields_, i)
        {
            if (mesh.foundObject<volScalarField>(fields_[i]))
            {
                nScalar++;
            }
            else if (mesh.foundObject<volVectorField>(fields_[i]))
            {
                nVector++;
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find field " << fields_[i] << nl
                    << "Valid scalar fields are:"
                    << mesh.names(volScalarField::typeName) << nl
                    << "Valid vector fields are:"
                    << mesh.names(volVectorField::typeName)
                    << exit(FatalError);
            }
        }
        vsInterp.setSize(nScalar);
        nScalar = 0;
        vvInterp.setSize(nVector);
        nVector = 0;

        forAll(fields_, i)
        {
            if (mesh.foundObject<volScalarField>(fields_[i]))
            {
                const volScalarField& f = mesh.lookupObject<volScalarField>
                (
                    fields_[i]
                );
                vsInterp.set
                (
                    nScalar++,
                    interpolation<scalar>::New
                    (
                        interpolationScheme_,
                        f
                    )
                );
            }
            else if (mesh.foundObject<volVectorField>(fields_[i]))
            {
                const volVectorField& f = mesh.lookupObject<volVectorField>
                (
                    fields_[i]
                );

                if (f.name() == UName_)
                {
                    UIndex = nVector;
                }

                vvInterp.set
                (
                    nVector++,
                    interpolation<vector>::New
                    (
                        interpolationScheme_,
                        f
                    )
                );
            }
        }
    }

    // Store the names
    scalarNames_.setSize(vsInterp.size());
    forAll(vsInterp, i)
    {
        scalarNames_[i] = vsInterp[i].psi().name();
    }
    vectorNames_.setSize(vvInterp.size());
    forAll(vvInterp, i)
    {
        vectorNames_[i] = vvInterp[i].psi().name();
    }

    // Check that we know the index of U in the interpolators.

    if (UIndex == -1)
    {
        FatalErrorInFunction
            << "Cannot find field to move particles with : " << UName_ << nl
            << "This field has to be present in the sampled fields " << fields_
            << " and in the objectRegistry."
            << exit(FatalError);
    }

    // Sampled data
    // ~~~~~~~~~~~~

    // Size to maximum expected sizes.
    allTracks_.clear();
    allTracks_.setCapacity(nSeeds);
    allScalars_.setSize(vsInterp.size());
    forAll(allScalars_, i)
    {
        allScalars_[i].clear();
        allScalars_[i].setCapacity(nSeeds);
    }
    allVectors_.setSize(vvInterp.size());
    forAll(allVectors_, i)
    {
        allVectors_[i].clear();
        allVectors_[i].setCapacity(nSeeds);
    }


    // Additional particle info
    streamLineParticle::trackingData td
    (
        particles,
        vsInterp,
        vvInterp,
        UIndex,         // index of U in vvInterp
        trackForward_,  // track in +u direction?
        nSubCycle_,     // automatic track control:step through cells in steps?
        trackLength_,   // fixed track length

        allTracks_,
        allScalars_,
        allVectors_
    );


    // Set very large dt. Note: cannot use GREAT since 1/GREAT is SMALL
    // which is a trigger value for the tracking...
    const scalar trackTime = Foam::sqrt(GREAT);

    // Track
    particles.move(td, trackTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamLine::streamLine
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    streamLineBase(name, obr, dict, loadFromFiles)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (setActive<fvMesh>())
    {
        read(dict_);
    }
    else
    {
        active_ = false;
        WarningInFunction
            << "No fvMesh available, deactivating."
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::streamLine::~streamLine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::streamLine::read(const dictionary& dict)
{
    if (active_)
    {
        Info<< type() << " " << name_ << ":" << nl;

        //dict_ = dict;
        dict.lookup("fields") >> fields_;
        if (dict.found("UName"))
        {
            dict.lookup("UName") >> UName_;
        }
        else
        {
            UName_ = "U";
            if (dict.found("U"))
            {
                IOWarningInFunction(dict)
                    << "Using deprecated entry \"U\"."
                    << " Please use \"UName\" instead."
                    << endl;
                dict.lookup("U") >> UName_;
            }
        }

        if (findIndex(fields_, UName_) == -1)
        {
            FatalIOErrorInFunction(dict)
                << "Velocity field for tracking " << UName_
                << " should be present in the list of fields " << fields_
                << exit(FatalIOError);
        }


        dict.lookup("trackForward") >> trackForward_;
        dict.lookup("lifeTime") >> lifeTime_;
        if (lifeTime_ < 1)
        {
            FatalErrorInFunction
                << "Illegal value " << lifeTime_ << " for lifeTime"
                << exit(FatalError);
        }


        bool subCycling = dict.found("nSubCycle");
        bool fixedLength = dict.found("trackLength");

        if (subCycling && fixedLength)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot both specify automatic time stepping (through '"
                << "nSubCycle' specification) and fixed track length (through '"
                << "trackLength')"
                << exit(FatalIOError);
        }

        nSubCycle_ = 1;
        if (dict.readIfPresent("nSubCycle", nSubCycle_))
        {
            trackLength_ = VGREAT;
            if (nSubCycle_ < 1)
            {
                nSubCycle_ = 1;
            }

            if (log_) Info
                << "    automatic track length specified through"
                << " number of sub cycles : " << nSubCycle_ << nl
                << endl;
        }
    }
}


// ************************************************************************* //
