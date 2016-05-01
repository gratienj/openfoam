/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "cloudInfo.H"
#include "dictionary.H"
#include "kinematicCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudInfo, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cloudInfo::writeFileHeader(const label i)
{
    writeHeader(os, "Cloud information");
    writeCommented(os, "Time");
    writeTabbed(os, "nParcels");
    writeTabbed(os, "mass");
    writeTabbed(os, "Dmax");
    writeTabbed(os, "D10");
    writeTabbed(os, "D32");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::cloudInfo
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFiles(obr, name),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    cloudNames_(),
    filePtrs_()
{
    read(dict);
}


Foam::autoPtr<Foam::functionObjects::cloudInfo>
Foam::functionObjects::cloudInfo::New
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
{
    return autoPtr<cloudInfo>
    (
        new cloudInfo(name, obr, dict, loadFromFiles)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::~cloudInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::cloudInfo::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFiles::resetNames(dict.lookup("clouds"));

        log_ = dict.lookupOrDefault<Switch>("log", true);
        dict.lookup("clouds") >> cloudNames_;

        if (log_)
        {
            Info<< type() << " " << name_ << ": ";

            if (cloudNames_.size())
            {
                Info<< "applying to clouds:" << nl;
                forAll(cloudNames_, i)
                {
                    Info<< "    " << cloudNames_[i] << nl;
                }
                Info<< endl;
            }
            else
            {
                Info<< "no clouds to be processed" << nl << endl;
            }
        }

        if (writeToFile())
        {
            filePtrs_.setSize(cloudNames_.size());
            forAll(filePtrs_, fileI)
            {
                const word& cloudName = cloudNames_[fileI];
                filePtrs_.set(fileI, createFile(cloudName));
                writeFileHeader(filePtrs_[fileI]);
            }
        }
    }
}


void Foam::functionObjects::cloudInfo::execute()
{}


void Foam::functionObjects::cloudInfo::end()
{}


void Foam::functionObjects::cloudInfo::timeSet()
{}


void Foam::functionObjects::cloudInfo::write()
{
    if (active_)
    {
        functionObjectFiles::write();

        forAll(names(), i)
        {
            const word& cloudName = cloudNames_[cloudI];

            const kinematicCloud& cloud =
                obr_.lookupObject<kinematicCloud>(cloudName);

            label nParcels = returnReduce(cloud.nParcels(), sumOp<label>());
            scalar massInSystem =
                returnReduce(cloud.massInSystem(), sumOp<scalar>());

            scalar Dmax = cloud.Dmax();
            scalar D10 = cloud.Dij(1, 0);
            scalar D32 = cloud.Dij(3, 2);

            if (Pstream::master() && writeToFile())
            {
                writeTime(file(i));
                file(i)
                    << token::TAB
                    << nParcels << token::TAB
                    << massInSystem << token::TAB
                    << Dmax << token::TAB
                    << D10 << token::TAB
                    << D32 << token::TAB
                    << endl;
            }

            if (log_) Info
                << type() << " " << name_ <<  " output:" << nl
                << "    number of parcels : " << nParcels << nl
                << "    mass in system    : " << massInSystem << nl
                << "    maximum diameter  : " << Dmax << nl
                << "    D10 diameter      : " << D10 << nl
                << "    D32 diameter      : " << D32 << nl
                << endl;
        }
    }
}


// ************************************************************************* //
