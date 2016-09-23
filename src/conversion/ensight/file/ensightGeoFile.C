/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "ensightGeoFile.H"
#include "foamVersion.H"

// Macros to stringify macro contents.
#define STRINGIFY(content)      #content
#define STRING_QUOTE(input)     STRINGIFY(input)


#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightGeoFile::initialize()
{
    #ifdef OPENFOAM_PLUS
    string desc2("Written by OpenFOAM+ " STRING_QUOTE(OPENFOAM_PLUS));
    #else
    string desc2("Written by OpenFOAM-" + string(Foam::FOAMversion));
    #endif

    writeBinaryHeader();
    write("Ensight Geometry File");  newline(); // description line 1
    write(desc2);                    newline(); // description line 2
    write("node id assign");         newline();
    write("element id assign");      newline();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightGeoFile::ensightGeoFile
(
    const fileName& pathname,
    IOstream::streamFormat format
)
:
    ensightFile(pathname, format)
{
    writeBinaryHeader();
    write("Ensight Geometry File");  newline();
    write("=====================");  newline();
    write("node id assign");         newline();
    write("element id assign");      newline();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightGeoFile::~ensightGeoFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::ensightGeoFile::writeKeyword(const string& key)
{
    // ensure we get ensightFile::write(const string&)
    write(static_cast<const string&>(key)); newline();

    return *this;
}


// ************************************************************************* //
