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

#include "forceCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffs, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffs, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffs::writeFileHeader(const label i)
{
    // Note: Only possible to create bin files after bins have been initialised

    if (writeToFile() && !coeffFilePtr_.valid())
    {
        coeffFilePtr_ = createFile("coefficient");
        writeIntegratedHeader("Coefficients", coeffFilePtr_());

        if (nBin_ > 1)
        {
            CmBinFilePtr_ = createFile("CmBin");
            writeBinHeader("Moment coefficient bins", CmBinFilePtr_());
            CdBinFilePtr_ = createFile("CdBin");
            writeBinHeader("Drag coefficient bins", CdBinFilePtr_());
            ClBinFilePtr_ = createFile("ClBin");
            writeBinHeader("Lift coefficient bins", ClBinFilePtr_());
        }
    }
}


void Foam::forceCoeffs::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, "Force coefficients");
    writeHeaderValue(os, "liftDir", liftDir_);
    writeHeaderValue(os, "dragDir", dragDir_);
    writeHeaderValue(os, "pitchAxis", pitchAxis_);
    writeHeaderValue(os, "magUInf", magUInf_);
    writeHeaderValue(os, "lRef", lRef_);
    writeHeaderValue(os, "Aref", Aref_);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "Cm");
    writeTabbed(os, "Cd");
    writeTabbed(os, "Cl");
    writeTabbed(os, "Cl(f)");
    writeTabbed(os, "Cl(r)");
    os  << endl;
}


void Foam::forceCoeffs::writeBinHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "bins", nBin_);
    writeHeaderValue(os, "start", binMin_);
    writeHeaderValue(os, "delta", binDx_);
    writeHeaderValue(os, "direction", binDir_);

    vectorField binPoints(nBin_);
    writeCommented(os, "x co-ords  :");
    forAll(binPoints, pointI)
    {
        binPoints[pointI] = (binMin_ + (pointI + 1)*binDx_)*binDir_;
        os << tab << binPoints[pointI].x();
    }
    os << nl;

    writeCommented(os, "y co-ords  :");
    forAll(binPoints, pointI)
    {
        os << tab << binPoints[pointI].y();
    }
    os << nl;

    writeCommented(os, "z co-ords  :");
    forAll(binPoints, pointI)
    {
        os << tab << binPoints[pointI].z();
    }
    os << nl;

    writeHeader(os, "");
    writeCommented(os, "Time");

        vectorField binPoints(nBin_);
        writeCommented(file(i), "x co-ords  :");
        forAll(binPoints, pointi)
        {
            binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
            file(i) << tab << binPoints[pointi].x();
        }
        file(i) << nl;

        writeCommented(file(i), "y co-ords  :");
        forAll(binPoints, pointi)
        {
            file(i) << tab << binPoints[pointi].y();
        }
    }

    os  << endl;
}

        writeCommented(file(i), "z co-ords  :");
        forAll(binPoints, pointi)
        {
            file(i) << tab << binPoints[pointi].z();
        }
        file(i) << nl;

void Foam::forceCoeffs::writeIntegratedData
(
    const word& title,
    const List<Field<scalar> >& coeff
) const
{
    scalar pressure = sum(coeff[0]);
    scalar viscous = sum(coeff[1]);
    scalar porous = sum(coeff[2]);
    scalar total = pressure + viscous + porous;

    if (log_)
    {
        Info<< "        " << title << "       : " << total << token::TAB
            << "("
            << "pressure: " << pressure << token::TAB
            << "viscous: " << viscous;

        if (porosity_)
        {
            Info<< token::TAB << "porous: " << porous;
        }

        Info<< ")" << endl;
    }
}


void Foam::forceCoeffs::writeBinData
(
    const List<Field<scalar> > coeffs,
    Ostream& os
) const
{
    os  << obr_.time().value();

    for (label binI = 0; binI < nBin_; binI++)
    {
        FatalErrorInFunction
            << "Unhandled file index: " << i
            << abort(FatalError);
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::forceCoeffs
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces(name, runTime, dict),
    liftDir_(Zero),
    dragDir_(Zero),
    pitchAxis_(Zero),
    magUInf_(0.0),
    lRef_(0.0),
    Aref_(0.0),
    coeffFilePtr_(),
    CmBinFilePtr_(),
    CdBinFilePtr_(),
    ClBinFilePtr_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffs::~forceCoeffs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffs::read(const dictionary& dict)
{
    forces::read(dict);

    // Directions for lift and drag forces, and pitch moment
    dict.lookup("liftDir") >> liftDir_;
    dict.lookup("dragDir") >> dragDir_;
    dict.lookup("pitchAxis") >> pitchAxis_;

    // Free stream velocity magnitude
    dict.lookup("magUInf") >> magUInf_;

    // Reference length and area scales
    dict.lookup("lRef") >> lRef_;
    dict.lookup("Aref") >> Aref_;

    return true;
}


bool Foam::functionObjects::forceCoeffs::execute(const bool postProcess)
{
    return true;
}


bool Foam::functionObjects::forceCoeffs::write(const bool postProcess)
{
    forces::calcForcesMoment();

    if (Pstream::master())
    {
        writeFiles::write();

    // Calculate coefficients
    scalar CmTot = 0;
    scalar CdTot = 0;
    scalar ClTot = 0;
    forAll(liftCoeffs, i)
    {
        momentCoeffs[i] = (moment_[i] & pitchAxis_)/(Aref_*pDyn*lRef_);
        dragCoeffs[i] = (force_[i] & dragDir_)/(Aref_*pDyn);
        liftCoeffs[i] = (force_[i] & liftDir_)/(Aref_*pDyn);

        CmTot += sum(momentCoeffs[i]);
        CdTot += sum(dragCoeffs[i]);
        ClTot += sum(liftCoeffs[i]);
    }

        List<Field<scalar>> coeffs(3);
        coeffs[0].setSize(nBin_);
        coeffs[1].setSize(nBin_);
        coeffs[2].setSize(nBin_);

    if (log_) Info
        << type() << " " << name_ << " output:" << nl
        << "    Coefficients" << nl;

    writeIntegratedData("Cm", momentCoeffs);
    writeIntegratedData("Cd", dragCoeffs);
    writeIntegratedData("Cl", liftCoeffs);

    if (log_) Info
        << "        Cl(f)    : " << ClfTot << nl
        << "        Cl(r)    : " << ClrTot << nl
        << endl;

        writeTime(file(0));
        file(0)
            << tab << Cm << tab  << Cd
            << tab << Cl << tab << Clf << tab << Clr << endl;

        Log << type() << " " << name() << " output:" << nl
            << "    Cm    = " << Cm << nl
            << "    Cd    = " << Cd << nl
            << "    Cl    = " << Cl << nl
            << "    Cl(f) = " << Clf << nl
            << "    Cl(r) = " << Clr << endl;

        if (nBin_ > 1)
        {
            if (binCumulative_)
            {
                forAll(liftCoeffs, i)
                {
                    for (label binI = 1; binI < nBin_; binI++)
                    {
                        liftCoeffs[i][binI] += liftCoeffs[i][binI-1];
                        dragCoeffs[i][binI] += dragCoeffs[i][binI-1];
                        momentCoeffs[i][binI] += momentCoeffs[i][binI-1];
                    }
                }
            }

            writeTime(file(1));

    // Write state/results information
    {
        setResult("Cm", CmTot);
        setResult("Cd", CdTot);
        setResult("Cl", ClTot);
        setResult("Cl(f)", ClfTot);
        setResult("Cl(r)", ClrTot);
    }

    if (writeFields_)
    {
        const volVectorField& force =
            obr_.lookupObject<volVectorField>(fieldName("force"));

        const volVectorField& moment =
            obr_.lookupObject<volVectorField>(fieldName("moment"));

        volVectorField& forceCoeff =
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>(fieldName("forceCoeff"))
            );

        volVectorField& momentCoeff =
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>(fieldName("momentCoeff"))
            );

        dimensionedScalar f0("f0", dimForce, Aref_*pDyn);
        dimensionedScalar m0("m0", dimForce*dimLength, Aref_*lRef_*pDyn);

        forceCoeff == force/f0;
        momentCoeff == moment/m0;
    }
}


void Foam::forceCoeffs::end()
{
    // Do nothing
}


void Foam::forceCoeffs::timeSet()
{
    // Do nothing
}


void Foam::forceCoeffs::write()
{
    if (!active_)
    {
        return;
    }

    if (writeFields_)
    {
        const volVectorField& forceCoeff =
            obr_.lookupObject<volVectorField>(fieldName("forceCoeff"));

        const volVectorField& momentCoeff =
            obr_.lookupObject<volVectorField>(fieldName("momentCoeff"));

        Log << endl;
    }

    return true;
}


// ************************************************************************* //
