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

\*---------------------------------------------------------------------------*/

#include "linearSpringDamper.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{
    defineTypeNameAndDebug(linearSpringDamper, 0);

    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        linearSpringDamper,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::linearSpringDamper
(
    const word& name,
    const dictionary& sDoFRBMRDict
)
:
    sixDoFRigidBodyMotionRestraint(name, sDoFRBMRDict)
{
    read(sDoFRBMRDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::~linearSpringDamper()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    if (anchor_.empty())
    {
        anchor_.reset
        (
            new TimeFunction1<vector>
            (
                motion.time(),
                "anchor",
                coeffDict()
            )
        );
    }

    scalar t = motion.time().value();

    restraintPosition = motion.transform(refAttachmentPt_);

    // Current axis of the spring
    vector r = restraintPosition - anchor_->value(t);
    vector rDir = r/(mag(r) + VSMALL);

    vector v = motion.velocity(restraintPosition);
    //vector a = motion.a();

    scalar m = motion.mass();

    restraintForce = f0_;
    restraintMoment = Zero;

    scalar dt =  motion.time().deltaTValue();

    if (mag(r) > restLength_)
    {
        oldError_ = error_;
        oldErrorIntegral_ = errorIntegral_;
        error_ = (mag(r) - restLength_)/restLength_;
        errorIntegral_ =
            oldErrorIntegral_ + 0.5*dt*(error_ + oldError_); //

        scalar errorDifferential = (error_ - oldError_)/dt;

        scalar factor =
            P_*error_ + I_*errorIntegral_ + D_*errorDifferential;

        wn_ *= (factor + 1);
        //wn_ = 3.14/(15*dt);

        //scalar wn =  3.14/C_;
        scalar damping = psi_*2*m*wn_;
        scalar stiffness = sqr(wn_)*m;

        if (error_ > 1e-3)
        {
            scalar magRest = mag(restraintForce)*(factor + 1);
            restraintForce = -magRest*rDir;
            //restraintForce =
            //    - damping*(rDir & v)*rDir
            //    - stiffness*(mag(r) - restLength_)*rDir;
        }

        //restraintMoment = (restraintPosition ^ restraintForce);

        if (motion.report())
        {
            //Info<< " damping :"  << damping << endl;
            //Info<< " stiffness :"  << stiffness << endl;
            //Info<< " damping force:"  <<  -damping*(rDir & v)*rDir << endl;
            //Info<< " stiffness force:"
            //    <<  -stiffness*(mag(r) - restLength_)*rDir << endl;

            //Info<< " factor :"  <<  factor << endl;
            //Info<< " error :"  <<  error_ << endl;
            //Info<< " errorIntegral :"  <<  errorIntegral_ << endl;
            //Info<< " errorDifferential :"  <<  errorDifferential << endl;

            Info<< " restraintPosition :"  <<  restraintPosition << endl;
            Info<< " anchor :"  <<  anchor_->value(t) << endl;
        }
    }
    else
    {
        restraintForce = Zero;
    }
    /*
    else if (mag(r) > 0.9*restLength_)
    {
        restraintForce =
              damping*mag(motion.v())*rDir
            - stiffness/10*(mag(r) - restLength_)*rDir;
    }
    */


    if (motion.report())
    {
        //Info<< " force :"  << restraintForce << endl;
        //Info<< " distance :"  << mag(r) - restLength_ << endl;
        //Info<< " velocity :"  << v << endl;
    }
}


bool Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::read
(
    const dictionary& sDoFRBMRDict
)
{
    sixDoFRigidBodyMotionRestraint::read(sDoFRBMRDict);

    sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    sDoFRBMRCoeffs_.readEntry("psi", psi_);
    sDoFRBMRCoeffs_.readEntry("C", C_);
    sDoFRBMRCoeffs_.readEntry("restLength", restLength_);

    sDoFRBMRCoeffs_.readEntry("P", P_);
    sDoFRBMRCoeffs_.readEntry("I", I_);
    sDoFRBMRCoeffs_.readEntry("D", D_);
    sDoFRBMRCoeffs_.readEntry("f0", f0_);

    wn_ = 3.14/C_;


    return true;
}


void Foam::sixDoFRigidBodyMotionRestraints::linearSpringDamper::write
(
    Ostream& os
) const
{
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("psi", psi_);
    os.writeEntry("C", C_);
    os.writeEntry("restLength", restLength_);
}


// ************************************************************************* //
