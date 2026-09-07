/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    A nd          |
  \\/     M anipulation     |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "consumptionSpeed.H"


// * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(consumptionSpeed, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::consumptionSpeed::consumptionSpeed
(
    const dictionary& dict
)
:
    omega0_(dict.get<scalar>("omega0")),
    eta_(dict.get<scalar>("eta")),
    sigmaExt_(dict.get<scalar>("sigmaExt")),
    omegaMin_(dict.get<scalar>("omegaMin"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::consumptionSpeed::~consumptionSpeed()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::consumptionSpeed::omega0Sigma
(
    scalar sigma,
    scalar a
) const
{
    /*
     * sigma is a strain-rate quantity.
     *
     * Only positive strain is used here.  The original implementation
     * contained:
     *
     *     1 - exp(eta*sigma)
     *
     * which rapidly becomes negative for positive sigma and therefore
     * collapses to omegaMin_.  We instead use a bounded exponential
     * decrease.
     */

    const scalar sigmaPos = max(sigma, scalar(0));

    if (sigmaPos >= sigmaExt_)
    {
        return omegaMin_;
    }

    scalar omega0 = a*omega0_*exp(-eta_*sigmaPos);

    /*
     * Protect against pathological values of a or omega0.
     */
    omega0 = max(omega0, omegaMin_);
    omega0 = min(omega0, a*omega0_);

    return omega0;
}


Foam::tmp<Foam::volScalarField>
Foam::consumptionSpeed::omega0Sigma
(
    const volScalarField& sigma
)
{
    auto tomega0 = volScalarField::New
    (
        "omega0",
        IOobject::NO_REGISTER,
        sigma.mesh(),
        dimensionedScalar
        (
            dimensionSet(1, -2, -1, 0, 0, 0, 0),
            Zero
        )
    );

    auto& omega0 = tomega0.ref();

    volScalarField::Internal& iomega0 = omega0;

    forAll(iomega0, celli)
    {
        iomega0[celli] =
            omega0Sigma
            (
                sigma[celli],
                1.0
            );
    }

    volScalarField::Boundary& bomega0 =
        omega0.boundaryFieldRef();

    forAll(bomega0, patchi)
    {
        forAll(bomega0[patchi], facei)
        {
            bomega0[patchi][facei] =
                omega0Sigma
                (
                    sigma.boundaryField()[patchi][facei],
                    1.0
                );
        }
    }

    return tomega0;
}


void Foam::consumptionSpeed::read
(
    const dictionary& dict
)
{
    dict.readEntry("omega0", omega0_);
    dict.readEntry("eta", eta_);
    dict.readEntry("sigmaExt", sigmaExt_);
    dict.readEntry("omegaMin", omegaMin_);
}


// ************************************************************************* //
