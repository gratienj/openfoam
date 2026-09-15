#ifdef REF
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].clamp_min(0);
}
#else
tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    Qdot = reaction->Qdot();

    /*
     * ================================================================
     * Molecular Schmidt numbers
     * ================================================================
     */

    const scalar ScH2  = 0.23;
    const scalar ScO2  = 0.80;
    const scalar ScH2O = 0.58;
    const scalar ScN2  = 0.91;


    /*
     * ================================================================
     * Identify species from their names
     * ================================================================
     */

    if
    (
        !composition.species().found("H2")
     || !composition.species().found("O2")
     || !composition.species().found("H2O")
     || !composition.species().found("N2")
    )
    {
        FatalErrorInFunction
            << "The reactingFoamH2Air solver requires exactly the "
            << "following species:"
            << nl
            << "    H2"
            << nl
            << "    O2"
            << nl
            << "    H2O"
            << nl
            << "    N2"
            << nl
            << "Available species: "
            << composition.species()
            << exit(FatalError);
    }

    const label H2Index =
        composition.species().find("H2");

    const label O2Index =
        composition.species().find("O2");

    const label H2OIndex =
        composition.species().find("H2O");

    const label N2Index =
        composition.species().find("N2");


    /*
     * ================================================================
     * Solve ALL FOUR species equations.
     *
     * N2 is no longer reconstructed from:
     *
     *     YN2 = 1 - YH2 - YO2 - YH2O
     *
     * This is necessary to give N2 its own Schmidt number.
     * ================================================================
     */

    forAll(Y, i)
    {
        if (composition.active(i))
        {
            volScalarField& Yi = Y[i];

            scalar ScYi;


            /*
             * Species-dependent Schmidt number
             */

            if (i == H2Index)
            {
                ScYi = ScH2;
            }
            else if (i == O2Index)
            {
                ScYi = ScO2;
            }
            else if (i == H2OIndex)
            {
                ScYi = ScH2O;
            }
            else if (i == N2Index)
            {
                ScYi = ScN2;
            }
            else
            {
                FatalErrorInFunction
                    << "No Schmidt number defined for species "
                    << composition.species()[i]
                    << exit(FatalError);
            }


            /*
             * ========================================================
             * Molecular diffusion coefficient
             *
             *     Gamma_i = rho*D_i
             *
             * with
             *
             *     D_i = nu/Sc_i
             *
             * and
             *
             *     nu = mu/rho
             *
             * therefore
             *
             *     Gamma_i = mu/Sc_i
             * ========================================================
             */
            const tmp<volScalarField> tDiffusionFactor =
                reaction->diffusionFactor();

            const volScalarField& diffusionFactor =
                tDiffusionFactor();

            volScalarField GammaYi
            (
                IOobject
                (
                    "Gamma_" + Yi.name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                diffusionFactor*
                turbulence->mu()/ScYi
            );


            Info<< "Species " << Yi.name()
                << " : Sc = " << ScYi
                << ", minmax(diffusion) = " << gMinMax(diffusionFactor)
                << ", min(Gamma) = " << gMin(GammaYi)
                << ", max(Gamma) = " << gMax(GammaYi)
                << endl;


            /*
             * ========================================================
             * Species equation
             * ========================================================
             */

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(GammaYi, Yi)
             ==
                reaction->R(Yi)
              + fvOptions(rho, Yi)
            );


            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            /*
             * Prevent negative mass fractions.
             */

            Yi.clamp_min(0);
        }
    }


    /*
     * ================================================================
     * Normalization of the four mass fractions
     *
     *     sum(Yi) = 1
     *
     * This replaces the old inert-species reconstruction.
     * ================================================================
     */

    volScalarField Ysum
    (
        IOobject
        (
            "Ysum",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 1.0)
    );

    Ysum = scalar(0);

    forAll(Y, i)
    {
        if (composition.active(i))
        {
            Ysum += Y[i];
        }
    }


    /*
     * Normalize all four species.
     */

    forAll(Y, i)
    {
        if (composition.active(i))
        {
            Y[i] /= max(Ysum, scalar(SMALL));
        }
    }

    /*
     * Check conservation of mass fractions
     */

    Info<< "Ysum min/max = "
        << gMin(Ysum)
        << " "
        << gMax(Ysum)
        << endl;
}
#endif
