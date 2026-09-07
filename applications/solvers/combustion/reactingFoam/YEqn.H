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
#ifdef TEST_REACTIVEFOAMDY
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

    /*
     * ================================================================
     * Species-dependent molecular Schmidt numbers
     *
     * Reference:
     * Male et al.:
     *
     *     Sc(H2)  = 0.23
     *     Sc(O2)  = 0.80
     *     Sc(H2O) = 0.58
     *     Sc(N2)  = 0.91
     *
     * Diagonal/Fickian approximation:
     *
     *     D_i = nu/Sc_i
     *
     * Since nu = mu/rho:
     *
     *     rho*D_i = mu/Sc_i
     *
     * Hence the conservative species equation contains:
     *
     *     div(mu/Sc_i grad(Y_i))
     *
     * ================================================================
     */

    const scalar ScH2  = 0.23;
    const scalar ScO2  = 0.80;
    const scalar ScH2O = 0.58;
    const scalar ScN2  = 0.91;


    /*
     * ---------------------------------------------------------------
     * Find species indices from their names.
     *
     * Do NOT assume:
     *
     *     H2  = 0
     *     O2  = 1
     *     H2O = 2
     *     N2  = 3
     *
     * because the CHEMKIN/OpenFOAM conversion may reorder species.
     * ---------------------------------------------------------------
     */

    const wordList requiredSpecies
    {
        "H2",
        "O2",
        "H2O",
        "N2"
    };

    forAll(requiredSpecies, s)
    {
        if (!composition.species().found(requiredSpecies[s]))
        {
            FatalErrorInFunction
                << "Required species "
                << requiredSpecies[s]
                << " not found in the thermodynamic mixture."
                << nl
                << "Available species: "
                << composition.species()
                << exit(FatalError);
        }
    }


    const label H2Index  = composition.species().find("H2");
    const label O2Index  = composition.species().find("O2");
    const label H2OIndex = composition.species().find("H2O");
    const label N2Index  = composition.species().find("N2");


    /*
     * ---------------------------------------------------------------
     * Check that N2 is the inert species.
     *
     * This version retains the original reactingFoam strategy:
     *
     *     solve H2, O2, H2O
     *     reconstruct N2 = 1 - H2 - O2 - H2O
     *
     * Therefore ScN2 is documented above but is not used in an
     * independent N2 equation in this version.
     * ---------------------------------------------------------------
     */

    if (inertIndex != N2Index)
    {
        FatalErrorInFunction
            << "This reactingFoamH2Air version expects N2 to be "
            << "the inert species."
            << nl
            << "inertSpecie = "
            << inertSpecie
            << nl
            << "N2 index = "
            << N2Index
            << nl
            << "inert index = "
            << inertIndex
            << exit(FatalError);
    }


    /*
     * ---------------------------------------------------------------
     * Solve the active species equations.
     * ---------------------------------------------------------------
     */

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            scalar ScYi;

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
            else
            {
                FatalErrorInFunction
                    << "Active species "
                    << composition.species()[i]
                    << " does not have a defined Schmidt number."
                    << nl
                    << "This solver is configured for:"
                    << nl
                    << "    H2   Sc = " << ScH2
                    << nl
                    << "    O2   Sc = " << ScO2
                    << nl
                    << "    H2O  Sc = " << ScH2O
                    << nl
                    << "    N2   Sc = " << ScN2
                    << exit(FatalError);
            }


            /*
             * ========================================================
             * Species-dependent molecular diffusion coefficient
             *
             *     Gamma_i = rho*D_i
             *              = mu/Sc_i
             *
             * NOTE:
             *
             * We deliberately use turbulence->mu(), NOT
             * turbulence->muEff().
             *
             * muEff() contains the turbulent contribution and would
             * therefore not correspond to the molecular Schmidt
             * numbers specified in the DNS reference.
             * ========================================================
             */

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
                turbulence->mu()/ScYi
            );


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

            Yi.clamp_min(0);

            Yt += Yi;
        }
    }


    /*
     * ---------------------------------------------------------------
     * Reconstruct inert N2.
     *
     * Sum(Y) = 1
     * ---------------------------------------------------------------
     */

    Y[N2Index] = scalar(1) - Yt;

    Y[N2Index].clamp_min(0);
}
#endif
