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
#include "HashSet.H"
#include "IOdictionary.H"

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
     * ================================================================
     */

    //- Default Schmidt number for species without a prescribed value
    const scalar ScDefault = 1.0;


    /*
     * ================================================================
     * Read the species Schmidt numbers from the "SchmidtNumber"
     * sub-dictionary of thermophysicalProperties (optional).
     *
     * The "importantSpecies" list is read from chemistryProperties:
     * it identifies the species for which the differential-diffusion
     * treatment is applied. The inert species (inertSpecie below) is
     * also treated individually. Any other active species uses the
     * default Schmidt number.
     * ================================================================
     */

    const dictionary chemistryDict
    (
        IOdictionary
        (
            IOobject
            (
                "chemistryProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );

    const wordHashSet importantSpecies
    (
        chemistryDict.found("importantSpecies")
      ? wordHashSet(chemistryDict.subDict("importantSpecies").toc())
      : wordHashSet()
    );

    //- Schmidt numbers for the important/inert species, if prescribed
    const dictionary* SchmidtNumberDictPtr = thermo.findDict("SchmidtNumber");


    /*
     * ================================================================
     * Solve ALL species equations.
     *
     * No species is reconstructed from unity (sum of mass fractions);
     * each active species is transported with its own Schmidt number.
     * ================================================================
     */

    forAll(Y, i)
    {
        if (composition.active(i))
        {
            const word& specieName = composition.species()[i];

            volScalarField& Yi = Y[i];


            /*
             * Species-dependent Schmidt number
             */

            const scalar ScYi
            (
                (
                    importantSpecies.found(specieName)
                 || specieName == inertSpecie
                )
             && SchmidtNumberDictPtr
              ? SchmidtNumberDictPtr->
                    getOrDefault<scalar>(specieName, ScDefault)
              : ScDefault
            );


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


            Info<< "Species " << Yi.name()
                << " : Sc = " << ScYi
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
     * Normalization of the mass fractions
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
     * Normalize all species.
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
