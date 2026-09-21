/*---------------------------------------------------------------------------*\
  initialiseTurbulentU

  OpenFOAM ESI v2606

  Initialise U with homogeneous isotropic synthetic turbulence.

  Mean velocity:
      Umean = (3.6 0 0) m/s

  Component RMS:
      u' = v' = w' = 0.36 m/s

  Integral/correlation length target:
      L = 0.002 m

  The field is generated from random Fourier modes with random
  wave-vector directions and random phases.

  The resulting fluctuations are:
      <u'> = <v'> = <w'> = 0
      RMS(u') = RMS(v') = RMS(w') = 0.36 m/s

  The fluctuations are approximately isotropic and have a
  characteristic spatial scale controlled by L.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"
#include "Pstream.H"
#include "mathematicalConstants.H"


using namespace Foam;


// ************************************************************************* //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Initialise U with homogeneous isotropic synthetic turbulence"
    );

    argList::addOption
    (
        "Umean",
        "vector",
        "Mean velocity, default = (3.6 0 0)"
    );

    argList::addOption
    (
        "Urms",
        "scalar",
        "RMS of each velocity fluctuation component, default = 0.36"
    );

    argList::addOption
    (
        "L",
        "scalar",
        "Target turbulence length scale, default = 0.002"
    );

    argList::addOption
    (
        "Nmodes",
        "label",
        "Number of Fourier modes, default = 64"
    );

    argList::addOption
    (
        "seed",
        "label",
        "Random seed, default = 123456"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    // ---------------------------------------------------------------------
    // User parameters
    // ---------------------------------------------------------------------

    vector Umean(3.6, 0.0, 0.0);
    scalar Urms = 0.36;
    scalar L = 0.002;
    label Nmodes = 64;
    label seed = 123456;


    if (args.optionFound("Umean"))
    {
        IStringStream is(args["Umean"]);
        is >> Umean;
    }

    if (args.optionFound("Urms"))
    {
        Urms = args.getOrDefault<scalar>("Urms", 0.36);
    }

    if (args.optionFound("L"))
    {
        L = args.getOrDefault<scalar>("L", 0.002);
    }

    if (args.optionFound("Nmodes"))
    {
        Nmodes = args.getOrDefault<label>("Nmodes", 64);
    }

    if (args.optionFound("seed"))
    {
        seed = args.getOrDefault<label>("seed", 123456);
    }


    if (L <= SMALL)
    {
        FatalErrorInFunction
            << "L must be greater than zero"
            << exit(FatalError);
    }

    if (Urms < 0)
    {
        FatalErrorInFunction
            << "Urms must be >= 0"
            << exit(FatalError);
    }

    if (Nmodes < 1)
    {
        FatalErrorInFunction
            << "Nmodes must be >= 1"
            << exit(FatalError);
    }


    Info<< nl
        << "============================================================"
        << nl
        << " initialiseTurbulentU"
        << nl
        << "============================================================"
        << nl
        << "Mean velocity Umean       = " << Umean << " m/s" << nl
        << "Component RMS             = " << Urms << " m/s" << nl
        << "Turbulence length L       = " << L << " m" << nl
        << "Number of Fourier modes   = " << Nmodes << nl
        << "Random seed               = " << seed << nl
        << "============================================================"
        << nl << endl;


    // ---------------------------------------------------------------------
    // Read U
    // ---------------------------------------------------------------------

    Info<< "Reading U..." << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // ---------------------------------------------------------------------
    // Random generator
    // ---------------------------------------------------------------------

    Random rnd(seed);


    // ---------------------------------------------------------------------
    // Fluctuation field
    //
    // Temporary internal field:
    //
    //     u'(x) = sum_m A_m cos(k_m.x + phi_m) e_m
    //
    // where e_m is a random unit vector.
    //
    // We use several modes with slightly different wave numbers
    // around k0 = 2*pi/L.
    // ---------------------------------------------------------------------

    volVectorField fluctuation
    (
        IOobject
        (
            "Ufluctuation",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "zero",
            dimVelocity,
            vector::zero
        )
    );


    const scalar pi = constant::mathematical::pi;
    const scalar k0 = 2.0*pi/L;


    Info<< "Generating synthetic turbulent fluctuations..." << endl;


    for (label m = 0; m < Nmodes; ++m)
    {
        // -------------------------------------------------------------
        // Random wave-vector direction
        // -------------------------------------------------------------

        vector kDir
        (
            rnd.GaussNormal<scalar>(),
            rnd.GaussNormal<scalar>(),
            rnd.GaussNormal<scalar>()
        );

        scalar kmagDir = mag(kDir);

        if (kmagDir < SMALL)
        {
            kDir = vector(1, 0, 0);
            kmagDir = 1.0;
        }

        kDir /= kmagDir;


        // -------------------------------------------------------------
        // Random wave-number magnitude
        //
        // Spread around k0 to avoid a single wavelength.
        // -------------------------------------------------------------

        scalar kFactor = 0.65 + 0.70*rnd.sample01<scalar>();

        vector k = kDir*(k0*kFactor);


        // -------------------------------------------------------------
        // Random phase
        // -------------------------------------------------------------

        scalar phase = 2.0*pi*rnd.sample01<scalar>();


        // -------------------------------------------------------------
        // Random polarization vector
        //
        // Build two vectors perpendicular to k.
        // The resulting polarization is transverse to k, giving
        // approximately solenoidal fluctuations.
        // -------------------------------------------------------------

        vector ref;

        if (mag(kDir.x()) < 0.9)
        {
            ref = vector(1, 0, 0);
        }
        else
        {
            ref = vector(0, 1, 0);
        }

        vector e1 = kDir ^ ref;

        scalar e1mag = mag(e1);

        if (e1mag < SMALL)
        {
            ref = vector(0, 0, 1);
            e1 = kDir ^ ref;
            e1mag = mag(e1);
        }

        e1 /= e1mag;

        vector e2 = kDir ^ e1;
        e2 /= mag(e2);


        // Random polarization angle

        scalar alpha = 2.0*pi*rnd.sample01<scalar>();

        vector polarization =
            Foam::cos(alpha)*e1
          + Foam::sin(alpha)*e2;


        // -------------------------------------------------------------
        // Mode amplitude
        //
        // Equal energy contribution per mode.
        // The final field is normalized exactly afterwards.
        // -------------------------------------------------------------

        const scalar amplitude = 1.0/Foam::sqrt(scalar(Nmodes));


        // -------------------------------------------------------------
        // Add mode to field
        // -------------------------------------------------------------

        vectorField& f = fluctuation.primitiveFieldRef();

        const vectorField& C = mesh.C();

        forAll(f, celli)
        {
            const scalar argument =
                (k & C[celli]) + phase;

            f[celli] +=
                amplitude
               *polarization
               *Foam::cos(argument);
        }
    }


    // ---------------------------------------------------------------------
    // Remove the global mean of each fluctuation component.
    //
    // Use cell-volume weighting.
    // ---------------------------------------------------------------------

    scalar localVolume = 0.0;
    vector localIntegral(vector::zero);

    const scalarField& V = mesh.V();

    const vectorField& f = fluctuation.internalField();

    forAll(f, celli)
    {
        localVolume += V[celli];
        localIntegral += V[celli]*f[celli];
    }


    Pstream::combineReduce
    (
        localVolume,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localIntegral,
        plusEqOp<vector>()
    );


    const vector globalMean =
        localIntegral/localVolume;


    Info<< "Initial fluctuation mean = "
        << globalMean << endl;


    // Remove mean

    vectorField& fRef = fluctuation.primitiveFieldRef();

    forAll(fRef, celli)
    {
        fRef[celli] -= globalMean;
    }


    // ---------------------------------------------------------------------
    // Calculate volume-weighted RMS of each component.
    // ---------------------------------------------------------------------

    scalar localU2 = 0.0;
    scalar localV2 = 0.0;
    scalar localW2 = 0.0;


    forAll(fRef, celli)
    {
        localU2 += V[celli]*sqr(fRef[celli].x());
        localV2 += V[celli]*sqr(fRef[celli].y());
        localW2 += V[celli]*sqr(fRef[celli].z());
    }


    Pstream::combineReduce
    (
        localU2,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localV2,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localW2,
        plusEqOp<scalar>()
    );


    const scalar rmsU =
        Foam::sqrt(localU2/localVolume);

    const scalar rmsV =
        Foam::sqrt(localV2/localVolume);

    const scalar rmsW =
        Foam::sqrt(localW2/localVolume);


    Info<< nl
        << "Before normalization:" << nl
        << "  RMS(u') = " << rmsU << " m/s" << nl
        << "  RMS(v') = " << rmsV << " m/s" << nl
        << "  RMS(w') = " << rmsW << " m/s" << nl
        << endl;


    // ---------------------------------------------------------------------
    // Normalize each component independently.
    //
    // This guarantees:
    //
    // RMS(u') = RMS(v') = RMS(w') = Urms
    // ---------------------------------------------------------------------

    const scalar scaleU =
        Urms/max(rmsU, SMALL);

    const scalar scaleV =
        Urms/max(rmsV, SMALL);

    const scalar scaleW =
        Urms/max(rmsW, SMALL);


    forAll(fRef, celli)
    {
        fRef[celli].x() *= scaleU;
        fRef[celli].y() *= scaleV;
        fRef[celli].z() *= scaleW;
    }


    // ---------------------------------------------------------------------
    // Add mean velocity
    // ---------------------------------------------------------------------

    vectorField& UInternal = U.primitiveFieldRef();

    forAll(UInternal, celli)
    {
        UInternal[celli] =
            Umean + fRef[celli];
    }


    // ---------------------------------------------------------------------
    // Recalculate statistics
    // ---------------------------------------------------------------------

    scalar localUmean = 0.0;
    scalar localVmean = 0.0;
    scalar localWmean = 0.0;

    localU2 = 0.0;
    localV2 = 0.0;
    localW2 = 0.0;


    forAll(UInternal, celli)
    {
        localUmean += V[celli]*(UInternal[celli].x() - Umean.x());
        localVmean += V[celli]*(UInternal[celli].y() - Umean.y());
        localWmean += V[celli]*(UInternal[celli].z() - Umean.z());

        localU2 +=
            V[celli]
           *sqr(UInternal[celli].x() - Umean.x());

        localV2 +=
            V[celli]
           *sqr(UInternal[celli].y() - Umean.y());

        localW2 +=
            V[celli]
           *sqr(UInternal[celli].z() - Umean.z());
    }


    Pstream::combineReduce
    (
        localUmean,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localVmean,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localWmean,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localU2,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localV2,
        plusEqOp<scalar>()
    );

    Pstream::combineReduce
    (
        localW2,
        plusEqOp<scalar>()
    );


    const vector fluctMean
    (
        localUmean/localVolume,
        localVmean/localVolume,
        localWmean/localVolume
    );


    const vector fluctRMS
    (
        Foam::sqrt(localU2/localVolume),
        Foam::sqrt(localV2/localVolume),
        Foam::sqrt(localW2/localVolume)
    );


    Info<< nl
        << "Final field statistics:" << nl
        << "  mean fluctuation = " << fluctMean << " m/s" << nl
        << "  RMS fluctuation  = " << fluctRMS << " m/s" << nl
        << "  target RMS       = " << Urms << " m/s" << nl
        << endl;


    // ---------------------------------------------------------------------
    // Correct boundary conditions.
    //
    // This preserves turbulentDFSEMInlet on coflowInletBottom.
    // ---------------------------------------------------------------------

    surfaceScalarField phi
    (
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
    );
    U.correctBoundaryConditions();


    // ---------------------------------------------------------------------
    // Write U
    // ---------------------------------------------------------------------

    Info<< "Writing initialized U..." << endl;

    U.write();


    Info<< nl
        << "============================================================"
        << nl
        << "Done."
        << nl
        << "U has been initialized with:"
        << nl
        << "  Umean = " << Umean << " m/s"
        << nl
        << "  RMS   = " << Urms << " m/s per component"
        << nl
        << "  L     = " << L << " m"
        << nl
        << "============================================================"
        << nl
        << endl;


    return 0;
}


// ************************************************************************* //
