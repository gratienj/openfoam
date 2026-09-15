/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    A nd          |
    \\/     M anipulation   |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include <limits>
#include "addToRunTimeSelectionTable.H"
#include "LESModel.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "FSDUnet.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
FSDUnet<ReactionThermo, ThermoType>::FSDUnet
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion<ReactionThermo, ThermoType>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    ftRaw_
    (
        IOobject
        (
            this->thermo().phasePropertyName("ftRaw"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    ft_
    (
        IOobject
        (
            this->thermo().phasePropertyName("ft"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    YFuelFuelStream_
    (
        dimensionedScalar
        (
            "YFuelStream",
            dimless,
            1.0
        )
    ),
    YO2OxiStream_
    (
        dimensionedScalar
        (
            "YOxiStream",
            dimless,
            0.232
        )
    ),
    Cv_
    (
        this->coeffs().getScalar("Cv")
    ),
    //C_(5.0),
    C_(1.0),
    ftMin_(0.0),
    ftMax_(1.0),
    ftDim_(300),
    ftVarMin_
    (
        this->coeffs().getScalar("ftVarMin")
    ),
    Pc_
    (
        IOobject
        (
            this->thermo().phasePropertyName("Pc"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    c_
    (
        IOobject
        (
            this->thermo().phasePropertyName("Pc"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    omegaFuelBar_
    (
        IOobject
        (
            this->thermo().phasePropertyName("omegaFuelBar"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
         ),
         this->mesh(),
         dimensionedScalar
         (
            dimensionSet(1, -2, -1, 0, 0, 0, 0),
            Zero
          )
    ),
    F0_
    (
        IOobject
        (
            this->thermo().phasePropertyName("F0"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    flameSensor_
    (
        IOobject
        (
            this->thermo().phasePropertyName("flameSensor"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    FTFM_
    (
        IOobject
        (
            this->thermo().phasePropertyName("FTFM"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, 1.0)
    ),
    deltaFlame_
    (
        IOobject
        (
            this->thermo().phasePropertyName("deltaFlame"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    cTilde_
    (
        IOobject
        (
            this->thermo().phasePropertyName("cTilde"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    phiTilde_
    (
        IOobject
        (
            this->thermo().phasePropertyName("phiTilde"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    deltaRatio_
    (
        IOobject
        (
            this->thermo().phasePropertyName("deltaRatio"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    ),
    deltaL_
    (
        this->coeffs().getScalar("deltaL")
    ),
    Ncell_
    (
        this->coeffs().getScalar("Ncell")
    ),
    Cs_
    (
        this->coeffs().getScalar("Cs")
    ),
    Fmax_
    (
        this->coeffs().getScalar("Fmax")
    ),
    EColin_
    (
        this->coeffs().template getOrDefault<scalar>("EColin", 1.0)
    ),
    TFM_
    (
        this->coeffs().template getOrDefault<bool>("TFM", false)
    )
{

    Info<< "FSDUnet::CONSTRUCTOR"<<endl ;
    const cellZoneMesh& zones = this->mesh().cellZones();

    Info<< "FSDUnet :: Nb zone : "<<zones.size() <<endl ;

    auto zone_info = [&](const std::string& zoneName)
    {
      Info<< "FSDUnet Zone Name     :" << zoneName << endl ;
      label zoneID = zones.findZoneID(zoneName);
      if(zoneID==-1)
      {
        Info << "FSDUnet Zone Name     :" << zoneName << "not found "<<endl;
      }
      else
      {
        const labelList& cells = zones[zoneID];
        Info << "FSDUnet Zone Name     :" << zoneName << nl
             << "             Id       :" << zoneID <<nl
             << "             Nb cells : " << cells.size() <<endl ;
      }
    } ;

    auto zone_xyz = [&](const std::string& zoneName) {
        label zoneID = zones.findZoneID(zoneName);
        std::vector<double> xyz;
        if(zoneID==-1)
        {
          Info << "FSDUnet Zone Name     :" << zoneName << "not found "<<endl;
        }
        else
        {
            const labelList& cells = zones[zoneID];
            xyz.reserve(3*cells.size());
            double xMin = std::numeric_limits<double>::max();
            double xMax = - std::numeric_limits<double>::max();
            double yMin = std::numeric_limits<double>::max();
            double yMax = - std::numeric_limits<double>::max();
            double zMin = std::numeric_limits<double>::max();
            double zMax = - std::numeric_limits<double>::max();

            forAll(cells, cellI)
            {
                label globalCellI = cells[cellI];

                const cell& c = this->mesh().cells()[globalCellI];

                double xC = 0. ;
                double yC = 0. ;
                double zC = 0. ;
                int countC = 0 ;
                forAll(c, faceI)
                {
                    label globalFaceI = c[faceI];

                    const face& f =  this->mesh().faces()[globalFaceI];
                    double xF = 0. ;
                    double yF = 0. ;
                    double zF = 0. ;
                    int countF = 0 ;
                    forAll(f, pointI)
                    {
                        label globalPointI = f[pointI];

                        const point& p =  this->mesh().points()[globalPointI];

                        //Info<< "point = " << p << endl;
                        xF += p[0] ;
                        yF += p[1] ;
                        zF += p[2] ;
                        ++countF;
                    }
                    xF /=countF ;
                    yF /=countF ;
                    zF /=countF ;

                    xC += xF ;
                    yC += yF ;
                    zC += zF ;

                    ++countC ;
                }
                xC /= countC ;
                yC /= countC ;
                zC /= countC ;

                xyz.push_back(xC) ;
                xyz.push_back(yC) ;
                xyz.push_back(zC) ;

                xMax = std::max(xC,xMax) ;
                xMin = std::min(xC,xMin) ;
                yMax = std::max(yC,yMax) ;
                yMin = std::min(yC,yMin) ;
                zMax = std::max(zC,zMax) ;
                zMin = std::min(zC,zMin) ;
            }
            Info << "FSDUnet Zone Name     :" << zoneName << nl
                 << "             Id       :" << zoneID <<nl
                 << "             Nb cells :" << cells.size() <<nl
                 << "              xMin    :" <<xMin <<nl
                 << "              xMax    :" <<xMax <<nl
                 << "              yMin    :" <<yMin <<nl
                 << "              yMax    :" <<yMax <<nl
                 << "              zMin    :" <<zMin <<nl
                 << "              zMax    :" <<zMax <<endl;
        }
        return xyz ;
    } ;


    if(zones.size() > 0)
    {
      zone_info("inletZone") ;
      auto inlet_xyz = zone_xyz("inletZone");
      zone_info("bottomCoflowZone") ;
      auto bottom_coflow_xyz = zone_xyz("bottomCoflowZone");
      zone_info("topCoflowZone") ;
      auto top_coflow_xyz = zone_xyz("topCoflowZone");
    }


    // Example: dump stoichiometric coefficients for each species (debug)
    const List<word>& names = this->thermo().composition().species();

    forAll(names, i)
    {
        const scalar nu = this->singleMixturePtr_->specieStoichCoeffs()[i];
        Info<< "FSDUnet: specie " << names[i] << " nu = " << nu << nl;
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
FSDUnet<ReactionThermo, ThermoType>::~FSDUnet()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void FSDUnet<ReactionThermo, ThermoType>::calculateSourceNorm()
{
    /*
     * Update the fresh mixture composition.
     */
    this->singleMixturePtr_->fresCorrect();

    const label fuelI =
        this->singleMixturePtr_->fuelIndex();

    const volScalarField& YFuel =
        this->thermo().composition().Y()[fuelI];

    const volScalarField& YO2 =
        this->thermo().composition().Y("O2");

    const dimensionedScalar s =
        this->singleMixturePtr_->s();


    // --------------------------------------------------------------------- //
    // Mixture fraction
    // --------------------------------------------------------------------- //

    //const volScalarField ftRaw
    ftRaw_ =
    (
        (
            s*YFuel
          - (YO2 - YO2OxiStream_)
        )
       /
        (
            s*YFuelFuelStream_
          + YO2OxiStream_
        )
    );

    /*
     * Physical mixture fraction:
     *
     *     0 <= ft <= 1
     *
     * This prevents unphysical values during the initial transient
     * from entering the beta-PDF calculation.
     */
    ft_ =
        min
        (
            max
            (
                ftRaw_,
                scalar(0)
            ),
            scalar(1)
        );


    // --------------------------------------------------------------------- //
    // Mixture-fraction gradient and flame-surface-density estimate
    // --------------------------------------------------------------------- //

    volVectorField nft
    (
        fvc::grad(ft_)
    );

    volScalarField mgft
    (
        mag(nft)
    );

    volScalarField cAux
    (
        scalar(1) - ft_
    );

    /*
     * Regularisation of |grad(ft)|.
     *
     * The regularisation is retained from the original FSDUnet model.
     */
    dimensionedScalar dMgft =
        1.0e-3
       *
        (
            ft_*cAux*mgft
        )().weightedAverage(this->mesh().V())
       /
        (
            (ft_*cAux)().weightedAverage(this->mesh().V())
          + SMALL
        )
      + dimensionedScalar
        (
            "ddMgft",
            mgft.dimensions(),
            SMALL
        );

    mgft += dMgft;

    mgft.max(SMALL);

    nft /= mgft;



    // --------------------------------------------------------------------- //
    // Stoichiometric mixture fraction
    // --------------------------------------------------------------------- //

    const scalar ftStoich =
        YO2OxiStream_.value()
       /
        (
            s.value()*YFuelFuelStream_.value()
          + YO2OxiStream_.value()
        );


    // --------------------------------------------------------------------- //
    // Local fields
    // --------------------------------------------------------------------- //

    auto& pc = Pc_;

    // --------------------------------------------------------------------- //
    // LES filter width
    // --------------------------------------------------------------------- //

    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
            turbulenceModel::propertiesName
        );

    const volScalarField& delta =
        lesModel.delta();

    F0_ = min
            (
                max
                (
                    Ncell_*delta
                  / dimensionedScalar
                    (
                        "deltaL",
                        dimLength,
                        deltaL_
                    ),
                    scalar(1)
                ),
                Fmax_
            ) ;

    Info<< "FSDUnet DEBUG:" << nl
        << "    delta min/max     = " << gMinMax(delta) << nl
        << "    F0   min/max      = " << gMinMax(F0_) << nl
        << endl;
    // --------------------------------------------------------------------- //
    // Sub-grid mixture-fraction variance
    // --------------------------------------------------------------------- //

    const volScalarField ftVarRaw
    (
        Cv_*sqr(delta)*sqr(mgft)
    );

    /*
     * For a bounded variable 0 <= ft <= 1:
     *
     *     Var(ft) <= 1/4
     *
     * The clipping prevents invalid beta distributions.
     */
    const volScalarField ftVar
    (
        min
        (
            max
            (
                ftVarRaw,
                scalar(0)
            ),
            scalar(0.25)
        )
    );


    // --------------------------------------------------------------------- //
    // Flame thickening factor
    // --------------------------------------------------------------------- //

    const volScalarField deltaF
    (
        delta
       /
        dimensionedScalar
        (
            "flame",
            dimLength,
            1.5e-3
        )
    );

    /*
     * Linear correlation between filter size and flame thickness.
     */
    const volScalarField omegaF
    (
        max
        (
            deltaF*(4.0/3.0) + (2.0/3.0),
            scalar(1)
        )
    );

    // --------------------------------------------------------------------- //
    // Numerical integration of the mixture fraction PDF
    // --------------------------------------------------------------------- //

    const scalar deltaFt =
        1.0/ftDim_;


    forAll(ft_, celli)
    {
        const scalar ftCell =
            ft_[celli];


        if
        (
            ftCell > ftMin_
         && ftCell < ftMax_
        )
        {
            if (ftVar[celli] > ftVarMin_)
            {
                /*
                 * Safe values for the beta-PDF parameters.
                 */
                const scalar ftSafe =
                    min
                    (
                        max(ftCell, 1e-6),
                        1.0 - 1e-6
                    );

                const scalar ftVarSafe =
                    max
                    (
                        ftVar[celli],
                        1e-8
                    );

                const scalar betaTerm =
                    ftSafe*(1.0 - ftSafe)/ftVarSafe - 1.0;

                const scalar a =
                    max
                    (
                        ftSafe*betaTerm,
                        1e-6
                    );

                const scalar b =
                    max
                    (
                        (1.0 - ftSafe)*betaTerm,
                        1e-6
                    );


                // --------------------------------------------------------- //
                // Beta PDF normalization
                // --------------------------------------------------------- //

                scalar pdfIntegral = 0.0;

                for
                (
                    label i = 1;
                    i < ftDim_;
                    ++i
                )
                {
                    const scalar ft =
                        i*deltaFt;

                    const scalar pdf =
                        pow(ft, a - 1.0)
                       *pow(1.0 - ft, b - 1.0);

                    pdfIntegral +=
                        pdf*deltaFt;
                }

                pc[celli] =
                    pdfIntegral;


                // --------------------------------------------------------- //
                // Filtered consumption speed
                // --------------------------------------------------------- //

                scalar omegaIntegral = 0.0;

                const scalar sigmaFt =
                    0.01*max
                    (
                            omegaF[celli],
                        scalar(1)
                    );

                for
                (
                    label i = 1;
                    i < ftDim_;
                    ++i
                )
                {
                    const scalar ft =
                        i*deltaFt;

                    const scalar pdf =
                        pow(ft, a - 1.0)
                       *pow(1.0 - ft, b - 1.0);

                    const scalar gaussian =
                        exp
                        (
                            -sqr(ft - ftStoich)
                           /(2.0*sqr(sigmaFt))
                        );
                    /*
                    omegaIntegral +=
                        omegaFuel[celli]
                       /max(    omegaF[celli], scalar(1))
                       *gaussian
                       *pdf
                       *deltaFt;*/
                }
            }
            else
            {
                /*
                 * No sub-grid PDF required.
                 */
                const scalar sigmaFt =
                    0.01*max
                    (
                            omegaF[celli],
                        scalar(1)
                    );
            }
        }
        else
        {
            pc[celli] =
                0.0;
        }
    }

    Info<< "FSDUnet DEBUG:" << nl
        << "    Pc min/max     = " << gMinMax(pc) << nl
        << endl;

    // --------------------------------------------------------------------- //
    // Progress variable probability
    // --------------------------------------------------------------------- //

    /*
     * For the H2/O2/H2O/N2 global reaction we expect one product:
     *
     *     H2O
     *
     * The original code used List<label>(2).  Using one product is safer
     * for the present four-species mechanism.
     */
    List<label> productsIndex(1, label(-1));

    {
        label i = 0;

        forAll
        (
            this->singleMixturePtr_->specieProd(),
            specieI
        )
        {
            if
            (
                this->singleMixturePtr_->specieProd()[specieI]
                < 0
            )
            {
                if (i < productsIndex.size())
                {
                    productsIndex[i] =
                        specieI;

                    ++i;
                }
            }
        }
    }


    // --------------------------------------------------------------------- //
    // Total product mass fraction from the fresh/flamelet state
    // --------------------------------------------------------------------- //

    scalar YprodTotal =
        0.0;

    forAll(productsIndex, j)
    {
        const label specieI =
            productsIndex[j];

        if (specieI >= 0)
        {
            YprodTotal +=
                this->singleMixturePtr_->Yprod0()[specieI];
        }
    }


    // --------------------------------------------------------------------- //
    // Flame probability from mixture fraction
    // --------------------------------------------------------------------- //

    /*
     * Pc is deliberately kept as a purely geometrical/mixing quantity.
     * It is NOT overwritten by the progress variable.
     *
     * Triangular flame-probability function:
     *
     *              ft
     * Pc = -----------------------------   for ft <= ftStoich
     *           ftStoich
     *
     *              1-ft
     * Pc = -----------------------------   for ft > ftStoich
     *           1-ftStoich
     *
     * Hence:
     *
     *     Pc(0) = 0
     *     Pc(ftStoich) = 1
     *     Pc(1) = 0
     */

    const scalar ftStoichSafe =
        min
        (
            max(ftStoich, 1e-6),
            1.0 - 1e-6
        );

    forAll(ft_, celli)
    {
        const scalar ftCell = ft_[celli];

        if (ftCell <= ftStoichSafe)
        {
            pc[celli] =
                ftCell/ftStoichSafe;
        }
        else
        {
            pc[celli] =
                (1.0 - ftCell)/(1.0 - ftStoichSafe);
        }

        pc[celli] =
            min
            (
                max(pc[celli], scalar(0)),
                scalar(1)
            );
    }

    // --------------------------------------------------------------------- //
    // U-net input 1/3: progress variable cTilde (Male et al. 2024, Eq. 2.3)
    //
    //     cTilde = (ft - YFuel) / (ft - max(0, (ft-ftStoich)/(1-ftStoich)))
    //
    // YFuel is Ytilde_H2, the resolved (Favre-filtered) H2 mass fraction.
    // The denominator is (ft - YFuelBurnt), the unburnt-minus-burnt H2
    // mass fraction at the same ft (BML-type relation); it is floored to
    // avoid division by zero far from the flame (ft -> 0).
    // --------------------------------------------------------------------- //

    const volScalarField YFuelBurnt
    (
        max
        (
            dimensionedScalar(dimless, Zero),
            (ft_ - ftStoichSafe)/(1.0 - ftStoichSafe)
        )
    );

    const dimensionedScalar smallDimless("small", dimless, SMALL);

    cTilde_ =
        min
        (
            max
            (
                (ft_ - YFuel)
               /max(ft_ - YFuelBurnt, smallDimless),
                scalar(0)
            ),
            scalar(1)
        );

    // --------------------------------------------------------------------- //
    // U-net input 2/3: local equivalence ratio phiTilde (Male et al. 2024,
    // Eq. 2.5), derived from the mixture fraction:
    //
    //     phiTilde = ft*(1-ftStoich) / (ftStoich*(1-ft))
    //
    // Floored denominator to avoid blow-up as ft -> 1 (pure fuel).
    // --------------------------------------------------------------------- //

    phiTilde_ =
        (ft_*(1.0 - ftStoichSafe))
       /max(ftStoichSafe*(1.0 - ft_), smallDimless);

    Info<< "FSDUnet DEBUG:" << nl
        << "    cTilde        min/max = " << gMinMax(cTilde_) << nl
        << "    phiTilde      min/max = " << gMinMax(phiTilde_) << nl
        << "    ftStoich      = " << ftStoichSafe << nl
        << "    ft            min/max = " << gMinMax(ft_) << nl
        << "    Pc triangular min/max = " << gMinMax(pc) << nl
        << endl;



    // --------------------------------------------------------------------- //
    // Actual products in the CFD field
    // --------------------------------------------------------------------- //

    auto tproducts =
        volScalarField::New
        (
            this->thermo().phasePropertyName("products"),
            IOobject::NO_REGISTER,
            YO2.mesh(),
            dimensionedScalar(dimless, Zero)
        );

    auto& products =
        tproducts.ref();


    forAll(productsIndex, j)
    {
        const label specieI =
            productsIndex[j];

        if (specieI >= 0)
        {
            products +=
                this->thermo().composition().Y()[specieI];
        }
    }


    // --------------------------------------------------------------------- //
    // Progress variable: diagnostic only
    // --------------------------------------------------------------------- //

    // --------------------------------------------------------------------- //
    // Normalised combustion progress
    // --------------------------------------------------------------------- //
    const scalar ftStoichValue =
        YO2OxiStream_.value()
       /
        (
            s.value()*YFuelFuelStream_.value()
          + YO2OxiStream_.value()
        );

    const scalar YprodBurnt =
        scalar(9.0)
       *ftStoichValue
       *YFuelFuelStream_.value();

    const scalar YprodBurntSafe =
        max(YprodBurnt, scalar(1e-12));

    auto& c = c_;
    c =
    (
        min
        (
            max
            (
                products/YprodBurntSafe,
                scalar(0)
            ),
            scalar(1)
        )
    );
    Info<< "FSDUnet DEBUG:" << nl
        << "    YFuelFuelStream = "
        << YFuelFuelStream_.value() << nl
        << "    YO2OxiStream    = "
        << YO2OxiStream_.value() << nl
        << "    s               = "
        << s.value() << nl
        << "    YprodBurnt(calc) = "
        << YprodBurnt
        << nl
        << endl;

    Info<< "FSDUnet DEBUG:" << nl
        << "    products      min/max = " << gMinMax(products) << nl
        << "    c diagnostic  min/max = " << gMinMax(c) << nl
        << endl;


    // --------------------------------------------------------------------- //
    // Flame sensor
    // --------------------------------------------------------------------- //

    /*
     * Pc is maximum at the stoichiometric mixture fraction and decreases
     * linearly away from it.  The previous sensor
     *
     *     tanh(Cs*Pc*(1-Pc))
     *
     * was zero at Pc=0 AND Pc=1 and therefore artificially suppressed the
     * centre of the flame, where Pc=1.
     *
     * The present sensor is a Gaussian centred on ftStoich:
     *
     *     S_f = exp[-(Cs*(ft-ftStoich))^2]
     *
     * With Cs=16, the sensor is close to one around the stoichiometric
     * mixture fraction and decays rapidly away from the flame region.
     *
     * The sensor is therefore independent of the CFD product field and does
     * not suffer from the circular Pc -> c -> Pc operation of the previous
     * implementation.
     */
    // --------------------------------------------------------------------- //
    // Flame sensor
    // --------------------------------------------------------------------- //

    flameSensor_ =
        min
        (
            max
            (
                scalar(4.0)*c*(scalar(1.0) - c)*pc,
                scalar(0)
            ),
            scalar(1)
        );
    Info<< "FSDUnet DEBUG:" << nl
        << "    flameSensor  min/max = " << gMinMax(flameSensor_) << nl
        << endl;


    // --------------------------------------------------------------------- //
    // Thickening factor
    // --------------------------------------------------------------------- //
    if(TFM_)
    {
        FTFM_ =
        (
            scalar(1)
          + flameSensor_*(F0_ - scalar(1))
        );


        deltaFlame_ =
        (
            FTFM_
          * dimensionedScalar
            (
                "deltaL",
                dimLength,
                deltaL_
            )
        );
    }
    else
    {

        // No artificial thickening: the resolved/filtered flame thickness
        // is the nominal laminar flame thickness (FTFM_ = 1).
        FTFM_ = dimensionedScalar(dimless, scalar(1));

        deltaFlame_ =
            dimensionedScalar
            (
                "deltaL",
                dimLength,
                deltaL_
            );

    }

    // --------------------------------------------------------------------- //
    // U-net input 3/3: deltaRatio = deltaL0/deltaL1 (Male et al. 2024,
    // Fig. 1), i.e. the ratio of the nominal (unfiltered) laminar flame
    // thickness to the locally resolved/thickened flame thickness.
    // With the TFM convention deltaFlame_ = FTFM_*deltaL_, this reduces to
    // 1/FTFM_ (equal to 1 when there is no thickening/filtering).
    // --------------------------------------------------------------------- //

    deltaRatio_ =
        dimensionedScalar("deltaL", dimLength, deltaL_)
       /max
        (
            deltaFlame_,
            dimensionedScalar("smallLength", dimLength, SMALL)
        );

    Info<< "FSDUnet DEBUG:" << nl
        << "    YprodBurnt = " << YprodBurnt << nl
        << "    products   min/max = " << gMinMax(products) << nl
        << "    c          min/max = " << gMinMax(c) << nl
        << "    Pc         min/max = " << gMinMax(pc) << nl
        << "    flameSensor min/max = " << gMinMax(flameSensor_) << nl
        << "    FTFM       min/max = " << gMinMax(FTFM_) << nl
        << "    deltaFlame    min/max = " << gMinMax(deltaFlame_) << nl
        << "    deltaRatio    min/max = " << gMinMax(deltaRatio_) << nl
        << endl;

    // --------------------------------------------------------------------- //
    // Final FSDUnet source
    // --------------------------------------------------------------------- //
    /*
    const volScalarField wFuelFSDUnet =
        max
        (
            mgft
           *max(pc, scalar(0))
           *max
            (
                omegaFuelBar,
                dimensionedScalar
                (
                    "zeroOmegaFuelBar",
                    omegaFuelBar.dimensions(),
                    0.0
                )
            ),
            dimensionedScalar
            (
                "zeroWFuel",
                this->wFuel_.dimensions(),
                0.0
            )
        );

    if (TFM_)
    {
        this->wFuel_ =  EColin_ / FTFM_ * wFuelFSDUnet;
    }
    else
    {
        this->wFuel_ = wFuelFSDUnet;
    }*/
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::FSDUnet<ReactionThermo, ThermoType>::diffusionFactor() const
{
    return EColin_*FTFM_;
}
// * * * * * * * * * * * * * * * * Correct * * * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void FSDUnet<ReactionThermo, ThermoType>::correct()
{
    Info<< "FSDUnet::CORRECT"<<endl ;
    this->wFuel_ = Zero;

    if (this->active())
    {
        calculateSourceNorm();
    }
}


// * * * * * * * * * * * * * * * * Read * * * * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
bool FSDUnet<ReactionThermo, ThermoType>::read()
{
    Info<< "FSDUnet::READ"<<endl ;
    if
    (
        singleStepCombustion
        <
            ReactionThermo,
            ThermoType
        >::read()
    )
    {
        this->coeffs().readEntry("Cv", Cv_);

        this->coeffs().readEntry
        (
            "ftVarMin",
            ftVarMin_
        );

        this->coeffs().readEntry("TFM", TFM_);
        /*
        reactionRateFlameArea_->read
        (
            this->coeffs()
        );*/

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// ************************************************************************* //
