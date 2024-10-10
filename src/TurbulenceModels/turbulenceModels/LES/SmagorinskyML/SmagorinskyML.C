/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "SmagorinskyML.H"
#include "fvOptions.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SmagorinskyML<BasicTurbulenceModel>::k
(
    const tmp<volTensorField>& gradU
) const
{
    std::cout<<"SmagorinskyML<BasicTurbulenceModel>::k()"<<typeid(*this).name()<<std::endl;
    volSymmTensorField D(symm(gradU));

    volScalarField a(this->Ce_/this->delta());
    volScalarField b((2.0/3.0)*tr(D));
    volScalarField c(2*Ck_*this->delta()*(dev(D) && D));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("k", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            sqr((-b + sqrt(sqr(b) + 4*a*c))/(2*a))
        )
    );
}


template<class BasicTurbulenceModel>
void SmagorinskyML<BasicTurbulenceModel>::correctNut()
{
    std::cout<<"SmagorinskyML<BasicTurbulenceModel>::correctNut()"<<typeid(*this).name()<<std::endl;
    volScalarField k(this->k(fvc::grad(this->U_)));

    this->nut_ = Ck_*this->delta()*sqrt(k);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
    std::cout<<"FIN SmagorinskyML<BasicTurbulenceModel>::correctNut()"<<std::endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SmagorinskyML<BasicTurbulenceModel>::SmagorinskyML
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Ck_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ck",
            this->coeffDict_,
            0.094
        )
    ),
    TSGS_
    (
        IOobject
        (
            IOobject::groupName("TSGS", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
#ifdef USE_ML4TURB
    ml4turb_enabled_ = this->LESDict_.template getOrDefault<Switch>("ml4turb",false);
    if(ml4turb_enabled_)
    {
      ml4trub_model_ = this->LESDict_.template getOrDefault<word>("MLModelConfig","ml4turb-config.json") ;
      clip_ = this->LESDict_.template getOrDefault<double>("clip",0.) ;
      std::cout<<"ML4TURB INFO : MODEL PATH "<<ml4trub_model_<<std::endl ;
      std::cout<<"               CLIP VALUE "<<clip_<<std::endl ;
    }
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SmagorinskyML<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Ck_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void SmagorinskyML<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    std::cout<<"SmagorinskyML::correct()"<<typeid(*this).name()<<std::endl;
    LESeddyViscosity<BasicTurbulenceModel>::correct();

    std::cout<<"SmagorinskyML::correctNut()"<<std::endl;
    correctNut();

#ifdef USE_ML4TURB
    if(!ml4turb_api_initialized_ && ml4turb_enabled_)
    {
      mlturb_api_ptr_.reset(new ml4turb::ML4Turb()) ;
      mlturb_api_ptr_->initFromConfigFile(ml4trub_model_) ;
      ml4turb_api_initialized_ = true;
    }
#endif

    std::cout<<"BUILD T SGS"<<std::endl ;
    volTensorField turb_t((this->alpha_*this->rho_*this->nut())*dev2(T(fvc::grad(this->U_)))) ;
    volSymmTensorField D(symm(fvc::grad(this->U_)));
    std::size_t icount = 0;
#ifdef USE_ML4TURB
    auto nb_cells = this->mesh_.cells().size() ;
    std::vector<bool> filter(nb_cells) ;
    double mag_d_max = 0. ;
    if(ml4turb_enabled_)
    {
      mlturb_api_ptr_->startCompute(nb_cells) ;
      forAll(TSGS_, celli)
      {
          const symmTensor& d = D[celli];
          double mag_d = mag(d) ;
          mag_d_max = std::max(mag_d_max,mag_d) ;
          if(mag_d>clip_)
          {
            std::vector<float> grad_vit = {  d.xx(), d.xy(), d.xz(),
                                             d.yx(), d.yy(), d.yz(),
                                             d.xz(), d.yz(), d.zz() };
            mlturb_api_ptr_->asynchCompute(grad_vit) ;
            filter[icount] = true ;
          }
          else
            filter[icount] = false ;
          ++icount ;
      }
      mlturb_api_ptr_->endCompute() ;
    }
    std::cout<<"MAG D MAX : "<<mag_d_max<<std::endl ;
    std::vector<float> tau_ij(6) ;
    icount = 0 ;
 #endif

    forAll(TSGS_, celli)
    {
      symmTensor& tsgs = TSGS_[celli];
#ifdef USE_ML4TURB
      if(filter[icount] && ml4turb_enabled_)
      {
        mlturb_api_ptr_->getNextResult(tau_ij) ;
        tsgs = symmTensor(tau_ij[0],tau_ij[1],tau_ij[2],
                                    tau_ij[3],tau_ij[4],
                                              tau_ij[5]) ;
      }
      else
#endif
      {
        auto& t = turb_t[celli] ;
        tsgs = symmTensor(t.xx(),t.xy(),t.xz(),
                                 t.yy(),t.yz(),
                                        t.zz()) ;
      }
      ++icount ;
    }
    std::cout<<"FIN SmagorinskyML::correct()"<<typeid(*this).name()<<std::endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
