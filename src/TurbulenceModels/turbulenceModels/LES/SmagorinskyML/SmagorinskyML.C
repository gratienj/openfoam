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

    std::cout<<"BUILD T SGS"<<std::endl ;
    volTensorField turb_t((this->alpha_*this->rho_*this->nut())*dev2(T(fvc::grad(this->U_)))) ;
    volSymmTensorField D(symm(fvc::grad(this->U_)));
    forAll(TSGS_, celli)
    {
        const symmTensor& d = D[celli];
        symmTensor& tsgs    = TSGS_[celli];
        auto& t             = turb_t[celli] ;
        //std::cout<<"||D["<<celli<<"]|| = "<<mag(d)<<std::endl ;
        //std::cout<<"turb_t["<<celli<<"]="<<typeid(t).name()<<" "<<t.xx()<<" "<<t.yy()<<" "<<t.zz()<<std::endl ;
        tsgs = symmTensor(t.xx(),t.xy(),t.xz(),
                                 t.yy(),t.yz(),
                                        t.zz()) ;
        //std::cout<<"Tensor D["<<celli<<"] = "<<d.xx()<<" "<<d.yy()<<" "<<d.zz()<<std::endl ;
        //std::cout<<"Tensor TSGS["<<celli<<"] = "<<tsgs.xx()<<" "<<tsgs.yy()<<" "<<tsgs.zz()<<std::endl ;
    }
    std::cout<<"FIN SmagorinskyML::correct()"<<typeid(*this).name()<<std::endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
