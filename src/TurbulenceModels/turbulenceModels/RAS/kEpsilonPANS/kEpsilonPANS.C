/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "kEpsilonPANS.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

/*
template<class BasicTurbulenceModel>
void kEpsilonPANS<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}

*/

template<class BasicTurbulenceModel>
void kEpsilonPANS<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmi_*sqr(kU_)/epsilonU_;
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonPANS<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonPANS<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonPANS<BasicTurbulenceModel>::kEpsilonPANS
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
    eddyViscosity<RASModel<BasicTurbulenceModel> >
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

    Cmi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmi",
            this->coeffDict_,
            0.22
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            -0.33
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    fEpsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fEpsilon",
            this->coeffDict_,
            1.0
        )
    ),

    cellVolume
    (
        IOobject
        (
            IOobject::groupName("cellVolume", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    fK_
    (
        IOobject
        (
            IOobject::groupName("fK", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    C2U_
    (
        IOobject
        (
            IOobject::groupName("C2U", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    kU_
    (
        IOobject
        (
            IOobject::groupName("kU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilonU_
    (
        IOobject
        (
            IOobject::groupName("epsilonU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )

{
    Info << "Defining kEpsilonPANS model" << endl;
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
        correctNut();
    }

    cellVolume.internalField() = this->mesh_.V();
    C2U_ = C1_ + (fK_/fEpsilon_)*(C2_-C1_);

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonPANS<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        //Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kEpsilonPANS<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(this->GName(), nut*(dev(twoSymm(tgradU())) && tgradU()));
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    
    // Unresolved Dissipation equation
    tmp<fvScalarMatrix> epsUEqn
    (
        fvm::ddt(alpha, rho, epsilonU_)
      + fvm::div(alphaRhoPhi, epsilonU_)
      - fvm::laplacian(alpha*rho*DepsilonUEff(), epsilonU_)
     ==
        C1_*alpha*rho*G*epsilonU_/kU_
      - fvm::SuSp(((2.0/3.0)*C1_ + C3_)*alpha*rho*divU, epsilonU_)
      - fvm::Sp(C2U_*alpha*rho*epsilonU_/kU_, epsilonU_)
      + epsilonSource()
    );

    epsUEqn().relax();
    epsUEqn().boundaryManipulate(epsilonU_.boundaryField());
    solve(epsUEqn);
    bound(epsilonU_, fEpsilon_ * this->epsilonMin_);

    // Unresolved Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkUEff(), kU_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, kU_)
      - fvm::Sp(alpha*rho*epsilonU_/kU_, kU_)
      + kSource()
    );

    kUEqn().relax();
    solve(kUEqn);
    bound(kU_, min(fK_)*this->kMin_);

    /*// Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha*rho*G*epsilon_/k_
      - fvm::SuSp(((2.0/3.0)*C1_ + C3_)*alpha*rho*divU, epsilon_)
      - fvm::Sp(C2_*alpha*rho*epsilon_/k_, epsilon_)
      + epsilonSource()
    );

    epsEqn().relax();
    epsEqn().boundaryManipulate(epsilon_.boundaryField());
    solve(epsEqn);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*epsilon_/k_, k_)
      + kSource()
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, this->kMin_);
*/

    // Calculation of Turbulent kinetic energy and Dissipation rate
    k_ = kU_/fK_;
    epsilon_ = epsilonU_/fEpsilon_;

    correctNut();

    // Re calculation of fK and C2U fo PANS fK is calculated cell by
    // cell and upper and lower limited to avoid unfisical solution
    for (int i=0; i<fK_.size(); i++)
    {
        fK_[i] = (1/sqrt(Cmi_.value()))*pow(pow(cellVolume[i],1.0/3.0)/
            (pow(k_[i],3.0/2.0)/epsilon_[i]),2.0/3.0);
        
        if (fK_[i]>1)       //upper limit
            fK_[i]=1;

        if (fK_[i]<0.1)     // lower limit
            fK_[i]=0.1;
    }

    C2U_ = C1_+ (fK_/fEpsilon_)*(C2_-C1_);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
