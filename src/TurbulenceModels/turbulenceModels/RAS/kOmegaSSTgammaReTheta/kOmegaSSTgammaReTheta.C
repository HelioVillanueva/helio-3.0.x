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

#include "kOmegaSSTgammaReTheta.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Return the DU/Ds for each cell

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::findGrad(const volVectorField& Vel) const
{
    volVectorField velgrad = fvc::grad(mag(Vel));
    
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "dUds",
            Vel.component(0)/(max(mag(Vel),minvel))*velgrad.component(0)
            +Vel.component(1)/max(mag(Vel),minvel)*velgrad.component(1)
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::kOmegaSSTgammaReTheta::F1(const volScalarField& CDkOmega) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    tmp<volScalarField> R_y = y_*sqrt(k_)*this->rho_/this->mu();

    tmp<volScalarField> F3_y = exp(-pow(R_y/120.0,8.0));

    return max(tanh(pow4(arg1)),F3_y);
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::kOmegaSSTgammaReTheta::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::kOmegaSSTgammaReTheta::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::kOmegaSSTgammaReTheta::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23() *= F3();
    }

    return f23;
}


template<class BasicTurbulenceModel>
void kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::correctNut(const volScalarField& S2)
{
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmegaSSTgammaReTheta<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTgammaReTheta<BasicTurbulenceModel>
::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kOmegaSSTgammaReTheta<BasicTurbulenceModel>::Qsas
(
    const volScalarField& S2,
    const volScalarField& gamma,
    const volScalarField& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTgammaReTheta<BasicTurbulenceModel>::kOmegaSSTgammaReTheta
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

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    freeStreamU
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "freeStreamU",
            this->coeffDict_,
            dimensionSet(0,1,-1,0,0,0,0),
            5.4
        )
    ),
    minvel
    (
        dimensionedScalar
        (
            "minvel",
            dimensionSet( 0, 1, -1, 0, 0, 0, 0),
            1.0e-06
        )
    ),

    mindis
    (
        dimensionedScalar
        (
            "mindis",
            dimensionSet( 0, 1, 0, 0, 0, 0, 0),
            1.0e-06
        )
    ),
    mintim
    (
        dimensionedScalar
        (
            "mintim",
            dimensionSet( 0, 0, -1, 0, 0, 0, 0),
            1.0e-06
        )
    ),
    F3_
    (
        Switch::lookupOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

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
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    intermittency_
    (
        IOobject
        (
            IOobject::groupName("intermittency", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    Ret_
    (
        IOobject
        (
            IOobject::groupName("Ret", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    Reo_
    (
        Ret_
    )

{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTgammaReTheta<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel> >::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());
        freeStreamU.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaSSTgammaReTheta<BasicTurbulenceModel>::correct()
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
    volScalarField nu = this->nu();

    eddyViscosity<RASModel<BasicTurbulenceModel> >::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField str(sqrt(S2));
    volScalarField vort(sqrt(2*magSqr(skew(tgradU()))));
    volScalarField GbyNu((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField G(this->GName(), nut*GbyNu);
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));

    {
        volScalarField gamma(this->gamma(F1));
        volScalarField beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            alpha*rho*gamma
           *min
            (
                GbyNu,
                (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, b1_*F23()*str)
            )
          - fvm::SuSp((2.0/3.0)*alpha*rho*gamma*divU, omega_)
          - fvm::Sp(alpha*rho*beta*omega_, omega_)
          - fvm::SuSp
            (
                alpha*rho*(F1 - scalar(1))*CDkOmega/omega_,
                omega_
            )
          + Qsas(S2, gamma, beta)
          + omegaSource()
        );

        omegaEqn().relax();

        omegaEqn().boundaryManipulate(omega_.boundaryField());

        solve(omegaEqn);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        min(alpha*rho*G, (c1_*betaStar_)*alpha*rho*k_*omega_)*intermittency_
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(min(max(intermittency_,0.1),1.0)*alpha*rho*betaStar_*omega_, k_)
      + kSource()
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

    correctNut(S2);


//   transition model

    volScalarField Rew = omega_*y_*y_/nu;
    volScalarField Rev = y_*y_*str/nu;

    volScalarField F_wake = exp(-pow(Rew/1e5,2.0));
    
    volScalarField delta = 50.0*vort*y_/max(mag(U),minvel)*15.0/2.0*nu*Ret_
                            /max(mag(U),minvel);
    volScalarField Rterm1 = F_wake*exp(-pow(y_/max(delta,mindis),4.0));
    volScalarField Rterm2 = 1.0-pow((intermittency_-1.0/50.0)/(1.0-1.0/50.0),2.0);
    
    volScalarField Ftheta = min(max(Rterm1,Rterm2),1.0);

    volScalarField Ptheta = 0.03*pow(mag(U),2.0)/(500.0*nu*(1.0-Ftheta,1.0));

    volScalarField F_turb = exp(-pow(k_/(omega_*nu*4.0),4.0));

    volScalarField F_reattach = exp(-pow(k_/(omega_*nu*20.0),4.0));

    //velocity gradient along a streamline
    volScalarField duds=findGrad(U);

    scalar Umag; // absolte value of the velocity
    scalar Tu;   // Turbulent intensity 
     
    rootFunction tf(1,2,3,4);  // creating the function object
    
    forAll(k_, cellI)
    {
        Umag = max(mag(U[cellI]),1.0e-8); // avoiding division by zero

        // bounding dUds for robustness
        duds[cellI] = max(min(duds[cellI],0.5),-0.5); 
        Tu = 100.0*sqrt(2.0/3.0*k_[cellI])/Umag;
        // modifying the object
        tf.modify(Tu,duds[cellI],nu[cellI],freeStreamU.value());
        // solving the nonlinear equation
        Reo_[cellI] = NewtonRoot<rootFunction>(tf, 1e-5).root(0.0);
        Reo_[cellI] = Reo_[cellI]*freeStreamU.value()/nu[cellI];
    }

    Reo_ = max(Reo_,20.0); //limiting Reo_ for robustness

    // Re_theta  equation
    tmp<fvScalarMatrix> ReEqn
    (
        fvm::ddt(alpha, rho, Ret_)
      + fvm::div(alphaRhoPhi, Ret_)
      - fvm::Sp(fvc::div(this->phi_), Ret_)
      - fvm::laplacian(DRetEff(), Ret_)
     ==
        Ptheta*Reo_
      - fvm::Sp(Ptheta, Ret_)
    );

    ReEqn().relax();
    solve(ReEqn);

    volScalarField Re_crit(Ret_);
    volScalarField F_length(Ret_);

    forAll(Ret_, cellI)
    {
     if (Ret_[cellI]<=1870.0) 
        { 
            Re_crit[cellI] = Ret_[cellI]-396.035*1.0e-02+120.656*1.0e-04
                            *Ret_[cellI]-868.230*1.0e-06*pow(Ret_[cellI],2.0)
                            +696.506*1.0e-09*pow(Ret_[cellI],3.0)-174.105
                            *1.0e-12*pow(Ret_[cellI],4.0);
        }  
     else
        { 
            Re_crit[cellI]=Ret_[cellI]-593.11-0.482*(Ret_[cellI]-1870.0);
        }

     if (Ret_[cellI]<400.0)
        {
            F_length[cellI] = 398.189*1.0e-01-119.270*1.0e-04*Ret_[cellI]
                            -132.567*1.0e-06*pow(Ret_[cellI],2.0);
        } 
     else if (Ret_[cellI]>=400.0 && Ret_[cellI]<596.0)
        {
            F_length[cellI] = 263.404-123.939*1.0e-02*Ret_[cellI]+194.548
                            *1.0e-05*pow(Ret_[cellI],2.0)-101.695*1.0e-08
                            *pow(Ret_[cellI],3.0);
        } 
     else if (Ret_[cellI]>=596.0 && Ret_[cellI]<1200.0)
        {
            F_length[cellI]=0.5-(Ret_[cellI]-596.0)*3.0*1.0e-04;
        } 
     else
        {
            F_length[cellI]=0.3188;
        }

    }

    //- Return the separation induced transition
    volScalarField separation_im = min(2.0*max(0.0,Rev/(3.235*Re_crit)-1.0)
                                    *F_reattach,2.0)*Ftheta;

    // For F_onset
    volScalarField F1onset = Rev/(2.193*max(Re_crit,1.0e-6));
    volScalarField F2onset = min(max(F1onset,pow(F1onset,4.0)),2.0);
    volScalarField F3onset = max(1.0-pow(k_/(omega_*nu*2.5),3.0),0.0);

    volScalarField F_onset = max(F2onset-F3onset,0.0);

    volScalarField F_sub = (exp(-pow(sqr(y_)*omega_/(500.0*nu*0.4),2.0)));
    F_length = F_length*(1.0-F_sub)+40.0*F_sub; 
    
    // ja tem const volScalarField vor(sqrt(2*magSqr(skew(fvc::grad(U)))));
    volScalarField Pr1(2.0*F_length*sqrt(intermittency_*F_onset)*str);
    volScalarField Pr2(0.06*vort*F_turb*intermittency_);
    
    // intermittency equation
    tmp<fvScalarMatrix> intermEqn
    (
        fvm::ddt(alpha, rho, intermittency_)
      + fvm::div(alphaRhoPhi, intermittency_)
      - fvm::Sp(fvc::div(this->phi_), intermittency_)
      - fvm::laplacian(DimEff(), intermittency_)
     == 
        
        Pr1
       +Pr2
      - fvm::Sp(Pr1,intermittency_)
      - fvm::Sp(Pr2*50.0,intermittency_)

    );

    intermEqn().relax();
    solve(intermEqn);
    
    //correction for separation induced transition 
    intermittency_ = max(intermittency_,separation_im);
    //bounding intermittency
    intermittency_ = min(max(intermittency_,0.00001),1.0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
