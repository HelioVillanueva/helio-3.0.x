/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

Class
    Foam::RASModels::kOmegaSSTSASnew

Group
    grpLESTurbulence

Description
    Scale-adaptive URAS model based on the k-omega-SST RAS model,
    updated 2009 version.

    References:
    \verbatim
        Egorov, Y., & Menter F.R. (2009).
        Formulation of of the Scale-Adaptive Simulation (SAS) Model during the
        DESIDER Project.
        DESider – A European Effort on Hybrid RANS-LES Modelling,
        Notes on Num. Fluid Mech. And Multidisciplinary Design,
        Volume 103, 19-67.
    \endverbatim

    The model coefficients are
    \verbatim
        kOmegaSSTSASnewCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;

            // Default SAS coefficients
            Cs          0.262;
            kappa       0.41;
            zeta2       1.47;
            sigmaPhi    2.0/3.0;
            C           2;

            // Delta must be specified for SAS e.g.
            delta cubeRootVol;

            cubeRootVolCoeffs
            {}
        }
    \endverbatim

SourceFiles
    kOmegaSSTSASnew.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTSASnew_H
#define kOmegaSSTSASnew_H

#include "kOmegaSST.H"
#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class kOmegaSSTSASnew Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTSASnew
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kOmegaSSTSASnew(const kOmegaSSTSASnew&);
        kOmegaSSTSASnew& operator=(const kOmegaSSTSASnew&);


protected:

    // Protected data

        // Model constants

            dimensionedScalar Cs_;
            dimensionedScalar kappa_;
            dimensionedScalar zeta2_;
            dimensionedScalar sigmaPhi_;
            dimensionedScalar C_;


        // Fields

            //- Run-time selectable delta model
            autoPtr<Foam::LESdelta> delta_;


    // Protected Member Functions

        //- SAS omega source
        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField& S2
        ) const;



    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

public:

    //- Runtime type information
    TypeName("kOmegaSSTSASnew");


    // Constructors

        //- Construct from components
        kOmegaSSTSASnew
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTSASnew()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kOmegaSSTSASnew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
