#include "TurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"
#include "eddyViscosity.H"
#include "linearViscousStress.H"
#include "RASModel.H"
#include "autoPtr.H"

/*namespace Foam
{
    typedef BasicTurbulenceModel<transportModel> 
        transportModelBasicTurbulenceModel;
    typedef RASModel<transportModelBasicTurbulenceModel>
        RAStransportModelBasicTurbulenceModel;
}*/

#include "kEpsilonPANS.H"

makeTemplatedTurbulenceModel
(
    turbulenceModel,
    RAS,
    kEpsilonPANS
);