// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);
std::cout<<"ASSEMBLE UEQ LINEARSYSTEM"<<std::endl ;
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
std::cout<<"AFTER ASSEMBLE UEQ LINEARSYSTEM"<<std::endl ;
fvVectorMatrix& UEqn = tUEqn.ref();

std::cout<<"UEQ RELAX"<<std::endl ;
UEqn.relax();

std::cout<<"FVOPTIONS CONSTRAINT"<<std::endl ;
fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
