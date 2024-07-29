// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  - fvm::laplacian(mixMu+turbulence->nut() * rho,U)
  - fvc::div(((mixMu+turbulence->nut()*rho)*dev2(fvc::grad(U)().T())))
 ==  
    fvModels.source(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
