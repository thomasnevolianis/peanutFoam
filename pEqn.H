
volScalarField rAU(1.0/UEqn.A());

surfaceScalarField rAUf = fvc::interpolate(rAU);

volVectorField HbyA("HbyA", U);

HbyA = rAU*UEqn.H();

surfaceScalarField phiHbyA
(
	"phiHbyA",
	(fvc::interpolate(HbyA) & mesh.Sf())
	+ fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
);

adjustPhi(phiHbyA,U,p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU);


// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
        // Pressure corrector

	fvScalarMatrix pEqn
	(
		fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );


        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

        if (piso.finalNonOrthogonalIter())
        {
         	phi = phiHbyA - pEqn.flux();
        }
            

	U = HbyA + rAU*fvc::reconstruct((- pEqn.flux())/rAUf);
        U.correctBoundaryConditions();
}
