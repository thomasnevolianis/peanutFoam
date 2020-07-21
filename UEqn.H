
// Solve the Momentum equation

fvVectorMatrix UEqn
(
    rho*fvm::ddt(U) 
  + rho*fvm::div(phi, U)
  - fvm::laplacian(mu, U)
);


//************* Surface tension scheme ***************//

int i = 0;
scalar area0=0;
		
volScalarField sqrGradC = magSqr(fvc::grad(C));
		
// Calculate Reference Values for Area

if (i == 0){
		
// Initial Surface Area
scalar intArea0 = 0;

			
forAll(C, celli){
	intArea0 +=((epsilon.value()/2)*sqrGradC[celli] + (1/(4*(epsilon.value())))*pow((pow(C[celli],2)-1),2)  )*C.mesh().V()[celli];
		}						

area0 = ( 3*pow(2,0.5)/4 )* intArea0;		
	   }
	

// Current Surface Area

scalar intArea = 0;
		

forAll(C,celli)
{
          intArea += ( (epsilon.value()/2)*sqrGradC[celli] + (1/(4*(epsilon.value())))*pow( (pow(C[celli],2)-1) ,2)  )*C.mesh().V()[celli];
}

        scalar area = ( 3*pow(2,0.5)/4 )*(intArea);




dimensionedScalar sigma = k_P*(area-area0);

//**************************************************************************//

solve(UEqn == 
	      - fvc::grad(p) 
	      - epsilon*fvc::div(sigma*fvc::grad(C)*fvc::grad(C)) 
	      + varEnergy*fvc::grad(C)
      );

