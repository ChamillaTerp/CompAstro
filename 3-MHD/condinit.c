#include "fargo3d.h"

void Init() {

  int i,j,k;

  //These are global structures
  real *rho = Density->field_cpu; //Centered
  real *cs  = Energy->field_cpu;  //Centered
  real *vx  = Vx->field_cpu;      //Staggered in x
  real *vy  = Vy->field_cpu;      //Staggered in y
  real *vz  = Vz->field_cpu;      //Staggered in z
  real *bx  = Bx->field_cpu;      //Staggered in x
  real *by  = By->field_cpu;      //Staggered in y
  real *bz  = Bz->field_cpu;      //Staggered in z
  
  real b0     = SOUNDSPEED*sqrt(2.0*MU0/BETA);
  real va     = b0/sqrt(MU0);
  real kmode  = 2*M_PI/(ZMAX-ZMIN);
  real ktilde = kmode*va/OMEGAFRAME;
  real sigma  = sqrt(sqrt(16*pow(ktilde,2)+1) - 2*pow(ktilde,2)-1)/sqrt(2);

  //We now compute the amplitudes:
  real dvy  =  sigma;
  real dvx  =  0.5*(pow(sigma,2) + pow(ktilde,2));
  real dby  =  ktilde;
  real dbx  = -2*ktilde*sigma/(pow(sigma,2) + pow(ktilde,2));
  
  real norm =  sqrt(pow(dvx,2) + pow(dvy,2) + pow(dbx,2) + pow(dby,2));

  dvx *= AMPLITUDE/norm;
  dvy *= AMPLITUDE/norm;
  dbx *= AMPLITUDE/norm;
  dby *= AMPLITUDE/norm;

  printf("\n==============================================\n");
  printf("You are computing the mode ktilde = %g\n", ktilde);
  printf("==============================================\n\n");
  
  if (ktilde > sqrt(3)) {
    printf("Error: k mode bigger than the cut-off mode.\n");
    printf("Hint: Use a smaller box or smaller a BETA.\n");
    exit(1);
  }
  
  //Loop over the mesh
  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      for (i=0; i<Nx+2*NGHX; i++) {

	//Note: index l is the cell index
	
	rho[l] =  1.0;
	cs[l]  =  SOUNDSPEED;
	
	vx[l]  = -1.5*OMEGAFRAME*ymed(j);
	vy[l]  =  0.0;
	vz[l]  =  0.0;
	bx[l]  =  0.0;
	by[l]  =  0.0;
	bz[l]  =  b0;
	
	//Now we add the perturbations using the eienvectors...
	vx[l] += dvx*sin(kmode*zmed(k));
	vy[l]  = dvy*sin(kmode*zmed(k));
	bx[l] += dbx*cos(kmode*zmed(k));
	by[l]  = dby*cos(kmode*zmed(k));

      }
    }
  }
}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   Init();
}
