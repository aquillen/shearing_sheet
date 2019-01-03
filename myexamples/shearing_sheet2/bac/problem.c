
/**
 * Shearing sheet
 *
 * This example simulates in
 * shearing sheet coordinates. If you have OpenGL enabled, 
 * you'll see one copy of the computational domain. Press `g` to see
 * the ghost boxes which are used to calculate gravity.
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void mkfilename();
void pos_output_ascii();
void drag_forces(struct reb_simulation* r);
void heartbeat(struct reb_simulation* const r);
void dowrap();
void vel_disp();
void add_massless();
void reset_massless();
double peak_sigma();
double QToom();
void compute_epi();

double wait_period;
double move_period;
double alpha1;
double alpha1_om;
double alpham; // massive partile damping rate units 1/t
double alpham_om; // to set alpham, units omega^-1
int  ntodamp;  // number of massive to damp
int nparticles_nomass;  // number of massless particles
double surfacedensity;
int *wrapx;  // for massless particles only
int *wrapy;
double *oldx;
double *oldy;
double rmassless;

int main(int argc, char* argv[]) {
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->opening_angle2	= .5;	// This determines the precission of the tree code gravity calculation.
	r->integrator			= REB_INTEGRATOR_SEI;
	r->boundary			= REB_BOUNDARY_SHEAR;
	r->gravity			= REB_GRAVITY_TREE;
	r->collision			= REB_COLLISION_NONE;
        r->additional_forces            = drag_forces;    
	r->heartbeat			= heartbeat;	// function pointer for heartbeat

        alpham                          = 0.0;  // will be damping rate for massive particles
        alpha1                          = 0.0;  // will be damping rate for massless  particles

/////////////////////stuff I want to set/input
        // Units: pc, km/s, Myr,  1km/s ~ pc/Myr, Msol
	r->G 				= 0.0045;
	double OMEGA 			= 0.030; // units 1/Myr
	double KAPPA_OMEGA 		= 1.4;	// unitless, is kappa/omega
        double OMEGAZ			= 1.8*OMEGA;
             // kappa/Omega = 1.0 for Keplerian
             //             = \sqrt{2} for flat rotation curve
	r->softening 			=  50.0;	// pc,   smoothing length
	double boxsize 			= 2000.0;	// pc, is half the box length
        nparticles_nomass           	= 1000;  // number of massless particles
        int nparticles                  = 50000;  // number of massive particles
        ntodamp 			= nparticles_nomass + nparticles;  
			// if less than sum then not all particles are damped

        double vel_std                  =  6.0;   // initial velocity dispersion 
                                                  // of massive particles
                                                  // in km/s
        double vel_stdz                 = 10.0;  // z velocity dispersion  massive particles
	surfacedensity 	        	= 10.0;	// surface density of massive particles 
                                                // in Msol/pc^2 
        rmassless			= 3.0;  // initial spatial disp of massless
						// in pc -- only if using a cluster
        alpham_om                       = 5e-3; // for setting alpham, damping rate for massive particles
        alpha1_om			= 0.0;  // massless
        // output file name for outputs?
        wait_period 			= 2.0; // number of periods for reset
        move_period 			= 5.0; // number of periods to run after reset
        double vel_stdz_massless        = 5.0; // z velocity dispersion of massless particles

// stuff set automatically
        double period = 2.0*M_PI/OMEGA;	// Myr
	r->dt 				= 1e-2*period;	// Myr
	r->ri_sei.OMEGA 		= OMEGA;
	r->ri_sei.KAPPA_OMEGA 		= KAPPA_OMEGA;
	r->ri_sei.OMEGAZ		= OMEGAZ;
	// This example uses two root boxes in the x and y direction. 
	// Although not necessary in this case, it allows for the parallelization using MPI. 
	// See Rein & Liu for a description of what a root box is in this context.
	reb_configure_box(r, boxsize, 2, 2, 1);
	r->nghostx = 2;
	r->nghosty = 2;
	r->nghostz = 0;
        double fourmkap3= (4.0 - KAPPA_OMEGA*KAPPA_OMEGA)/3.0;


/// stuff to print out
	double	particle_mass = surfacedensity*boxsize*boxsize*4/nparticles;
        printf("particle mass = %f\n",particle_mass);
        double kappa = KAPPA_OMEGA*OMEGA;
	printf("Toomre wavelength: %f\n",
           4.*M_PI*M_PI*surfacedensity*r->G/kappa/kappa);
	printf("Sigma: %.3f\n",surfacedensity);
        double vsoft = sqrt(r->G*particle_mass/r->softening);
        printf("vsoft = %.2e\n",vsoft);
	double total_mass = surfacedensity*boxsize*boxsize*4;
        printf("total mass simulated = %.2e msol\n",total_mass);
        printf("Amax = %.3f pc\n",vel_std/OMEGA/KAPPA_OMEGA);
        printf("box length = %.0f pc\n",2*boxsize);
        printf("box boxsize.x=%.2f boxsize.y=%.2f boxsize=%.2f \n",r->boxsize.x, r->boxsize.y, boxsize);


/////////////////////
// add massless particles to simulation 
// these are first!!!!!!!!!!
        add_massless( r, nparticles_nomass, rmassless,vel_stdz_massless);

///////////////////////////////////////
	// Add massive particles to simulation
        int ip0 = r->N;
	for(int ip=ip0;ip<ip0+nparticles;ip++){ 
		struct reb_particle pt;
              	// randomly generate current guiding centers
                double xg  = reb_random_uniform(-r->boxsize.x/2.,r->boxsize.x/2.);
                double yg  = reb_random_uniform(-r->boxsize.x/2.,r->boxsize.x/2.);
        	// randomly generate epicyclic amplitude and phase in plane
                double phi = reb_random_uniform(0.0, 2.0*M_PI); // epicyclic angle
                double A = 2*vel_std/OMEGA/KAPPA_OMEGA*reb_random_uniform(0,1.0);  
                        // epicyclic amplitude
                double x = xg + A*cos(phi); 
		double y = yg - 2.0/KAPPA_OMEGA * A* sin(phi); 
                if (x< -r->boxsize.x/2){ 
                    x +=  r->boxsize.x;  
                    xg+=  r->boxsize.x;  
                }
                if (x> r->boxsize.x/2){   
                    x -=  r->boxsize.x;  
                    xg-=  r->boxsize.x;  
                }
                if (y< -r->boxsize.x/2)  y +=  r->boxsize.x;  
                if (y> r->boxsize.x/2)   y -=  r->boxsize.x;  
                
		pt.x 		= x;
		pt.y 		= y;
		pt.vx 		= -OMEGA*KAPPA_OMEGA*A*sin(phi);
		pt.vy 		=  -3.0/2.0*fourmkap3*OMEGA*xg
                                         - 2.0*OMEGA*A*cos(phi);
              
                double Az = 2*vel_stdz/OMEGAZ *reb_random_uniform(0,1.0);  
                        // epicyclic amplitude
                double phiz = reb_random_uniform(0.0, 2.0*M_PI); // epicyclic angle
		double	z	= Az*cos(phiz); 
                if (z > r->boxsize.z ) z -= r->boxsize.z;
                if (z <-r->boxsize.z ) z += r->boxsize.z;
                pt.z            =z;
		pt.vz 		= Az*OMEGAZ*sin(phiz);
              

		pt.ax= 0; pt.ay= 0; pt.az= 0; 
		pt.r 		= 0.0;		
		pt.m 		= particle_mass;
		pt.hash         = ip;
		reb_add(r, pt);
	}

// more stuff to print!
        double ovx,ovy,ovz;
        vel_disp(r,&ovx,&ovy,&ovz);
        double Q_Toomre = QToom(r, surfacedensity); 
	printf("Toomre Q: %f\n",Q_Toomre);
        printf("vel-dispersion (%.3f %.3f %.3f) tot=%.3f\n",ovx,ovy,ovz, 
              sqrt(ovx*ovx + ovy*ovy + ovz*ovz));


       double tmax = (wait_period+move_period)*period + 0.1;  // numbers of orbital periods
// integrate!
	reb_integrate(r, tmax);
}


void heartbeat(struct reb_simulation* const r){
        static int first =0;
        double period  = 2.0*M_PI/r->ri_sei.OMEGA;
	if (reb_output_check(r, 1e-3*period)){
	   // double KAPPA_OMEGA   = r->ri_sei.KAPPA_OMEGA;
           double Q_Toomre = QToom(r, surfacedensity); 
           printf(" Q=%.2f",Q_Toomre);
           // double ovx,ovy,ovz;
           // vel_disp(r,&ovx,&ovy,&ovz); 
           // printf(" vdisp=(%.3f %.3f %.3f)",ovx,ovy,ovz);

           double sigmax = peak_sigma(r,250.0,0.0);
           printf(" Sgm=%.1f",sigmax);

           // timing for output checks depends on OMEGA
	   reb_output_timing(r, 0);
	
           if (first==0){
              if (Q_Toomre>1.0) { // turn on damping
	          // double OMEGA 	= r->ri_sei.OMEGA;
                  alpham = alpham_om*r->ri_sei.OMEGA;
                  alpha1 = alpha1_om*r->ri_sei.OMEGA;
                  first=1;
              }
           }
	}
        // wrapping of massless particles
        dowrap(r,nparticles_nomass);

        // printouts
        char filename[30];
        static const char froot[10] = "ss";
        static int id=0;
        double fperiods = 0.5;  // print every half period
	if (reb_output_check(r, fperiods*period)){
            mkfilename(filename,froot, id);
            pos_output_ascii(filename,r); // write an ascii file of positions and vels
		//reb_output_ascii("position.txt");
            printf("\n printing out %d\n",id);
            id++;
	}

        static int massless_reset = 0;
        if ((r->t > 0.1) && (massless_reset==0)) {
	   if (reb_output_check(r, wait_period*period))  { // reseting massless particles
              massless_reset=1;
              printf("reseting massless \n");
              reset_massless(r,nparticles_nomass,rmassless);
           }
        }
}

// drag coef is alpha in units of MYR
// if you are doing things different for different particles
// it should be done with hash
void drag_forces(struct reb_simulation* r){
        struct reb_particle* particles = r->particles;
	// double OMEGA 	= r->ri_sei.OMEGA;
	// double KAPPA_OMEGA = r->ri_sei.KAPPA_OMEGA;
        // double fourmkap3= (4.0 - KAPPA_OMEGA*KAPPA_OMEGA)/3.0;

        for(int i=0;i < r->N;i++){
             const struct reb_particle p = particles[i];   
             double dvx = p.vx;   
             // double dvy = p.vy + 3.0/2.0*fourmkap3*p.x*OMEGA;  
                // velocity difference from shearing sheet
             // double dvz = p.vz;
             int ihash = p.hash;
	     if (ihash < nparticles_nomass){ // massless particles
                 particles[i].ax -= alpha1*dvx;
                 // particles[i].ay -= alpha1*dvy;
                 // particles[i].az -= alpha1*dvz; 
             }
             else { // massive particles
                 if (ihash < ntodamp){
                     particles[i].ax -= alpham*dvx;
                     // particles[i].ay -= alpham*dvy;
                     // particles[i].az -= alpham*dvz; 
                 }
             }

        }
}

// give a filename ending with .dat and with file root passed
void mkfilename(char *fname,char *froot, int id){
   char junks[20];
   // strcpy(fname,"");
   strcpy(fname,froot);
   if (id<10) strcat(fname,"0");
   if (id<100) strcat(fname,"0");
   if (id<1000) strcat(fname,"0");
   sprintf(junks,"%d",id);
   strcat(fname,junks);
   strcat(fname,".txt");
   printf("%s\n",fname);
}

// output all particle positions and velocities via hash #
void pos_output_ascii(char *filename,struct reb_simulation* const r){
   FILE *fp;
   fp = fopen(filename,"w");
   fprintf(fp,"# t= %.3f\n",r->t);
   double Q_Toomre = QToom(r, surfacedensity); 
   fprintf(fp,"# Q=%.2f \n", Q_Toomre);
   double sigmax = peak_sigma(r,250.0,0.0);
   fprintf(fp,"# Surface density =%.1f Msol/pc2\n",surfacedensity);
   fprintf(fp,"# Sigmax=%.1f msol/pc2 \n",sigmax);
   double OMEGA       = r->ri_sei.OMEGA;
   double KAPPA_OMEGA = r->ri_sei.KAPPA_OMEGA;
   double kappa = KAPPA_OMEGA*OMEGA;
   fprintf(fp, "# Toomre wavelength: %.1f (pc)\n",
           4.*M_PI*M_PI*surfacedensity*r->G/kappa/kappa);
   fprintf(fp, "# softening: %.1f (pc)\n",r->softening); // pc,   smoothing length
   fprintf(fp, "# N=%d  Nmassive=%d  Nmassless=%d\n",
         r->N, r->N-nparticles_nomass, nparticles_nomass);
   fprintf(fp, "# boxlength %.1f pc \n", r->boxsize.x);

   fprintf(fp,"# i xyz vxvyvz m wxi wyi xg yg A phi \n");
   for(int i=0;i< r->N;i++){
      struct reb_particle* pthash = reb_get_particle_by_hash(r, i); 
      fprintf(fp,"%d ",i);
      fprintf(fp,"%8.3f ",pthash->x);
      fprintf(fp,"%8.3f ",pthash->y);
      fprintf(fp,"%8.3f ",pthash->z);
      fprintf(fp,"%8.3f ",pthash->vx);
      fprintf(fp,"%8.3f ",pthash->vy);
      fprintf(fp,"%8.3f ",pthash->vz);
      fprintf(fp,"%8.3f ",pthash->m);
      if (i < nparticles_nomass){
        fprintf(fp,"%d ",wrapx[i]);
        fprintf(fp,"%d ",wrapy[i]);
      }
      else 
        fprintf(fp,"0  0 ");
      double xg,yg,A,phi;  // compute epicyclic and guiding radii
      compute_epi(r,i,&xg,&yg,&A,&phi);
      fprintf(fp,"%8.3f %8.3f ",xg,yg);
      fprintf(fp,"%8.3f %8.3f ",A,phi);
      fprintf(fp,"\n");
   }
   fclose(fp);
}


// compute velocity dispersion of massive particles
void vel_disp(struct reb_simulation* const r, double *ovx, double *ovy, double *ovz){
     struct reb_particle* particles = r->particles;
     double OMEGA 	= r->ri_sei.OMEGA;
     double KAPPA_OMEGA = r->ri_sei.KAPPA_OMEGA;
     double fourmkap3= (4.0 - KAPPA_OMEGA*KAPPA_OMEGA)/3.0;
     int nmassive=0;
     double vx_sum = 0.0; double vx2_sum = 0.0;
     double vy_sum = 0.0; double vy2_sum = 0.0;
     double vz_sum = 0.0; double vz2_sum = 0.0;
     for(int i=0;i< r->N;i++){
         const struct reb_particle p = particles[i];   
         if (p.m >0){
            vx_sum +=  p.vx;
            double dvy = p.vy + 3.0/2.0*fourmkap3*p.x*OMEGA;  
            vy_sum +=  dvy;
            vz_sum +=  p.vz;
            vx2_sum +=  p.vx*p.vx;
            vy2_sum +=  dvy*dvy;
            vz2_sum +=  p.vz*p.vz;
            nmassive++;
         }
     }
     double vxmean = vx_sum/nmassive;
     double vymean = vy_sum/nmassive;
     double vzmean = vz_sum/nmassive;
     double fac = 1.0*nmassive/(nmassive-1.0);
     double vx_var = vx2_sum/(nmassive-1.) - vxmean*vxmean*fac;
     double vy_var = vy2_sum/(nmassive-1.) - vymean*vymean*fac;
     double vz_var = vz2_sum/(nmassive-1.) - vzmean*vzmean*fac;
     *ovx = sqrt(vx_var); // standard deviations
     *ovy = sqrt(vy_var);
     *ovz = sqrt(vz_var);
}

// Add massless particles to simulation
// these are first
#define CLUSTER 0  // do a cluster setup with rmassless around origin
void add_massless(struct reb_simulation* const r, int nmassless, double rmassless,
         double vel_stdz){
     double OMEGA 	= r->ri_sei.OMEGA;
     double KAPPA_OMEGA = r->ri_sei.KAPPA_OMEGA;
     double fourmkap3= (4.0 - KAPPA_OMEGA*KAPPA_OMEGA)/3.0;
     int ip0 = 0; // massless are first
     // allocate global arrays for keeping track of wrapping
     wrapx = malloc(nmassless*sizeof(int));
     wrapy = malloc(nmassless*sizeof(int));
     oldx  = malloc(nmassless*sizeof(double));
     oldy  = malloc(nmassless*sizeof(double));
     for(int ip=ip0;ip<ip0+nparticles_nomass;ip++){ 
		struct reb_particle pt;
                if (CLUSTER){
		   pt.x 		=  reb_random_normal(rmassless);
		   pt.y 		=  reb_random_normal(rmassless);
		   pt.z 		=  reb_random_normal(rmassless); 
                }
                else{ // all in the line
                   pt.z = 0.0;
                   pt.x = 0.0;
                   pt.y = reb_random_uniform(-r->boxsize.x/2.,r->boxsize.x/2.);
                }
		pt.vx 		= 0.0;
		pt.vy 		= -3.0/2.0*fourmkap3*pt.x*OMEGA;  // shear
		pt.vz 		= 0.0;
		pt.ax= 0; pt.ay= 0; pt.az= 0; 
        	pt.r= 0.0;		
		pt.m 		= 0.0;
		pt.hash         = ip;

                if (vel_stdz > 0.0){
                   double OMEGAZ = r->ri_sei.OMEGAZ;
                   double Az = 2*vel_stdz/OMEGAZ *reb_random_uniform(0,1.0);
                        // epicyclic amplitude
                   double phiz = reb_random_uniform(0.0, 2.0*M_PI); // epicyclic angle
                   double  z       = Az*cos(phiz);
                   if (z > r->boxsize.z ) z -= r->boxsize.z;
                   if (z <-r->boxsize.z ) z += r->boxsize.z;
                   pt.z            =z;
                   pt.vz           = Az*OMEGAZ*sin(phiz);
                }

		reb_add(r, pt);
                oldx[ip] = pt.x;
                oldy[ip] = pt.y;
                wrapx[ip] = 0;
                wrapy[ip] = 0;
    }
}

void reset_massless(struct reb_simulation* const r, int nmassless, double rmassless){
     double OMEGA 	= r->ri_sei.OMEGA;
     double KAPPA_OMEGA = r->ri_sei.KAPPA_OMEGA;
     double fourmkap3= (4.0 - KAPPA_OMEGA*KAPPA_OMEGA)/3.0;
     for(int i=0;i< nmassless;i++){ // these are first!
        struct reb_particle* pthash = reb_get_particle_by_hash(r, i); 
        if (CLUSTER){
		   pthash->x 		=  reb_random_normal(rmassless);
		   pthash->y 		=  reb_random_normal(rmassless);
		   pthash->z 		=  reb_random_normal(rmassless); 
        }
        else{ // all in the line
                   pthash->z = 0.0;
                   pthash->x = 0.0;
                   pthash->y = reb_random_uniform(-r->boxsize.x/2.,r->boxsize.x/2.);
        }
	pthash->vx 		= 0.0;
	pthash->vy 		= -3.0/2.0*fourmkap3*pthash->x*OMEGA;  // shear
	pthash->vz 		= 0.0;
	pthash->ax= 0.0; 
        pthash->ay= 0.0; 
	pthash->az= 0.0; 
        wrapx[i]=0;
        wrapy[i]=0;
        oldx[i] = pthash->x;
        oldy[i] = pthash->y;
     }
}

// keep track of wrapping of massless particles
// these are first
void dowrap(struct reb_simulation* const r, int nmassless){
    for(int i=0;i< nmassless;i++){ // these are first!
       struct reb_particle* pthash = reb_get_particle_by_hash(r, i); 
       double x= pthash->x;
       double y= pthash->y;
       if (x-oldx[i]  >  r->boxsize.x/2) wrapx[i]--;
       if (oldx[i]-x  >  r->boxsize.x/2) wrapx[i]++;
       if (y-oldy[i]  >  r->boxsize.y/2) wrapy[i]--;
       if (oldy[i]-y  >  r->boxsize.y/2) wrapy[i]++;
       oldx[i] = x;
       oldy[i] = y;
    }
}


// estimate a peak surface density use boxes with width x,y width h
double peak_sigma(struct reb_simulation* const r, double h, double y0){
   int nh = (int)(r->boxsize.x /h);
   double *sumarr;
   sumarr = malloc(nh*sizeof(double));
   for(int i=0;i<nh;i++) sumarr[i] = 0.0;
   for(int i=0;i<r->N;i++){   
     double x =  r->particles[i].x;
     double y =  r->particles[i].y;
     double m =  r->particles[i].m;
     if (fabs(y - y0)<h/2){ // add up all mass with |y-y0|<h/2
       int xi = (int)((x + r->boxsize.x/2) /h); // but in box specified by h
       if (xi < nh)
          sumarr[xi] += m;
     }  
   }
   // now find the maximum
   double maxm = 0.0;
   for(int i=0;i<nh;i++){   
      if (sumarr[i] > maxm){
          maxm = sumarr[i];
      }
   }
   double sig = maxm/(h*h); // h*h is area
   return sig;
}


// return Toomre Q parameter right now computed with PI
double QToom(struct reb_simulation* const r, double surfaceden){
   double ovx,ovy,ovz;
   vel_disp(r,&ovx,&ovy,&ovz);
   double OMEGA         = r->ri_sei.OMEGA;
   double KAPPA_OMEGA   = r->ri_sei.KAPPA_OMEGA;
   double Q_Toomre = ovx*KAPPA_OMEGA*OMEGA/(surfaceden * r->G*M_PI);
   // printf(" Q=%.2f  vdisp=(%.3f %.3f %.3f)",Q_Toomre,ovx,ovy,ovz);
   return Q_Toomre;
}


// compute guiding center and epicyclic amplitude and phi for the shearing sheet
// for particle with hash number hashi
// does not take into account wrapping
void compute_epi(struct reb_simulation* const r, int hashi,
  double *xg, double *yg, double *A, double *phi){

   double OMEGA         = r->ri_sei.OMEGA;
   double KAPPA_OMEGA   = r->ri_sei.KAPPA_OMEGA;
   double OMEGA_KAPPA   = 1.0/KAPPA_OMEGA;
   double KAPPA         = KAPPA_OMEGA*OMEGA;
   struct reb_particle* pthash = reb_get_particle_by_hash(r, hashi); 
       // could add in wraps for x, y here 
   double x = pthash->x + wrapx[hashi]*r->boxsize.x;
   double y = pthash->y;
   double vx = pthash->vx;
   double vy = pthash->vy;
   double x_guide = OMEGA_KAPPA*OMEGA_KAPPA *(4.0*x + 2*vy/OMEGA);
   double y_guide = y - OMEGA_KAPPA*OMEGA_KAPPA*2.0*vx/OMEGA;
   double A2 = pow(x - x_guide,2.0) + pow(vx/KAPPA,2.0);
   double phi_ep = atan2(-vx/KAPPA,x-x_guide);
   *xg = x_guide; *yg = y_guide;
   *A = sqrt(A2);  *phi = phi_ep;
}

