/* ilmsfxd.c: Quartic potential, Fixed boundary conditions */
/* This program can be continued from where it left off */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define chainlngth 100
#define nsprngs (chainlngth+1)
#define nmode 100  /* Ignore; not used here */
#define nprntstps 10001 /* Number of output lines */
#define dt 0.001 /* Time step */
#define beta 0.7 /* Beware! beta is the nonlinear coefficient! */
		/* Usually alpha and beta appaears interchanged in literature */
#define alpha .16 /* alpha is the coefficient of the linear term! */

void accel(double *, double *);

int main() {
	struct timeval start, end;
	double delta;
	printf("Alpha is:  %lf \n", alpha);
	printf("Beta is :  %lf \n", beta);
	
	#pragma omp parallel
		printf("Hello, world.\n");
	
	FILE *fp, *fp1, *fp2, *fp3, *fp5, *fp6, *fp7, *fp8;
	int k, n1, prncnt, prntstps;
	double v[chainlngth], x[chainlngth], tke, tpe, te;
	double acc[chainlngth], ke[chainlngth], pe[nsprngs], fac, y[chainlngth];
	double dx, hdt, hdt2, alphaby4, cmass, cmom;
	
	prntstps = (int) (1.0 / dt);
	
	// Start Timing
	gettimeofday(&start, NULL);
	
	//char *buf=(char *)malloc(sizeof(char)*10000000);
	//setvbuf(fp1, buf, _IOFBF, sizeof(buf));
	
	alphaby4 = beta / 4.0;
	dx = tke = tpe = te = 0.0;
	
	// Open Files for Writing
	fp = fopen("toten.dat","w");
	fp1 = fopen("strsh.dat","w");
	fp3 = fopen("velsh.dat","w");
	fp5 = fopen("cmass.dat","w");
	fp6 = fopen("ke.dat","w");
	fp7 = fopen("acce.dat","w");
	fp8 = fopen("pe.dat","w");
	
	/* Initialize the position, velocity, acceleration arrays */
	for (int a = 0; a < chainlngth; a++) { 
		x[a] = 0.0;
		acc[a] = 0.0;
		v[a] = 0.0;
	}
	
	/* Initial perturbation at the center of the chain and it is a double
	particle perturbation (EVEN Parity) */
	x[50] = -0.98;
	x[51] = +0.98; // Even Parity
	
	
	for (int a = 0; a < chainlngth; a++) { 
		ke[a] = 0.0;
		fprintf(fp6,"%.10f\t", ke[a]);
		if (a == 0) {
			dx = x[a];
		} else {
			dx = x[a] - x[a - 1]; 
		}
		fac = dx * dx;
		pe[a] = alpha * 0.5 * fac + alphaby4 * fac * fac;
		fprintf(fp8,"%.10f\t", pe[a]);
		tpe += pe[a];
	}
	fprintf(fp6,"\n");
	
	dx = -x[chainlngth - 1];
	fac = dx * dx;
	pe[chainlngth - 1] = alpha * 0.5 * fac + alphaby4 * fac * fac;
	fprintf(fp8,"%.10f\n", pe[chainlngth - 1]);
	tpe += pe[chainlngth - 1];
	te = tpe + tke;
	
	fprintf(fp,"%d\t%.10f\n", 0, te);
	
	for (int c = 0; c < chainlngth; c++) {
		fprintf(fp1,"%.10f\t", x[c]);
		fprintf(fp3,"%.10f\t", v[c]);
		fprintf(fp7,"%.10f\t", acc[c]);
	}
	fprintf(fp1,"\n"); 
	fprintf(fp3,"\n");
	fprintf(fp7,"\n"); 
	
	///////////////////////////////////////
	
	hdt = 0.5 * dt;
	hdt2 = dt * hdt;
	n1 = 1;
	for (int n = 1; n < nprntstps;) {
		
		
		/* new positions and mid-velocities; velocity-Verlet algorithm  */
		
		for (int b = 0; b < chainlngth; b++) {
			x[b] += dt * v[b] + hdt2 * acc[b];
			v[b] += hdt * acc[b];
		}
		
		/* new accelerations */
		accel(x, acc);
		
		/* new final velocities; Ignore the variables cmom and cmass */
		cmom = 0.0;
		for (int b = 0; b < chainlngth; b++) {
			v[b] += hdt * acc[b];
			cmom += v[b];
		}
		cmom /= chainlngth;
		
		/* Kinetic energies */
		prncnt = n1 / prntstps;  //percent completion
		if (prncnt == 1) { //if 100% done, print all
			tke = tpe = te = dx = 0.0;  //reset all variables
			cmass = 0.0;
			for (int b = 0; b < chainlngth; b++) {
				ke[b] = 0.5 * v[b] * v[b]; 
				fprintf(fp6,"%.10f\t",ke[b]);
				tke += ke[b];
				if (b == 0) {
					dx = x[b];
				} else {
					dx = x[b] - x[b - 1];
				}
				fac = dx * dx;
				double temp = 0;
				temp = alpha * 0.5 * fac + alphaby4 * fac * fac;
				fprintf(fp8,"%.10f\t", temp);
				tpe += temp;
				cmass += x[b];
			}
			fprintf(fp6,"\n");
			
			dx = -x[chainlngth - 1];
			fac = dx * dx;
			
			double temp2 = 0;
			temp2 = alpha * 0.5 * fac + alphaby4 * fac * fac;
			tpe += temp2;
			
			fprintf(fp8,"%.10f\n", temp2);
			
			fprintf(fp5, "%d\t%.10f\n", 0, cmass);
			
			cmass /= chainlngth;
			te = tpe + tke;
			
			fprintf(fp,"%d\t%.10f\n", n, te);
			
			for (k=0; k < chainlngth; k++) {
				y[k] = x[k] - cmass;
				fprintf(fp1,"%.10f\t", y[k]); 
				fprintf(fp3,"%.10f\t", v[k]); 
				fprintf(fp7,"%.10f\t", acc[k]); 
			}
			fprintf(fp1,"\n"); 
			fprintf(fp3,"\n"); 
			fprintf(fp7,"\n"); 
			
			n1 = 1;
			n++;
		}
		n1++;
	}
	
	// Close Files
	fclose(fp);
	fclose(fp1);
	fclose(fp3);
	fclose(fp5);
	fclose(fp6);
	fclose(fp7);
	fclose(fp8);
	
	fp2 = fopen("restart.dat","w");
	fprintf(fp2, "%d\t%d\n", -1, nprntstps);
	for (int b = 0; b < chainlngth; b++) {
		fprintf(fp2, "%.15f\t%.15f\t%.15f\n", x[b], v[b], acc[b]);
	}
	fclose(fp2);
	
	// Time how long the operation took
	gettimeofday(&end, NULL);
	delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
	end.tv_usec - start.tv_usec) / 1.e6;
	printf("Total time=%f seconds\n", delta);
}

void accel(double *x, double *acc) {
	for (int a = 0; a < chainlngth; a++) {
		int j = a - 1;
		int k = a + 1;
		double dximn1 = 0.0;
		double dxipl1 = 0.0;
		if (j == -1) {
			dximn1 = 0.0;
		} else {
			dximn1 = x[j];
		}
		if (k == chainlngth) {
			dxipl1 = 0.0;
		} else {
			dxipl1 = x[k];
		}
		double dxi = x[a];
		double twodxi = 2.0 * dxi;
		double fac = dxipl1 + dximn1 - twodxi;
		double fac1 = dxipl1 - dxi;
		double fac2 = dxi - dximn1;
		double fac13 = fac1 * fac1 * fac1;
		double fac23 = fac2 * fac2 * fac2;
		acc[a] = alpha * fac + beta * (fac13 - fac23);
	}
}

