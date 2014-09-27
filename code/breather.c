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
	double a = alpha;
	double b = beta;
	printf("Alpha is:  %lf \n", a);
	printf("Beta is :  %lf \n", b);
	
	FILE *fp, *fp1, *fp2, *fp3, *fp5, *fp6, *fp7, *fp8;
	int i, j, k, n, n1, prncnt, prntstps;
	double v[chainlngth], x[chainlngth], tke, tpe, te;
	double acc[chainlngth], ke[chainlngth], pe[nsprngs], fac, y[chainlngth];
	double dx, hdt, hdt2, alphaby4, cmass, cmom;
	
	prntstps = (int) (1.0 / dt);
	
	// Start Timing
	gettimeofday(&start, NULL);
	
	//char *buf=(char *)malloc(sizeof(char)*10000000);
	//setvbuf(fp1, buf, _IOFBF, sizeof(buf));
	alphaby4 = beta / 4.0;
	fp = fopen("toten.dat","w");
	fp1 = fopen("strsh.dat","w");
	fp3 = fopen("velsh.dat","w");
	fp5 = fopen("cmass.dat","w");
	fp2 = fopen("restart.dat","w");
	fp6 = fopen("ke.dat","w");
	fp7 = fopen("acce.dat","w");
	fp8 = fopen("pe.dat","w");
	
	/* Initialize the position, velocity, acceleration arrays */
	for (int a = 0; a < chainlngth; a++) { 
		x[a] = 0.0;
		acc[a] = 0.0;
		v[a] = 0.0;
	}
	
	/*
	Initial perturbation at the center of the chain and it is a double particle
	perturbation (EVEN Parity)
	*/
	x[50] = -0.98;
	x[51] = +0.98; // Even Parity
	
	dx = tke = tpe = te = 0.0;
	
	for (int a = 0; a < chainlngth; a++) { 
		ke[a] = 0.5 * v[a] * v[a];
		fprintf(fp6,"%.10f\t", ke[a]);
		tke += ke[a];
		j = a - 1;
		if (j == -1) {
			dx = x[a];
		} else {
			dx = x[a] - x[j]; 
		}
		fac = dx * dx;
		pe[a] = alpha * 0.5 * fac + alphaby4 * fac * fac;
		fprintf(fp8,"%.10f\t", pe[a]);
		tpe += pe[a];
	}
	dx = -x[chainlngth - 1];
	fac = dx * dx;
	pe[chainlngth - 1] = alpha * 0.5 * fac + alphaby4 * fac * fac;
	fprintf(fp8,"%.10f\n", pe[chainlngth - 1]);
	tpe += pe[chainlngth - 1];
	fprintf(fp6,"\n");
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
	i = 0;
	
	hdt = 0.5 * dt;
	hdt2 = dt * hdt;
	n = 1; n1 = 1;
	while (n < nprntstps) {
		
		/* new positions and mid-velocities; velocity-Verlet algorithm  */
		
		for (j = 0; j < chainlngth; j++) {
			x[j] += dt * v[j] + hdt2 * acc[j];
			v[j] += hdt * acc[j];
		}
		
		/* new accelerations */
		accel(x, acc);
		
		/* new final velocities; Ignore the variables cmom and cmass */
		cmom = 0.0;
		for (j = 0; j < chainlngth; j++) {
			v[j] += hdt * acc[j];
			cmom += v[j];
		}
		cmom /= chainlngth;
		
		/* Kinetic energies */
		prncnt = n1 / prntstps;  //percent completion
		if (prncnt == 1) {  //if 100% done, print all
			tke = tpe = te = dx = 0.0;  //reset all variables
			cmass = 0.0;
			for (j = 0; j < chainlngth; j++) {
				ke[j] = 0.5 * v[j] * v[j]; 
				fprintf(fp6,"%.10f\t",ke[j]);
				tke += ke[j];
				k = j-1;
				if (k == -1) {
					dx = x[j];
				} else {
					dx = x[j] - x[k];
				}
				fac = dx * dx;
				pe[i] = alpha * 0.5 * fac + alphaby4 * fac * fac;
				fprintf(fp8,"%.10f\t", pe[i]);
				tpe += pe[i];
				cmass += x[j];
			}
			dx = -x[chainlngth - 1];
			fac = dx * dx;
			pe[i] = alpha * 0.5 * fac + alphaby4 * fac * fac;
		
			fprintf(fp8,"%.10f\n", pe[i]);
			tpe += pe[i];
			fprintf(fp5, "%d\t%.10f\n", i, cmass);
			cmass /= chainlngth;
			fprintf(fp6,"\n");
			te = tpe + tke;
			fprintf(fp,"%d\t%.10f\n", n, te);
			for (k=0; k < chainlngth; k++) {
				y[k] = x[k] - cmass;
				fprintf(fp1,"%.10f\t", y[k]); 
			}
			fprintf(fp1,"\n"); 
			for (k=0; k < chainlngth; k++) {
				fprintf(fp3,"%.10f\t", v[k]); 
			}
			fprintf(fp3,"\n"); 
			for (k=0; k < chainlngth; k++) {
				fprintf(fp7,"%.10f\t", acc[k]); 
			}
			fprintf(fp7,"\n"); 
			n++; n1 = 1;
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
	
	fprintf(fp2, "%d\t%d\n", i-1, n-1);
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


