#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// ---------- Dimitriou Eleftherios ---------- //

// ---------- Astrodynamics : Project 1 - a) ---------- //

double GM_E = 398600.433;   //  km^3/s^2
double R_E = 6378.0;        //  km
double h, t;
	
void RK(double *X, double *Xn);                  // Runge-Kutta 4
void derivatives(double *Y, double *dotX);       // equations of motion
void orbital_elements(double *X, double *P);     // state vector to orbital elements


int main() {
	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	fp = fopen("H.dat", "w");
	fp1 = fopen("E.dat", "w");
	fp2 = fopen("xy_vs_t_1.dat", "w");
	
	double a;                              // km
    double H = 400.0;                      // km
	a = R_E + H;                           // circular orbit
	double n = sqrt(GM_E/(a*a*a));         // rad/sec  
	double T = 2.0*M_PI/n;                 // sec
	
	printf("The sattelite's period is T = %lf sec\n", T);

	double X[4];
	double Xn[4];
	
	double P[2];
	
    // ---------- initial conditions ---------- //

	X[0] = 0.0;                        // x_0
	X[1] = a;                          // y_0
	X[2] = -sqrt(GM_E/a);              // ux_0
	X[3] = 0.0;                        // uy_0
	
	P[0] = a;
	P[1] = 0.0;

	Xn[0] = 0.0;
	Xn[1] = 0.0;
	Xn[2] = 0.0;
	Xn[3] = 0.0;
	
    // ---------------------------------------- //
    
    double R_s_i, E_s_i, R_s, E_s;
	R_s_i = sqrt(X[0]*X[0] + X[1]*X[1]);                    // initial distance from planet's center
	E_s_i = 0.5*(X[2]*X[2] + X[3]*X[3]) - GM_E/R_s_i;       // initial energy
	
	printf("Initial values (t=0) : a = %.15lf and R_s = %.15lf \n", a, R_s_i);
	printf("E_initial = %.15lf \n", E_s_i);
	fprintf(fp,"%lf %.15lf \n", 0.0, fabs(a - R_E));
	fprintf(fp1,"%lf %.15lf \n", 0.0, E_s_i);
	
	int m, p, max;
	int s = 0;

	h = 0.5;                             // time step = 1 sec
	max = lround(6000.0*T/h);                    
	printf("max steps = %d\n", max);     // total integration time (total steps)
	for(m=1; m <= max; m++) {

		t = m*h;
		RK(X, Xn);

		for(p=0; p<4; p++) {
			X[p] = Xn[p];
		}
		
		orbital_elements(X,P);

		R_s = sqrt(X[0]*X[0] + X[1]*X[1]);                        // new values from integration
		E_s = 0.5*(X[2]*X[2] + X[3]*X[3]) - GM_E/R_s;             // new values from integration
		if(m%100000 == 0)
		{
			fprintf(fp,"%lf %.15lf \n", t, fabs(P[0]-R_E));        // for H vs t plot
		    fprintf(fp1,"%lf %.15lf \n", t, E_s);                  // for energy vs time plot
		    fprintf(fp2,"%lf %.15lf %.15lf \n", t, X[0], X[1]);    // for x,y vs t plot
		}
	}
	
	printf("Energy: E_final = %.15lf \n", E_s);
	
	double H_f;
	H_f = P[0] - R_E;
	
	printf("The error of H is %.15lf km\n", 100.0*fabs(H_f - H)/H);
	
	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	
	return 0;

}


void RK(double* X, double* Xn) {
	int l;
	double k1[4];
	double k2[4];
	double k3[4];
	double k4[4];
	double dotX[4];
	double Y[4];

	for(l=0; l<4; l++) {
		Y[l] = X[l];
	}

	for(l=0; l<4; l++) {
		derivatives(Y, dotX);
		k1[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Y[l] = X[l] + 0.5*k1[l];
	}

	for(l=0; l<4; l++) {
		derivatives(Y, dotX);
		k2[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Y[l] = X[l] + 0.5*k2[l];
	}

	for(l=0; l<4; l++) {
		derivatives(Y, dotX);
		k3[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Y[l] = X[l] + k3[l];
	}

	for(l=0; l<4; l++) {
		derivatives(Y, dotX);
		k4[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Xn[l] = X[l] + (k1[l] + 2.0*k2[l] + 2.0*k3[l] + k4[l])/6.0;
	}

}

void derivatives(double *Y, double *dotX) {
	double rad,rad3;
	rad=sqrt(Y[0]*Y[0]+Y[1]*Y[1]);
	rad3=rad*rad*rad;
	dotX[2] = -GM_E*Y[0]/rad3;              // ux_dot
	dotX[3] = -GM_E*Y[1]/rad3;              // uy_dot
	dotX[0] = Y[2];                         // ux
	dotX[1] = Y[3];                         // uy
	
}

void orbital_elements(double *X, double *P)
{
	double r, E, hx, hy, hz, h, p;
	r = sqrt(X[0]*X[0] + X[1]*X[1]);
	E = 0.5*(X[2]*X[2] + X[3]*X[3]) - GM_E/r;
	hx = 0.0;
	hy = 0.0;
	hz = X[0]*X[3] - X[1]*X[2];
	h = sqrt(hx*hx + hy*hy + hz*hz);
	p = h*h/GM_E;

	P[0] = - GM_E/(2.0*E);	                                                             // (h*h/GM_E)*(1.0/(1.0-P[1]*P[1]));        // a
	P[1] = sqrt(1.0 - p/P[0]);                       //sqrt(1.0 + (h*h/(GM_E*GM_E))*(X[2]*X[2] + X[3]*X[3] - 2.0*GM_E/r));           // e
}
