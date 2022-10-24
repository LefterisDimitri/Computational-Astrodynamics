Copyright 2022, Dimitriou Eleftherios, All rights reserved.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// ---------- Dimitriou Eleftherios - AEM : 4399 ---------- //

// ---------- Astrodynamics : Project 1 - b) ---------- //

double GM_E = 398600.433;   //  km^3/s^2
double R_E = 6378.0;        //  km
double h, t;

void RK(double *X, double *Xn);                 // Runge-Kutta 4
void derivatives(double *Y, double *dotX);      // equations of motion
void orbital_elements(double *X, double *P);    // state vector to orbital elements
double Kepler_solve(double n, double t);        // solving the Kepler's equation


int main() {
	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	fp = fopen("orbit.dat", "w");
	fp1 = fopen("xy_vs_t.dat", "w");
	fp2 = fopen("gauss.dat", "w");
	fp3 = fopen("a_e_vs_t.dat", "w");

	double a;
    double H = 400.0;                      // km
	a = R_E + H;                           // initial circular orbit
	double n = sqrt(GM_E/(a*a*a));         // rad / sec  
	double T = 2.0*M_PI/n;                 // sec
	
	printf("The sattelite's period is T = %lf and n = %.15lf\n", T, n);

	double X[6];
	double Xn[6];
	double P[2];
	
	// ---------- initial conditions ---------- //

	X[0] = 0.0;                        // x_0
	X[1] = a;                          // y_0
	X[2] = -sqrt(GM_E/a);              // ux_0
	X[3] = 0.0;                        // uy_0
	
	X[4] = a;                          // initial a for gauss
	X[5] = pow(10.0,-15);              // initial e for gauss
	
	P[0] = a;                          // a
	P[1] = pow(10.0,-15);              // e

	Xn[0] = 0.0;
	Xn[1] = 0.0;
	Xn[2] = 0.0;
	Xn[3] = 0.0;
	Xn[4] = 0.0;
	Xn[5] = 0.0;
	
	
    // ---------------------------------------- //


    double R_s_i, E_s_i, R_s, E_s;
	R_s_i = sqrt(X[0]*X[0] + X[1]*X[1]);                    // initial distance from planet's center
	E_s_i = 0.5*(X[2]*X[2] + X[3]*X[3]) - GM_E/R_s_i;       // initial energy
	
	printf("E_initial = %.15lf \n", E_s_i);
	
	int m, p, max;

	h = 0.5;     
	max = lround(6000.0*T/h);                    
	printf("max steps = %d\n", max);
	for(m=1; m <= max; m++) {
		
		t = m*h;
		
		RK(X, Xn);
		
		for(p=0; p<6; p++) {
			X[p] = Xn[p];
		}
		 
		orbital_elements(X,P);
		
		fprintf(fp,"%.15lf %.15lf \n", X[0], X[1]);                // for orbit
		if(m%100000 == 0)
		{
			fprintf(fp1,"%lf %.15lf %.15lf \n", t, X[0], X[1]);    // for x,y vs t plot
			fprintf(fp2,"%lf %.15lf %.15lf \n", t, X[4], X[5]);    // a,e from gauss equations integration
			fprintf(fp3,"%lf %.15lf %.15lf \n", t, P[0], P[1]);    // a,e from integration
		}
		
		R_s = sqrt(X[0]*X[0] + X[1]*X[1]);
		E_s = 0.5*(X[2]*X[2] + X[3]*X[3]) - GM_E/R_s;
		
		if(R_s <= R_E) {
			printf("Distance from Earth's surface is d = %.15lf \n", R_s);
			printf("Time until impact: t = %lf sec", t);
			return 0;
		}

	}

	printf("E_final = %.15lf \n", E_s);
	
	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	return 0;

}


void RK(double *X, double *Xn) {
	int l;
	double k1[6];
	double k2[6];
	double k3[6];
	double k4[6];
	double dotX[6];
	double Y[6];
	
	for(l=0; l<6; l++) {
		Y[l] = X[l];
	}

	for(l=0; l<6; l++) {
		derivatives(Y, dotX);
		k1[l] = h*dotX[l];
	}

	for(l=0; l<6; l++) {
		Y[l] = X[l] + 0.5*k1[l];
	}

	for(l=0; l<6; l++) {
		derivatives(Y, dotX);
		k2[l] = h*dotX[l];
	}

	for(l=0; l<6; l++) {
		Y[l] = X[l] + 0.5*k2[l];
	}

	for(l=0; l<6; l++) {
		derivatives(Y, dotX);
		k3[l] = h*dotX[l];
	}

	for(l=0; l<6; l++) {
		Y[l] = X[l] + k3[l];
	}

	for(l=0; l<6; l++) {
		derivatives(Y, dotX);
		k4[l] = h*dotX[l];
	}
	
	for(l=0; l<6; l++) {
		Xn[l] = X[l] + (k1[l] + 2.0*k2[l] + 2.0*k3[l] + k4[l])/6.0;
	}
	
}

void derivatives(double *Y, double *dotX) {
	double H;
	double rho, C, S, fx, fy, f, u, ms;
	double rad,rad3;
	rad=sqrt(Y[0]*Y[0]+Y[1]*Y[1]);
	rad3=rad*rad*rad;
	u=sqrt(Y[2]*Y[2] + Y[3]*Y[3]);
	
	double E, a, theta, n_der;
	E = 0.5*(Y[2]*Y[2] + Y[3]*Y[3]) - GM_E/rad;
	a = -GM_E/(2.0*E);
	n_der = sqrt(GM_E/(a*a*a));

	/*
	ex = (1.0/GM_E)*(u*u - GM_E/rad)*Y[0] - (Y[0]*Y[2] + Y[1]*Y[3])*Y[3]/GM_E;
	ey = (1.0/GM_E)*(u*u - GM_E/rad)*Y[1] - (Y[0]*Y[2] + Y[1]*Y[3])*Y[4]/GM_E;
	e = sqrt(ex*ex + ey*ey);
	
	if(ur >= 0.0)
	{
		theta = acos((ex*Y[0] + ey*Y[1])/(e*rad));
	}
	else
	{
		theta = 2.0*M_PI - acos((ex*Y[0] + ey*Y[1])/(e*rad));
	}
	*/
	
	H = rad - R_E;                

	if (H >= 15.0 && H <= 400.0) {
		rho = pow(0.1*H,-7.5)*pow(10.0, 9);
	} else if (H > 0.0 && H < 15.0) {
		rho = 0.1*pow(10.0, 9);
	}

	C = 1.5;
	S = 8000.0*pow(10.0, -6);    // km^(-2)
	ms = 400000.0;               // kg

	f = 0.5*rho*u*u*C*S/ms;
	fx = f*Y[2]/u;
	fy = f*Y[3]/u;
	
	dotX[2] = -GM_E*Y[0]/rad3 - fx;         // ux_dot
	dotX[3] = -GM_E*Y[1]/rad3 - fy;         // uy_dot
	dotX[0] = Y[2];                               // ux
	dotX[1] = Y[3];                               // uy
	
	// ---------- Gauss ODEs for e = 0 ---------- //
	double p_R, p_S;
	
	theta = Kepler_solve(n_der,t);
	p_R = fx*cos(theta) + fy*sin(theta);
	p_S = -fx*sin(theta) + fy*cos(theta);
	dotX[4] = -rho*sqrt(GM_E*a)*C*S/ms;                                   // da/dt
	dotX[5] = (1.0/(n_der*a))*(sin(theta)*p_R + 2.0*cos(theta)*p_S);      // de/dt
	// ------------------------------------------ //
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

double Kepler_solve(double n, double t)
{
	double f;
	
	f = n*t;  // only for circular orbit (e=0)

	return f;
}
