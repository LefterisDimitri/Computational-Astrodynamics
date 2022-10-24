Copyright 2022, Dimitriou Eleftherios, All rights reserved.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// ---------- Dimitriou Eleftherios - AEM : 4399 ---------- //

// ---------- Astrodynamics : Project 1 - c) ---------- //

double GM_E = 398600.433;   //  km^3/s^2
double R_E = 6378.0;        //  km
double h;


void RK(double t, double *X, double *Xn);              // Runge-Kutta 4
void derivatives(double t, double *Y, double *dotX);   // equations of motion
void orbital_elements(double *Xn, double *P);          // state vector to orbital elements


int main() {
	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	fp = fopen("a_vs_t.dat","w");
	fp1 = fopen("e_vs_t.dat","w");
	fp2 = fopen("xy_vs_t_3.dat", "w");

	double a;
    double H = 400.0;                      // km
	a = R_E + H;                           // circular orbit
	double n = sqrt(GM_E/(a*a*a));         // rad / sec  
	double T = 2.0*M_PI/n;                 // sec
	
	printf("The sattelite's period is T = %lf and n = %.15lf\n", T, n);
	
	double X[4];
	double Xn[4];
	double P[2];

	// ---------- initial conditions ---------- //

	X[0] = 0.0;                        // x_0
	X[1] = a;                          // y_0
	X[2] = -sqrt(GM_E/a);              // ux_0
	X[3] = 0.0;                        // uy_0
	
	P[0] = a;                          // a
	P[1] = 0.0;                        // e

	Xn[0] = 0.0;
	Xn[1] = 0.0;
	Xn[2] = 0.0;
	Xn[3] = 0.0;

    // ---------------------------------------- //
	
	fprintf(fp,"%lf %lf \n", 0.0, P[0]);
	fprintf(fp1,"%lf %lf \n", 0.0, P[1]);
	
	double R_s_i, E_s_i, R_s, E_s;
	R_s_i = sqrt(X[0]*X[0] + X[1]*X[1]);                    // initial distance from planet's center
	E_s_i = 0.5*(X[2]*X[2] + X[3]*X[3]) - GM_E/R_s_i;       // initial energy
	
	printf("E_initial = %.15lf \n", E_s_i);
	
	double s = 0.0;
	int m, p, max;
	double t;
	h = 1.0;      
	max = lround(6000.0*T/h); 
	printf("max = %d\n", max);
	for(m=1; m <= max; m++) {
		
		t = m*h;
		
		RK(t, X, Xn);
		
		for(p=0; p<4; p++) {
			X[p] = Xn[p];
		}
		
		orbital_elements(X,P);
		
		//fprintf(fp,"%lf %.15lf \n", t, P[0]);
	    //fprintf(fp1,"%lf %.15lf \n", t, P[1]);
	    
	    if(m%10000 == 0)
	    {
	    	s++;
	    	fprintf(fp2,"%lf %.15lf %.15lf \n", t, X[0], X[1]);  
	    	fprintf(fp,"%lf %.15lf \n", t, P[0]);
	        fprintf(fp1,"%lf %.15lf \n", t, P[1]);
		}
		
		R_s = sqrt(X[0]*X[0] + X[1]*X[1]);		
	}
	
	printf("s = %lf\n", s);
	
	E_s = 0.5*(X[2]*X[2] + X[3]*X[3]) - GM_E/R_s_i;

	printf("E_final = %.15lf \n", E_s);
	
	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	
	return 0;

}


void RK(double t, double* X, double* Xn) {
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
		derivatives(t, Y, dotX);
		k1[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Y[l] = X[l] + 0.5*k1[l];
	}

	for(l=0; l<4; l++) {
		derivatives(t, Y, dotX);
		k2[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Y[l] = X[l] + 0.5*k2[l];
	}

	for(l=0; l<4; l++) {
		derivatives(t, Y, dotX);
		k3[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Y[l] = X[l] + k3[l];
	}

	for(l=0; l<4; l++) {
		derivatives(t, Y, dotX);
		k4[l] = h*dotX[l];
	}

	for(l=0; l<4; l++) {
		Xn[l] = X[l] + (k1[l] + 2.0*k2[l] + 2.0*k3[l] + k4[l])/6.0;
	}

}

void derivatives(double t, double *Y, double *dotX) {
	double peekaboo, AB;
	double xs, ys, uxs, uys, urel, urel_x, urel_y, u_sun, theta, angular_vel_sun, xrel, yrel, drel;
	double H;
	double rho, C, S, fx, fy, u, ms;
	double rad,rad3;
	rad=sqrt(Y[0]*Y[0]+Y[1]*Y[1]);
	rad3=rad*rad*rad;
	u=sqrt(Y[2]*Y[2] + Y[3]*Y[3]);
	
	angular_vel_sun = 0.0174532925/86400.0;   // rad/sec
	//u_sun = angular_vel_sun*149597871.0;    
	
	// 1 AU -> km
	theta = angular_vel_sun*t;
	xs = 149597871.0*cos(theta);     // km 
	ys = 149597871.0*sin(theta);     // km
	//uxs = u_sun*sin(theta);
	//uys = u_sun*cos(theta);
	//urel_x = Y[2]-uxs;
	//urel_y = Y[3]-uys;
	//urel=sqrt(urel_x*urel_x + urel_y*urel_y);
	
	xrel = Y[0] - xs;
	yrel = Y[1] - ys;
	drel = sqrt(xrel*xrel + yrel*yrel);
	
	ms = 400000.0;                                   // kg
	S = 0.008;  // km^2
	C = 1.5;
	
	AB = sqrt((Y[0]*Y[0] + Y[1]*Y[1]) - ((Y[0]*xs + Y[1]*ys)/sqrt(xs*xs + ys*ys))*((Y[0]*xs + Y[1]*ys)/sqrt(xs*xs + ys*ys)));
	if((Y[0]*xs + Y[1]*ys) < 0.0 && AB < R_E)
	{
		peekaboo = 0.0;
	}
	else
	{
		peekaboo = 1.0;
	}
	
	double phi_over_c = 4.56/(10.0*10.0*10.0);

	fx = - (peekaboo*phi_over_c*C*S/ms)*xrel/drel;
	fy = - (peekaboo*phi_over_c*C*S/ms)*yrel/drel;

	dotX[2] = -GM_E*Y[0]/rad3 + fx;    // ux_dot
	dotX[3] = -GM_E*Y[1]/rad3 + fy;    // uy_dot
	dotX[0] = Y[2];                    // ux
	dotX[1] = Y[3];                    // uy
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

	P[0] = -GM_E/(2.0*E);	                                                             // (h*h/GM_E)*(1.0/(1.0-P[1]*P[1]));        // a
	P[1] = sqrt(1.0 - p/P[0]);                       //sqrt(1.0 + (h*h/(GM_E*GM_E))*(X[2]*X[2] + X[3]*X[3] - 2.0*GM_E/r));           // e
	
}
