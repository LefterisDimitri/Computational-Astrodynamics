Copyright 2022, Dimitriou Eleftherios, All rights reserved.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// ---------- Astrodynamics : Project 2 ---------- //

double GM_E = 398600.433;   //  km^3/s^2
double R_E = 6378.0;        //  km
double h, t;

double J2 = 0.00108263;

void orbital_elements_to_state_vec(double *P, double *X);               // stoixeia troxias ---> theseis kai taxuthtes
void state_vec_to_orbital_elements(double *X, double *P);               // theseis kai taxuthtes ---> stoixeia troxias

void RK(double *X, double *Xn);
void derivatives(double *Y, double *dotX);
void orbital_elements(double *X, double *P);


int main() {
	FILE *file_LEO[20];
	for (int pointer = 0; pointer < 20; pointer++) {
		char filename[100];
		sprintf(filename, "orbits_LEO_%d.dat", pointer);
		file_LEO[pointer] = fopen(filename, "w");

	}

	double H_LEO = 400.0;
	double aL = R_E + H_LEO;
	double nL = sqrt(GM_E/(aL*aL*aL));
	double TL = 2.0*M_PI/nL;

	printf("LEO: The sattelite's period is T = %lf and n = %.15lf\n", TL, nL);
	
	double XL[6];
	double PL[6];
	double XLn[6];

	int i;
	for(i=0; i<6; i++) {
		XL[i] = 0.0;
		PL[i] = 0.0;
		XLn[i] = 0.0;
	}

	// ---------- LEO orbit: initial conditions ---------- //

	PL[0] = aL;
	PL[1] = pow(10.0,-10);
	PL[2] = 0.0;
	PL[3] = 0.0;
	PL[4] = 0.0;
	PL[5] = 0.0;
	
	// ---------------------------------------- //
	

	double i_matrix_LEO[20];

	for(i=0; i<20; i++) {
		i_matrix_LEO[i] = (85.0 + i*1.0)*M_PI/180.0;
	}

	for(i=0; i<20; i++) {
		printf("i = %lf \n", i_matrix_LEO[i]);
	}

	double dOmegaL_dt, domega_tinyL_dt;

	int m, p, maxL, j;

	h = 1.0;
	maxL = lround(1000.0*TL/h);
	printf("maxL = %d\n", maxL);
	for(j=0; j<20; j++) {
		
		// ---------- LEO orbit: initial conditions ---------- //

		PL[0] = aL;                       // a
		PL[1] = pow(10.0,-10);            // e
		PL[2] = i_matrix_LEO[j];          // i
		PL[3] = 0.0;                      // Omega
		PL[4] = 0.0;                      // omega_tiny
		PL[5] = 0.0;                       // f

		// ---------------------------------------- //

		//PL[2] = i_matrix_LEO[j];
	
		orbital_elements_to_state_vec(PL,XL);
		
		for(m=1; m <= maxL; m++) 
		{
			
			t = m*h;
			
			RK(XL, XLn);
			
			for(p=0; p<6; p++) 
			{
				XL[p] = XLn[p];
			}

			state_vec_to_orbital_elements(XL,PL);

			dOmegaL_dt = -(3.0*J2*R_E*R_E*nL/(2.0*PL[0]*PL[0]*(1.0 - PL[1]*PL[1])*(1.0 - PL[1]*PL[1])))*cos(PL[2]);

			domega_tinyL_dt = (3.0*J2*R_E*R_E*nL/(2.0*PL[0]*PL[0]*(1.0 - PL[1]*PL[1])*(1.0 - PL[1]*PL[1])))*(2.0 - 5.0*sin(PL[2])*sin(PL[2])/2.0);

			if(m%10000 == 0)
			{
				fprintf(file_LEO[j],"%lf %lf %lf %lf %lf %lf %.15lf %.15lf\n", t, PL[0], PL[1], PL[2], PL[3], PL[4], dOmegaL_dt, domega_tinyL_dt);
			}
		}

	}

	/*int p_i;
	for(p_i = 0; p_i < 20; p_i ++) {
		fclose(file_LEO[p_i]);
	}*/

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
	double constant, px, py, pz;
	double rad, rad2, rad3, rad4;
	rad=sqrt(Y[0]*Y[0] + Y[1]*Y[1] + Y[2]*Y[2]);
	rad2 = Y[0]*Y[0] + Y[1]*Y[1] + Y[2]*Y[2];
	rad3=rad*rad*rad;
	rad4=rad3*rad;

	constant = (3.0/2.0)*J2*GM_E*R_E*R_E/rad4;
	px = constant*(Y[0]/rad)*(5.0*Y[2]*Y[2]/rad2 - 1.0);
	py = constant*(Y[1]/rad)*(5.0*Y[2]*Y[2]/rad2 - 1.0);
	pz = constant*(Y[2]/rad)*(5.0*Y[2]*Y[2]/rad2 - 3.0);


	dotX[3] = -GM_E*Y[0]/rad3 + px;               // ux_dot
	dotX[4] = -GM_E*Y[1]/rad3 + py;               // uy_dot
	dotX[5] = -GM_E*Y[2]/rad3 + pz;               // uz_dot
	dotX[0] = Y[3];                               // ux
	dotX[1] = Y[4];                               // uy
	dotX[2] = Y[5];                               // uz

}

void state_vec_to_orbital_elements(double *X, double *P) {
	double r_svtoe, u_svtoe, E_svtoe;
	double hx, hy, hz, h, p_svtoe, N, nx, ny, nz, ur, ex, ey, ez, e;
	double n_svtoe;

	r_svtoe = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
	u_svtoe = sqrt(X[3]*X[3] + X[4]*X[4] + X[5]*X[5]);
	E_svtoe = 0.5*(X[3]*X[3] + X[4]*X[4] + X[5]*X[5]) - GM_E/r_svtoe;

	hx = X[1]*X[5] - X[2]*X[4];
	hy = X[2]*X[3] - X[0]*X[5];
	hz = X[0]*X[4] - X[1]*X[3];
	h = sqrt(hx*hx + hy*hy + hz*hz);
	p_svtoe = h*h/GM_E;
	
	ex = (1.0/GM_E)*(u_svtoe*u_svtoe - GM_E/r_svtoe)*X[0] - (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])*X[3]/GM_E;
	ey = (1.0/GM_E)*(u_svtoe*u_svtoe - GM_E/r_svtoe)*X[1] - (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])*X[4]/GM_E;
	ez = (1.0/GM_E)*(u_svtoe*u_svtoe - GM_E/r_svtoe)*X[2] - (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])*X[5]/GM_E;
	e = sqrt(ex*ex + ey*ey + ez*ez);

	nx = -hy;
	ny = hx;
	nz = 0.0;
	N = sqrt(nx*nx + ny*ny + nz*nz);

	ur = (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])/r_svtoe;

	P[0] = -GM_E/(2.0*E_svtoe);
	
	n_svtoe = sqrt(GM_E/(P[0]*P[0]*P[0]));
	
	P[1] = sqrt(1.0 - p_svtoe/P[0]);
	
	P[2] = acos(hz/h);

	if(ny >= 0.0) {
		P[3] = acos(nx/N);
	} else {
		P[3] = 2.0*M_PI - acos(nx/N);
	}

	if(ez >= 0.0) {
		P[4] = acos((nx*ex + ny*ey + nz*ez)/(N*e));
	} else {
		P[4] = 2.0*M_PI - acos((nx*ex + ny*ey + nz*ez)/(N*e));
	}
	
	//P[5] = n_svtoe*t;

	if(ur >= 0.0) {
		P[5] = n_svtoe*t;
	} else {
		P[5] = 2.0*M_PI - n_svtoe*t;
	}
}

void orbital_elements_to_state_vec(double *P, double *X) {
	int i, j, k;
	double r_oetsv;
	double r_PQW[3];
	double u_PQW[3];
	double r_ECI[3];
	double u_ECI[3];
	double temp_matrix[3][3];
	double T_PQW_to_ECI[3][3];

	r_oetsv = P[0]*(1.0-P[1]*P[1])/(1.0 + P[1]*cos(P[5]));

	double R3_tr[3][3] = {
		{cos(P[3]), -sin(P[3]), 0.0},
		{sin(P[3]), cos(P[3]), 0.0},
		{0.0, 0.0, 1.0}
	};

	double R1_tr[3][3] = {
		{1.0, 0.0, 0.0},
		{0.0, cos(P[2]), -sin(P[2])},
		{0.0, sin(P[2]), cos(P[2])}
	};

	double R2_tr[3][3] = {
		{cos(P[4]), -sin(P[4]), 0.0},
		{sin(P[4]), cos(P[4]), 0.0},
		{0.0, 0.0, 1.0}
	};

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			temp_matrix[i][j] = 0.0;
			T_PQW_to_ECI[i][j] = 0.0;
		}
		r_ECI[i] = 0.0;
		u_ECI[i] = 0.0;
	}

	// ---------- PQW ---------- //
	r_PQW[0] = r_oetsv*cos(P[5]);
	r_PQW[1] = r_oetsv*sin(P[5]);
	r_PQW[2] = 0.0;

	u_PQW[0] = sqrt(GM_E/(P[0]*(1.0 - P[1]*P[1])))*(-sin(P[5]));
	u_PQW[1] = sqrt(GM_E/(P[0]*(1.0 - P[1]*P[1])))*(P[1] + cos(P[5]));
	u_PQW[2] = 0.0;
	// ------------------------ //

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				temp_matrix[i][j] += R1_tr[i][k]*R2_tr[k][j];
			}
		}
	}

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			for(k=0; k<3; k++) {
				T_PQW_to_ECI[i][j] += R3_tr[i][k]*temp_matrix[k][j];
			}
		}
	}

	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			r_ECI[i] += T_PQW_to_ECI[i][j]*r_PQW[j];
			u_ECI[i] += T_PQW_to_ECI[i][j]*u_PQW[j];
		}
	}

	for(i=0; i<3; i++) {
		X[i] = r_ECI[i];
		X[i+3] = u_ECI[i];
	}
}

