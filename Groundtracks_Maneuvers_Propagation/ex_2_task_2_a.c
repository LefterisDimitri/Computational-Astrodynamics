#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Draw the Keplerian ellipse and compute the distance from Moon's surface 

double h_rk, t_rk;
double GM_moon;

void orbital_elements_to_state_vec(double *P, double *X);    // apo stoixeia troxias se theseis kai taxuthtes
void RK(double* X, double* Xn);                              // Runge-Kutta 4
void derivatives(double* Y, double* dotX);                   // eksiswseis kinhshs
void draw_the_ellipses(double *P, double *X, FILE *fp1);     

int main()
{
	FILE *fp;
	FILE *fp1;
	fp = fopen("R_s_vs_t.dat","w");
	fp1 = fopen("draw_the_ellipse.dat","w");
	
	double R_moon = 1738.1;
	GM_moon = 4902.8;
	
	double R_s;

	double X[6];
	double Xn[6];

	X[0] = 0.0;
	X[1] = 0.0;
	X[2] = 0.0;
	X[3] = 0.0;
	X[4] = 0.0;
	X[5] = 0.0;
	
	Xn[0] = 0.0;
	Xn[1] = 0.0;
	Xn[2] = 0.0;
	Xn[3] = 0.0;
	Xn[4] = 0.0;
	Xn[5] = 0.0;
	
	double Xi[6];
	double Pi[6];
	
	Xi[0] = 0.0;
	Xi[1] = 0.0;
	Xi[2] = 0.0;
	Xi[3] = 0.0;
	Xi[4] = 0.0;
	Xi[5] = 0.0;
	
	Pi[0] = 5737.4;         // a
	Pi[1] = 0.61;           // e
	Pi[2] = 1.0091493735;   // i
	Pi[3] = 0.0;            // Omega
	Pi[4] = 1.5707963268;   // omega_tiny
	Pi[5] = 0.0;            // f
	
	double n = sqrt(GM_moon/(Pi[0]*Pi[0]*Pi[0]));
	double T = 2.0*M_PI/n;
	
	
	orbital_elements_to_state_vec(Pi,Xi);
	draw_the_ellipses(Pi,Xi,fp1);

	R_s = sqrt(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2]);
	fprintf(fp,"%lf %.15lf \n", 0.0, R_s);
	
	int m, p, max;
    double r_max = 1.0;
    double t_max;
	h_rk = 1.0;                   // time step dt = 1 sec
	max = lround(T/h_rk);
	for(m=1; m <= max; m++) {
		t_rk = m*h_rk;
		RK(Xi, Xn);
		for(p=0; p<6; p++) {
			Xi[p] = Xn[p];
		}
		
		R_s = sqrt(Xi[0]*Xi[0] + Xi[1]*Xi[1] + Xi[2]*Xi[2]);   // distance from Moon's center
		fprintf(fp,"%lf %.15lf \n", t_rk, R_s-R_moon);         // print the distance from Moon's surface
		
		if(R_s - R_moon > r_max)
		{
			r_max = R_s - R_moon;
			t_max = t_rk;
		}
		
	}
	
	// print the max distance as dmax = a(1+e) - R_moon
	printf("The max distance is (theoretically) d_max = a*(1 + e) - R_moon = %lf \nand (computationally) r_max = %lf\n\n", Pi[0]+Pi[0]*Pi[1] - R_moon, r_max); 
	printf("For the max distance we have t = %lf sec\n\n", t_max);
	
	// ---------- Verification of results via Kepler's equation and Newton-Raphson's method ---------- //
	
	double f, E_end, E_start, M_start, error;
	
	M_start = n*t_max;       // input = t_max , in Kepler's equation in order to find the true anomaly at the time of the max distance from Moon's surface
	E_start = M_start;
	while(error > 0.0001) {
		E_end = E_start - (E_start - Pi[1]*sin(E_end) - M_start)/(1.0-Pi[1]*cos(E_start));
		error = 100.0*abs((E_end-E_start)/E_end);
		E_start = E_end;
	}
	
	if(tan(E_end/2.0) < 0.0) {
		f = 2.0*(M_PI - atan(-sqrt((1.0+Pi[1])/(1.0-Pi[1]))*tan(E_start/2.0)));
	} else {
		f = 2.0*atan(sqrt((1.0+Pi[1])/(1.0-Pi[1]))*tan(E_start/2.0));
	}
	
	double r_max_exp;
	
	r_max_exp = Pi[0]*(1.0 - Pi[1]*Pi[1])/(1.0 + Pi[1]*cos(f));                   // the max distance from the Moon's center
	
	printf("VERIFICATION : The max distance is d = %lf\n", r_max_exp - R_moon);   // print the max distance from Moon's surface
	
	// ----------------------------------------------------------------------------------------------- //
	
	fclose(fp);
	fclose(fp1);
	
	return 0;
}


void orbital_elements_to_state_vec(double *P, double *X)
{
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
	
	u_PQW[0] = sqrt(GM_moon/(P[0]*(1.0 - P[1]*P[1])))*(-sin(P[5]));
	u_PQW[1] = sqrt(GM_moon/(P[0]*(1.0 - P[1]*P[1])))*(P[1] + cos(P[5]));
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
	
	for(i=0; i<3; i++)
	{
		X[i] = r_ECI[i];
		X[i+3] = u_ECI[i];
	}
	
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
		k1[l] = h_rk*dotX[l];
	}

	for(l=0; l<6; l++) {
		Y[l] = X[l] + 0.5*k1[l];
	}

	for(l=0; l<6; l++) {
		derivatives(Y, dotX);
		k2[l] = h_rk*dotX[l];
	}

	for(l=0; l<6; l++) {
		Y[l] = X[l] + 0.5*k2[l];
	}

	for(l=0; l<6; l++) {
		derivatives(Y, dotX);
		k3[l] = h_rk*dotX[l];
	}

	for(l=0; l<6; l++) {
		Y[l] = X[l] + k3[l];
	}

	for(l=0; l<6; l++) {
		derivatives(Y, dotX);
		k4[l] = h_rk*dotX[l];
	}

	for(l=0; l<6; l++) {
		Xn[l] = X[l] + (k1[l] + 2.0*k2[l] + 2.0*k3[l] + k4[l])/6.0;
	}

}

void derivatives(double *Y, double *dotX) 
{
	double rad,rad3;
	rad=sqrt(Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2]);
	rad3=rad*rad*rad;

	dotX[3] = -GM_moon*Y[0]/rad3;    // ux_dot
	dotX[4] = -GM_moon*Y[1]/rad3;    // uy_dot
	dotX[5] = -GM_moon*Y[2]/rad3;    // uz_dot
	dotX[0] = Y[3];               // ux
	dotX[1] = Y[4];               // uy
	dotX[2] = Y[5];               // uz
}

void draw_the_ellipses(double *P, double *X, FILE *fp1)
{
	int i_dte;
	double X_dte[6];
	double P_dte[6];
	for(i_dte=0; i_dte<6; i_dte++)
	{
		X_dte[i_dte] = 0.0;
		P_dte[i_dte] = 0.0;
	}
	
	for(i_dte=0; i_dte<6; i_dte++)
	{
		X_dte[i_dte] = X[i_dte];
		P_dte[i_dte] = P[i_dte];
	}
		
	double f_start, f_end, E_start, E_end, M_start, M_end, error, E_end_new;
	
	double t_dte, T_dte, n_dte;
	
	n_dte = sqrt(GM_moon/(P_dte[0]*P_dte[0]*P_dte[0]));
	T_dte = 2.0*M_PI/n_dte;
		
	f_start = P_dte[5];

	if(tan(f_start/2.0) >= 0.0)
	{
		E_start = 2.0*atan(sqrt((1.0-P_dte[1])/(1.0+P_dte[1]))*tan(f_start/2.0));
	}
	else
	{
		E_start = 2.0*(M_PI - atan(-sqrt((1.0-P_dte[1])/(1.0+P_dte[1]))*tan(f_start/2.0)));
	}
	
	M_start = E_start - P_dte[1]*sin(E_start);
	
	for(t_dte=0.0; t_dte<=T_dte; t_dte = t_dte + 1.0)
	{
		M_end = M_start + n_dte*t_dte;
		
		E_end = M_end;

		while(error > 0.0001) {
			E_end_new = E_end - (E_end - P_dte[1]*sin(E_end) - M_end)/(1.0-P_dte[1]*cos(E_end));
			error = 100.0*abs((E_end_new-E_end)/E_end_new);
			E_end = E_end_new;
		}
		
		if(tan(E_end/2.0) < 0.0) {
			f_end = 2.0*(M_PI - atan(-sqrt((1.0+P_dte[1])/(1.0-P_dte[1]))*tan(E_end/2.0)));
		} else {
			f_end = 2.0*atan(sqrt((1.0+P_dte[1])/(1.0-P_dte[1]))*tan(E_end/2.0));
		}
		
		P_dte[5] = f_end;
		orbital_elements_to_state_vec(P_dte, X_dte);
		fprintf(fp1,"%lf %lf %lf \n", X_dte[0], X_dte[1], X_dte[2]);
		M_start = M_end;
	}
}
