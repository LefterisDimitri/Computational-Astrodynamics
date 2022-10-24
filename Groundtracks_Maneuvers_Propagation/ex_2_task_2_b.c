#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Orbital elements after 30 days via numerical integration

double h_rk, t_rk;
double GM_moon;

void orbital_elements_to_state_vec(double *P, double *X);       // apo stoixeia troxias se theseis kai taxuthtes
void state_vec_to_orbital_elements(double *X, double *P);       // apo theseis kai taxuthtes se stoixeia troxias  

void RK(double* X, double* Xn);                                 // Runge-Kutta 4
void derivatives(double* Y, double* dotX);                      // eksiswseis kinhshs 

int main()
{
	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	FILE *fp5;
	fp = fopen("a.dat","w");
	fp1 = fopen("e.dat","w");
	fp2 = fopen("i.dat","w");
	fp3 = fopen("Omega.dat","w");
	fp4 = fopen("omega_tiny.dat","w");
	fp5 = fopen("f.dat","w");
	
	int i;
	double R_moon = 1738.1;	
	GM_moon = 4902.8;
	
	double R_s;

	double X[6];
	double P[6];
	
	double Xn[6];

	double Xi[6];
	double Pi[6];
	
	for(i=0; i<6; i++)
	{
		X[i] = 0.0;
		P[i] = 0.0;
		Xi[i] = 0.0;
		Pi[i] = 0.0;
	}
		
	Pi[0] = 5737.4;
	Pi[1] = 0.61;
	Pi[2] = 1.0091493735; 
	Pi[3] = 0.0;
	Pi[4] = 1.5707963268;
	Pi[5] = 0.0;
	
	double n = sqrt(GM_moon/(Pi[0]*Pi[0]*Pi[0]));
	double T = 2.0*M_PI/n;
	
	orbital_elements_to_state_vec(Pi,Xi);            
	
	int m, p, max;

	h_rk = 1.0;                  // time step = 1 sec
	max = 2592000;               // 30 days = 2592000 sec
	for(m=1; m <= max; m++) {
		t_rk = m*h_rk;
		RK(Xi, Xn);
		for(p=0; p<6; p++) {
			Xi[p] = Xn[p];
		}
		
		fprintf(fp,"%lf %.10lf\n", t_rk, Pi[0]);
		fprintf(fp1,"%lf %.10lf\n", t_rk, Pi[1]);
		fprintf(fp2,"%lf %.10lf\n", t_rk, Pi[2]);
		fprintf(fp3,"%lf %.10lf\n", t_rk, Pi[3]);
		fprintf(fp4,"%lf %.10lf\n", t_rk, Pi[4]);
		fprintf(fp5,"%lf %.10lf\n", t_rk, Pi[5]);
	}
	
	state_vec_to_orbital_elements(Xi,Pi);                  
	
	printf("after 30 days the orbital elemens are: \n");
	for(i=0; i<6; i++)
	{
		printf("\nPi[%d] = %.15lf \n", i, Pi[i]);
	}
	
	
	return 0;
}


void state_vec_to_orbital_elements(double *X, double *P)
{
	double r_svtoe, u_svtoe, E_svtoe;
	double hx, hy, hz, h, p_svtoe, N, nx, ny, nz, ex, ey, ez, e, ur;
	
	r_svtoe = sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
	u_svtoe = sqrt(X[3]*X[3] + X[4]*X[4] + X[5]*X[5]);
	E_svtoe = 0.5*(X[3]*X[3] + X[4]*X[4] + X[5]*X[5]) - GM_moon/r_svtoe;
		
	hx = X[1]*X[5] - X[2]*X[4];
	hy = X[2]*X[3] - X[0]*X[5];
	hz = X[0]*X[4] - X[1]*X[3];
	h = sqrt(hx*hx + hy*hy + hz*hz);
	p_svtoe = h*h/GM_moon;
	
	ex = (1.0/GM_moon)*(u_svtoe*u_svtoe - GM_moon/r_svtoe)*X[0] - (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])*X[3]/GM_moon;
	ey = (1.0/GM_moon)*(u_svtoe*u_svtoe - GM_moon/r_svtoe)*X[1] - (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])*X[4]/GM_moon;
	ez = (1.0/GM_moon)*(u_svtoe*u_svtoe - GM_moon/r_svtoe)*X[2] - (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])*X[5]/GM_moon;
	e = sqrt(ex*ex + ey*ey + ez*ez);
	
	
	nx = -hy;
	ny = hx;
	nz = 0.0;
	N = sqrt(nx*nx + ny*ny + nz*nz);
	
	ur = (X[0]*X[3] + X[1]*X[4] + X[2]*X[5])/r_svtoe;
	
	P[0] = - GM_moon/(2.0*E_svtoe);      
	P[1] = e;                                      //sqrt(1.0 - p/P[0]);          
	P[2] = acos(hz/h);
	
	if(ny >= 0.0)
	{
		P[3] = acos(nx/N);
	}
	else
	{
		P[3] = 2.0*M_PI - acos(nx/N);
	} 
	
	if(ez >= 0.0)
	{
		P[4] = acos((nx*ex + ny*ey)/(N*e));
	}
	else
	{
		P[4] = 2.0*M_PI - acos((nx*ex + ny*ey)/(N*e));
	}
	
	if(ur >= 0.0)
	{
		P[5] = acos((ex*X[0] + ey*X[1] + ez*X[2])/(e*r_svtoe));
	}
	else
	{
		P[5] = 2.0*M_PI - acos((ex*X[0] + ey*X[1] + ez*X[2])/(e*r_svtoe));
	}
	
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
