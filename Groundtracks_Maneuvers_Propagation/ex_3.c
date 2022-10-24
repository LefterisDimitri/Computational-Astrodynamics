#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Earth satellite maneuvering 

double GM_E = 398600.433;
double h_rk;


void orbital_elements_to_state_vec(double *P, double *X);               // stoixeia troxias ---> theseis kai taxuthtes
void draw_the_ellipses(int orbit, double *P, double *X, FILE *fp);      // zwgrafikh ths troxias
void state_vec_to_orbital_elements(double *X, double *P);               // theseis kai taxuthtes ---> stoixeia troxias

void RK(double *X, double *Xn);                                         // Runge-Kutta 4th order
void derivatives(double *Y, double *dotX);                              // oi eksiswseis kinhshs


void for_RK(double *TOF, double theta, double *X, double *P, FILE *fp3) // oloklhrwsh twn eksiswsewn kinhshs apo thn arxikh sthn telikh troxia
{
	//printf("\nf = %.10lf \n", theta);
	int i_temp;
	double X_temp[6];
	double Xn_temp[6];
	
	for(i_temp=0; i_temp<6; i_temp++)
	{
		X_temp[i_temp] = 0.0;
		Xn_temp[i_temp] = 0.0;
	}
	
	for(i_temp=0; i_temp<6; i_temp++)
	{
		X_temp[i_temp] = X[i_temp];     
	}
	
	int m, max, p_temp;
	double r_temp, u_temp, ur_temp, ex_temp, ey_temp, ez_temp, e_temp;
	double T_temp, n_temp, t_temp;
	
	n_temp = sqrt(GM_E/(P[0]*P[0]*P[0]));
	T_temp = 2.0*M_PI/n_temp;
	
	h_rk = 0.1;                                                                  // time step // best precision at h = 0.0001 
	max = lround(T_temp/h_rk);
	for(m = 0; m <= max; m++)
	{
		t_temp = m*h_rk;
		RK(X_temp,Xn_temp);                                                       // RK4
		for(p_temp=0; p_temp<6; p_temp++) 
		{
			X_temp[p_temp] = Xn_temp[p_temp];
		}
		
		fprintf(fp3,"%.10lf %.10lf %.10lf \n", X_temp[0], X_temp[1], X_temp[2]);
		
		r_temp = sqrt(X_temp[0]*X_temp[0] + X_temp[1]*X_temp[1] + X_temp[2]*X_temp[2]);
		u_temp = sqrt(X_temp[3]*X_temp[3] + X_temp[4]*X_temp[4] + X_temp[5]*X_temp[5]);
		
		ur_temp = (X_temp[0]*X_temp[3] + X_temp[1]*X_temp[4] + X_temp[2]*X_temp[5])/r_temp;
		
		ex_temp = (1.0/GM_E)*(u_temp*u_temp - GM_E/r_temp)*X_temp[0] - (X_temp[0]*X_temp[3] + X_temp[1]*X_temp[4] + X_temp[2]*X_temp[5])*X_temp[3]/GM_E;
	    ey_temp = (1.0/GM_E)*(u_temp*u_temp - GM_E/r_temp)*X_temp[1] - (X_temp[0]*X_temp[3] + X_temp[1]*X_temp[4] + X_temp[2]*X_temp[5])*X_temp[4]/GM_E;
	    ez_temp = (1.0/GM_E)*(u_temp*u_temp - GM_E/r_temp)*X_temp[2] - (X_temp[0]*X_temp[3] + X_temp[1]*X_temp[4] + X_temp[2]*X_temp[5])*X_temp[5]/GM_E;
	    e_temp = sqrt(ex_temp*ex_temp + ey_temp*ey_temp + ez_temp*ez_temp);
	    
	    // ---------- true anomaly ---------- //
		if(ur_temp >= 0.0)
	    {
		    P[5] = acos((ex_temp*X_temp[0] + ey_temp*X_temp[1] + ez_temp*X_temp[2])/(e_temp*r_temp));
	    }
	    else
	    {
		    P[5] = 2.0*M_PI - acos((ex_temp*X_temp[0] + ey_temp*X_temp[1] + ez_temp*X_temp[2])/(e_temp*r_temp));
	    }
	    // ---------------------------------- //
	    
	    
	    if(fabs(theta - P[5]) <= 0.0001)    // best precision at epsilon = 0.0000001
	    {
	    	printf(" ----------------------------- ORBITAL CHANGE ----------------------------\n");
	    	printf("theta input = %.10lf and output = %.10lf\n", theta, P[5]);
	    	printf("\nsugklish: diafora = %.10lf \n", fabs(theta - P[5]));
	    	printf(" -------------------------------------------------------------------------\n");
	    	//m = max;
	    	TOF[0] = TOF[0] + t_temp;
	    	break;
		}
	}
}
	
int main()
{	
	FILE *fp;
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	FILE *fp5;

	fp = fopen("maneuvers_final_orbit.dat", "w");
	fp1 = fopen("maneuvers_orbit_1_initial.dat", "w");
	fp2 = fopen("maneuvers_orbit_2.dat", "w");
	fp3 = fopen("maneuvers_orbit_3.dat", "w");
	fp4 = fopen("maneuvers_transfer_orbit.dat", "w");
	fp5 = fopen("orbital_propagation.dat", "w");

	
	int p_i;
	double X[6];       // state vector
	double P[6];       // elements vector
	
	double Pi[6];      // state vector for t = 0
	double Xi[6];      // elements vector for t = 0
	
	double Pf[6];      // state vector for final orbit
	double Xf[6];      // elements vector for final orbit
	
	for(p_i=0; p_i<6; p_i++)
	{
		X[p_i] = 0.0;
		P[p_i] = 0.0;
		Xi[p_i] = 0.0;
		Pi[p_i] = 0.0;
		Xf[p_i] = 0.0;
		Pf[p_i] = 0.0;
	}
	
	
	// ----- time of flight ----- //
	double TOF[1];
	TOF[0] = 0.0;
	// -------------------------- //
	
	
	// ---------- ECI ---------- //
	Xi[0] = -9141.878;                  //   x(t=0)
	Xi[1] = -1648.0758;                 //   y(t=0)
	Xi[2] = 4141.679;                   //   z(t=0)
	Xi[3] = -1.153;                     //   ux_t=))
	Xi[4] = -5.31;                      //   uy(t=0)
	Xi[5] = -2.898;                     //   uz(t=0)
	// ------------------------ //
	
	
	// ---------- final orbital elements ---------- //
	Pf[0] = 12940.0;                                      // a
	Pf[1] = 0.2173;                                       // e
	Pf[2] = 0.8692;                                       // i
	Pf[3] = 1.448;                                        // Omega
	Pf[4] = 2.721;                                        // omega_tiny
	Pf[5] = 2.827;                                        // f
	// -------------------------------------------- //
	
	
	// ---------- initial orbit ---------- //
	
	state_vec_to_orbital_elements(Xi,Pi);              // ECI state vector --> orbital elements for t = 0 
	draw_the_ellipses(1,Pi,Xi,fp1);                    // draw the 1st orbit (initial)
	
	printf("---------- INITIAL ORBIT ----------\n");
	for(p_i = 0; p_i<6; p_i++)
	{
		printf("Pi[%d] = %lf \n", p_i, Pi[p_i]);      // print the initial orbital elements (t=0)
	}
	printf("-----------------------------------\n");
	
	// ----------------------------------- //
	
	printf("\n");	
		
	// ---------- final orbit ---------- //
	
	orbital_elements_to_state_vec(Pf,Xf);               // orbital elements --> ECI state vector (final orbit)
	draw_the_ellipses(5,Pf,Xf,fp);	                    // drawing the final orbit
	
	printf("---------- FINAL ORBIT ----------\n");
	for(p_i = 0; p_i<6; p_i++)
	{
		printf("Xf[%d] = %lf \n", p_i, Xf[p_i]);      // print the state vector of final orbit
	}
	printf("-----------------------------------\n");
	
	// --------------------------------- //
			
		
	// ---------- change orbital plane ---------- //
	
	double D_Omega, D_i;                                           
	double a, cosa, sinu1, sinu2;
	double u1, u2, theta1, theta2, omega2;
	double u_theta, Du1;
	
	D_Omega = Pf[3] - Pi[3];                                            // Omega_2 - Omega_1 > 0                                           
	D_i = Pf[2] - Pi[2];                                                // i_2 - i_1 > 0                                              
		
	cosa = cos(Pi[2])*cos(Pf[2]) + sin(Pi[2])*sin(Pf[2])*cos(D_Omega);  // cosa > 0
	
	a = acos(cosa);                                                     
	printf("a=%lf", a);
	                  	
	sinu1 = sin(Pf[2])*sin(D_Omega)/sin(a);                             // sinu1 > 0
	sinu2 = sin(Pi[2])*sin(D_Omega)/sin(a);                             // sinu2 > 0
	
	u1 = asin(sinu1);                                                 
	u2 = asin(sinu2);                                                   

	theta1 = u1 - Pi[4];                                                // theta1 = u1 - omega_tiny_1
	theta2 = theta1;
	omega2 = u2 - u1 + Pi[4]; 
	
	printf("\ntheoretical : theta1 = %.10lf\n", theta1);
	
	// ---------- PROPAGATION ----------- //
	for_RK(TOF,theta1,Xi,Pi,fp5);                                    // theta_1 convergence via integration
	// ---------------------------------- //         
	                                             	
	printf("\napproximation : theta1 = %.10lf \n", Pi[5]);
	printf("1st propagate TOF = %lf \n\n", TOF[0]);
	
	u_theta = sqrt(GM_E/(Pi[0]*(1.0-Pi[1]*Pi[1])))*(1.0 + Pi[1]*cos(theta1));
	Du1 = 2.0*u_theta*sin(a/2.0);                                                 // Du : plane change
	
	// ---------- orbit 2 ---------- //
	
	double omega3, D_omega;
	double theta2_a, theta3_a;
	double Du2;
	
	double P2[6];
	double X2[6];
	
	for(p_i=0; p_i<6; p_i++)
	{
		X2[p_i] = 0.0;
		P2[p_i] = 0.0;
	}
	
	P2[0] = Pi[0];
	P2[1] = Pi[1];
	P2[2] = Pf[2];
	P2[3] = Pf[3];
	P2[4] = omega2;
	P2[5] = theta2;
	
	orbital_elements_to_state_vec(P2,X2);
	draw_the_ellipses(2,P2,X2,fp2);                                //draw the 2nd orbit
	
	omega3 = Pf[4];                                                // omega_tiny_3 = omega_tiny_final
	D_omega = omega3 - omega2;
	
	theta2_a = D_omega/2.0;
	theta3_a = 2.0*M_PI - D_omega/2.0;
	
	printf("\ntheoretical : theta2_a = %.10lf \n", theta2_a);
	printf("\ntheoretical : theta3_a = %.10lf \n", theta3_a);
	
	// ---------- PROPAGATION ----------- //
	for_RK(TOF,theta2_a,X2,P2,fp5);                              // theta_2_a convergence via integration
	// ---------------------------------- //  
	              
	printf("2nd propagate TOF = %lf \n", TOF[0]);
	printf("\napproximation : theta2_a_ = %.10lf \n\n", P2[5]);
		
	// ---------- end of orbit 2 ---------- //
	
	
	// ---------- orbit 3 ---------- //
	double P3[6];
	double X3[6];
	
	for(p_i=0; p_i<6; p_i++)
	{
		X3[p_i] = 0.0;
		P3[p_i] = 0.0;
	}
		
	P3[0] = Pi[0];
	P3[1] = Pi[1];
	P3[2] = Pf[2];
	P3[3] = Pf[3];
	P3[4] = omega3;                                                              //omega_tiny_3 == Pf[4] == omega_tiny_final
	P3[5] = theta3_a;
	
	Du2 = 2.0*sqrt(GM_E/(P3[0]*(1.0-P3[1]*P3[1])))*P3[1]*sin(D_omega/2.0);       // Du : change periapsis
	
	orbital_elements_to_state_vec(P3,X3);
	draw_the_ellipses(3,P3,X3,fp3);                                               // draw the 3rd orbit   
	
	// ---------- PROPAGATION ----------- //
	for_RK(TOF,0.0,X3,P3,fp5);                                                  // 0.0 convergence via integration
	// ---------------------------------- //                              
	
	printf("3rd propagate TOF = %lf \n", TOF[0]);
	printf("\napproximation of theta = 0 at 1st transfer point :  %.10lf \n", P3[5]);
	
	// ---------- end of orbit 3 ---------- //
	
	// ---------- transfer orbit ---------- //
	double P_transfer[6];	
	double X_transfer[6];
	for(p_i=0; p_i<6; p_i++)
	{
		X_transfer[p_i] = 0.0;
		P_transfer[p_i] = 0.0; 
	}
	
	P_transfer[0] = (P3[0]*(1.0-P3[1]) + Pf[0]*(1.0+Pf[1]))/2.0;                                          
	P_transfer[1] = (Pf[0]*(1.0+Pf[1]) - P3[0]*(1.0-P3[1]))/(Pf[0]*(1.0+Pf[1]) + P3[0]*(1.0-P3[1]));      
	P_transfer[2] = Pf[2];
	P_transfer[3] = Pf[3];
	P_transfer[4] = Pf[4];
	P_transfer[5] = 0.0; 
	
	printf("for transfer orbit : a = %lf \n", P_transfer[0]);
	printf("for transfer orbit : e = %lf \n", P_transfer[1]);
	
	// ---------- for Dv ---------- //
	double Du1_tr, Du2_tr, uP3, uPt, uAf, uAt, Du_total, Dv_total;
	
	uPt = sqrt((GM_E/P_transfer[0])*(1.0+P_transfer[1])/(1.0-P_transfer[1]));
	uAt = sqrt((GM_E/P_transfer[0])*(1.0-P_transfer[1])/(1.0+P_transfer[1]));
	
	uP3 = sqrt((GM_E/P3[0])*(1.0+P3[1])/(1.0-P3[1]));
	uAf = sqrt((GM_E/Pf[0])*(1.0-Pf[1])/(1.0+Pf[1]));
	
	Du1_tr = uPt - uP3;
	Du2_tr = uAf - uAt;
	
	Du_total = Du1_tr + Du2_tr;                                                   // Du : total impulse
	
	orbital_elements_to_state_vec(P_transfer,X_transfer);
	draw_the_ellipses(4,P_transfer,X_transfer,fp4);                               // draw the transfer orbit
	
	// ---------- PROPAGATION ----------- //
	for_RK(TOF,M_PI,X_transfer,P_transfer,fp5);                                   // pi convergence via integration
	// ---------------------------------- //
	
	printf("4th propagate TOF = %lf \n", TOF[0]);
	printf("\napproximation of theta = pi at 2nd transfer point = %.10lf \n", P_transfer[5]);
	
	// ---------- end of transfer orbit ---------- //
	
	// ---------- final orbit ----------//
	Pf[0] = Pf[0];
	Pf[1] = Pf[1];
	Pf[2] = Pf[2];
	Pf[3] = Pf[3];
	Pf[4] = Pf[4];
	Pf[5] = M_PI;
	
	orbital_elements_to_state_vec(Pf,Xf);
	
	// ---------- PROPAGATION ----------- //
	for_RK(TOF,2.827,Xf,Pf,fp5);                                                   // theta_final convergence via integration
	// ---------------------------------- //
	
	printf("final propagate TOF = %lf \n", TOF[0]);
	printf("\napproximation : theta_final = %.10lf \n", Pf[5]);
	printf("\nORBIT NO 5 (final ellipse) period: T = %lf sec\n", 2.0*M_PI/sqrt(GM_E/(Pf[0]*Pf[0]*Pf[0])));
	
	// ---------- end of final orbit ---------- //
	Dv_total = Du1 + Du2 + Du_total;
	
	printf("\n\n---------- RESULTS ----------\n\n");
	printf("Time Of Flight (TOF) = %lf \n", TOF[0]);
	printf("Total cost : Dv = %lf\n", Dv_total);
		
	return 0;
}

void state_vec_to_orbital_elements(double *X, double *P)
{
	double r_svtoe, u_svtoe, E_svtoe;
	double hx, hy, hz, h, p_svtoe, N, nx, ny, nz, ex, ey, ez, e, ur;
	
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
	P[1] = e;                                             // sqrt(1.0 - p_svtoe/P[0]);          
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
		P[4] = acos((nx*ex + ny*ey + nz*ez)/(N*e));
	}
	else
	{
		P[4] = 2.0*M_PI - acos((nx*ex + ny*ey + nz*ez)/(N*e));
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
	
	for(i=0; i<3; i++)
	{
		X[i] = r_ECI[i];
		X[i+3] = u_ECI[i];
	}
	
}

void draw_the_ellipses(int orbit, double *P, double *X, FILE *fp)
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
	
	n_dte = sqrt(GM_E/(P_dte[0]*P_dte[0]*P_dte[0]));
	T_dte = 2.0*M_PI/n_dte;
	
	if(orbit == 1)
	{
		printf("ORBIT NO 1 (initial) period: T = %lf sec\n\n", T_dte);
	}
	else if(orbit == 2)
	{
		printf("ORBIT NO 2 period: T = %lf sec\n\n", T_dte);
	}
	else if(orbit == 3)
	{
		printf("ORBIT NO 3 period: T = %lf sec\n\n", T_dte);
	}
	else if(orbit == 4)
	{
		printf("ORBIT NO 4 (transfer ellipse) period: T = %lf sec\n\n", T_dte);
	}
	else if(orbit == 5)
	{
		printf("ORBIT NO 5 (final ellipse) period: T = %lf sec\n\n", T_dte);
	}
		
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
		fprintf(fp,"%lf %lf %lf \n", X_dte[0], X_dte[1], X_dte[2]);
		M_start = M_end;
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

	dotX[3] = -GM_E*Y[0]/rad3;    // ux_dot
	dotX[4] = -GM_E*Y[1]/rad3;    // uy_dot
	dotX[5] = -GM_E*Y[2]/rad3;    // uz_dot
	dotX[0] = Y[3];               // ux
	dotX[1] = Y[4];               // uy
	dotX[2] = Y[5];               // uz
}
