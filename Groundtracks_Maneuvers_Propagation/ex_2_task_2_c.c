#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Groundtracks 

double GM;

//void orbital_elements_to_state_vec(double *P, double *X);       

int main() {
	FILE* fp;
	fp = fopen("groundtracks3.dat", "w");

	double t;
	int i, j, k;
	double p, r, n, theta, lon, lat;
	double E_start, E_end, M_start, M_end, f_end, E_end_new, error;
	double P[6];
	double r_PQW[3];
	double u_PQW[3];
	double r_ECI[3];
	double u_ECI[3];
	double r_ECEF[3];
	double Rz[3][3];
	double temp_matrix[3][3];
	double T_PQW_to_ECI[3][3];


	GM = 4902.8;

	P[0] = 5737.4;          // a
	P[1] = 0.61;            // e
	P[2] = 1.0091493735;    // i
	P[3] = 0.0;             // Omega
	P[4] = 1.5707963268;    // omega_tiny
	P[5] = 0.0;             // f

	n = sqrt(GM/(P[0]*P[0]*P[0]));
	p = P[0]*(1.0 - P[1]*P[1]);
	r = p/(1.0 + P[1]*cos(P[5]));

	E_start = 2.0*atan(sqrt((1.0-P[1])/(1.0+P[1]))*tan(P[5]/2.0));
	M_start = E_start;

	r_PQW[0] = r*cos(P[5]);
	r_PQW[1] = r*sin(P[5]);
	r_PQW[2] = 0.0;
	
	/*double Xi[6];
	Xi[0] = 0.0;
	Xi[1] = 0.0;
	Xi[2] = 0.0;
	Xi[3] = 0.0;
	Xi[4] = 0.0;
	Xi[5] = 0.0;*/

	//printf("%lf\n%lf\n%lf\n\n\n", r_PQW[0], r_PQW[1], r_PQW[2]);

	/*u_PQW[0] = sqrt(GM_E/p)*(-sin(P[5]));
	u_PQW[1] = sqrt(GM_E/(P[0]*(1.0 - P[1]*P[1])))*(cos(P[5]));
	u_PQW[2] = 0.0; */

	/*double R3[3][3] = {
	{cos(P[3]), sin(P[3]), 0.0},
	{-sin(P[3]), cos(P[3]), 0.0},
	{0.0, 0.0, 1.0}
	};

	double R1[3][3] = {
	{1.0, 0.0, 0.0},
	{0.0, cos(P[2]), sin(P[2])},
	{0.0, -sin(P[2]), cos(P[2])}
	};

	double R2[3][3] = {
	{cos(P[4]), sin(P[4]), 0.0},
	{-sin(P[4]), cos(P[4]), 0.0},
	{0.0, 0.0, 1.0}
	};*/

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
			Rz[i][j] = 0.0;
		}
		r_ECI[i] = 0.0;
		r_ECEF[i] = 0.0;
	}

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
			//u_ECI[i] = T_PQW_to_ECI[i][j]*u_PQW[j];
		}
	}
		
	
	//orbital_elements_to_state_vec(P,Xi);
	
	for(t=0.0; t<=259200.0; t = t + 0.1) {
		/*for(i=0; i<3; i++)
		{
			r_ECEF[i] = 0.0;
		}*/
		//printf("%lf \n", t);

		theta = 0.00000266186*t;

		Rz[0][0] = cos(theta);
		Rz[0][1] = sin(theta);
		Rz[0][2] = 0.0;
		Rz[1][0] = -sin(theta);
		Rz[1][1] = cos(theta);
		Rz[1][2] = 0.0;
		Rz[2][0] = 0.0;
		Rz[2][1] = 0.0;
		Rz[2][2] = 1.0;

		for(i=0; i<3; i++) {
			r_ECEF[i] = 0.0;
			for(j=0; j<3; j++) {
				r_ECEF[i] += Rz[i][j]*r_ECI[j];
			}
		}

		lat = atan2(r_ECEF[2],sqrt(r_ECEF[0]*r_ECEF[0] + r_ECEF[1]*r_ECEF[1]));
		lon = atan2(r_ECEF[1],r_ECEF[0]);

		/*if(r_ECEF[1]/r_ECEF[0] < 0.0) {
			lon = M_PI - atan(-r_ECEF[1]/r_ECEF[0]);
		} else {
			lon = atan(r_ECEF[1]/r_ECEF[0]);
		}
		
		if(r_ECEF[2]/sqrt(r_ECEF[0]*r_ECEF[0] + r_ECEF[1]*r_ECEF[1]) < 0.0)
		{
			lat = M_PI - atan(-r_ECEF[2]/sqrt(r_ECEF[0]*r_ECEF[0] + r_ECEF[1]*r_ECEF[1]));
		}
		else
		{
			lat = atan(r_ECEF[2]/sqrt(r_ECEF[0]*r_ECEF[0] + r_ECEF[1]*r_ECEF[1]));
		}*/


		fprintf(fp,"%lf %lf \n", lon, lat);


		M_end = n*t;

		E_end = M_end;

		while(error > 0.0001) {
			E_end_new = E_end - (E_end - P[1]*sin(E_end) - M_end)/(1.0-P[1]*cos(E_end));
			error = 100.0*abs((E_end_new-E_end)/E_end_new);
			E_end = E_end_new;
		}
		
		

		if(tan(E_end) < 0.0) {
			f_end = 2.0*(M_PI - atan(-sqrt((1.0+P[1])/(1.0-P[1]))*tan(E_end/2.0)));
		} else {
			f_end = 2.0*atan(sqrt((1.0+P[1])/(1.0-P[1]))*tan(E_end/2.0));
		}
		
		//f_end = 2.0*atan(sqrt((1.0+P[1])/(1.0-P[1]))*tan(E_end/2.0)); // ok both
		
			
		P[5] = f_end;
		
		r = p/(1.0 + P[1]*cos(P[5]));

		r_PQW[0] = r*cos(P[5]);
		r_PQW[1] = r*sin(P[5]);
		r_PQW[2] = 0.0;

		/*for(i=0; i<3; i++) {
			for(j=0; j<3; j++) {
				temp_matrix[i][j] = 0.0;
				for(k=0; k<3; k++) {
					temp_matrix[i][j] += R1_tr[i][k]*R2_tr[k][j];
				}
			}
		}

		for(i=0; i<3; i++) {
			for(j=0; j<3; j++) {
				T_PQW_to_ECI[i][j] = 0.0;
				for(k=0; k<3; k++) {
					T_PQW_to_ECI[i][j] += R3_tr[i][k]*temp_matrix[k][j];
				}
			}
		}*/

		for(i=0; i<3; i++) {
			r_ECI[i] = 0.0;
			for(j=0; j<3; j++) {
				r_ECI[i] += T_PQW_to_ECI[i][j]*r_PQW[j];
			}
		}
	}

	fclose(fp);


	return 0;
}

/*void orbital_elements_to_state_vec(double *P, double *X)
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
	
	u_PQW[0] = sqrt(GM/(P[0]*(1.0 - P[1]*P[1])))*(-sin(P[5]));
	u_PQW[1] = sqrt(GM/(P[0]*(1.0 - P[1]*P[1])))*(P[1] + cos(P[5]));
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
	
} */

