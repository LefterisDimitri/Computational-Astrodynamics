#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// ---------- Time in Earth's shadow ---------- //

int main()
{
	double GM_E, R_E, a, e, n, x, r, cosf, f1, f2, E1, E2, M1, M2, t;
	
	GM_E = 398600.433;
	R_E = 6378.0;
	a = 32000.0;
	e = 0.2;
	n = sqrt(GM_E/(a*a*a));	
	
	x = a*sqrt(1.0 - R_E*R_E/(a*a*(1.0 - e*e)));
	
	printf("x = %lf \n", x);
	
	r = sqrt((x + a*e)*(x + a*e) + R_E*R_E);
	
	printf("r = %lf \n", r);
	
	cosf = (1.0/e)*(a*(1.0 - e*e)/r - 1.0);
	
	printf("cosf = %lf \n", cosf);
	
	f1 = M_PI - acos(-cosf);                                        // cosf < 0
	f2 = 2.0*M_PI - f1;
	
	printf("f1 = %lf \n", f1);
	printf("f2 = %lf \n", f2);
	
	printf("tan(f1/2.0) = %lf \n", tan(f1/2.0));
	printf("tan(f2/2.0) = %lf \n", tan(f2/2.0));
	
	E1 = 2.0*atan(sqrt((1.0-e)/(1+e))*tan(f1/2.0));
	E2 = 2.0*(M_PI - atan(-sqrt((1.0-e)/(1+e))*tan(f2/2.0)));         //tan(f2/2) < 0
	
	printf("E1 = %lf \n", E1);
	printf("E2 = %lf \n", E2);
	
	M1 = E1 - e*sin(E1);
	M2 = E2 - e*sin(E2);
	
	printf("M1 = %lf \n", M1);
	printf("M2 = %lf \n", M2);
	
	t = (M2 - M1)/n;                                                  // [t] = sec	
	printf("Total time in Earth's shadow: Dt = %.15lf", t);
	
	return 0;
}
