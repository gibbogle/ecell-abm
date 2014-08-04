/* 
The principal change in Version3.0 is to update the conjugate 
gradient routine to CG_DESCENT 6.0. Note that the objective 
structure also contains a pointer to the array of free indices */

#include <math.h>
#include "asa_user.h" // needed by the program which calls asa_cg 

double a1, b1, *c1, *o1, a2, b2, *c2, *o2;
double r1, r2;

// prototypes for the function and gradient evaluation routines 
double myvalue
(
    asa_objective *asa
) ;

void mygrad
(
    asa_objective *asa
) ;

double myvalgrad
(
    asa_objective *asa
) ;

// We will want to call this from Fortran, passing all the required info about the two ellipsoids:
// For each ellipsoid:
//     double a
//     double b
//     double *centre
//     double *orient
// We will solve for s1 = x[0] and s2 = x[1]
// The centreline point P(s) has bounds: 
//     -(a^2-b^2)/a < x0 < (a^2-b^2)/a
// which need to be translated into equivalent s values (0 - 1):
//     (-(a^2-b^2)/a^2 + 1)/2 < s < ((a^2-b^2)/a + 1)/2
// =>
//     (1/2).b^2/a^2 <= s <= 1 - (1/2).b^2/a^2

void __declspec(dllexport) min_dist(double aval1, double bval1, double *centre1, double*orient1, 
	double aval2, double bval2, double *centre2, double*orient2, 
	double *s1, double *s2, double *rad1, double *rad2, double *d, int *res)
{
	double x[2], lo[2], hi[2];
	double work[1000];
	int iwork[1000];
	INT n;
	int i;
	asa_stat Stat; // structure with statistics (can be NULL) 

	a1 = aval1;
	b1 = bval1;
	c1 = centre1;
	o1 = orient1;
	a2 = aval2;
	b2 = bval2;
	c2 = centre2;
	o2 = orient2;

	n = 2;

	// bounds
	lo[0] = 0.5*b1*b1/(a1*a1);
	hi[0] = 1 - lo[0];
	lo[1] = 0.5*b2*b2/(a2*a2);
	hi[1] = 1 - lo[1];

	// initial guesses
	x[0] = *s1;
	x[1] = *s2;
	for (i=0; i<2; i++) {
		if (x[i] < lo[i]) x[i] = lo[i];
		if (x[i] > hi[i]) x[i] = hi[i];
	}
	*res = asa_cg (x, lo, hi, n, &Stat, NULL, NULL, 1.e-6, myvalue, mygrad, myvalgrad, NULL, NULL) ;
//	*res = asa_cg (x, lo, hi, n, &Stat, NULL, NULL, 1.e-6, myvalue, mygrad, myvalgrad, work, iwork) ;
	*s1 = x[0];
	*s2 = x[1];
	*rad1 = r1;
	*rad2 = r2;
	*d = Stat.f;
	if (!QUIET) printf("bounds: lo: %f hi: %f\n",lo[0],hi[0]);
	return;
}

double myvalue // evaluate the objective function 
(
    asa_objective *asa
)
{
    double f, *x ;
	double p1, p2, D, s1, s2, x0;
    INT n ;
	int i;
    x = asa->x ;
    n = asa->n ;
	s1 = x[0];
	s2 = x[1];
	D = 0;
	for (i=0; i<3; i++) {
		p1 = c1[i] + a1*(2*s1-1)*o1[i];
		p2 = c2[i] + a2*(2*s2-1)*o2[i];
		D += (p1-p2)*(p1-p2);
	}
	D = sqrt(D);
	// ellipsoid 1, x0: -a1 - a1
	x0 = a1*(2*s1-1);
	r1 = b1*sqrt(1-x0*x0/(a1*a1-b1*b1));
	// ellipsoid 2, x0: -a2 - a2
	x0 = a2*(2*s2-1);
	r2 = b2*sqrt(1-x0*x0/(a2*a2-b2*b2));
	f = D - r1 - r2;
//	printf("myvalue: D,r1,r2: %f %f %f\n",D,r1,r2);
	return f;
}

void mygrad // evaluate the gradient of the objective function 
(
    asa_objective *asa
)
{
    double *g, *x ;
	double p1[3], p2[3], dx[3], D, s1, s2, x0, dD_ds1, dD_ds2, dr1_ds1, dr2_ds2;
    INT n;
	int i;

    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
	s1 = x[0];
	s2 = x[1];
	D = 0;
	for (i=0; i<3; i++) {
		p1[i] = c1[i] + a1*(2*s1-1)*o1[i];
		p2[i] = c2[i] + a2*(2*s2-1)*o2[i];
		dx[i] = p1[i] - p2[i];
		D += dx[i]*dx[i];
	}
	D = sqrt(D);
	dD_ds1 = 0;
	dD_ds2 = 0;
	for (i=0; i<3; i++) {
//		dD_dx[0] = dx[0]/D;
//		dx_ds1[0] = o1[0];
//		dD_ds1 += dD_dx[0]*dx_ds1[0];
		dD_ds1 += (dx[i]/D)*2*a1*o1[i];
		dD_ds2 += -(dx[i]/D)*2*a2*o2[i];
	}
	// ellipsoid 1, x0: -a1 - a1
	x0 = a1*(2*s1-1);
	r1 = b1*sqrt(1-x0*x0/(a1*a1-b1*b1));
	dr1_ds1 = -2*b1*a1*a1*(2*s1-1)/(r1*(a1*a1-b1*b1));
	// ellipsoid 2, x0: -a2 - a2
	x0 = a2*(2*s2-1);
	r2 = b2*sqrt(1-x0*x0/(a2*a2-b2*b2));
	dr2_ds2 = -2*b2*a2*a2*(2*s2-1)/(r2*(a2*a2-b2*b2));
	g[0] = dD_ds1 - dr1_ds1;
	g[1] = dD_ds2 - dr2_ds2;
	return;
}

double myvalgrad // evaluate the value and gradient of the objective function 
(
    asa_objective *asa
)
{
    double *g, *x ;
	double p1[3], p2[3], dx[3], D, s1, s2, x0, dD_ds1, dD_ds2, dr1_ds1, dr2_ds2, f;
    INT n;
	int i;

    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
	s1 = x[0];
	s2 = x[1];
	D = 0;
	for (i=0; i<3; i++) {
		p1[i] = c1[i] + a1*(2*s1-1)*o1[i];
		p2[i] = c2[i] + a2*(2*s2-1)*o2[i];
		dx[i] = p1[i] - p2[i];
//		printf("dx: %d %f\n",i,dx[i]);
		D += dx[i]*dx[i];
	}
	D = sqrt(D);
	dD_ds1 = 0;
	dD_ds2 = 0;
	for (i=0; i<3; i++) {
//		dD_dx[0] = dx[0]/d;
//		dx_ds1[0] = o1[0];
//		dD_ds1 += dD_dx[0]*dx_ds1[0];
		dD_ds1 += (dx[i]/D)*2*a1*o1[i];
		dD_ds2 += -(dx[i]/D)*2*a2*o2[i];
	}
	// ellipsoid 1, x0: -a1 - a1
	x0 = a1*(2*s1-1);
	r1 = b1*sqrt(1-x0*x0/(a1*a1-b1*b1));
	dr1_ds1 = -2*b1*a1*a1*(2*s1-1)/(r1*(a1*a1-b1*b1));
	// ellipsoid 2, x0: -a2 - a2
	x0 = a2*(2*s2-1);
	r2 = b2*sqrt(1-x0*x0/(a2*a2-b2*b2));
	dr2_ds2 = -2*b2*a2*a2*(2*s2-1)/(r2*(a2*a2-b2*b2));
	g[0] = dD_ds1 - dr1_ds1;
	g[1] = dD_ds2 - dr2_ds2;
	f = D - r1 - r2;
	if (!QUIET) printf("\nmyvalgrad: s1, s2: %f %f\n  p2: %f %f %f\nD,r1,r2,f: %f %f %f %f\n",s1,s2,p2[0],p2[1],p2[2],D,r1,r2,f);
	return f;
}

/*
int main (void)
{
	double aval1, bval1, centre1[3], orient1[3];
	double aval2, bval2, centre2[3], orient2[3];
	double s1, s2, d;

	aval1 = 5;
	bval1 = 3;
	centre1[0] = 0;
	centre1[1] = 0;
	centre1[2] = 0;
	orient1[0] = 1;
	orient1[1] = 0;
	orient1[2] = 0;
	aval2 = 5;
	bval2 = 3;
	centre2[0] = 0;
	centre2[1] = 10;
	centre2[2] = 0;
	orient2[0] = 0;
	orient2[1] = 1;
	orient2[2] = 0;
	s1 = 0.5;
	s2 = 0.5;
	min_dist(aval1,bval1,centre1,orient1,aval2,bval2,centre2,orient2,&s1,&s2,&d);
	printf("\ns1, s2: %f %f  d: %f\n",s1,s2,d);
	return 0;
}


int main (void)
{
    double *x, *lo, *hi ;
    INT i, n ;

    // if you want to change parameter value, you need the following: 
    asacg_parm cgParm ;
    asa_parm asaParm ;

    // allocate arrays for problem solution and bounds 
    n = 100 ; // problem dimension 
    x  = (double *) malloc (n*sizeof (double)) ;
    lo = (double *) malloc (n*sizeof (double)) ;
    hi = (double *) malloc (n*sizeof (double)) ;
    for (i = 0; i < n; i++) lo [i] = (double) 0 ;
    for (i = 0; i < n; i++) hi [i] = (double) 1 ;

    // if you want to change parameter value, initialize strucs with default 
    asa_cg_default (&cgParm) ;
    asa_default (&asaParm) ;

    // if you want to change parameters, change them here: 
    cgParm.PrintParms = TRUE ;
    cgParm.PrintLevel = 0 ;
    asaParm.PrintParms = TRUE ;
    asaParm.PrintLevel = 0 ;

    // starting guess 
    for (i = 0; i < n; i++) x [i] = 1 ;

    // run the code 
    asa_cg (x, lo, hi, n, NULL, &cgParm, &asaParm,
                     1.e-8, myvalue, mygrad, myvalgrad, NULL, NULL) ;

    // if no change in parameters, you could replace Parm arguments by NULL
    for (i = 0; i < n; i++) x [i] = 1 ; // starting guess 
    asa_cg (x, lo, hi, n, NULL, NULL, NULL,
                     1.e-8, myvalue, mygrad, myvalgrad, NULL, NULL) ;

    // with some loss of efficiency, you could omit the valgrad routine 
    for (i = 0; i < n; i++) x [i] = 1 ; // starting guess 
    asa_cg (x, lo, hi, n, NULL, NULL, NULL, 1.e-8, myvalue, mygrad, NULL,
            NULL, NULL);

    free (x) ;
    free (lo) ;
    free (hi) ;
}


double myvalue // evaluate the objective function 
(
    asa_objective *asa
)
{
    double f, t, *x ;
    INT i, n ;
    x = asa->x ;
    n = asa->n ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        f += exp (x [i]) - t*x [i] ;
    }
    return (f) ;
}

void mygrad // evaluate the gradient of the objective function 
(
    asa_objective *asa
)
{
    double t, *g, *x ;
    INT i, n ;
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        g [i] = exp (x [i]) -  t ;
    }
    return ;
}

double myvalgrad // value and gradient of the objective function 
(
    asa_objective *asa
)
{
    double f, xi, t, *g, *x ;
    INT i, n ;
    x = asa->x ;
    g = asa->g ;
    n = asa->n ;
    f = 0 ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        xi = x [i] ;
        f += exp (xi) - t*xi ;
        g [i] = exp (xi) -  t ;
    }
    return (f) ;
}

*/