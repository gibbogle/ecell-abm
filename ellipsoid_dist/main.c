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
