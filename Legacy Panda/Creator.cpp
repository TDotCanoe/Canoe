#include "Creator.h"


int counter;

int bowpower;			// these shouldn't be global variables but I'm tired and it's 2:15 AM
int sternpower;

int main() {
	bowpower = 4;
	sternpower = 4;
	counter = 0;
	ifstream c4table;
	c4table.open("c4table.txt");
	ifstream weighttable;
	weighttable.open("weighttable.txt");
	ReadC4Table(c4table);
	ReadWeightsTable(weighttable);


	cout << "*****************\n";
	cout << "* CANOE CREATOR *\n";
	cout << "*****************\n\n";
	cout << "Version 1.02\n\n\n\n\n";//Created for the University of Toronto Concrete Canoe Team by David Ruggiero

	int flag = 0;
	while (flag != 1) {
		flag = UserInterface();
	}

	/*


//	ofstream meshout;
//	meshout.open("Mesh.txt");
//	c1.output.open("Output.txt");

	ofstream resultsout;
	resultsout.open("Results.txt");

	ifstream inputsetup;
	inputsetup.open("inputsetup.txt");

	ofstream writeinput;
	writeinput.open("inputtable.txt");

	CreateInputs(inputsetup,writeinput);
	writeinput.close();
	
	ifstream inputtable;
	inputtable.open("inputtable.txt");

	int numcanoes;
	char buffer[128];
	inputtable >> numcanoes;			// first line of file
	int i;
	for (i=0; i < 15; i++) {
		inputtable >> buffer;
	}

	ofstream meshout;
	meshout.open("Mesh.txt");

	Canoe c;
	int flag;
	for (i=0; i < numcanoes; i++) {
		flag = 0;
		double canoenum,L,Lp,Ld,Lf,W,t1,t2,d,h,b,s,f,n,dens;
		inputtable >> canoenum >> L >> Lp >> Ld >> Lf >> W >> t1 >> t2 >> d >> h >> b >> s >> f >> n >> dens;
		c.InitializeCanoe(L,Lp,Ld,Lf,W,t1,t2,d,h,b,s,f,n,dens);
		flag = c.AnalyzeAll();
		resultsout << canoenum;

		if (flag == 1) {
			resultsout << '\t' << "Analysis Failed\n";
			cout << "Canoe " << canoenum << " analysis failed.\n";
		} else {
			for (int j=0; j<17; j++) {
				resultsout << '\t' << c.outputs[j];
			}
			resultsout << '\t' << Score(c.outputs,targets[0],targets[1],targets[2],17) << endl;
			cout << "Canoe " << canoenum << " successfully analyzed.\n";
	//		cout << canoenum << '\t' << Score(c.outputs,targets[0],targets[1],targets[2],17) << endl;
		}

		c.OutputMesh(meshout);

		c.Destruct();
	}


	//c.OutputMesh(meshout);

//	c1.InitializeCanoe(6,2,3,2,0.71,1,1,0.5,0,0.05,0.05,0,0.7,10);

//	ofstream canoeout;
//	canoeout.open("Canoe Output.txt");
//	c1.AnalyzeAll();
//	c1.OutputAll(canoeout);

	cout << counter << endl;


	cin >> i;


  */

	return 0;
}


double SplineArea(double P0y, double P0z, double P1y, double P1z,
						double P2y, double P2z, double P3y, double P3z,
						double ti, double tf) {
	
	double A1 = P0y - 3*P1y + 3*P2y - P3y;
	double B1 = 3*P1y - 6*P2y + 3*P3y;
	double C1 = 3*P2y - 3*P3y;
	double D1 = P3y;

	double A2 = 0;
	double B2 = 3*P0z - 9*P1z + 9*P2z - 3*P3z;
	double C2 = 6*P1z - 12*P2z + 6*P3z;
	double D2 = 3*P2z - 3*P3z;

	double E[7];
	E[6] = A1*A2;
	E[5] = A1*B2 + B1*A2;
	E[4] = A1*C2 + B1*B2 + C1*A2;
	E[3] = A1*D2 + B1*C2 + C1*B2 + D1*A2;
	E[2] = B1*D2 + C1*C2 + D1*B2;
	E[1] = C1*D2 + D1*C2;
	E[0] = D1*D2;

	double area = 0;
	int i;
	for (i = 0; i < 7; i++) {
		area += E[i]/(i+1) * (pow(tf,(i+1)) - pow(ti,(i+1)));
	}
	return area;
}

int SplineMoment(double P0y, double P0z, double P1y, double P1z,
							 double P2y, double P2z, double P3y, double P3z,
							 double ti, double tf, double& momentareay, double& momentareaz) {

	double A1 = P0y - 3*P1y + 3*P2y - P3y;
	double B1 = 3*P1y - 6*P2y + 3*P3y;
	double C1 = 3*P2y - 3*P3y;
	double D1 = P3y;

	double A2 = 0;
	double B2 = 3*P0z - 9*P1z + 9*P2z - 3*P3z;
	double C2 = 6*P1z - 12*P2z + 6*P3z;
	double D2 = 3*P2z - 3*P3z;

	double A3 = P0z - 3*P1z + 3*P2z - P3z;
	double B3 = 3*P1z - 6*P2z + 3*P3z;
	double C3 = 3*P2z - 3*P3z;
	double D3 = P3z;

	double A4 = 0;
	double B4 = 3*P0y - 9*P1y + 9*P2y - 3*P3y;
	double C4 = 6*P1y - 12*P2y + 6*P3y;
	double D4 = 3*P2y - 3*P3y;

	double E[7];
	E[6] = A1*A2;
	E[5] = A1*B2 + B1*A2;
	E[4] = A1*C2 + B1*B2 + C1*A2;
	E[3] = A1*D2 + B1*C2 + C1*B2 + D1*A2;
	E[2] = B1*D2 + C1*C2 + D1*B2;
	E[1] = C1*D2 + D1*C2;
	E[0] = D1*D2;

	double L[10];
	L[9] = A3*E[6];
	L[8] = B3*E[6] + A3*E[5];
	L[7] = C3*E[6] + B3*E[5] + A3*E[4];
	L[6] = D3*E[6] + C3*E[5] + B3*E[4] + A3*E[3];
	L[5] = D3*E[5] + C3*E[4] + B3*E[3] + A3*E[2];
	L[4] = D3*E[4] + C3*E[3] + B3*E[2] + A3*E[1];
	L[3] = D3*E[3] + C3*E[2] + B3*E[1] + A3*E[0];
	L[2] = D3*E[2] + C3*E[1] + B3*E[0];
	L[1] = D3*E[1] + C3*E[0];
	L[0] = D3*E[0];

	double Ey[7];
	Ey[6] = A3*A4;
	Ey[5] = A3*B4 + B3*A4;
	Ey[4] = A3*C4 + B3*B4 + C3*A4;
	Ey[3] = A3*D4 + B3*C4 + C3*B4 + D3*A4;
	Ey[2] = B3*D4 + C3*C4 + D3*B4;
	Ey[1] = C3*D4 + D3*C4;
	Ey[0] = D3*D4;

	double Ly[10];
	Ly[9] = A1*Ey[6];
	Ly[8] = B1*Ey[6] + A1*Ey[5];
	Ly[7] = C1*Ey[6] + B1*Ey[5] + A1*Ey[4];
	Ly[6] = D1*Ey[6] + C1*Ey[5] + B1*Ey[4] + A1*Ey[3];
	Ly[5] = D1*Ey[5] + C1*Ey[4] + B1*Ey[3] + A1*Ey[2];
	Ly[4] = D1*Ey[4] + C1*Ey[3] + B1*Ey[2] + A1*Ey[1];
	Ly[3] = D1*Ey[3] + C1*Ey[2] + B1*Ey[1] + A1*Ey[0];
	Ly[2] = D1*Ey[2] + C1*Ey[1] + B1*Ey[0];
	Ly[1] = D1*Ey[1] + C1*Ey[0];
	Ly[0] = D1*Ey[0];

	momentareaz = 0;
	momentareay = 0;
	int i;
	for (i = 0; i < 10; i++) {
		momentareaz +=  L[i]/(i+1) * (pow(tf,(i+1)) - pow(ti,(i+1)));
		momentareay += Ly[i]/(i+1) * (pow(tf,(i+1)) - pow(ti,(i+1)));
	}

	double Piy = Spline(A1,B1,C1,D1,ti);
	double Piz = Spline(A3,B3,C3,D3,ti);
	double Pfy = Spline(A1,B1,C1,D1,tf);
	double Pfz = Spline(A3,B3,C3,D3,tf);

	momentareay += Piy*Piy*Piz/2 - Pfy*Pfy*Pfz/2;	// extra rectangular regions

	momentareaz = momentareaz;		// I think the negatives should be here...
	momentareay = -momentareay;			// But not here?

	return 0;
}

double SolveCubic(double a, double b, double c, double d) {
/* Only works if there's 1 real root
	double q,r,s,t;
	
	q = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a);
	r = sqrt(pow(((3*a*c - b*b)/(9*a*a)),3) + q*q);
	s = pow((q+r),1./3.);
	if (q>r)
		t = pow((q-r),1./3.);
	else
		t = -pow((r-q),1./3.);

	return (s + t - b/(3.*a));
*/

	double ap,bp,cp,p,q,f;
	double rootfr,rootfi,gr,gi,magg,argg,magu1,argu1,u1r,u1i,u2r,u2i,u3r,u3i,x1r,x1i,x2r,x2i,x3r,x3i;

	ap = b/a;
	bp = c/a;
	cp = d/a;

	p = bp - ap*ap/3.;
	q = cp + (2*ap*ap*ap - 9*ap*bp)/27.;

	if (p == 0 && q == 0)
		return -ap/3.;
	if (p == 0) {
		magg = q;
		argg = 0;
	}
	else {
		f = q*q/4. + p*p*p/27.;
		if (f < 0) {
			rootfr = 0;
			rootfi = sqrt(-f);
		}
		else {
			rootfr = sqrt(f);
			rootfi = 0;
		}

		gr = rootfr + q/2.;
		gi = rootfi;

		magg = sqrt(gr*gr+gi*gi);
		argg = atan2(gi,gr) + (rootfi < 0 ? 2*PI : 0);
	}

	magu1 = pow(magg,1./3.);
	argu1 = argg/3.;

	u1r = magu1*cos(argu1);
	u1i = magu1*sin(argu1);

	u2r = -0.5*u1r - sqrt(3)/2.*u1i;
	u2i = -0.5*u1i + sqrt(3)/2.*u1r;

	u3r = -0.5*u1r + sqrt(3)/2.*u1i;
	u3i = -0.5*u1i - sqrt(3)/2.*u1r;

	x1r = p/3.*u1r/(u1r*u1r+u1i*u1i) - u1r - ap/3.;
	x1i = -p/3.*u1i/(u1r*u1r+u1i*u1i) - u1i;

	x2r = p/3.*u2r/(u2r*u2r+u2i*u2i) - u2r - ap/3.;
	x2i = -p/3.*u2i/(u2r*u2r+u2i*u2i) - u2i;

	x3r = p/3.*u3r/(u3r*u3r+u3i*u3i) - u3r - ap/3.;
	x3i = -p/3.*u3i/(u3r*u3r+u3i*u3i) - u3i;

//	cout << x1r << " + " << x1i << "i\n";
//	cout << x2r << " + " << x2i << "i\n";
//	cout << x3r << " + " << x3i << "i\n";
	
	// determine which real root falls between 0 and 1 - if none do, return -1
	double threshold = 1e-10;		// zero threshold
	double ROUNDOFF = threshold*1000;	// this one is killing me
	if (x1r <= 1+ROUNDOFF && x1r >= 0-ROUNDOFF && fabs(x1i)<threshold)
		return x1r;
	else if (x2r <= 1+ROUNDOFF && x2r >= 0-ROUNDOFF && fabs(x2i)<threshold)
		return x2r;
	else if (x3r <= 1+ROUNDOFF && x3r >= 0-ROUNDOFF && fabs(x3i)<threshold)
		return x3r;

	return -1;

}

int SectionResultant(double P0y, double P0z, double P1y, double P1z,
						double P2y, double P2z, double P3y, double P3z,
						double d, double theta,												// d is freeboard
						double& magnitude, double& yloc, double& zloc) {
	
	if (d > -P3z) {
		magnitude = 0;
		yloc = 0;
		zloc = 0;
		return 0;
	}

	int flag = 0;

	double A1 = P0y - 3*P1y + 3*P2y - P3y;
	double B1 = 3*P1y - 6*P2y + 3*P3y;
	double C1 = 3*P2y - 3*P3y;
	double D1 = P3y;

	double A3 = P0z - 3*P1z + 3*P2z - P3z;
	double B3 = 3*P1z - 6*P2z + 3*P3z;
	double C3 = 3*P2z - 3*P3z;
	double D3 = P3z;

	double ti, tfr, tfl;
	if (d == 0 && theta == 0) {
		ti = 0;
		tfr = 1;
		tfl = 1;
	} else {
		ti = 0;								// lower bound
		tfr = SolveCubic(A3-tan(theta)*A1, B3-tan(theta)*B1, C3-tan(theta)*C1, D3-tan(theta)*D1+d-P0z);	// upper bound right side
		if (tfr < 0)
			flag = -1;				// leak point
		tfl = SolveCubic(A3+tan(theta)*A1, B3+tan(theta)*B1, C3+tan(theta)*C1, D3+tan(theta)*D1+d-P0z);	// upper bound left side
		if (tfl < 0)
			flag = -1;
	}

	double arearspline = SplineArea(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfr);
	double arealspline = SplineArea(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfl);

	double momyrspline, momzrspline, momylspline, momzlspline;
	SplineMoment(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfr,momyrspline,momzrspline);
	SplineMoment(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfl,momylspline,momzlspline);

	momylspline = - momylspline;

	double ytfr = Spline(A1,B1,C1,D1,tfr);
	double ztfr = Spline(A3,B3,C3,D3,tfr);
	double ytfl = -Spline(A1,B1,C1,D1,tfl);
	double ztfl = Spline(A3,B3,C3,D3,tfl);

	double areartri = -ytfr*((P0z-d)-ztfr)/2;
	double arealtri = -ytfl*((P0z-d)-ztfl)/2;

	double arear = arearspline - areartri;		// subtract triangular portion - there's an explicit negative here
	double areal = arealspline + arealtri;		// add triangular portion

	magnitude = arear + areal;

	double ymoment = momyrspline - 1./3.*ytfr*areartri + momylspline + 1./3.*ytfl*arealtri;
	double zmoment = momzrspline - (ztfr - 1./3.*(ztfr - (P0z-d)))*areartri
						  + momzlspline + (P0z-d - 2./3.*((P0z-d) - ztfl))*arealtri;

	yloc = (magnitude == 0 ? 0 : ymoment/magnitude);	
	zloc = (magnitude == 0 ? 0 : zmoment/magnitude);

	if (flag == -1)
		return -1;
	return 0;
}

double Spline(double A, double B, double C, double D, double t) {
	return (A*t*t*t + B*t*t + C*t + D);
}

int Canoe::ControlPoints(int station, double& P0y, double& P0z, double& P1y, double& P1z,
				  double& P2y, double& P2z, double& P3y, double& P3z) {
//	double x = station*increment;
	double w = wlvalues[station];
	double d = klvalues[station];
	double flareangle;
	double fmax = asin(w/(2.*d));
	if (fmax < flare)
		flareangle = fmax;
	else
		flareangle = flare;

	double u1max = d/cos(flareangle);
	double u2max = w - d*tan(flareangle);
	double u1 = u1max*(0.25+0.75*shapeparam);
	double u2 = u2max*(0.5+0.5*shapeparam);

	P0y = w;
	P0z = 0;
	P1y = P0y - u1*sin(flareangle);
	P1z = P0z - u1*cos(flareangle);
	P3y = 0;
	P3z = -d;
	P2y = P3y + u2;
	P2z = P3z;

	return 0;
}

// OVERLOADED FUNCTION
// This version works for any x value, not just at predefined stations. This one is more computationally intensive,
// so it should be used only sparingly
int Canoe::ControlPoints(double x, double& P0y, double& P0z, double& P1y, double& P1z,
				  double& P2y, double& P2z, double& P3y, double& P3z) {
	double w = Waterline(x);
	double d = Keelline(x);
	double flareangle;
	double fmax = acos(w/(2.*d));
	if (fmax < flare)
		flareangle = fmax;
	else
		flareangle = flare;

	double u1max = d/cos(flareangle);
	double u2max = w - d*tan(flareangle);
	double u1 = u1max*(0.25+0.75*shapeparam);
	double u2 = u2max*(0.5+0.5*shapeparam);

	P0y = w;
	P0z = 0;
	P1y = P0y - u1*sin(flareangle);
	P1z = P0z - u1*cos(flareangle);
	P3y = 0;
	P3z = -d;
	P2y = P3y + u2;
	P2z = P3z;

	return 0;

}

double Ramp(double x, double t) {
	return (t/2. * log(0.5*exp(x/t) + 0.5*exp(-x/t)) + x/2. + t*log(2.)/2.);
}

double Canoe::Waterline(double x) {
	double m1 = width/(2.*lfirst);
	double m2 = width/(2.*(length - lpaddler - lfirst));
	double rbow = width/2. - m1*Ramp(width/(2.*m1),smooth1) - m2*Ramp(-length + width/(2.*m2),smooth2);
	double rstern = width/2. - m1*Ramp(width/(2.*m1)-length,smooth1) - m2*Ramp(width/(2.*m2),smooth2);

	return (width/2. - m1*Ramp(width/(2.*m1)-x,smooth1) - m2*Ramp(x-length+width/(2.*m2),smooth2) - (rstern-rbow)*x/length - rbow);
}

double Canoe::Keelline(double x) {
	if (x < ldeepest) {
		return (depth - brocker*fabs(pow((1.-x/ldeepest),bowpower)));
	}
	else
		return (depth - srocker*fabs((pow((ldeepest-x)/(length-ldeepest),sternpower))));
	return 0;
}

int Canoe::InitializeCanoe(double L, double Lp, double Ld, double Lf, double W, double t1, double t2, double d, double h, double b,
						  double s, double f, double n, double density) {
	numstations = 101;
	increment = L/(numstations-1);

	// Default tolerance values
	wltol = 1e-4;	
	tippingtol = 1e-2;
	bowtol = 1e-4;
	wlmaxit = 50;
	tippingmaxit = 50;

	wlvalues = new double[numstations];
	klvalues = new double[numstations];

	length = L;
	lpaddler = Lp;
	ldeepest = Ld;
	lfirst = Lf;
	width = W;
	smooth1 = t1;
	smooth2 = t2;
	depth = d;
	hflange = h;
	brocker = b;
	srocker = s;
	flare = f;
	shapeparam = n;
	this->density = density;

	int i;
	for (i=0; i<numstations; i++) {
		wlvalues[i] = Waterline(i*increment);
		klvalues[i] = Keelline(i*increment);
	}
	
	// 3 men load case - turning is important
	loadcase[0].numpaddlers = 3;
	loadcase[0].paddlerweight = 85;			// about 185 pounds
	loadcase[0].paddlercm = .40;		// 40 cm (rough estimate)

	// 2+2 load case - speed is important
	loadcase[1].numpaddlers = 4;
	loadcase[1].paddlerweight = 72.5;			// about 135 pounds for women
	loadcase[1].paddlercm = .35;			// about 30 cm for women

	return 0;
}

int Canoe::Destruct() {
	delete[] wlvalues;
	delete[] klvalues;
	return 0;
}

Canoe::Canoe(double L, double Lp, double Ld, double Lf, double W, double t1, double t2, double d, double h, double b,
						  double s, double f, double n) {
	length = L;
	lpaddler = Lp;
	ldeepest = Ld;
	lfirst = Lf;
	width = W;
	smooth1 = t1;
	smooth2 = t2;
	depth = d;
	hflange = h;
	brocker = b;
	srocker = s;
	flare = f;
	shapeparam = n;
}

Canoe::Canoe() {

}

int Canoe::Analyze(double d, double theta, double trim, double& volume, double& cenx, double& ceny, double& cenz) {

	double* areas = new double[numstations];
	double* centroidx = new double[numstations];
	double* centroidy = new double[numstations];
	double* centroidz = new double[numstations];

	double x=0;
	int i;

	double P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z;
	int flag = 0;
	
	for (i=0; i < numstations; i++) {
		ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z);
		// freeboard at distance x is given by d - trim/length*x
		if (SectionResultant(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z, (trim == 0 ? d : d-trim/length*x),theta, areas[i],centroidy[i],centroidz[i]) == -1) {
			//cout << "LEAK\n";			// let it pass for now
			flag = 1;
			//return 1;
		}
		centroidx[i] = x;

		x += increment;
	}

	volume = 0;
	for (i=1; i < numstations; i++) {							// numerically integrate using parallelograms
		volume += (areas[i]+areas[i-1])/2;
	}
	volume *= increment;

	cenx = 0;
	ceny = 0;
	cenz = 0;

	for (i=0; i < numstations; i++) {
		cenx += centroidx[i]*areas[i]*increment;			// this is not 100% correct... but it's not too far wrong
		ceny += centroidy[i]*areas[i]*increment;
		cenz += centroidz[i]*areas[i]*increment;
	}
	cenx /= volume;
	ceny /= volume;
	cenz /= volume;

	delete[] areas;
	delete[] centroidx;
	delete[] centroidy;
	delete[] centroidz;

	if (flag == 1)
		return 1;

	return 0;
}

double Canoe::FindWLine(double vtarget, double theta, double trim) {
	double v,cx,cy,cz;
	int MAXITERATIONS = wlmaxit;
	double tolerance = wltol;

	double testd = depth/2;
	double width = depth/2;

	int flag;
	// Check if it's possible to attain equilibrium (only valid for the zero trim case)
	if (trim == 0) {
		double width = GetBWL(0)/2.;
		double dcheck = width*tan(theta);
		Analyze(dcheck,theta,trim,v,cx,cy,cz);
		if (v < vtarget)
			return -1;
	}

	// binary-type algorithm
	int i;
	for (i=0; i < MAXITERATIONS; i++) {
		flag = Analyze(testd,theta,trim,v,cx,cy,cz);
		if (fabs(v-vtarget) < tolerance)
			return testd;
//		if (v < vtarget && flag == 1)
//			return -1;								// equilibrium not possible -- leak
		if (v > vtarget || flag == 1)
			testd += width/2.;
		else 
			testd -= width/2.;
		width /= 2.;
	}

//	cout << "WARNING: LOSS OF PRECISION\n";
	return -1;

}

double Moment(double y1, double z1, double y2, double z2, double weight, double theta) {
	return -(cos(theta)*(y2-y1) - sin(theta)*(z2-z1))*weight;
}


// NOTE - cmheight is the centre of mass height taken ABOVE THE TOP OF THE GUNWALE
// The tipping angle calculation assumes the leak happens at the widest point (which is always true for zero trim)
int Canoe::CrossCurve(double weight, double cmheight, double& tangent, double& tipping) {
	double theta = 0;
	double tinc = 0.01;
	double v,cx,cy,cz;
	
	// roughly get tangent	
	double theta1 = theta;
	double theta2 = theta+tinc;
	double d1 = FindWLine(weight, theta1, 0);
	Analyze(d1,theta1,0,v,cx,cy,cz);
	double moment1 = Moment(cy,cz,0,cmheight,weight,theta);
	double d2 = FindWLine(weight, theta2, 0);	
	Analyze(d2,theta2,0,v,cx,cy,cz);
	double moment2 = Moment(cy,cz,0,cmheight,weight,theta);
	tangent = (moment2-moment1)/(theta2-theta1);

	// find critical tipping angle
	double tolerance = tippingtol;
	double width = GetBWL(0)/2.;
	double dtest = d1;					// start off with d for 0 heel
	double dmax = depth;					// initial maximum value of d is depth of canoe
	double dmin = 0.;						// initial minimum value of d is 0


	theta = atan2(dtest,width);
	Analyze(dtest,theta,0,v,cx,cy,cz);
	int MAXITERATIONS = tippingmaxit;
	int count = 0;
	int flag = 0;
	while (fabs(weight-v) > tolerance) {
		if (count >= MAXITERATIONS) {
			flag = 1;
			break;
		}
		count++;
		if (weight < v) {
			dmin = dtest;
			dtest = (dtest+dmax)/2.;
		} else {
			dmax = dtest;
			dtest = (dtest+dmin)/2.;
		}
		theta = atan2(dtest,width);
		Analyze(dtest,theta,0,v,cx,cy,cz);
	}

	tipping = theta;

	return flag;
}



/* OLD algorithm - it works too but it is very slow
// NOTE - cmheight is the centre of mass height taken ABOVE THE TOP OF THE GUNWALE
int Canoe::CrossCurve(double weight, double cmheight, double& tangent, double& tipping) {
	double theta = 0;
	double tinc = 0.01;
	double v,cx,cy,cz;
	
	// roughly get tangent	
	double d = FindWLine(weight, theta, 0);
	double theta1 = theta;
	double theta2 = theta+tinc;
	Analyze(d,theta1,0,v,cx,cy,cz);
	double moment1 = Moment(cy,cz,0,cmheight,weight,theta);
	Analyze(d,theta2,0,v,cx,cy,cz);
	double moment2 = Moment(cy,cz,0,cmheight,weight,theta);
	tangent = (moment2-moment1)/(theta2-theta1);

	while (d >= 0) {
		Analyze(d,theta,0,v,cx,cy,cz);
		
		output << theta << '\t' << Moment(cy,cz,0,cmheight,weight,theta) << endl;

		theta += tinc;

		d = FindWLine(weight, theta, 0);
	}

	tipping = theta - tinc;

	return 0;
}
*/

int Canoe::DisplacementCurve() {
	double d = depth;
	double dinc = 0.01;
	double v,cx,cy,cz;
	while (d >= 0) {
		Analyze(d,0,0,v,cx,cy,cz);

		output << d << '\t' << v << endl;

		d -= dinc;
	}

	return 0;
}

double SplineLength(double P0y, double P0z, double P1y, double P1z,
						double P2y, double P2z, double P3y, double P3z,
						double d, double& cmz) {
	if (d > -P3z) {
		return 0;
	}

	double A1 = P0y - 3*P1y + 3*P2y - P3y;
	double B1 = 3*P1y - 6*P2y + 3*P3y;
	double C1 = 3*P2y - 3*P3y;
	double D1 = P3y;

	double A3 = P0z - 3*P1z + 3*P2z - P3z;
	double B3 = 3*P1z - 6*P2z + 3*P3z;
	double C3 = 3*P2z - 3*P3z;
	double D3 = P3z;

	double ti = 0;								// lower bound
	double tf = SolveCubic(A3, B3, C3, D3+d-P0z);	// upper bound
	if (tf < 0)
		tf = 0;			// if there's a leak we ignore the section - presumably there should not be any leaks in this routine
	if (d == 0)
		tf = 1;

	double t;
	double tinc = 0.01;	// UNDEFINED CONSTANT
	double y1,z1,y2,z2;
	double l = 0;

	y1 = Spline(A1,B1,C1,D1,ti);
	z1 = Spline(A3,B3,C3,D3,ti);

	double zmoment = 0;

	for (t = ti+tinc; t <= tf; t+= tinc) {
		y2 = Spline(A1,B1,C1,D1,t);
		z2 = Spline(A3,B3,C3,D3,t);
		double segment = sqrt((y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
		l += segment;
		zmoment += segment*(z1+z2)/2.;		// moment = length * average height

		y1 = y2;
		z1 = z2;
	}

	cmz = zmoment/l;		// z centre of mass

	return l;
}

double Canoe::SurfaceArea(double d, double& cx, double& cz) {
	double area = 0;
	double areamoment = 0;
	double areazmoment = 0;
	int i;
	double* lengths = new double[numstations];
	double* cenzs = new double[numstations];
	double P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z;
//	double x=0;
	double cmz;
	for (i=0; i < numstations; i++) {
		ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z);
		lengths[i] = SplineLength(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,d,cmz);
		//x += increment;
		cenzs[i] = cmz;			// 
	}
	for (i=1; i < numstations; i++) {							// numerically integrate using parallelograms
		area += (lengths[i]+lengths[i-1])/2;
		areamoment += (lengths[i]*i*increment);
		areazmoment += cenzs[i]*lengths[i];
	}
	cx = areamoment/area;
	cz = areazmoment/area;

	area *= increment;
	area *= 2.;

	delete[] lengths;
	delete[] cenzs;

	return area;
}


// John Winters' KAPER formula - disp is in long tonnes, lengths in feet, speed in knots, LCB is a fraction of length
double Canoe::Kaper(double BWL, double EWL, double WS, double Cp, double Cv, double LCB, double disp, double le, double v) {
	double vtol = v/sqrt(EWL);

	double c1 = 0.002*sqrt(BWL/EWL)*pow((4*vtol),4);
	double c2 = 0.005*(sin(le*360/2./PI))*(4*vtol)*(4*vtol);
	double c3 = (vtol < 1.5 ? 1 : 0.7)*0.8*cos(3.65*vtol+0.07) + (vtol < 1.5 ? 1 : 0.96)*1.2;
	double c4 = GetC4(Cv,Cp,vtol);
	double c5 = pow((0.5/LCB),0.35);
	double c6 = 1.0;

	double Rr = (4*c5*disp*pow(vtol,4) + c1 + c2)*(c3*c4*c5);
	if (v >= 1.6 && v < 3)
		Rr *= 0.75;
	if (v >= 1.4 && v < 1.6)
		Rr *= 0.85;

	double reynold = v*1.6889*EWL/1.2791*100000;
	double Cf = 0.075/pow((log10(reynold)-2),2);
	double Rf = 0.99525*Cf*WS*(v*1.6889)*(v*1.6889);		// speed in feet/sec

	return (Rf+Rr);
}


int ReadC4Table(ifstream& c4) {
	int i,j,k;
	for (i=0; i<3; i++) {
		for (j=0; j<17; j++) {
			for (k=0; k<14; k++) {
				c4 >> C4TABLE[i][j][k];
			}
		}
	}
	return 0;
}

// Look up c4 in the table
double GetC4(double Cv, double Cp, double vtol) {
	int Cvind1 = (int)((Cv - .001)/.0005);
	int Cvind2 = (int)((Cv - .001)/.0005) + 1;
	if (Cvind1 < 0) Cvind1 = 0;
	if (Cvind1 > 2) Cvind1 = 2;
	if (Cvind2 < 0) Cvind2 = 0;
	if (Cvind2 > 2) Cvind2 = 2;

	int Cpind1 = (int)((Cp - 0.48)/0.01);
	int Cpind2 = (int)((Cp - 0.48)/0.01) + 1;
	if (Cpind1 < 0) Cpind1 = 0;
	if (Cpind1 > 16) Cpind1 = 16;
	if (Cpind2 < 0) Cpind2 = 0;
	if (Cpind2 > 16) Cpind2 = 16;

	int vtolind1 = (int)((vtol - 0.4)/0.1);
	int vtolind2 = (int)((vtol - 0.4)/0.1) + 1;
	if (vtolind1 < 0) vtolind1 = 0;
	if (vtolind1 > 13) vtolind1 = 13;
	if (vtolind2 < 0) vtolind2 = 0;
	if (vtolind2 > 13) vtolind2 = 13;

	double Cvratio, Cpratio, vtolratio;
	if (Cvind1 == Cvind2)
		Cvratio = 0;
	else
		Cvratio = (Cv - (Cvind1*.0005 + .001))/(Cvind2*.0005 - Cvind1*.0005);
	if (Cpind1 == Cpind2)
		Cpratio = 0;
	else
		Cpratio = (Cp - (Cpind1*.01 + .48))/(Cpind2*.01 - Cpind1*.01);
	if (vtolind1 == vtolind2)
		vtolratio = 0;
	else
		vtolratio = (vtol - (vtolind1*.1 + .4))/(vtolind2*.1 - vtolind1*.1);

	// 8 corners of the prism to be interpolated
	double corners[2][2][2];
	int i,j,k;

	for (i=Cvind1; i<=Cvind2; i++) {
		for (j=Cpind1; j<=Cpind2; j++) {
			for (k=vtolind1; k<=vtolind2; k++) {
				corners[i-Cvind1][j-Cpind1][k-vtolind1] = C4TABLE[i][j][k];
			}
		}
	}

	double actualCv[2][2];
	for (i=0; i<=1; i++) {
		for (j=0; j<=1; j++) {
			actualCv[i][j] = Cvratio*corners[1][i][j] + (1-Cvratio)*corners[0][i][j];
		}
	}

	double actualCp[2];
	for (i=0; i<=1; i++) {
		actualCp[i] = Cpratio*actualCv[1][i] + (1-Cpratio)*actualCv[0][i];
	}

	double c4 = vtolratio*actualCp[1] + (1-vtolratio)*actualCp[0];
/*
	cout << "Cv from " << Cvind1 << " to " << Cvind2 << endl;
	cout << "Cp from " << Cpind1 << " to " << Cpind2 << endl;
	cout << "VtoL from " << vtolind1 << " to " << vtolind2 << endl;

	cout << "Cv ratio: " << Cvratio << endl;
	cout << "Cp ratio: " << Cpratio << endl;
	cout << "VtoL ratio: " << vtolratio << endl;

	cout << "C4 = " << c4 << endl;
*/
	return c4;
}

// Get "Effective waterline length" (metres)
double Canoe::GetEWL(double d, double theta, double trim) {
	double* areas = new double[numstations];

	double x=0;
	int i;
	double dummy1,dummy2;
	double areamax = 0;
	int xmax = 0;

	double P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z;
	
	for (i=0; i < numstations; i++) {
		ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z);
		// freeboard at distance x is given by d - trim/length*x
		if (SectionResultant(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z, d-trim/length*x,theta, areas[i],dummy1,dummy2) == -1) {
			cout << "LEAK\n";
			delete[] areas;
			return 1;
		}
		if (areas[i] > areamax) {
			areamax = areas[i];
			xmax = i;
		}
		x += increment;
	}

/*
	output << "\n\nCurve of Areas\n";
	for (i=0; i<numstations; i++) {
		output << i*increment << '\t' << areas[i] << endl;
	}
*/

	int xstart,xend;
	bool flag;

	for (xstart=0; xstart<xmax; xstart++) {
		flag = false;
		for (i=xstart; i<xmax; i++) {
			if ((i-xstart)/(xmax-xstart)*areamax > areas[i]) {			// Find (approx) tangent to the curve
				flag = true;
				break;
			}
		}
		if (!flag)
			break;	// We now have xstart for EWL
	}

	for (xend=numstations-1; xend>xmax; xend--) {
		flag = false;
		for (i=xend; i>xmax; i--) {
			if ((xend-i)/(xend-xmax)*areamax > areas[i]) {			// Find (approx) tangent to the curve
				flag = true;
				break;
			}
		}
		if (!flag)
			break;	// We now have end for EWL
	}

	double EWL = (xend - xstart)*increment;
//	cout << "EWL = " << EWL << endl;
	delete[] areas;

	return EWL;
}

// Get waterline length (metres)
double Canoe::GetLWL(double d) {
	int i;
	bool flag = false;
	int starti,endi;

	for (i=0; i<numstations; i++) {
		if (klvalues[i] > d) {
			flag = true;
			starti = i;
			break;
		}
	}

	if (!flag)
		return length;

	flag = false;
	for (i=numstations-1; i>=0; i--) {
		if (klvalues[i] > d) {
			flag = true;
			endi = i;
			break;
		}
	}

	return ((endi-starti)*increment);
}

// Get waterline beam (metres)
double Canoe::GetBWL(double d) {
	double BWL = 0;
	double P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z;

	double width;

	int i;
	for (i=0; i<numstations; i++) {
		ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z);
		double A3 = P0z - 3*P1z + 3*P2z - P3z;
		double B3 = 3*P1z - 6*P2z + 3*P3z;
		double C3 = 3*P2z - 3*P3z;
		double D3 = P3z;		
		double tf = SolveCubic(A3, B3, C3, D3+d-P0z);

		double A1 = P0y - 3*P1y + 3*P2y - P3y;
		double B1 = 3*P1y - 6*P2y + 3*P3y;
		double C1 = 3*P2y - 3*P3y;
		double D1 = P3y;
		
		width = Spline(A1,B1,C1,D1,tf);
		if (width > BWL) {
			BWL = width;
		}
	}
	return BWL*2;
}

// Get waterline bow half angle (degrees)
double Canoe::Getle(double d) {
	double xentrance;
	double TOLERANCE = bowtol;
	double xinterval = ldeepest;
	double xpos = ldeepest/2.;
	if (depth - brocker > d) {
		xentrance = 0;
	} else {
		bool flag = false;
		while (!flag) {
			if (fabs(Keelline(xpos)-d) < TOLERANCE) {
				flag = true;
			} else {
				if (Keelline(xpos) > d) {
					xpos = xpos - xinterval/2.;
				}
				else {
					xpos = xpos + xinterval/2.;
				}
				xinterval /= 2.;
			}
		}
		xentrance = xpos;
	}

	double P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z;

	// have to use the expensive ControlPoints function
	ControlPoints(xentrance, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z);
	double A3 = P0z - 3*P1z + 3*P2z - P3z;
	double B3 = 3*P1z - 6*P2z + 3*P3z;
	double C3 = 3*P2z - 3*P3z;
	double D3 = P3z;		
	double tf = SolveCubic(A3, B3, C3, D3+d-P0z);
	double A1 = P0y - 3*P1y + 3*P2y - P3y;
	double B1 = 3*P1y - 6*P2y + 3*P3y;
	double C1 = 3*P2y - 3*P3y;
	double D1 = P3y;
	double w1 = Spline(A1,B1,C1,D1,tf);

	double dx = length/1000.;													// a sufficiently small dx
	ControlPoints(xentrance+dx, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z);			
	A3 = P0z - 3*P1z + 3*P2z - P3z;
	B3 = 3*P1z - 6*P2z + 3*P3z;
	C3 = 3*P2z - 3*P3z;
	D3 = P3z;		
	tf = SolveCubic(A3, B3, C3, D3+d-P0z);
	A1 = P0y - 3*P1y + 3*P2y - P3y;
	B1 = 3*P1y - 6*P2y + 3*P3y;
	C1 = 3*P2y - 3*P3y;
	D1 = P3y;
	double w2 = Spline(A1,B1,C1,D1,tf);

	double angle = atan2(w2-w1,dx);

	angle = angle/PI*180;

	return angle;
}


// Get Cv
double Canoe::GetCv(double disp, double LWL) {
	return (disp/(LWL*LWL*LWL));
}


// Get Cp
double Canoe::GetCp(double d, double LWL) {
	double* areas = new double[numstations];

//	double x=0;
	int i;
	double dummy1,dummy2;
	double areamax = 0;

	double P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z;
	
	for (i=0; i < numstations; i++) {
		ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z);
		// freeboard at distance x is given by d - trim/length*x
		if (SectionResultant(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z, d,0, areas[i],dummy1,dummy2) == -1) {
			cout << "LEAK\n";
			delete[] areas;
			return 1;
		}
		if (areas[i] > areamax) {
			areamax = areas[i];
		}
//		x += increment;
	}

	double v,cx,cy,cz;
	Analyze(d,0,0,v,cx,cy,cz);

	delete[] areas;

	return (v/LWL/areamax);
}


// Get Friction Force (Kaper wrapper function)
// Input v in m/s
double Canoe::GetFriction(double d, double v) {
	double volume,cx,cy,cz,ccenx,ccenz;
	if (Analyze(d,0,0,volume,cx,cy,cz)==1) {
		cout << "LEAK - Terminating\n";
		return -1;
	}
	double LWL = GetLWL(d);
	double BWL = GetBWL(d);
	double EWL = GetEWL(d,0,0);
	double WS = SurfaceArea(d,ccenx,ccenz);
	double Cp = GetCp(d,LWL);
	double Cv = GetCv(volume,LWL);
	double LCB = cx/length;
	double disp = volume;
	double le = Getle(d);

	// Convert units just to be sure...
	LWL *= 3.2808399;			// metres to feet
	BWL *= 3.2808399;			// metres to feet
	EWL *= 3.2808399;			// metres to feet
	WS *= 10.7639104;			// square metres to square feet
	disp *= 0.984206528;		// cubic metres of water to long tons
	v *= 1.94384449;			// m/s to knots

	double force = Kaper(BWL,EWL,WS,Cp,Cv,LCB,disp,le,v);

	return force;
}


int Canoe::OutputMesh (ofstream& meshout) {
	int i,j;
	double n = 20;
	double x,y,z,t;
	double P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z;
	meshout << "0\n\n";
	double A1;
	double B1;
	double C1;
	double D1;

	double A3;
	double B3;
	double C3;
	double D3;
	for (i=0;i<numstations;i++) {
		x = i*increment;			
		ControlPoints(i,P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z);																																																																	//Chen was here... and so was Kilroy
		A1 = P0y - 3*P1y + 3*P2y - P3y;
		B1 = 3*P1y - 6*P2y + 3*P3y;
		C1 = 3*P2y - 3*P3y;
		D1 = P3y;

		A3 = P0z - 3*P1z + 3*P2z - P3z;
		B3 = 3*P1z - 6*P2z + 3*P3z;
		C3 = 3*P2z - 3*P3z;
		D3 = P3z;
		for (j=0;j<=n;j++) {
			t = j/(double)n;
			y = Spline(A1,B1,C1,D1,t);
			z = Spline(A3,B3,C3,D3,t);
			meshout << length - x << "\t" << y << "\t" << z << "\n"; // x values need to be reversed for freeship
		}
		meshout << endl;
	}
	meshout << "EOF\n";
	return 0;
}


int Canoe::Waterplane(double d, double& wpcen,double& wparea, double& wpmom2) {
	double* depths = new double[numstations];

	double keel;
	double submerged;

	int i;
	for (i=0; i < numstations; i++) {
		keel = klvalues[i];
		submerged = keel - d;
		if (submerged < 0) submerged = 0;
		depths[i] = submerged;
	}	

	double area = 0;
	double firstmoment = 0;
	double segcentroid;
	double segarea;
	double centroid;
	for (i=1; i < numstations; i++) {
		segarea = (depths[i]+depths[i-1])/2.;
		area += segarea;				// parallelogram
		segcentroid = increment*(i  -  (2*depths[i-1]+depths[i])/3./(depths[i-1]+depths[i]));		// centroid of parallelogram
		firstmoment += segarea*segcentroid;
	}
	area *= increment;
	centroid = firstmoment/area*increment;	

	double secondmoment = 0;
	for (i=1; i < numstations; i++) {			// find second moments
		segarea = (depths[i]+depths[i-1])/2.;
		segcentroid = increment*(i  -  (2*depths[i-1]+depths[i])/3./(depths[i-1]+depths[i]));
		secondmoment += (segcentroid-centroid)*(segcentroid-centroid)*segarea;			// parallel axis theorem
	}
	secondmoment *= increment;

	wpcen = centroid;
	wparea = area;
	wpmom2 = secondmoment;

	delete[] depths;

	return 0;
}

int Canoe::AnalyzeAll() {
	double cx,cz;
	surfacearea = SurfaceArea(0,cx,cz);
	CMx = cx;
	CMz = cz;
	weight = surfacearea*density;
	actualbeam = GetBWL(0);

	int i;
	bool error = false;
	for (i=0; i<2; i++) {
		loadcase[i].totweight = (loadcase[i].numpaddlers*loadcase[i].paddlerweight + weight)/1000.;
		loadcase[i].totcm = (loadcase[i].numpaddlers*loadcase[i].paddlerweight*(loadcase[i].paddlercm - depth) + weight*CMz)/loadcase[i].totweight;
		
		double d = FindWLine(loadcase[i].totweight*1.0, 0, 0);		// density of water = 1.0
		loadcase[i].freeboard = d;
		if (d < 0) {
			error = true;
			break;
		}

		double LCB,ceny,cenz,volume;
		Analyze(d,0,0,volume,LCB,ceny,cenz);
		// report paddler centre as percentage of paddler box length away from the "ideal" paddler centre
		double paddlercentreabs = (LCB*volume - CMx*weight/1000.)/(loadcase[i].numpaddlers*loadcase[i].paddlerweight/1000.);
		loadcase[i].paddlercentre = (paddlercentreabs - (lfirst + 0.5*lpaddler))/lpaddler;				// I THINK THIS EQUATION IS FLAWED

		double LWL = GetLWL(d);
		loadcase[i].Cp = GetCp(d,LWL);
		double friction = GetFriction(d,2.57222222);	// 5 knots = about 10 km/h
		if (friction < 0) {
			error = true;
			break;
		}
		loadcase[i].friction5 = friction;		
	
		double tangent,tipping;
		if (CrossCurve(loadcase[i].totweight*1.0,loadcase[i].totcm,tangent,tipping) == 1) {
			error = true;
			break;
		}
		loadcase[i].initstability = tangent;
		loadcase[i].tippingangle = tipping;

		double wpcen,wparea,wpmom2;
		Waterplane(d,wpcen,wparea,wpmom2);
		loadcase[i].waterplane2 = wpmom2;
		loadcase[i].waterplanecentroid = wpcen/length;		// give it as a percentage of length
	}

	outputs[0] = weight;
	outputs[1] = loadcase[0].Cp;
	outputs[2] = loadcase[0].freeboard;
	outputs[3] = loadcase[0].friction5;
	outputs[4] = loadcase[0].initstability;
	outputs[5] = loadcase[0].tippingangle;
	outputs[6] = loadcase[0].waterplane2;
	outputs[7] = loadcase[0].waterplanecentroid;
	outputs[8] = loadcase[0].paddlercentre;
	outputs[9] = loadcase[1].Cp;
	outputs[10] = loadcase[1].freeboard;
	outputs[11] = loadcase[1].friction5;
	outputs[12] = loadcase[1].initstability;
	outputs[13] = loadcase[1].tippingangle;
	outputs[14] = loadcase[1].waterplane2;
	outputs[15] = loadcase[1].waterplanecentroid;
	outputs[16] = loadcase[1].paddlercentre;

	return (error ? 1 : 0);
}


// Run this only after AnalyzeAll has been run
int Canoe::OutputAll(ofstream& canoeout) {
	canoeout << "Canoe Data\n\n";

	canoeout << "Input Values\n\n";
	canoeout << "Length\tNominal Beam\tDepth\tNominal Flare Angle\tLength to Deepest Point\tLength to First Paddler\tLength of Paddlers' Box\tShape Parameter\tBow Rocker\tStern Rocker\tArea Density\n"
		<< length << '\t' << width << '\t' << depth << '\t' << flare << '\t' << ldeepest << '\t' << lfirst << '\t' << lpaddler << '\t' << shapeparam << '\t' << brocker << '\t' << srocker << '\t' << density << "\n\n";

	canoeout << "Output Values\n\n";
	canoeout << "Surface Area\tWeight\tLCM\tActual Beam\n"
		<< surfacearea << '\t' << weight << '\t' << CMx << '\t' << actualbeam << "\n\n";

	canoeout << "3 Men Load Case\n";
	int i = 0;
	canoeout << "Cp\tFreeboard\tDrag at 5 knots\tInitial Stability\tLeak Angle\tSecond Moment of Waterplane Area\tWaterplane Centroid\tPaddler Centre\n"
		<< loadcase[i].Cp << '\t' << loadcase[i].freeboard << '\t' << loadcase[i].friction5 << '\t' << loadcase[i].initstability << '\t' << loadcase[i].tippingangle << '\t' << loadcase[i].waterplane2 << '\t' << loadcase[i].waterplanecentroid << '\t' << loadcase[i].paddlercentre << "\n\n";

	i = 1;
	canoeout << "4 Mixed Load Case\n";
	canoeout << "Cp\tFreeboard\tDrag at 5 knots\tInitial Stability\tLeak Angle\tSecond Moment of Waterplane Area\tWaterplane Centroid\n"
		<< loadcase[i].Cp << '\t' << loadcase[i].freeboard << '\t' << loadcase[i].friction5 << '\t' << loadcase[i].initstability << '\t' << loadcase[i].tippingangle << '\t' << loadcase[i].waterplane2 << '\t' << loadcase[i].waterplanecentroid << '\t' << loadcase[i].paddlercentre << "\n\n";

	return 0;
}

double Score(double* values, double* targetmeans, double* targetstds, double* weights, int n) {
	int i;
	double score = 0;
	for (i=0; i < n; i++) {
		score += weights[i]*Gaussian(targetmeans[i],targetstds[i],values[i]);
	}
	return score;
}

double Gaussian(double mean, double std, double x) {
	return 1.*exp(-(x-mean)*(x-mean)/2/std/std);
}

int ReadWeightsTable(ifstream& targettable) {
	int i,j;
	char buffer[128];
	for (i=0; i<17; i++) {
		targettable >> buffer;
		for (j=0; j<3; j++) {
			targettable >> targets[j][i];
		}
	}
	return 0;
}

int CreateInputs(ifstream& inputsetup, ofstream& writeinput) {

	int numinputs;
	inputsetup >> numinputs;

	char buffer[128];
	double* mins = new double[numinputs];
	double* maxes = new double[numinputs];
	int* nums = new int[numinputs];
	int* ranges = new int[numinputs];

	int i=0;
	int product = 1;
	for (i=0; i < numinputs; i++) {
		inputsetup >> buffer;
		inputsetup >> mins[i];
		inputsetup >> maxes[i];
		inputsetup >> nums[i];
		ranges[i] = 1;
		product *= nums[i];
	}

	writeinput << product << endl;
	writeinput << "Canoe\tLength\tLp\tLd\tLf\tW\tt1\tt2\td\th\tb\ts\tf\tn\tdensity\n";

	i = 0;
	int count = 1;
	writeinput << count++ << '\t';
	OneInput(maxes,mins,nums,ranges,numinputs,writeinput);
	while(i<numinputs) {
		if (ranges[i] < nums[i]) {
			ranges[i]++;
			writeinput << count++ << '\t';
			OneInput(maxes,mins,nums,ranges,numinputs,writeinput);
			i = 0;
		} else if (ranges[i] == nums[i]){
			ranges[i] = 1;
			i++;
		}
	}

	delete[] mins;
	delete[] maxes;
	delete[] nums;
	delete[] ranges;

	return 0;
}


int OneInput(double* maxes, double* mins, int* nums, int* ranges, int numinputs, ofstream& writeinput) {
	int i;

	for (i=0; i < numinputs; i++) {
		writeinput << (nums[i]==1? mins[i] : (maxes[i]-mins[i])*(1.0*ranges[i]-1.0)/(1.0*nums[i]-1.0)+mins[i]) << '\t';
	}
	writeinput << endl;
	


	return 0;
}

int UserInterface() {

	// First level menu
	int numeric;

	cout << "(1) Analyze a single canoe design\n(2) Bulk analyze several canoes\n(3) Quit\n(4) Set curvature parameters\n\n> ";
	cin >> numeric;
	cout << "\n\n";

	int flag = 0;

	switch (numeric) {
	case 0:
		cout << "ABOUT\n\nCreated for the Universtiy of Toronto Concrete Canoe Team by David Ruggiero with input from Chen Chen, " <<
			"Lyle Gordon, Justin Shum, Jonathan Ho, Michael Ferri, Eva Chau, Timothy Reyes and others. Based on a formulation " <<
			"by Nicolas Lee, Adam Steinberg, Daniel Zaide, Cameron SR Fraser, and David Ruggiero.\n2007\n\n\n";
		break;
	case 1:
		UISingle();
		break;
	case 2:
		UIBulk();
		break;
	case 3:
		return 1;
	case 4:
		cout << "Input bow polynomial power: ";
		cin >> bowpower;
		cout << "input stern polynomial power: ";
		cin >> sternpower;
		cout << "\n\n";
		break;
	default:
		return 1;
	}




	return 0;
}

int UISingle() {
	cout << "(1) Create from input file\n(2) Create from manual input\n\n";
	cout << "> ";
	int numeric;
	cin >> numeric;
	cout << "\n\n";

	Canoe c;
	double L,Lp,Ld,Lf,W,t1,t2,d,h,b,s,f,n,dens;
	ifstream input;

	switch (numeric) {
	case 1:
		char buffer[128];
		cout << "Input filename: ";
		cin >> buffer;
		cout << "\n\n";
		input.open(buffer);
		if (input.fail()) {
			cout << "File not found\n\n\n";
			return 1;
		}

		input >> L >> Lp >> Ld >> Lf >> W >> t1 >> t2 >> d >> h >> b >> s >> f >> n >> dens;
		break;
	case 2:
		cout << "Length: ";	cin >> L;
		cout << "Length of paddlers' box: "; cin >> Lp;
		cout << "Length to deepest point: "; cin >> Ld;
		cout << "Length to first paddler: "; cin >> Lf;
		cout << "Maximum Width: "; cin >> W;
		cout << "Bow smoothing parameter: "; cin >> t1;
		cout << "Stern smoothing parameter: "; cin >> t2;
		cout << "Maximum Depth: "; cin >> d;
		h = 0;
		cout << "Bow rocker: "; cin >> b;
		cout << "Stern rocker: "; cin >> s;
		cout << "Flare angle: "; cin >> f;
		cout << "Shape parameter: "; cin >> n;
		cout << "Surface area density: "; cin >> dens;
		break;
	default:
		return 1;
	}

	cout << endl;
	c.InitializeCanoe(L,Lp,Ld,Lf,W,t1,t2,d,h,b,s,f,n,dens);

	cout << "Canoe created\n\n\n";
	int flag = 0;

	while (flag == 0) {
		cout << "(1) Analyze all aspects\n(2) Output mesh\n(3) Calculate resistance\n(4) Change default paddler data\n(5) Change default program settings\n(6) Back to main menu\n\n";

		cout << "> ";
		cin >> numeric;
		cout << "\n\n";

		ofstream output;
		int i, leakflag, loadcasenum;
		double speed,freeboard;

		char* outputs[8] = {"Cp","Freeboard","Drag at 5 knots","Initial Stability","Leak Angle","Second Moment of Waterplane Area","Waterplane Centroid","Paddler Centre"};

		switch (numeric) {
		case 1:
			leakflag = c.AnalyzeAll();
			if (leakflag == 1) {
				cout << "Canoe has at least one leak point\n\n\n";
				return 1;
			}

			cout << "Canoe Data\n\n";

			cout << "Input Values\n\n";
			cout << "Length: " << c.length << "\nNominal Beam: " << c.width << "\nDepth: " << c.depth << "\nNominal Flare Angle: " << c.flare <<
				"\nLength to Deepest Point: " << c.ldeepest << "\nLength to First Paddler: " << c.lfirst << "\nLength of Paddlers' Box: " << c.lpaddler <<
				"\nShape Parameter: " << c.shapeparam << "\nBow Rocker: " << c.brocker << "\nStern Rocker: " << c.srocker << "\nArea Density: " << c.density << "\n\n";

			cout << "Load Cases\n\t\t1\t\t2\n";
			cout << "Number of paddlers\n";
			cout << "\t\t" << c.loadcase[0].numpaddlers << "\t\t" << c.loadcase[1].numpaddlers << endl;
			cout << "Individual paddler weight\n";
			cout << "\t\t" << c.loadcase[0].paddlerweight << "\t\t" << c.loadcase[1].paddlerweight << endl;
			cout << "Paddler centroid (above knees)\n";
			cout << "\t\t" << c.loadcase[0].paddlercm << "\t\t" << c.loadcase[1].paddlercm << endl;

			cout << endl;


			cout << "Output Values\n\n";
			cout << "Surface Area: " << c.surfacearea << "\nWeight: " << c.weight << "\nLCM: " << c.CMx << "\nActual Beam: " << c.actualbeam << "\n\n";

			cout << "Load Case\n";
			cout << "\t\t1\t\t2\n";

			int i,j;
			for (j = 0; j < 8; j++) {
				cout << outputs[j] << endl << "\t\t";
				cout << c.outputs[j+1] << "\t\t" << c.outputs[8+j+1] << iostream::left << endl;
			}
			cout << "\n\n";
			break;
		case 2:
			char buffer[128];
			cout << "Input mesh filename: ";
			cin >> buffer;
			cout << "\n\n\n";
			output.open(buffer);
			c.OutputMesh(output);
			break;
		case 3:
			cout << "Input velocity: ";
			cin >> speed;
			cout << "Input freeboard: ";
			cin >> freeboard;
			cout << "\n\n";
			cout << "Resistance: " << c.GetFriction(freeboard,speed) << "\n\n\n";
			break;
		case 4:
			cout << "Input loadcase: ";
			cin >> loadcasenum;
			if (loadcasenum != 0 && loadcasenum != 1) {
				cout << "\n\n";
				break;
			}
			cout << "Input number of paddlers: ";
			cin >> c.loadcase[loadcasenum].numpaddlers;
			cout << "Input individual paddler weight: ";
			cin >> c.loadcase[loadcasenum].paddlerweight;
			cout << "Input paddler centre of mass: ";
			cin >> c.loadcase[loadcasenum].paddlercm;
			cout << "\n\n";
			break;
		case 5:
			cout << "Input tolerance for waterline calculations [" << c.wltol << "]: ";
			cin >> c.wltol;
			cout << "Input tolerance for tipping angle calculations [" << c.tippingtol << "]: ";
			cin >> c.tippingtol;
			cout << "Input tolerance for bow entrance angle calculations [" << c.bowtol << "]: ";
			cin >> c.bowtol;
			cout << "Input maximum number of iterations for waterline calculations [" << c.wlmaxit << "]: ";
			cin >> c.wlmaxit;
			cout << "Input maximum number of iterations for tipping angle calculations [" << c.tippingmaxit << "]: ";
			cin >> c.tippingmaxit;
			cout << "\n\n";
			break;
		default:
			c.Destruct();
			flag = 1;
			break;
		}
	}

	return 0;
}

int UIBulk() {
	cout << "(1) Analyze from an existing input setup file\n(2) Create a new input setup file and analyze from it\n(3) Back to main menu\n\n";

	int numeric;
	cout << "> ";
	cin >> numeric;
	cout << "\n\n";
	ifstream input;
	ofstream writeinput;
	double n1,n2,n3;
	ofstream inputfile;
	char buffer[128];


	switch(numeric) {
	case 1:		
		cout << "Input filename: ";
		cin >> buffer;
		cout << "\n";
		input.open(buffer);
		if (input.fail()) {
			cout << "File not found\n\n\n";
			return 1;
		}
		break;
	case 2:
		cout << "Input input setup filename: ";
		cin >> buffer;
		cout << "\n\n";
		writeinput.open(buffer);

		cout << "Input all fields in the format [Minimum Value] [Maximum Value] [Number of Values]\n";

		writeinput << 14 << endl;
		cout << "Length: ";	cin >> n1 >> n2 >> n3;
		writeinput << "Length\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "Lp: ";	cin >> n1 >> n2 >> n3;
		writeinput << "Lp\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "Ld: ";	cin >> n1 >> n2 >> n3;
		writeinput << "Ld\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "Lf: ";	cin >> n1 >> n2 >> n3;
		writeinput << "Lf\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "W: ";	cin >> n1 >> n2 >> n3;
		writeinput << "W\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "t1: ";	cin >> n1 >> n2 >> n3;
		writeinput << "t1\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "t2: ";	cin >> n1 >> n2 >> n3;
		writeinput << "t2\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "d: ";	cin >> n1 >> n2 >> n3;
		writeinput << "d\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		writeinput << "h\t" << 0 << '\t' << 0 << '\t' << 1 << endl;
		cout << "b: ";	cin >> n1 >> n2 >> n3;
		writeinput << "b\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "s: ";	cin >> n1 >> n2 >> n3;
		writeinput << "s\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "f: ";	cin >> n1 >> n2 >> n3;
		writeinput << "f\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "n: ";	cin >> n1 >> n2 >> n3;
		writeinput << "n\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		cout << "Density: ";	cin >> n1 >> n2 >> n3;
		writeinput << "Density\t" << n1 << '\t' << n2 << '\t' << n3 << endl;
		writeinput.close();
		input.open(buffer);
		break;
	default:
		return 1;
	}

	// Do the bulk analysis

	cout << "\nInput input table filename: ";
	cin >> buffer;
	cout << "\n";
	inputfile.open(buffer);

	CreateInputs(input,inputfile);
	inputfile.close();
	
	ifstream inputtable;
	inputtable.open(buffer);

	cout  << "Input output table filename: ";
	cin >> buffer;
	cout << "\n";
	ofstream resultsout;
	resultsout.open(buffer);

	int numcanoes;
	inputtable >> numcanoes;			// first line of file
	int i;
	for (i=0; i < 15; i++) {
		inputtable >> buffer;
	}

	Canoe c;
	int flag;
	for (i=0; i < numcanoes; i++) {
		flag = 0;
		double canoenum,L,Lp,Ld,Lf,W,t1,t2,d,h,b,s,f,n,dens;
		inputtable >> canoenum >> L >> Lp >> Ld >> Lf >> W >> t1 >> t2 >> d >> h >> b >> s >> f >> n >> dens;
		c.InitializeCanoe(L,Lp,Ld,Lf,W,t1,t2,d,h,b,s,f,n,dens);
		flag = c.AnalyzeAll();
		resultsout << canoenum;

		if (flag == 1) {
			resultsout << '\t' << "Analysis Failed\n";
			cout << "Canoe " << canoenum << " analysis failed.\n";
		} else {
			for (int j=0; j<17; j++) {
				resultsout << '\t' << c.outputs[j];
			}
			resultsout << '\t' << Score(c.outputs,targets[0],targets[1],targets[2],17) << endl;
			cout << "Canoe " << canoenum << " successfully analyzed.\n";
	//		cout << canoenum << '\t' << Score(c.outputs,targets[0],targets[1],targets[2],17) << endl;
		}

		c.Destruct();
	}

	cout << "\n\n";


	return 0;
}

