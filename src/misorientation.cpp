#include "misorientation.h"

MisorientationHdl::MisorientationHdl(void) {
	firstMisAll = NULL;
	lastMisAll = NULL;
	firstMisDist = NULL;
	c_value = 8; //Settings::C_VALUE
	a_value = 3; //Settings::A_VALUE
	idxMatrix = NULL;
	vectorMatrix = NULL;
	idxMatrixBravais = NULL;
	createIndexMatrixBravais();
}

MisorientationHdl::~MisorientationHdl(void) {
	free(idxMatrix);
	free(vectorMatrix);
	free(idxMatrixBravais);
}

Misori::Misori(misHdlP own, Real xm, Real ym, Real angle, int h = 0, int k = 0,
		int i = 0, int l = 0, int H = 0, int K = 0, int L = 0, Real U = 0,
		Real V = 0, Real W = 0) {
	theta = angle;
	axis[0] = h;
	axis[1] = k;
	axis[2] = i;
	axis[3] = l;
	miller[0] = H;
	miller[1] = K;
	miller[2] = L;
	u = U;
	v = V;
	w = W;
	x = xm;
	y = ym;
	owner = own;

	if (own->firstMisAll)
		own->firstMisAll->prev = this;

	next = own->firstMisAll;
	prev = NULL;
	own->firstMisAll = this;

}

Misori::Misori(void) {
	theta = 0;
	axis[0] = 0;
	axis[1] = 0;
	axis[2] = 0;
	axis[3] = 0;
}

int MisorientationHdl::readData(char *filename, int *XCells, int *YCells) {
	int N;
	Real stepsiz;
	Real X, Y;
	Real max1, max2, max3;

	max1 = max2 = max3 = 0;

	FILE *data = fopen(filename, "r");
	char *buf = (char *) calloc(BUFSIZ, sizeof(char));

	if (!data)
		exitus("Cannot open the input file");

	fgets(buf, BUFSIZ, data);
	sscanf(buf, "%d%lf%lf%lf\n", &N, &stepsiz, &X, &Y);

	this->AllOri = (OriP) calloc(N, sizeof(Ori));
	if (!AllOri)
		exitus(
				"Cannot allocate enough memory for your project; dynamic allocation might be necessary");

	this->stepsize = stepsiz;

	*XCells = ((int) (X / stepsiz + 1.1));
	*YCells = ((int) (Y / stepsiz + 1.1));

	assert(*XCells * *YCells == N);

	int index = 0;

	while (fgets(buf, BUFSIZ, data) != 0) {
		int idx;

		Real x, y, val1, val2, val3;
		euler angles;

		if (!(sscanf(buf, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", &idx, &x, &y, &val1,
				&val2, &val3)))
			exitus("Something wrong with the input file");

		AllOri[index].p1 = val1 * _PI_ / 180;
		AllOri[index].P = val2 * _PI_ / 180;
		AllOri[index].p2 = val3 * _PI_ / 180;
		AllOri[index].x = x + 0.5 * stepsiz;
		AllOri[index].y = y + 0.5 * stepsiz;
		if (val1 > max1)
			max1 = val1;
		if (val2 > max2)
			max2 = val2;
		if (val3 > max3)
			max3 = val3;

		if (index != idx - 1)
			exitus("Indices are different");
		index++;

	}
	fclose(data);
	printf("p1:%lf P:%lf p2:%lf\n", max1, max2, max3);

	return N;
}
double MisorientationHdl::calculateMisorientation_hexagonal(double* p,
		double* q) {

	/*double p1 = oria.p1*_PI_/180;   double p12 = orib.p1*_PI_/180;     //Degrees
	 double t = oria.P*_PI_/180;     double t2 = orib.P*_PI_/180;
	 double p2 = oria.p2*_PI_/180;   double p22 = orib.p2*_PI_/180;*/

	// Quaternions from Euler angles
	//double p[4] = {co1*cos((p1+p2)/2),s1*cos((p1-p2)/2),s1*sin((p1-p2)/2),co1*sin((p1+p2)/2)};
	//double q[4] = {co2*cos((p12+p22)/2),s2*cos((p12-p22)/2),s2*sin((p12-p22)/2),co2*sin((p12+p22)/2)};

	double qm1[4]; //Inverse of quaternion q

	for (int i = 0; i < 4; i++)
		qm1[i] = q[i];

	qm1[0] *= -1; //Inverting unit quaternion; yes the inverse of a unit quaternion is a simple operation

	double r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	r[0] = p[0] * qm1[0] - p[1] * qm1[1] - p[2] * qm1[2] - p[3] * qm1[3];
	r[1] = p[1] * qm1[0] + p[0] * qm1[1] - p[2] * qm1[3] + p[3] * qm1[2];
	r[2] = p[2] * qm1[0] + p[0] * qm1[2] - p[3] * qm1[1] + p[1] * qm1[3];
	r[3] = p[3] * qm1[0] + p[0] * qm1[3] - p[1] * qm1[2] + p[2] * qm1[1];

	//Now, we have to determine the smallest angle.

	double r0[6][2]; //There are 12 possible angles

	double a, b, c, d;
	double rt3 = sqrt(3.0);

	a = r[0];
	b = r[1];
	c = r[2];
	d = r[3];

	r0[0][0] = a;
	r0[0][1] = d;
	r0[1][0] = b;
	r0[1][1] = c;
	r0[2][0] = 0.5 * (a + rt3 * d);
	r0[2][1] = 0.5 * (d - rt3 * a);
	r0[3][0] = 0.5 * (b - rt3 * c);
	r0[3][1] = 0.5 * (c + rt3 * b);
	r0[4][0] = 0.5 * (a - rt3 * d);
	r0[4][1] = 0.5 * (d + rt3 * a);
	r0[5][0] = 0.5 * (b + rt3 * c);
	r0[5][1] = 0.5 * (c - rt3 * b);

	double omega = 0.0;

	int gp = -1; //General position
	int sp = -1; //specific position

	// The component with the maximal value is determined
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 2; j++)
			if (fabs(r0[i][j]) > omega) {
				omega = fabs(r0[i][j]);
				gp = i;
				sp = j;
			}

	if (gp < 0 || sp < 0 || gp > 5 || sp > 1) {
		printf("Couldn't find the misorientation angle\n");
	}

	double eqq[6][4];

	// Once this component is known, it is possible to determine the 6 possible quaternions which
	// corresponds to the Disorientation

	int ip = gp + 1;
	int iq = gp + 3;
	int ir = gp + 5;
	int isp; //inverse specific position

	if (sp == 0)
		isp = 1;
	else
		isp = 0;
	if (ip >= 6)
		ip -= 6;
	if (iq >= 6)
		iq -= 6;
	if (ir >= 6)
		ir -= 6;

	//These quaternions are here defined.

	eqq[0][0] = r0[gp][sp];
	eqq[0][1] = r0[ip][0];
	eqq[0][2] = r0[ip][1];
	eqq[0][3] = r0[gp][isp];
	eqq[1][0] = r0[gp][sp];
	eqq[1][1] = r0[ip][1];
	eqq[1][2] = r0[ip][0];
	eqq[1][3] = r0[gp][isp];
	eqq[2][0] = r0[gp][sp];
	eqq[2][1] = r0[iq][0];
	eqq[2][2] = r0[iq][1];
	eqq[2][3] = r0[gp][isp];
	eqq[3][0] = r0[gp][sp];
	eqq[3][1] = r0[iq][1];
	eqq[3][2] = r0[iq][0];
	eqq[3][3] = r0[gp][isp];
	eqq[4][0] = r0[gp][sp];
	eqq[4][1] = r0[ir][0];
	eqq[4][2] = r0[ir][1];
	eqq[4][3] = r0[gp][isp];
	eqq[5][0] = r0[gp][sp];
	eqq[5][1] = r0[ir][1];
	eqq[5][2] = r0[ir][0];
	eqq[5][3] = r0[gp][isp];

	int isst = -1; //index of the quaternion in the standard stereographic triangle

	for (int i = 0; i <= 5; i++) {
		// The quaternion in the SST must meet these requirements.
		a = fabs(eqq[i][0]);
		b = fabs(eqq[i][1]);
		c = fabs(eqq[i][2]);
		d = fabs(eqq[i][3]);
		if (a >= b && b >= rt3 * c && rt3 * c >= 0 && a >= 0.5 * (rt3 * b + c)
				&& a >= (2 + rt3) * d && (2 + rt3) * d >= 0) {
			int ib, id;

			ib = sp;
			if (i % 2)
				id = 1;
			else
				id = 0;

			if (ib == id)
				isst = i; //The angle is determined by comparison with a Databank of Indices
		}
	}

	Misori nmis;

	if (isst < 0) {
		//printf("Couldn't find an axis in the SST, exiting program\n");
		//restart=0;
		//break;

		//Maybe useful for debugging. In the actual programm if a desorientation cannot be calculated
		//it is set to 0 and therefore ignored

	} else {
		//MisoriP nmis = new Misori( 0,0,0,0,0 );
		determineAngleAxis(eqq[isst], &nmis);
	}

	//std::cout << "misorientation: " << nmis.theta << "\n";

	return nmis.theta;

}

void MisorientationHdl::calculateMisorientation(Ori oria, Ori orib, MisoriP mis) {

	double p1 = oria.p1;
	double p12 = orib.p1; //Radians
	double t = oria.P;
	double t2 = orib.P;
	double p2 = oria.p2;
	double p22 = orib.p2;

	/*double p1 = oria.p1*_PI_/180;   double p12 = orib.p1*_PI_/180;     //Degrees
	 double t = oria.P*_PI_/180;     double t2 = orib.P*_PI_/180;
	 double p2 = oria.p2*_PI_/180;   double p22 = orib.p2*_PI_/180;*/

	double co1 = cos(t / 2);
	double co2 = cos(t2 / 2);
	double s1 = sin(t / 2);
	double s2 = sin(t2 / 2);

	// Quaternions from Euler angles
	double p[4] = { co1 * cos((p1 + p2) / 2), s1 * cos((p1 - p2) / 2), s1
			* sin((p1 - p2) / 2), co1 * sin((p1 + p2) / 2) };
	double q[4] = { co2 * cos((p12 + p22) / 2), s2 * cos((p12 - p22) / 2), s2
			* sin((p12 - p22) / 2), co2 * sin((p12 + p22) / 2) };

	double qm1[4]; //Inverse of quaternion q

	for (int i = 0; i < 4; i++)
		qm1[i] = q[i];

	qm1[0] *= -1; //Inverting unit quaternion; yes the inverse of a unit quaternion is a simple operation

	double r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	r[0] = p[0] * qm1[0] - p[1] * qm1[1] - p[2] * qm1[2] - p[3] * qm1[3];
	r[1] = p[1] * qm1[0] + p[0] * qm1[1] - p[2] * qm1[3] + p[3] * qm1[2];
	r[2] = p[2] * qm1[0] + p[0] * qm1[2] - p[3] * qm1[1] + p[1] * qm1[3];
	r[3] = p[3] * qm1[0] + p[0] * qm1[3] - p[1] * qm1[2] + p[2] * qm1[1];

	//Now, we have to determine the smallest angle.

	double r0[6][2]; //There are 12 possible angles

	double a, b, c, d;
	double rt3 = sqrt(3.0);

	a = r[0];
	b = r[1];
	c = r[2];
	d = r[3];

	r0[0][0] = a;
	r0[0][1] = d;
	r0[1][0] = b;
	r0[1][1] = c;
	r0[2][0] = 0.5 * (a + rt3 * d);
	r0[2][1] = 0.5 * (d - rt3 * a);
	r0[3][0] = 0.5 * (b - rt3 * c);
	r0[3][1] = 0.5 * (c + rt3 * b);
	r0[4][0] = 0.5 * (a - rt3 * d);
	r0[4][1] = 0.5 * (d + rt3 * a);
	r0[5][0] = 0.5 * (b + rt3 * c);
	r0[5][1] = 0.5 * (c - rt3 * b);

	double omega = 0.0;

	int gp = -1; //General position
	int sp = -1; //specific position

	// The component with the maximal value is determined
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 2; j++)
			if (fabs(r0[i][j]) > omega) {
				omega = fabs(r0[i][j]);
				gp = i;
				sp = j;
			}

	if (gp < 0 || sp < 0 || gp > 5 || sp > 1) {
		printf("Couldn't find the misorientation angle\n");
	}

	double eqq[6][4];

	// Once this component is known, it is possible to determine the 6 possible quaternions which
	// corresponds to the Disorientation

	int ip = gp + 1;
	int iq = gp + 3;
	int ir = gp + 5;
	int isp; //inverse specific position

	if (sp == 0)
		isp = 1;
	else
		isp = 0;
	if (ip >= 6)
		ip -= 6;
	if (iq >= 6)
		iq -= 6;
	if (ir >= 6)
		ir -= 6;

	//These quaternions are here defined.

	eqq[0][0] = r0[gp][sp];
	eqq[0][1] = r0[ip][0];
	eqq[0][2] = r0[ip][1];
	eqq[0][3] = r0[gp][isp];
	eqq[1][0] = r0[gp][sp];
	eqq[1][1] = r0[ip][1];
	eqq[1][2] = r0[ip][0];
	eqq[1][3] = r0[gp][isp];
	eqq[2][0] = r0[gp][sp];
	eqq[2][1] = r0[iq][0];
	eqq[2][2] = r0[iq][1];
	eqq[2][3] = r0[gp][isp];
	eqq[3][0] = r0[gp][sp];
	eqq[3][1] = r0[iq][1];
	eqq[3][2] = r0[iq][0];
	eqq[3][3] = r0[gp][isp];
	eqq[4][0] = r0[gp][sp];
	eqq[4][1] = r0[ir][0];
	eqq[4][2] = r0[ir][1];
	eqq[4][3] = r0[gp][isp];
	eqq[5][0] = r0[gp][sp];
	eqq[5][1] = r0[ir][1];
	eqq[5][2] = r0[ir][0];
	eqq[5][3] = r0[gp][isp];

	int isst = -1; //index of the quaternion in the standard stereographic triangle

	for (int i = 0; i <= 5; i++) {
		// The quaternion in the SST must meet these requirements.
		a = fabs(eqq[i][0]);
		b = fabs(eqq[i][1]);
		c = fabs(eqq[i][2]);
		d = fabs(eqq[i][3]);
		if (a >= b && b >= rt3 * c && rt3 * c >= 0 && a >= 0.5 * (rt3 * b + c)
				&& a >= (2 + rt3) * d && (2 + rt3) * d >= 0) {
			int ib, id;

			ib = sp;
			if (i % 2)
				id = 1;
			else
				id = 0;

			if (ib == id)
				isst = i; //The angle is determined by comparison with a Databank of Indices
		}
	}

	Misori nmis;

	if (isst < 0) {
		//printf("Couldn't find an axis in the SST, exiting program\n");
		//restart=0;
		//break;

		//Maybe useful for debugging. In the actual programm if a desorientation cannot be calculated
		//it is set to 0 and therefore ignored

	} else {
		//MisoriP nmis = new Misori( 0,0,0,0,0 );
		determineAngleAxis(eqq[isst], &nmis);
	}

	mis->theta = nmis.theta;
	mis->axis[0] = nmis.axis[0];
	mis->axis[1] = nmis.axis[1];
	mis->axis[2] = nmis.axis[2];
	mis->axis[3] = nmis.axis[3];
	mis->miller[0] = nmis.miller[0];
	mis->miller[1] = nmis.miller[1];
	mis->miller[2] = nmis.miller[2];
	mis->u = nmis.u;
	mis->v = nmis.v;
	mis->w = nmis.w;

}

void MisorientationHdl::determineAngleAxis(Real *quaternion, MisoriP misori) {
	double _sr3 = 1 / SQRT3;
	double c_ah = c_value / a_value;

	double a = fabs(quaternion[0]);
	double b = fabs(quaternion[1]);
	double c = fabs(quaternion[2]);
	double d = fabs(quaternion[3]);
	double omega = a;

	double maxC = 0.0; //********I think this is not necessary but cannot remember; need verification
	if (b >= c && b >= d)
		maxC = b;
	if (c >= b && c >= d)
		maxC = c;
	if (d >= b && d >= c)
		maxC = d;

	double u = b;
	double v = c;
	double w = d;

	if (maxC == 0)
		maxC = 1.0;

	//assert( maxC>0 );

	u = u / maxC;
	v = v / maxC;
	w = w / maxC;

	if (a > 0.9999) // for same neighboring orientation
	{
		a = 0.9999;
	}

	//*******Axis in orthogonal coordinate system is here first normalized

	double fact = 1 / sqrt((1 - SQR(a)));
	double Uc = fact * b;
	double Vc = fact * c;
	double Wc = fact * d;

	//******Transformation to hexagonal coordinate system

	double Uh = Uc + Vc * _sr3;
	double Vh = 2 * _sr3 * Vc;
	double Wh = Wc / c_ah;
	/*double UVWh = SQR(Uh)-Uh*Vh+SQR(Vh)+SQR(c_ah*Wh);
	 assert( SQR(UVWh-1)<TOL );*/

	//******Multiplying by reciprocal basis of hexagonal lattice

	double Uhr = Uh * 2 * _sr3 / a_value;
	double Vhr = Vh * 2 * _sr3 / a_value;
	double Whr = Wh / c_value;

	double _mag = 1 / sqrt(SQR(Uhr) + SQR(Vhr) + SQR(Whr));

	Uhr *= _mag; //Normalizing
	Vhr *= _mag;
	Whr *= _mag;

	/*	double UU = ( 2 * Uhr - Vhr ) / 3.0; //Vector form of Miller-Bravais; they are not used; here only for Debugging
	 double VV = ( 2 * Vhr - Uhr ) / 3.0;
	 double LL = -( Uhr + Vhr ) / 3.0;
	 double WW = Whr;

	 double _magh = 1/sqrt(SQR(UU)+SQR(VV)+SQR(WW));

	 UU *= _magh;
	 VV *= _magh;
	 LL *= _magh;
	 WW *= _magh;     */

	int miller = 0;
	int mBravais = 0;
	double maxs = 0.0;

	//*********Calculating Index with smallest deviation from real vector from the set generated previously
	for (int i = indexes - 3; i >= 0; i -= 3) {
		double sinus = Uhr * vectorMatrix[i] + Vhr * vectorMatrix[i + 1] + Whr
				* vectorMatrix[i + 2];
		if (sinus >= maxs) {
			maxs = sinus;
			miller = i;
			mBravais = miller * 4 / 3;
		}
	}

	assert(miller % 3 == 0);

	int h = idxMatrix[miller];
	int k = idxMatrix[miller + 1];
	int l = idxMatrix[miller + 2];

	int hhh = idxMatrixBravais[mBravais];
	int kkk = idxMatrixBravais[mBravais + 1];
	int iii = idxMatrixBravais[mBravais + 2];
	int lll = idxMatrixBravais[mBravais + 3];

	/*	double hhex = h + k * _sr3;
	 double khex = 2 * _sr3 * k;
	 double lhex = l / c_value;     */

	assert(omega <= 1.000001);

	if (a > 1.0)
		a = 1.0;
	if (a < -1.0)
		a = -1.0;

	omega = 2 * acos(a);

	misori->theta = omega;
	misori->axis[0] = hhh;
	misori->axis[1] = kkk;
	misori->axis[2] = iii;
	misori->axis[3] = lll;
	misori->miller[0] = h;
	misori->miller[1] = k;
	misori->miller[2] = l;
	misori->u = vectorMatrix[miller];
	misori->v = vectorMatrix[miller + 1];
	misori->w = vectorMatrix[miller + 2];
}

void MisorientationHdl::createIndexMatrixBravais(void) {

	Real _sr3 = 1 / SQRT3;

	//Real c_a = c_value/a_value;
	int maxIndex = 3;

	idxMatrixBravais = (int *) calloc(4 * CUBE(218), sizeof(int));
	vectorMatrix = (double *) calloc(3 * CUBE(218), sizeof(double));
	idxMatrix = (int *) calloc(3 * CUBE(218), sizeof(int));

	if (!idxMatrix || !vectorMatrix || !idxMatrixBravais)
		exitus(
				"Cannot allocate enough memory\nPlease select another maximal index\n");

	//Indices are supossed to be already in hexagonal coordinate system. Therefore a transformation is not neccessary
	//However, they have to be multiplied by the reciprocal hexagonal basis and normalized.

	char flag = 0;
	int l = 0;
	int m, i, j, k;
	int n = 0;

	for (i = -maxIndex; i <= maxIndex; i++) //-maxIndex
		for (j = -maxIndex; j <= maxIndex; j++)
			for (k = -maxIndex; k <= maxIndex; k++) {
				flag = 0;

				if (fabs((double) -(i + j)) > 3)
					continue;

				for (m = 2; m <= maxIndex; m++) {
					if ((i % m == 0 && j % m == 0 && (-(i + j)) % m == 0 && k
							% m == 0)) {
						flag = 1;
					}
				}

				if (!flag) {

					double x = (i - j) * 2 * _sr3 / a_value; // Hexagonal Basis
					double y = -(i + 2 * j) * 2 * _sr3 / a_value;
					double z = k / c_value;

					double _mag = 1 / sqrt((SQR(x) + SQR(y) + SQR(z))); //Magnitude

					idxMatrixBravais[n] = i;
					n++; //Miller-Bravais
					idxMatrixBravais[n] = -(i + j);
					n++; //K is defined by H and L; makes the programming a bit easier
					idxMatrixBravais[n] = j;
					n++; //I=-(H+K) and restricted to maxIndex=3
					idxMatrixBravais[n] = k;
					n++;

					idxMatrix[l] = i - j; // Miller
					vectorMatrix[l] = x * _mag; //Normalizing
					l++;

					idxMatrix[l] = -(i + 2 * j);
					vectorMatrix[l] = y * _mag;
					l++;

					idxMatrix[l] = k;
					vectorMatrix[l] = z * _mag;
					l++;
				}
			}
	this->indexes = l;
}

void MisorientationHdl::determineGrainBoundaries(int N, int XCells, int YCells) {
	Ori oria, orib, oric;

	for (int j = 0; j < YCells; j++)
		for (int i = 0; i < XCells; i++) {
			int idxCell = i + XCells * j;

			oria.p1 = AllOri[idxCell].p1;
			oria.P = AllOri[idxCell].P;
			oria.p2 = AllOri[idxCell].p2;
			oria.x = AllOri[idxCell].x;
			oria.y = AllOri[idxCell].y;

			if (i < XCells - 1) {
				orib.p1 = AllOri[idxCell + 1].p1;
				orib.P = AllOri[idxCell + 1].P;
				orib.p2 = AllOri[idxCell + 1].p2;
				orib.x = AllOri[idxCell + 1].x;
				orib.y = AllOri[idxCell + 1].y;

				Real xm = 0.5 * (oria.x + orib.x);
				Real ym = 0.5 * (oria.y + orib.y);

				double dp1_2 = SQR( oria.p1 - orib.p1 );
				double dP_2 = SQR( oria.P - orib.P );
				double dp2_2 = SQR( oria.p2 - orib.p2 );

				Misori misCalc;

				if (dp1_2 > TOLANGLE2 && dP_2 > TOLANGLE2 && dp2_2 > TOLANGLE2)
					calculateMisorientation(oria, orib, &misCalc);

				if (misCalc.theta > MINMIS) {
					MisoriP misp = new Misori(this, xm, ym, misCalc.theta,
							misCalc.axis[0], misCalc.axis[1], misCalc.axis[2],
							misCalc.axis[3], misCalc.miller[0],
							misCalc.miller[1], misCalc.miller[2], misCalc.u,
							misCalc.v, misCalc.w);
					if (!misp)
						exitus("Cannot allocate more memory");
				}

			}

			if (j >= 1) {
				oric.p1 = AllOri[idxCell - XCells].p1;
				oric.P = AllOri[idxCell - XCells].P;
				oric.p2 = AllOri[idxCell - XCells].p2;
				oric.x = AllOri[idxCell - XCells].x;
				oric.y = AllOri[idxCell - XCells].y;

				Real xm = 0.5 * (oria.x + orib.x);
				Real ym = 0.5 * (oria.y + orib.y);

				double dp1_2 = SQR( oria.p1 - oric.p1 );
				double dP_2 = SQR( oria.P - oric.P );
				double dp2_2 = SQR( oria.p2 - oric.p2 );

				Misori misCalc;

				if (dp1_2 > TOLANGLE2 && dP_2 > TOLANGLE2 && dp2_2 > TOLANGLE2)
					calculateMisorientation(oria, oric, &misCalc);

				if (misCalc.theta > MINMIS) {
					MisoriP misp = new Misori(this, xm, ym, misCalc.theta,
							misCalc.axis[0], misCalc.axis[1], misCalc.axis[2],
							misCalc.axis[3], misCalc.miller[0],
							misCalc.miller[1], misCalc.miller[2], misCalc.u,
							misCalc.v, misCalc.w);
					if (!misp)
						exitus("Cannot allocate more memory");
				}

			}

		}

}

void exitus(char* s) {
	if (s)
		printf("%s\n", s);
	printf("Press enter to exit...\n");
	getchar();
	getchar();
	exit(1);
}

