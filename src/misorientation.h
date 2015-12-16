
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<iostream>



#define _PI_ 3.1415926535897932384626433832795
#define SQR(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define SQRT3 1.7320508075688772935274463415059
#define TOL 3e-3      //3e-4
#define TOLANGLE2 3e-5  //3e-6
#define MINMIS (0.0349) //0.052

typedef double Real;
typedef long windowRef;
typedef long fileRef;

class Misori;
class orientation;
class MisorientationHdl;

typedef MisorientationHdl misHdl;
typedef misHdl* misHdlP;
typedef Misori* MisoriP;
typedef orientation Ori;
typedef orientation* OriP;
typedef struct eulerAngles euler;
typedef euler eulerP;
typedef struct misOriAngleAxis misOriA_A;
typedef misOriA_A* misOriA_AP;

//extern Real a_value, c_value;
extern int *idxMatrix;

extern Real deviation;

void exitus( char* hint );

struct eulerAngles{
        Real p1;
        Real p2;
        Real P;
};

struct misOriAngleAxis{
        Real theta;
        int h;
        int k;
        int i;
        int l;
};

class MisorientationHdl{
        public:
                MisoriP firstMisAll;
                MisoriP lastMisAll;
                MisoriP firstMisDist;
                euler EulerAngles;

                int *idxMatrix;
                double *vectorMatrix;
                int *idxMatrixBravais;
                double c_value;
                double a_value;


                int indexes;
                OriP AllOri;

                OriP firstOri;
                OriP lastOri;

                MisorientationHdl( void );
                ~MisorientationHdl( void );

                int readData( char *filename, int *XCells, int *YCells );
                void determineGrainBoundaries( int N, int XCells, int YCells );
                void createIndexMatrixBravais( void );
                double calculateMisorientation_hexagonal( double* p , double* q );
                void calculateMisorientation( Ori oria, Ori orib, MisoriP mis );
                void determineAngleAxis(Real *quaternion, MisoriP misori);
                void determineMisorientationAccordingAxis( int hu, int ku, int iu, int lu );
                Real calculateMisoriFraction( misOriA_A a );

        private:

                long filePtr;
                long winPtr;
                Real stepsize;

};

struct orientation{

        Real p1;
        Real p2;
        Real P;
        Real x;
        Real y;
};

class Misori{

        public:
                Misori( misHdlP own, Real xm, Real ym, Real angle, int h, int k, int i, int l, int H, int K, int L, Real u, Real v, Real w );
                Misori( void );
                ~Misori( void ){};
                Real theta;
                int axis[4];
                Real miller[3];
                Real u;  //Vectors defining the axis already in hexagonal coordinate system
                Real v;
                Real w;
                MisoriP next;
                MisoriP prev;
                misHdlP owner;

                Real x;
                Real y;
};



