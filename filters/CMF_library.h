double getAntipodalPoint(double p);
double dist(double a, double b);
double getMB(double* y, int i, int j, int M, int N);
double getMBInd(int j, int N);

void checkDims(int M, int N, int R, int T);

void medfiltCirc2D(double* u, double* y, int M, int N, int r, int t);
void medfiltCirc2DQuant(double* u, double* y, int M, int N, int r, int t, double* v, int S);
void medfiltGeoRN2D(double* u, double* y, int M, int N, int K, int RD, int TD, int nIter, double stopTol);




