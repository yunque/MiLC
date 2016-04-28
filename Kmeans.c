/* K-MEANS */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define DistanceDefault 1.7e308
#define ErrorDefault 1.7e308
#define ErrorThreshold 0.000000000001
#define MaxIterations 10000

// Global variables
int nK = 5;
int nFeatures = 2;
int nData = 10000;
//double errorOld=1000000;
double errorNew=999999;//((errorOld)-1);
double errorSum = 0;

// Global pointers
double ** coor; // ptrData
double * points; // ptrFt
double ** centroids;
int * klabel;
double ** distance;


void readInput()
{
	const char* filenameInput="input.csv";
	FILE *fin = fopen(filenameInput,"r");
	int line=0;
	double xCSV,yCSV;
	for (line=0;line<nData;line++)
	{
		fscanf(fin, "%lf, %lf", &xCSV, &yCSV);
		printf("Here 2 \n");
		coor[line]=(double *) malloc(nFeatures*sizeof(double)); //coor[line]==1-D nData
		coor[line][0] = xCSV;									//coor[line][ft]==2-D
		coor[line][1] = yCSV;
	}
	fclose(fin);
}


// INITIALIZE CENTROIDS
void initCentroids()
{
	int k=0;
	for (k=0;k<nK;k++)
	{
		centroids[k] = (double*)malloc(nFeatures*sizeof(double));
	}

	centroids[0][0]=-0.357;	//A
	centroids[0][1]=-0.253;
	centroids[1][0]=-0.055; //B
	centroids[1][1]=4.392;
	centroids[2][0]=2.674; 	//C
	centroids[2][1]=-0.001;
	centroids[3][0]=1.044; 	//D
	centroids[3][1]=-1.251;
	centroids[4][0]=-1.495;	//E
	centroids[4][1]=-0.090;
}


// E-STEP
/* Assign points to their nearest cluster */
double Estep()
{
	double errorOld=errorNew;
	errorNew=0;
	int line=0;
	for (line=0;line<nData;line++)	// for every point...
	{
		double distanceMin = DistanceDefault;
		distance[line] = (double*)malloc(nK*sizeof(double));
		int k=0;
		for (k=0;k<5;k++)			// ...& every cluster
		{
			distance[line][k]=0;
			int f=0;
			for (f=0;f<nFeatures;f++)
			{
				distance[line][k] += (coor[line][f]-centroids[k][f])*(coor[line][f]-centroids[k][f]);
			}

			if (distance[line][k]<distanceMin)
			{
				distanceMin = distance[line][k];
				klabel[line] = k;		// keep only minimum distance
			}
		}
		errorNew += sqrt(distanceMin);
	}
	return errorOld;
}


// M-STEP
/* Recalibrate cluster centroids */
void Mstep()
{
	double * DataSum = (double*)malloc(nFeatures*sizeof(double));
	int k=0;
	for (k=0;k<nK;k++)					// for each cluster...
	{
		DataSum[0]=0;
		DataSum[1]=0;
		int nClusterPoints = 0;
		int line=0;
		for (line=0;line<nData;line++)	// ...sum all points in that cluster
		{
			if (klabel[line]==k)
			{
				nClusterPoints++;
				int f=0;
				for (f=0;f<nFeatures;f++)
				{
					DataSum[f] += coor[line][f];
				}
			}
		}
		if (nClusterPoints>0)				// find new cluster centroids...
		{
			int f=0;
			for (f=0;f<nFeatures;f++)
			{
				centroids[k][f] = DataSum[f]/nClusterPoints;
				printf("DataSum[f]=%f, nData=%d\n",DataSum[f],nData);
			}
		}
		else								// ...or randomize if cluster is empty
		{
			int f=0;
			for (f=0;f<nFeatures;f++)
			{
				centroids[k][f] =rand()%5 -2.5;
			}
		}
	}
}


// CONVERGENCE?
int convergence(double errorOld)
{
	if ((errorOld-errorNew) > ErrorThreshold)
		{return 0;}
	else 					//error converges
		{return 1;}
}


// OUTPUT
void writeOutput(double error)
{

	const char *filenameOutput = "OUTPUT.TXT";
	const char *filenameComments = "OUTPUT.TXT";
	FILE *fout = fopen(filenameOutput,"w");

	//char *Kname = (char*)malloc(nK*sizeof(char*));
	const char *Kname[] = {"Adam","Bob","Charley","David","Edward"};

	fprintf(fout, "error = %.3f\n", error);
	int line=0;
	for (line=0;line<nData;line++)
		{fprintf(fout, "%s \n", Kname[klabel[line]]);}
	fclose(fout);
}


// MAIN
int main()
{

	int nIterations = 0;
	double errorOld=errorNew+1;
	srand(time(0));
	coor = (double**)malloc(nData*sizeof(double*));
	points = (double*)malloc(nFeatures*sizeof(double));
	distance = (double**)malloc(nData*sizeof(double*));
	centroids = (double**)malloc(nK*sizeof(double*));
	klabel = (int*)malloc(nData*sizeof(int));


	readInput();
 	printf("Here \n");
	initCentroids();

	//while (nIterations<MaxIterations)
	while (!convergence(errorOld) && (nIterations<MaxIterations))
	{
		nIterations++;
		printf("Iteration #%d\n",nIterations);
		errorOld = Estep();
		Mstep();
		printf("errorOld = %f\n",errorOld);
		printf("errorNew = %f\n",errorNew);
		printf("errorOld-errorNew = %f\n",(errorOld-errorNew));
	}
	if (nIterations>MaxIterations)
	{
		printf("Error: Maximum number of iterations has been reached.\n");
		printf("       Algorithm does not converge fast enough.\n");
	}
	writeOutput(errorOld);
}
