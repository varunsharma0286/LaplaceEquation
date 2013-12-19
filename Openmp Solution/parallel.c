#include<stdio.h>
#include<math.h>
#include "omp.h"

#define M 1
#define INTERVAL 0.01
#define INIT_GUESS 0.0//0.4
#define PI 3.14159265
#define E 2.71828
#define CONV .00001
#define TRUE 1
#define DEBUG 0

void checkSolution(double **buffer);

int main(int argc,char **argv)
{

	double input;
        if(argc < 2)
        {
                printf("\n Incorrect usage");
                printf("\n Correnct usage: <input>\n");
                return 1;
        }
        sscanf(argv[1],"%lf",&input);
        printf("\n gridsize = %lf",input);

	int i,j;
	int xDim = M/input;
	int yDim = M/input;
	long int steps = 0;
	float x = 0.0, y=0.0; 
	int reps = 0;
	double start, end;
	double exp = pow(E,-PI);
	double totalDiff = 0.0;
	double tempDiff = 0.0;
	double **matrix = (double **)malloc(sizeof(double*)*(yDim+1));
	double **matrix1 = (double **)malloc(sizeof(double*)*(yDim+1));
	#pragma omp parallel
	{	
		#pragma omp master
		{	
			printf("\nnum of threads = %d\n",omp_get_num_threads());
		}
	}
	//This loop will initialize all to an initial guess value.
	for(i=0;i<=xDim;i++)
	{
		matrix[i] = (double *)malloc(sizeof(double)*(xDim+1));
		matrix1[i] = (double *)malloc(sizeof(double)*(xDim+1));
		for(j=0;j<=xDim;j++)
		{
			if(y==0.0)
			{				
				matrix[i][j] =	sin(x*PI);
				matrix1[i][j] = sin(x*PI);
			}
			else if(y>=1.0)
			{
				matrix[i][j] = sin(x*PI)*exp;
				matrix1[i][j] = sin(x*PI)*exp;
			}
			else if(x==0.0||x>=1.0)
			{
				matrix[i][j] = 0;
				matrix[i][j] = 0;
			}	
			else
			{
				matrix[i][j] = INIT_GUESS;
				matrix1[i][j] = INIT_GUESS;
			}
			x = x+(float)input;
			if(x>=1.0)
				x = 1.0;
		}
		x=0;
		y = y+(float)input;
	}
	start = omp_get_wtime();
	while(TRUE)
	{
		#pragma omp parallel for reduction(+:totalDiff) private(j,tempDiff)
		for(i=1;i<xDim;i++)
		{
			for(j=1;j<xDim;j++)
			{
				matrix1[i][j] = (matrix[i-1][j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1])/4.0;
				tempDiff = (matrix1[i][j] - matrix[i][j]);
				totalDiff =totalDiff + (tempDiff* tempDiff);
			}
		}
		totalDiff = sqrt(totalDiff);

		if(totalDiff <= CONV)
		{
			printf("\n number of steps = %ld, grid size = %d",steps, xDim);
			break;
		}
		//printf("\n steps = %ld, totalDiff = %f",steps, totalDiff);
		totalDiff = 0.0;
		steps++;

		#pragma omp parallel for private(j)
                for(i=1;i<xDim;i++)
                {
			/*for(j=1;j<xDim;j++)
			{
				matrix[i][j] = matrix1[i][j];
			}*/
			memcpy(matrix[i], matrix1[i],sizeof(double)*xDim);
                }
	}
	end = omp_get_wtime();
	printf("\ntime taken = %f\n",(end-start));
	if(DEBUG)
		checkSolution(matrix);

        for(i=0;i<=xDim;i++)
        {
                free(matrix[i]);
                free(matrix1[i]);
        }
        free(matrix);
        free(matrix1);
	return 0;
}


void checkSolution(double **buffer)
{
	int x, y;
	while(1)
	{
		printf("\n Enter row = ");
		scanf("%d",&y);
                printf("\n Enter column = ");
                scanf("%d",&x);
		printf("matrix[%d][%d] = %lf",y,x,buffer[y][x]);
	}
}
