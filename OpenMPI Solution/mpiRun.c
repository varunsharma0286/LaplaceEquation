/*Solution Strategy
 * First create the initial boundary conditions
 * Then iteratively compute the value of each element using the formula
 * Need to maintain the image of the entire matrix in the root process
 * Gather the sum from all the individual processes using the mpi_allReduce;
 * Root will check for the convergence
 */

#include<stdio.h>
#include<math.h>
#include "mpi.h"

#define M 1
#define INTERVAL 0.1
#define INIT_GUESS 0.0
#define PI 3.14159265
#define E 2.71828
#define CONV .00001
#define TRUE 1
#define DEBUG 0

#define SEND_UP 1
#define SEND_LOW 1
#define RECV_UP 1
#define RECV_LOW 1

#define FIRSTP 0

void checkSolution(double **buffer);
void printInitCondition(double **buffer, int row, int column, int rank);

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

	int ierr, rank, size;
        MPI_Status status;
        MPI_Init(&argc,&argv);

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

	int i,j;
	int Dim = M/input;
	int last = 0;
	int RowPerProcessor;
	int column = Dim;
	double **matrix;
	double **matrix1;
	float x = 0.0, y=0.0;
	double *upRow;
	double *lowRow;
	double localDiff = 0.0;
	//long int totalSteps;
	long int localSteps = 0;
	long int globalSteps = 0;
	int startIdx = 0, endIdx;
	double startTime=0.0,endTime=0.0;
 	int LASTP = size -1; //Total number of processes

	double start, end;
	double exp = pow(E,-PI);
	double totalDiff = 0.0;
	double tempDiff = 0.0;

        if(rank == FIRSTP)
                printf("\n gridsize = %lf",input);

        if(Dim%size == 0)
        {
                RowPerProcessor = Dim/size;
                last = RowPerProcessor;
        }
        else
        {
                RowPerProcessor = Dim/size;
                last = (RowPerProcessor) + (Dim - (RowPerProcessor * size));
        }
	
	if(rank == LASTP)
	{
		RowPerProcessor = last;
	}

	if(DEBUG && (rank == 0))
	{
		printf("\n Total Rows = %d \n RowPerProcessor = %d\n last = %d\n",Dim,RowPerProcessor,last);
	}		

        matrix = (double **)malloc(sizeof(double*)*(RowPerProcessor));
	matrix1 = (double **)malloc(sizeof(double*)*(RowPerProcessor));

       	for(i=0;i<RowPerProcessor;i++)
       	{
               	matrix[i] = (double *)malloc(sizeof(double)*(column+1));
                matrix1[i] = (double *)malloc(sizeof(double)*(column+1));
       	        if(rank == FIRSTP)
               	{
               		for(j=0;j<=column;j++)
                	{
      	                	if(i==0)
               	        	{
                       	        	matrix[i][j] =  sin(x*PI);
                               		matrix1[i][j] = sin(x*PI);
                       		}
                        	else if(j==0||j==(column))
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
                       	x=0.0;
               	}
		else if(rank == LASTP)
		{
			for(j=0;j<=column;j++)
			{
	         		if(i==RowPerProcessor-1)
       	                       	{
               	               		matrix[i][j] = sin(x*PI)*exp;
                    	                matrix1[i][j] = sin(x*PI)*exp;
                      	       	}
                      		else if(j==0 || j == column)
                               	{
                                        matrix[i][j] = 0.0;
       	                                matrix1[i][j] = 0.0;
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
			x=0.0;
		}
               	else
               	{
               		for(j=0;j<=column;j++)
               		{
               			if(j==0 || j==column)
                		{
       	        			matrix[i][j] = 0.0;
               				matrix1[i][j] = 0.0;
               			}
               			else
                		{
       	            			matrix[i][j] = INIT_GUESS;
        	                	matrix1[i][j] = INIT_GUESS;
       	        	        }	
                	}	
       	        }
       	}

	if(DEBUG)
	{
		printInitCondition(matrix,RowPerProcessor,column,rank);
	}

	//Barrier here to make sure that all the processes are done with their initial condition
	MPI_Barrier(MPI_COMM_WORLD);

	//Now we have the initial condition done so start with the usual stuff.
	upRow = (double *)malloc(sizeof(double)*(column+1));
	lowRow = (double *)malloc(sizeof(double)*(column+1));

	if(rank == LASTP)
	{
		startIdx = 0;
		endIdx = RowPerProcessor-1;
	}
	else if(rank == FIRSTP)
	{
		startIdx = 1;
		endIdx = RowPerProcessor;	
	}
	else
	{
		startIdx =0;
		endIdx = RowPerProcessor; 
	}
	
	startTime = MPI_Wtime();
	while(1)
	{
		if(rank == FIRSTP)
		{
			
			MPI_Send(matrix[RowPerProcessor-1],column+1,MPI_DOUBLE,1,SEND_UP,MPI_COMM_WORLD);
			MPI_Recv(upRow,column+1,MPI_DOUBLE,1,RECV_UP,MPI_COMM_WORLD,&status);
		}
		else if(rank == LASTP)
		{
                        MPI_Send(matrix[0],column+1,MPI_DOUBLE,size-2,SEND_LOW,MPI_COMM_WORLD);
                        MPI_Recv(lowRow,column+1,MPI_DOUBLE,size-2,RECV_LOW,MPI_COMM_WORLD,&status);			
		}
		else
		{
                        MPI_Send(matrix[0],column+1,MPI_DOUBLE,rank-1,SEND_LOW,MPI_COMM_WORLD);
                        MPI_Send(matrix[RowPerProcessor-1],column+1,MPI_DOUBLE,rank+1,SEND_UP,MPI_COMM_WORLD);

                        MPI_Recv(lowRow,column+1,MPI_DOUBLE,rank-1,RECV_LOW,MPI_COMM_WORLD,&status);
                        MPI_Recv(upRow,column+1,MPI_DOUBLE,rank+1,RECV_UP,MPI_COMM_WORLD,&status);
		}
		//MPI_Barrier(MPI_COMM_WORLD);
		
                for(i=startIdx;i<endIdx;i++)
                {
                        for(j=1;j<column;j++)
                        {
				if(rank == FIRSTP && i == endIdx-1) //First process
				{
                                	matrix1[i][j] = (matrix[i-1][j] + upRow[j] + matrix[i][j-1] + matrix[i][j+1])/4.0;
				}
				else if(rank == LASTP && (i==startIdx)) //Last process
				{
					matrix1[i][j] = (matrix[i+1][j] + lowRow[j] + matrix[i][j-1] + matrix[i][j+1])/4.0;
				}
				else if((i==startIdx || i==endIdx-1) && (rank>FIRSTP && rank<LASTP))
				{
					if(i == startIdx)
					{
                     		        	matrix1[i][j] = (lowRow[j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1])/4.0;
					}
					else if(i==endIdx-1)
					{
						matrix1[i][j] = (matrix[i-1][j] + upRow[j] + matrix[i][j-1] + matrix[i][j+1])/4.0;
					}					
				}
				else	
					matrix1[i][j] = (matrix[i-1][j] + matrix[i+1][j] + matrix[i][j-1] + matrix[i][j+1])/4.0;

                                tempDiff = (matrix1[i][j] - matrix[i][j]);
                                localDiff =localDiff + (tempDiff* tempDiff);
                        }
                }	 
		//localSteps++;
		MPI_Allreduce(&localDiff,&totalDiff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		totalDiff = sqrt(totalDiff);
                if(DEBUG && rank == 0)
                        printf("\n 0-> steps = %ld, 0->totalDiff = %f", localSteps,totalDiff);

		if(totalDiff <= CONV)
		{
			endTime = MPI_Wtime();
			if(rank == 0)
			{
				printf("\n number of steps = %ld,totalDiff = %f, Time = %f",localSteps,totalDiff,(endTime-startTime));
			}
			break;
		}
		else
		{
	                totalDiff = 0.0;
        	        localDiff = 0.0;
			localSteps++;
		}

                for(i=0;i<RowPerProcessor;i++)
                {
			//printf("\n copy->Process = %d",rank);
			memcpy(matrix[i], matrix1[i],sizeof(double)*column+1);
                }
		//MPI_Barrier(MPI_COMM_WORLD);
	}

	if(DEBUG)
	{
		printf("\n rank = %d\n",rank);
		for(j=0;j<=column;j++)
               	{
                       	printf("%f ",lowRow[j]);
               	}
		printf("\n");
	}

	MPI_Finalize();
	return 0;
}

void printInitCondition(double **buffer, int row, int column, int rank)
{	
	int i,j;
	printf("\n Process number = %d\n", rank);
	for(i=0;i<row;i++)
	{
		for(j=0;j<=column;j++)
		{
			printf("%f ",buffer[i][j]);
		}
		printf("\n");
	}	
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
		printf("\n correct = %f", sin(.01*x*PI));//*pow(E,(-0.01*y*PI)));
        }
}

