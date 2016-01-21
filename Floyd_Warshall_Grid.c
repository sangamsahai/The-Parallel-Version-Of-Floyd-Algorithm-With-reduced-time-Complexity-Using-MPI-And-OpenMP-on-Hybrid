/* This is Floyd Algorithm with time complexity of O(N log_2 P / sqrt(P)) 
This parallel Algorithm breaks the bigger adjacency matrix into smaller grids.
The grids are of size n/sqrt(p)*n/sqrt(p) .
Where n*n is the size of the bigger adjacency matrix and p is the total number of processes.
And finally , the computation is done on grids parallely  by different processes. 
*/
#include <stdio.h>
#include <mpi.h>
#include <limits.h>
#include <math.h>
#include <omp.h>

 MPI_Status status;
 int ROWS_IN_ORIGINAL_ARRAY=4;
 int COLS_IN_ORIGINAL_ARRAY=4;
 
 int min(int,int);
 int add(int,int);

int main (int argc, char **argv)
{

 int rank, size;
 MPI_Init (&argc, &argv); /* starts MPI */
 
 MPI_Comm_rank (MPI_COMM_WORLD, &rank); /* get current process id */
 MPI_Comm_size (MPI_COMM_WORLD, &size); /* get number of processes */
 
 //size = p =total number of processes
 //ROWS_IN_ORIGINAL_ARRAY = n = total number of rows in the adjacency matrix
 
 if(rank == 0)//Master Process code starts
 {
 //------Process to break the big matrix into grid chunks STARTS
 
 
 int size_of_grid_chunk=ROWS_IN_ORIGINAL_ARRAY / (sqrt((int)size));  //  g = n/sqrt(p)
 
 
//Initializing the Adjacency Matrix in the Master Process
int originalArray_W0[ROWS_IN_ORIGINAL_ARRAY][COLS_IN_ORIGINAL_ARRAY];  //Initializing the original Array

//Populating the original Adjacency Array in the Master Process

originalArray_W0[0][0]=0;
originalArray_W0[0][1]=INT_MAX;
originalArray_W0[0][2]=3;
originalArray_W0[0][3]=0;
originalArray_W0[1][0]=-2;
originalArray_W0[1][1]=0;
originalArray_W0[1][2]=INT_MAX;
originalArray_W0[1][3]=1;
originalArray_W0[2][0]=INT_MAX;
originalArray_W0[2][1]=INT_MAX;
originalArray_W0[2][2]=0;
originalArray_W0[2][3]=5;
originalArray_W0[3][0]=INT_MAX;
originalArray_W0[3][1]=4;
originalArray_W0[3][2]=INT_MAX;
originalArray_W0[3][3]=0;

int counter=0;
int i=0;
int j=0;
int k=0;
int g=size_of_grid_chunk; //g is the size of grid chunk

//sending the array chunks(grids) to the slave processes
 for( i=0;i<g;i++)
 {
 for( j=0;j<g;j++)
 {
 for( k=0;k<g;k++)
 {
 MPI_Send(&originalArray_W0[(i*g)+k][j*g],g,MPI_INT,counter,k,MPI_COMM_WORLD);
 }//end of k loop
 counter=counter+1;
 }//end of j loop
 }//end of i loop
 
 k=0;
 
 
//receiving the Master's chunk(grid)
int original_chunk[g][g]; 
for(i=0;i<g;i++)
{
MPI_Status status;
MPI_Recv(&original_chunk[i][0],2,MPI_INT,0,i,MPI_COMM_WORLD,&status);
}


//-----Process to break the big matrix into grid chunks ENDS
 
 //MAIN LOOP (in master) which runs k times where k is the number of vertices in the graph
 for(k=0;k<ROWS_IN_ORIGINAL_ARRAY;k++)  
 {
 
  //Sending the i rows 
  for(i=0;i<size;i++)
  {
  int i_by_g=i/g;
   for(j=i_by_g*g;j<(i_by_g*g)+g;j++)
		{
        MPI_Send(&originalArray_W0[j][0],COLS_IN_ORIGINAL_ARRAY,MPI_INT,i,(j-(i_by_g*g)),MPI_COMM_WORLD);//This will run g times for each process 		
		}   //end of j loop
  }//end of i loop


	//Receiving the i rows
    int i_rows_chunk[g][COLS_IN_ORIGINAL_ARRAY];   
    i=0;
    for(i=0;i<g;i++)
	{
	MPI_Status status;
	MPI_Recv(&i_rows_chunk[i][0],COLS_IN_ORIGINAL_ARRAY,MPI_INT,0,i,MPI_COMM_WORLD,&status);
	}
	
	
	//Sending the k row
	for(i=0;i<size;i++)
    {
    MPI_Send(&originalArray_W0[k][0],COLS_IN_ORIGINAL_ARRAY,MPI_INT,i,k,MPI_COMM_WORLD);//tag is k	
	}
    //Receiving k the row 
	int k_row[COLS_IN_ORIGINAL_ARRAY];
	MPI_Status status1;
	MPI_Recv(&k_row,COLS_IN_ORIGINAL_ARRAY,MPI_INT,0,k,MPI_COMM_WORLD,&status1);
	
    
	//Computation starts using MULTI THREADING
	omp_set_num_threads(g);
#pragma omp parallel
	{
	int ID = omp_get_thread_num();j=0;
	for(j=0;j<g;j++)
	{
	int absolute_i= (((int)rank/g)*g)+ID;
	int absolute_j= ((rank%g)*g)+j;
	original_chunk[ID][j]=min(original_chunk[ID][j], add(i_rows_chunk[ID][k],k_row[absolute_j]) ) ;
	}
	}
	
	
	//Sending the results to master to update
    for(i=0;i<g;i++)
	{
	MPI_Send(&original_chunk[i][0],g,MPI_INT,0,i,MPI_COMM_WORLD);
    }//end of i loop	
 
	//Receiving from slave processes and updating the original array
	for(i=0;i<size;i++)//i represents the process from which we receive
	{
	int start_i= ((int)i/g)*g;
	int start_j= (i%g)*g;
	for(j=start_i;j<start_i+g;j++)
	{
	MPI_Status status;
	MPI_Recv(&originalArray_W0[j][start_j],g,MPI_INT,i,(start_i-j),MPI_COMM_WORLD,&status);
	}//end of j loop
	}//end of i loop
 
 
 } //MAIN k loop ends
 
 
 printf("\n After the computation, the values of grid in process %d are  - \n",rank);
 
    for(i=0;i<g;i++)
	{
	for(j=0;j<g;j++)
	{
	printf(" %d ",original_chunk[i][j]);
	}
	printf("\n");
	}
    printf("\n\n");
 
 
 
 
 }//end of Master Process
 
 //******************************************************************************************************************
 
 
 else //Slave Process
 {
 
 int size_of_grid_chunk=ROWS_IN_ORIGINAL_ARRAY / (sqrt((int)size));  //  g = n/sqrt(p)
 int g=size_of_grid_chunk;

 //receiving the grid chunk from master
int original_chunk[g][g]; 
int i=0;int j=0;
for(i=0;i<g;i++)
{
MPI_Status status;
MPI_Recv(&original_chunk[i][0],2,MPI_INT,0,i,MPI_COMM_WORLD,&status);
}


int k=0;
//MAIN LOOP (in slave) which runs k times where k is the number of vertices in the graph
for(k=0;k<ROWS_IN_ORIGINAL_ARRAY;k++)  
 {
 
  //Receiving the i rows from master
  int i_rows_chunk[g][COLS_IN_ORIGINAL_ARRAY];   
  i=0;
  for(i=0;i<g;i++)
	{
	MPI_Status status;
	MPI_Recv(&i_rows_chunk[i][0],COLS_IN_ORIGINAL_ARRAY,MPI_INT,0,i,MPI_COMM_WORLD,&status);
	}

	//Receiving k the row 
	int k_row[COLS_IN_ORIGINAL_ARRAY];
	MPI_Status status1;
	MPI_Recv(&k_row,COLS_IN_ORIGINAL_ARRAY,MPI_INT,0,k,MPI_COMM_WORLD,&status1);
	
	

	//Computation starts using MULTI THREADING
	omp_set_num_threads(g);
#pragma omp parallel
	{
	int ID = omp_get_thread_num();j=0;
	for(j=0;j<g;j++)
	{
	int absolute_i= (((int)rank/g)*g)+ID;
	int absolute_j= ((rank%g)*g)+j;
	original_chunk[ID][j]=min(original_chunk[ID][j], add(i_rows_chunk[ID][k],k_row[absolute_j]) ) ;
	}
	}
	
	
	
	
	//Sending the results to master to update
    for(i=0;i<g;i++)
	{
	MPI_Send(&original_chunk[i][0],g,MPI_INT,0,i,MPI_COMM_WORLD);
    }//end of i loop
	
  
 } //MAIN k loop ends

 
 
 printf("\n After the computation, the values of grid in process %d are  - \n",rank);
 
    for(i=0;i<g;i++)
	{
	for(j=0;j<g;j++)
	{
	printf(" %d ",original_chunk[i][j]);
	}
	printf("\n");
	}
    printf("\n\n");

 }//End of Slave Process
 
 
 MPI_Finalize();
 return 0;
}

//This function returns the minimum of 2 numbers
int min(int a,int b)
{
if(a<b) return a;
else return b;
}

//This function returns the sum of 2 numbers. This also handles the cases where the input is infinity.
int add(int a,int b)
{
if(a==INT_MAX && b==INT_MAX) return INT_MAX;
if(a==INT_MAX || b==INT_MAX) return INT_MAX;

return a+b;
}