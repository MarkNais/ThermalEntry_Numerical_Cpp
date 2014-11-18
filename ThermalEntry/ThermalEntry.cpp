/***********************************************************************************
* Mech 7171 Engineering Programing
* LAB #4
*
* Purpose: The of this lab is to solve the heat conduction equation in 2D
*          using finite difference (F-D) methods.  The solution will be compared
*          to an analytic solution and solution convergence will be tracked.
*          This lab will exercise a facets of C taught in Mech 7171.
*
* Author: Corbin Turner, Cyrus Ang, Mark Naismith
* Student ID: A00780890, A00781218, A00819714
*
* Declaration:
* I, Corbin Turner, Cyrus Ang, Mark Naismith, declare that the following 
* program was written by us.
*
* Date Created: Dec.6/2013
* Revised: Mark Naismith, 2/25/2014 - Edited for Fluid Flow Project
* Revised: Mark Naismith, 9/20/2014 - Edited for Koorosh's Fin Lab 01
***********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <new.h>

#pragma warning(disable:4996)

#define PI                       (4.0*atan(1.0))
#define MAX_ITER                 100000            // maximum iterations for F-D
#define DATA_FILENAME			 "data.txt"
#define MAX_LINE_LENGTH			 1024

//  structure for finite-difference solution
typedef struct
{
	double r;   // grid node unit-less position along the length of the plate
	double Temp;   // grid node temperature (finite-difference solution)
	double u;        //Fully developed flow velocity
}
PLATEPOINT;


typedef struct
{
	int scase;		//Case counter (which simulation is it?)
	int Ny;			//Node count in y (int?)
	double Pe;		//Biot number
	double Tfinal;	//When to stop the simulation - The Theta that determines when the simulation ends.
	double dx;		//the unitless node spacing along pipe
	double dy;		//the unitless node spacing through pipe
	int iter;		//Tracking the iteration count

	//Interior limits of the simulation, ex: of a 9 by 5 grid, only the inside 7 by 3 grid is calculated.
	int NyInter;
	int run;		//Run the simulatuon? 1=run 0=dont run

	PLATEPOINT *pp;		//The pointer to the first  array of plate points. (t) "current"
	PLATEPOINT *pp2;	//The pointer to the second array of plate points. (t+1) "future"

	double *tri_hold[4];//Components of Koorosh's Tri-Diagonal solver. Should not be passed to the function

	double *tri_a;		//a-Component of Koorosh's Tri-Diagonal solver
	double *tri_b;		//b-Component of Koorosh's Tri-Diagonal solver
	double *tri_c;		//c-Component of Koorosh's Tri-Diagonal solver
	double *tri_y;		//c-Component of Koorosh's Tri-Diagonal solver

}
PROGRAMDATA;

//typedef struct
//{
//	int iter;
//	double MaxRes;
//	//double sumResSq;
//	//double RMS;
//
//	//Interior limits of the simulation, ex: of a 9 by 5 grid, only the inside 7 by 3 grid is calculated.
//	int NxInter, NyInter;
//	double dx, dy;
//}
//SIMULATIONDATA;

//------------------------- FUNCTION PROTOTYPES -----------------------------------------

// Get data from text file
PROGRAMDATA GetProgramData(FILE*);

//Initialize some of the variables inside of the pd struct.
PROGRAMDATA pdInit(PROGRAMDATA);

// set boundary functions for each case
void boundary_set(PROGRAMDATA);  //constant temp.

// Function to contain the specific order of regions to be simulated
void num_simulation(PROGRAMDATA);
// mathematic function for simulating the interior of the body (and the boundaries)
PROGRAMDATA num_sim_body(PROGRAMDATA);

// functions for each simulation
void simulate(PROGRAMDATA);

//Copy the contents of tri_hold
PROGRAMDATA triCopy(PROGRAMDATA);

// allocating and free allocation fucntions
PROGRAMDATA allocate(PROGRAMDATA);
void freepp(PROGRAMDATA);

// reading and writing functions
FILE *filewrite(PROGRAMDATA);
//Dear Dave. No, I shan't call it fileopen, as we're prepping to write to a file 
//with this function, not just open one, which might imply only reading.
//-Mark N.
void printLabPlates(PROGRAMDATA);

//rounding function
int nint(double);

//Koorosh's TriDiagonal solver
void Tridiagonal( int, double*, double*, double*, double*);
//------------------------- END OF FUNCTION PROTOTYPES ----------------------------------


int main()

{  
	int i;  //loop counter


	//header
	printf("MECH 7171 Engineering programming                %c\n",179); 
	printf("Lab #4                                           %c\n",179); 
	printf("Mark Naismith, A00819714, Set 5B                 %c\n",179);
	printf("Cyrus Ang,     A00781218, Set 5B                 %c\n",179);
	printf("Corbin Turner, A00780890, Set 5B                 %c\n",179);
	printf("December 6, 2013                                 %c\n",179);
	printf("This program reads data from a data.txt file and %c\n",179);
	printf("stores the data in a structure within the program%c\n",179);
	printf("that is transfered between modular functions.    %c\n",179);
	printf("The tasks that this program completes are the    %c\n",179); 
	printf("analytical solution to 2D conductive             %c\n",179); 
	printf("heat transfer probrlems as well as               %c\n",179); 
	printf("the numerical solutions. The end result is a     %c\n",179); 
	printf("given heat transfer problem                      %c\n",179); 
	printf("that is solved numerically.                      %c\n",179);

	//this loop "closes the box" of the header. Run program for details.
	for(i=0;i<49;i++)printf("%c",196); 
	printf("%c\n",217); 

	//creation of the structure for simulations. Holds all relevant data
	PROGRAMDATA pd;

	//File stream pointer for the input file. This will be passed to the GetProgramData
	//function for every simulation
	FILE *f;

	//Opens stream to the input file. Read only
	f=fopen(DATA_FILENAME, "r");

	//Check to see if the file exists, if not, exits function, main() complains
	if(f==NULL)
	{
		printf("Input file \"%s\" doesn't exist! Exiting...",DATA_FILENAME);
		getchar();
		return 5;   
	}

	i = 0;
	//The proceding blocks are all very similar. Only the first will be explained.
	do{
		i++;
		pd = GetProgramData(f);
		if (pd.run == 1)//Checks to see if user wants to run 
		{
			printf("\nRunning Simulation %d :\n", i);
			simulate(pd);
		}
	} while (pd.run != -1);

	// closes data file once all simulations are complete 
	//(could be closed before Given problem)
	fclose(f);

	//End Program
	printf("Press Enter to end the program....");
	fflush(stdin);//needed if the last file needed to be overwriten.
	getchar();
	return 0;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                          FUNCTIONS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/********************************************************
* GetProgramData *
* *
* Purpose: Retrieves only the appropriate amount of data from the input file.
* Note: it does not error check for bad lines or an incomplete file. 
* *
* Parameters *
* i - loop counter for retrieving data , 
* j - static int to determine case for filewrite function
* *pdp[] - array of pointers to addresses of varibales in PROGRAMDATA structure
* *    ^good idea Dave^
* *
* Returns *
* pd to PROGRAMDATA function  *
********************************************************/
PROGRAMDATA GetProgramData(FILE *f)
{
	//create a new program data structure to overwrite the old one.
	PROGRAMDATA pd;

	char buff[MAX_LINE_LENGTH];   //Receiver of lines from the input file
	int i;                    //loop counter  
	static int j=1;           //case counter (how many times 
	//has this function been called?)
	//An array of pointers where the read data is to be stored.
	double *pdp[] = { &pd.Pe, &pd.Tfinal, &pd.dx};
	int *pdpi[]={&pd.run,&pd.Ny};

   //assume the simulations should not run, Minor protection from read error.
   pd.run=0;

	// retrieves 3 more lines of 'int' data, and points them to their destination
	for(i=0;i<3;i++)
	{
      fgets(buff,MAX_LINE_LENGTH,f);
		if(feof(f)!=0){
         pd.run=-1;
         return pd;
      }
		*pdpi[i]=atoi(buff);                         
	}//End of for

	// retrieves 2 more lines of 'double' data, and points them to their destination
	for(i=0;i<2;i++)
	{
		fgets(buff,MAX_LINE_LENGTH,f);
		if(feof(f)!=0){
         pd.run=-1;
         return pd;
      }
		*pdp[i]=atof(buff);                         
	}//End of for

   // if the end of the file is reached, stop trying to run simulations
   

	pd.scase=j; //sets the case for filewrite later
	j++;        //increments j
	fgets(buff,MAX_LINE_LENGTH,f); // skips a line down in data file

	// since we have not used fclose, the file is still open and the line that this
	// function is looking at remains the same
	return pd;     //return pd to main()   
}//End of GetProgramData

/********************************************************
* allocate *
* *
* Purpose: allocate an array to be used for simulations
* *
* Parameters *
* w,h - loop counters for nodes horizontal, vertical
* W,H - horizontal and vertical nodes calculated 
*       from input data *
* *
* Returns *
* Initialized allocated array into PLATEPOINT structure nested in PROGRAMDATA *
********************************************************/
PROGRAMDATA allocate(PROGRAMDATA pd)
{
	int i;
	double *pda[]={pd.tri_hold[0], pd.tri_hold[1], pd.tri_hold[2], pd.tri_hold[3], pd.tri_a, pd.tri_b, pd.tri_c, pd.tri_y};

	if(pd.Ny<3)
	{
		printf("Sorry, you must use a positive integer of 3 or larger");
		printf("for Node Count\n");
		printf("\nPress Enter to end this program...");
		getchar();
		exit(0);
	}

	// allocate memory for a horizontal array of pointers
	pd.pp = (PLATEPOINT*)malloc(pd.Ny*sizeof(PLATEPOINT));
	pd.pp2 = (PLATEPOINT*)malloc(pd.Ny*sizeof(PLATEPOINT));

	// check allocation for array of pointers
	if (pd.pp == NULL || pd.pp2 == NULL)
	{
		printf("Cannot allocate pd.pp or pd.pp2, exiting program...\n");
		getchar();
		exit(0);
	}

	//introduce allocation for the new freepp additions.
	for(i=0;i<8;i++){
		pda[i] = (double*)malloc(pd.Ny*sizeof(double));
		if (pda[i] == NULL)
		{
			printf("Cannot allocate matrix-array number: %d, exiting program...\n",i);
			getchar();
			exit(0);
		}
	}
	return pd;
}

/********************************************************
* simulate1 *
* *
* Purpose: Run const. temp coarse simulation
* *
********************************************************/
void simulate(PROGRAMDATA pd)
{
	/*
	Function calls are in the order of:
	allocate array
	initialize the other variables inside of the structure
	initialize borders
	numerical simulation
	print to files 
	free allocation 
	*/
	pd=allocate(pd);
	pd=pdInit(pd);
	boundary_set(pd);
	num_simulation(pd);
	printLabPlates(pd);
	freepp(pd);
}

/********************************************************
* pdInit *
* *
* Purpose: Set some of the pd's variables that do not come from the data.txt file
* *
* Parameters *
* int i - for loop index
* *
* Returns *
* None (void) *
********************************************************/
PROGRAMDATA pdInit(PROGRAMDATA pd)
{
	pd.iter=0;
   pd.NyInter = pd.Ny - 1;
	pd.dy = 1 / pd.NyInter;
	return pd;
}

/********************************************************
* boundary_set *
* *
* Purpose: Set boundaries for analytical or numerical for constant temperature
* *
* Parameters *
* w - loop counter for horizontal nodes , 
* pd.Nx,pd.Ny - horizontal and vertical nodes calculated 
*       from input data *
* *
* Returns *
* None (void) *
********************************************************/
void boundary_set(PROGRAMDATA pd)
{
   int y;
   double a = 4.0/(pd.Pe*pow(pd.dy,2)), bc = -2.0 / (pd.Pe * pow(pd.dy,2)) ;

	// initialize PLATEPOINT variables in allocated array 
	for(y=0;y<pd.Ny;y++)
	{
		pd.pp[y].Temp = 0.0;
		pd.pp2[y].Temp = pd.pp[y].Temp;

		pd.pp[y].r = (double)y*pd.dy;
		pd.pp2[y].r = pd.pp[y].r;

		pd.pp[y].u = 2*(1-pow(pd.pp[y].r,2));
		pd.pp2[y].u = pd.pp[y].u;
      
      pd.tri_hold[3][y] = pd.pp[y].u / pd.dx;
      pd.tri_hold[0][y] = pd.tri_hold[3][y] + a;
      pd.tri_hold[1][y] = bc * (1 + pd.dy/(2*pd.pp[y].r));
      pd.tri_hold[2][y] = bc * (1 - pd.dy/(2*pd.pp[y].r));

	}
   pd.pp[pd.NyInter].Temp = 1.0;
   pd.pp2[pd.NyInter].Temp = 1.0;
   
   pd.tri_hold[0][0]= 3.0/2.0; //R1a
   pd.tri_hold[1][0]= -2.0; //R1b
   pd.tri_hold[2][0]= 1.0/2.0; //R1c
   pd.tri_hold[3][0]= 0.0; //R1y
   pd.tri_hold[0][pd.NyInter] = 1.0; //a
   pd.tri_hold[1][pd.NyInter] = 0.0; //b
   pd.tri_hold[2][pd.NyInter] = 0.0; //c
   //tri_hold[3] at NyInter is the same from the loop above.
}

/********************************************************
* num_simulation *
* *
* Purpose: Simulate numerical solution for Thermal Entry
* *
* Parameters *
* FILE *f - file pointer
* i,j,iter - loop counters for nodes horizontal, vertical, 
*          and iteration number *
* W,H - horizontal and vertical nodes calculated 
*       from input data *
* *
* Returns *
* None (void) *
********************************************************/
void num_simulation(PROGRAMDATA pd)
{
	FILE *f;
	PLATEPOINT *Hold;

	f=filewrite(pd);

	// error checking
	if(f==NULL)
	{
		printf("File can be read but not writen to.\n");
		getchar();
		exit (1);
	}

	fprintf(f,"iter,log(MaxRes)\n");

	//This loop will look through every core temperature node
	/*The loop will continue as long as the Maximum residual of the
	current iteration is larger than the defined maximum residual*/
	do
	{
		pd.iter++;
		//Need to reset or else MaxRes will always > abs(Res) after the first iteration
		pd.MaxRes=0.0;

		//Simulates the innerbody, then the boundaries. Writes to pd.pp2
		pd=num_sim_body(pd);

		//swap the addresses
		//The addresses will be swapped so that the latest timestep (pp2) will take the place of the old data (pp),
		//the old data will be overwriten with even newer data.
		//pp is only ever written to durring initialization. During simulation, it is only read.
		Hold = pd.pp;
		pd.pp = pd.pp2;
		pd.pp2 = Hold;

		//Print off the maximum residual and root mean squared for each iteration
		fprintf(f,"%d,%le\n",pd.iter,log10(pd.MaxRes));
	}
	while(pd.pp[pd.NxInter].Temp>pd.Tfinal && pd.iter<=MAX_ITER ); //|| iter<=MAX_ITER

	fclose(f);
}

/********************************************************
* triCopy *
* *
* Purpose: Copy the contents from pd.tri_hold to the component vectors
* *
* Parameters *
* i - for loop index
* *
* Returns *
* pd *
********************************************************/
PROGRAMDATA triCopy(PROGRAMDATA pd){
   int i;
   for(i=0; i<pd.Ny; i++){
      pd.tri_a[i]=pd.tri_hold[0][i];
      pd.tri_b[i]=pd.tri_hold[1][i];
      pd.tri_c[i]=pd.tri_hold[2][i];
      pd.tri_y[i]=pd.tri_hold[3][i] * pd.pp[i].u / pd.dx;
   }
   return pd;
}

/********************************************************
* num_sim_body *
* *
* Purpose: Simulate numerical solution for the interior region
* *
* Parameters *
* i, and j. Index trackers
* *
* Returns *
* None (void) *
********************************************************/
PROGRAMDATA num_sim_body(PROGRAMDATA pd)
{


	return pd;
}

/********************************************************
* filewrite *
* *
* Purpose: Creates file names, checks for file existance, and streams files for writing.
* *
* Parameters *
* f - File stream, also the return
* resp - yes or no response
* filen - The filename to be composed procedurally
* j - static int acting as a switch. See below.
* ocase - copy of pd.scase, and comparator to the last time filewrite was called.
* *
* Returns *
* Initialized allocated array into PLATEPOINT structure nested in PROGRAMDATA *
********************************************************/
FILE *filewrite(PROGRAMDATA pd)
{

	FILE *f;    //Filestream pointer - the return
	char resp, filen[MAX_LINE_LENGTH]; //response char and filename composer

	//This static int will be used to first differentiate
	//between the RMS file and the plate simulation results file
	static int j = 0;

	//This static int is used to check to see if a new 
	//simulation has started since the last time the function was called.
	static int ocase = pd.scase;

	//If the simulation has changed, reset j to 0
	if (ocase != pd.scase) j = 0;
	ocase = pd.scase;//set ocase to scase, whether or not scase is differemt or not.

	//Write to first part of the file name.
	sprintf(filen, "Time-Flow Simulation Nx-");

	//write the second part of the file name

	switch (ocase)
	{
	case 1:
		strcat(filen, "11");
		break;
	case 2:
		strcat(filen, "21");
		break;
	case 3:
		strcat(filen, "41");
		break;
	case 4:
		strcat(filen, "81");
		break;
	case 5:
		strcat(filen, "161");
		break;

	default:
		break;
	}

	//If it's the first time it's been called in a simulation, it's an RMS file,
	//Next it's the plate file. If the function is called again,
	//it'll store the data in an extra txt file, warns the user.
	if (j == 0) strcat(filen, "RRRMS.csv");
	else if (j == 1) strcat(filen, " contour.dat");
	else
	{
		strcat(filen, "EXTRA.txt");
		printf("Function \"filewrite\" called too many times!\n");
	}


	//attempts to open the file in read mode for checking
	f = fopen(filen, "r");


	if (f != NULL) //If the file exists...
	{

		printf("File \"%s\" exists. Ok to overwrite? (y/n): ", filen);

		fflush(stdin);
		scanf("%c", &resp);        // get response from user

		// if response is "no", tell user the file will not be overwritten
		if (resp == 'n' || resp == 'N')
		{
			printf("The existing file will not be overwritten.\n");
			fclose(f);     //close the file stream
			j++;
			return NULL;   //exit the function with no stream.
		}
	}//end overwrite check

	printf("Opening \"%s\" for writing\n\n", filen);

	//opens the file for writing
	f = fopen(filen, "w");

	j++;//increment j

	// if f returns a NULL send the user a message and return a 
	//NULL to the appropriate function that called filewrite 
	if (f == NULL)
	{
		printf("File can be read but not writen to!\n");
		getchar();
		return NULL;
	}

	return f;
}

/********************************************************
* printLabPlates *
* *
* Purpose: prints the simulation results to the appropriate file
* *
* Parameters *
* f - File stream, from filewrite
* W,H - horizontal and vertical nodes calculated 
*       from input data *
* i,j - horizontal and vertical nodes calculated 
*       from input data *
* *
********************************************************/
void printLabPlates(PROGRAMDATA pd)
{
	FILE *f;
	int i;

	f = filewrite(pd);//get file stream
	if (f == NULL) return;//if something went wrong, exit and continue

	//header
	fprintf(f, "VARIABLES = ""X"", ""T""\n");
	fprintf(f, "ZONE I=%hd, F=POINT\n", pd.Nx);

	//loop heights, starting from the bottom
	
		for (i = 0; i<pd.Nx; i++)//loop width, starting from the left
		{
			//All if the values, in order, in scientific notation. See header for order
			fprintf(f, "%+12.7lg %+12.7lg\n"
				, pd.pp[i].x, 1-pd.pp[i].Temp);
		}
	fclose(f);//close the file stream
	return;
}

/********************************************************
* nint *
* *
* Purpose: floors resulting double to nearest int
* *
* Returns *
* int for W and H calculations*
********************************************************/
int nint(double d)
{
	//Rounds to the nearest int
	return (int)floor(d+0.5);
}

/********************************************************
* freepp *
* *
* Purpose: free allocated array
********************************************************/
void freepp(PROGRAMDATA pd)
{
	free(pd.pp);
	free(pd.pp2);
	free(pd.tri_hold);
	free(pd.tri_a);
	free(pd.tri_b);
	free(pd.tri_c);
}


//********************************************************************
//
// Tridiagonal
//
// m is the size of the matrix =(m*m)
// Care to the values vs. pointers
// The array destroys the variables, thus give copies
// 
// Function to solve a tridiagonal system of equations.
// The non zero terms are:
//    a32 b-2 c1/2
//    c   a   b
//       c   a   b
//          c   a   b
//          b0  c0  a1
void Tridiagonal( int m, double *c, double *a, double *b, double *y )
{
   int i;

   // Process first row.
   b[0] /= a[0];
   c[0] /= a[0];
   // Adjust b[1]
   b[1] -= c[1] * c[0];
    
   // Process rows 2 to n-1.
   for ( i = 1; i < (m-1); i++ )
      {
      a[i] -= c[i] * b[i-1];
      b[i] /= a[i];
      }

   // Last row.
   c[m-1] -= b[m-1] * b[m-3];
   a[m-1] -= c[m-1] * b[m-2];


   //
   // Forward solution step.
   //
   y[0] /= a[0];
   for ( i = 1; i < (m-1); i++ )
      y[i] = (y[i] - (c[i] * y[i-1])) / a[i];
   y[m-1] = (y[m-1] - (c[m-1] * y[m-2]) - (b[m-1] * y[m-3])) / a[m-1];

   //
   // Backward solution step.
   //
   for ( i = m-2; i > 0; i-- )
      y[i] -= b[i] * y[i+1];
   y[0] -= (b[0] * y[1]) + (c[0] * y[2]);

   return;
}