#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ------ 2-D HEAT CONDUCTION SOLVER -------*/
// explicit method
// Unsteady 2D conduction

//v4 includes a matrix sum for de div matrix that saves each point divergence
/*
  	if div=1    point converged
  	if div=0    no converggence

  	if sum div[n][m] == n*m  divergence = 0 that means convergence of all mesh elements
*/

void matprop(double **prope);
void material(double *x, double *y, double xa, double ya, double yb, int n, int m,  int **mat, double Lx, double Ly,int info1);
void linspace(double *xy, double L, int nodes, int info6);
void printfullmatrix(int n, int m, double **A, int info);
void printinteriormatrix(int n, int m, double **A, int info);
void temperaturemap(int n, int m, double **T, double Tini, double Tw1, double Tw2,double Tw3,double Tw4, int info2);
int position(int j, int i, int n, int m);
void stabcriterion(double ae, double an, double aw, double as, double ap0);
void heatconduction(double *x, double *y, double **T, double **prop, int *divergence, double theta, int n, int m,double dx, double dy,double dt,int **mat, int **div, double time,int explicit,int implicit);
void intmatrixsum(int n, int m, int **A,int *divergence);
int inposition(int j, int i, int n, int m);
void temperaturebound(int n, int m, double **T, double Tw1, double Tw2,double Tw3,double Tw4);


int main(int argc, char const *argv[])
{

	/* ---- PROBLEM STATEMENT ---------*/
	/* 	->  Centered nodes
		->  Extra dimension for bouundary nodes for better computation
		->  Computation only done for interior nodes
	*/

	/*

							Wall 1
			________________________________________
			|c	-----------Posi 1------------      c|	c=Corner									
			|-									   -|	
			|-									   -|
			|P									   P|
	Wall 2	|o									   o|  Wall 4    
			|s			Interior 				   s|
			|i			nodes    				   i|
			|2									   4|
			|-									   -|
			|-									   -|
			|c	--------------------------------   c|
			 _______________________________________
							Wall 3
			
			Interioor nodes  -> Posi 0   ------>Grid formulation [n][m]
			
			Posi 1 - Twall1
			Posi 2 - Twall2     ----->  They correspond to nodes placed at walls(boundaries)
			Posi 3 - Twall3
			Posi 4 - Twall4

			Interior nodes + Boundary nodes    ----> Grid formulation [n+2][m+2]
	*/

	//Phisical Dimensions 2-D plate
	/*----------------------------*/
	double Lx = 1.1;		
	double Ly = 0.8;
	
	//In percent
	double xa = 0.4555555*Lx;
	double ya = 0.125*Ly;
	double yb = 0.5*Ly;	
	

	//Discretization
	/*-----------------------------*/

	double n = 200;		// X Axis nodes [1,2,3,4,5,..., n]
	double m = 200;		// Y Axis nodes [1,2,3,4,5,..., m]

	double dx;
	double dy;

	double *x;			// X Axis geom nodes [0, 1,2,3,4,5,...,n, n+1]
	double *y;			// Y Axis geom nodes [0, 1,2,3,4,5,...,m, m+1]

	dx = Lx/n;
	dy = Ly/m;

	//Properties
	/*-----------------------------*/
	int nprop=4,nmat=4;  //nprop=number of properties  // nmat = number of materials
	double **prop;
	int 	**mat;
	double 	**T;

	//Boundary and interior temperature
	/*---------------------------------*/
	//Right:	Isoterm T=25ºC
	//Left: 	Uniform temperature T=10+0.005t ºC
	//Top: 		Contact with fluid Tg=35ºC  Heat transfer coeff = 9 W/mK
	//Bottom    Uniform Qflow=9 W/m

	//Boundaries and interior   // Temperature [K]
	double Tini	  = 8 + 273;	
	double Twall1 = Tini;			//Q flow over boundary	
	double Twall2 = 33+273;		
	double Twall3 = 23+273;
	double Twall4 = 8+273;
	







	//Solver
	/*-----------------------------*/
	double 	error;
	double 	theta=0.00001;
	int 	divergence=1;
	int 	iterations=0;

	double t;
	double time=5000;//Seconds
	double dt=1;	//Seconds

	int implicit = 1;		//Implicit = 1 = Gauss-Seidel     
							//Implicit = 2 = Line-By-Line
	int explicit =0; 

	if(implicit==explicit) exit(20000);

	//Information
	int info1=0;	//Matherial
	int info2=0;	//Temperaturemap
	int info3=0;	//Tp explicit 
	int info4=0;	//Tp implicit total matrix
	int info5=1;    //Tp implicit inner matrix
	int info6=0; 	//Linspace information
	int info7=1;    //Iterative parameters
	int info10=1;	//Print position matrix
	int info11=0;   //Print div matrix

	int **div;		//Auxiliar parameter

	//Dinamically Memory Allocation
	/*-----------------------------*/
	mat =	(int **)malloc(sizeof(int*)*(n+2));
	T 	=	(double **)malloc(sizeof(double*)*(n+2));
	x 	=	(double *)malloc(sizeof(double)*(n+2));
	y 	= 	(double *)malloc(sizeof(double)*(m+2));
	prop=	(double **)malloc(sizeof(double*)*nprop);
	div =	(int **)malloc(sizeof(int*)*(n)); //Only for interior nodes 

	if (prop == NULL){
		fprintf(stderr, "Out of memory");
		exit(1000);
	}	
	if (T == NULL){
		fprintf(stderr, "Out of memory");
		exit(100);
	}
	if (mat == NULL){
		fprintf(stderr, "Out of memory");
		exit(101);
	}
	if (x == NULL){
		fprintf(stderr, "Out of memory");
		exit(106);
	}
	if (y == NULL){
		fprintf(stderr, "Out of memory");
		exit(107);
	}
	if (div == NULL){
		fprintf(stderr, "Out of memory");
		exit(10700);
	}

	int i,j;
	for(i=0; i<n+2; i++){
		mat[i]	=	(int*)malloc(sizeof(int)*(m+2));
		T[i]	=	(double *)malloc(sizeof(double)*(m+2));
		div[i]	=	(int*)malloc(sizeof(int)*(m));
		if (mat[i] == NULL)
		{
			fprintf(stderr, "Out of memory");
			exit(102);
		}
		if (T[i] == NULL)
		{
			fprintf(stderr, "Out of memory");
			exit(103);
		}
		if (div[i] == NULL)
		{
			fprintf(stderr, "Out of memory");
			exit(100003);
		}
	}

	for(i=0; i<nprop; i++){
		prop[i]=(double *)malloc(sizeof(double)*nmat);
		if (prop[i] == NULL)
		{
			fprintf(stderr, "Out of memory");
			exit(1002);
		}
	}

	

/*----------------------------------------------CODE CORE-------------------------------------------------------*/

	matprop(prop);
	linspace(x, Lx, n,info6);
	linspace(y, Ly, m,info6);
	material(x, y,  xa,  ya,  yb,  n,  m,  mat, Lx, Ly,info1 );
	temperaturemap(n,m,T,Tini,Twall1,Twall2,Twall3,Twall4,info2);
	

	if(explicit){
		for(t=0.000001; t<=time+0.000001; t=t+dt){

			iterations=iterations+1;
			if(info3)printf("Time: %f  \n",t );
			heatconduction(x,y,T,prop,&divergence,theta,n,m,dx,dy,dt,mat,div,time,explicit,implicit);
			printfullmatrix( n,  m,  T,  info3);

		}
	}
	
	if(implicit==1){

		for(t=dt; t<=time; t=t+dt){
		
					
			while(divergence==1){

				iterations=iterations+1;
				
				heatconduction(x,y,T,prop,&divergence,theta,n,m,dx,dy,dt,mat,div,time,explicit,implicit);		

			}

			divergence=1;	
			iterations=0;
			if(info4)printf("Time: %f Iterations = %d\n",t,iterations );

			Twall4 = Twall4 + 0.005*t;
			temperaturebound( n,  m,  T,  Twall1,  Twall2, Twall3, Twall4);

		}

		printinteriormatrix( n,  m,  T,  info5);
		printfullmatrix( n,  m,  T,  info4);
			
			
	}	

	
	/*----------------AUXILIAR INFORMATION ------------------------------*/
	int posi;
	if(info10){
		for (j=0; j<m+2; j++ ){
			for(i=0; i<n+2; i++){

				posi=inposition( j,  i,  n, m);
				
				if(i==0) printf("row %d \t:  ",j);
				printf("%d ", posi);
					
				
			} 
			printf("\n");
		}
	}


	/*-----------------FREE MEMORY SPACE---------------------------------*/

	free(mat);
	free(T);
	free(x);
	free(y);
	free(prop);
	free(div);

	return 0;
}

void matprop(double **prope){

	/*-----   For saving the needed  -------*/
	/*-----   material properties    -------*/


	/*
	int prop[5][4] = {
						{1500,1700,1900,2600},	//Prop 0 - Density (kg/m^3)
						{750,770,810,930},		//Prop 1 - Cp (J/Kg*K), the specific heat measured at constant pressure
						{170,140,200,140},		//Prop 2 - Thermal Conductivity k (W/m*K)
						{1,2,3,4}				//Prop 3
					};*/

	//Prop 1 - Density (kg/m^3)
	prope[0][0] = 1500;
	prope[0][1] = 1700;
	prope[0][2] = 1900;
	prope[0][3] = 2600;

	//Prop 2 - Cp (J/Kg*K), the specific heat measured at constant pressure
	prope[1][0] = 750;
	prope[1][1] = 770;
	prope[1][2] = 810;
	prope[1][3] = 930;

	//Prop 3 - Thermal Conductivity k (W/m*K)
	prope[2][0] = 170;
	prope[2][1] = 140;
	prope[2][2] = 200;
	prope[2][3] = 140;

	//Prop 4
	prope[3][0] = 1111;
	prope[3][1] = 2222;
	prope[3][2] = 3333;
	prope[3][3] = 4444;
}

void material(double *x, double *y, double xa, double ya, double yb, int n, int m,  int **mat, double Lx, double Ly,int info1){

	/* ---    Finds Matherial over the 2D domain  ---- */

	//For saving the data into a .txt
	FILE * fPointer;
	fPointer = fopen("matherial.txt", "w" );

	//For finding the material 
	int i,j;
	for(j=0; j<m+2; j++){
		for(i=0; i<n+2; i++){


			//Conditions for 4 material selection
			if(x[i] <= xa && y[j] >= yb) mat[i][j] = 1;
			if(x[i] >= xa && y[j] <= ya) mat[i][j] = 3; 
			if(x[i] < xa && y[j] < yb) mat[i][j] = 0;
			if(x[i] > xa && y[j] > ya) mat[i][j] = 2; 
			
			if(info1){
				if(i==0) printf("row %d \t:  ",j);
				printf("%d ", mat[i][j]);				
			}
		
			fprintf(fPointer, "%d", mat[i][j] );
			fprintf(fPointer, " ");

			//printf("Material[%d][%d] \t %p \t %d \n", i,j , &mat[i][j], mat[i][j]);
		}
		if(info1)printf("\n");
		fprintf(fPointer, "\n" );		
	}
	fclose(fPointer);
}

void linspace(double *xy, double L, int nodes, int info6){
	/*--- Generates Uniform grid spacing	-------*/
	// \Delta x = Lx/n
	// \Delta y = Ly/m
	int i;
	for(i=0; i<nodes+2; i++){
		
		if(i==0){
			xy[i]=0;
		}
		if(i==1){
			xy[i]=L/(nodes*2);
		}
		if(i>1 && i<nodes+1){
			xy[i]=xy[i-1]+L/nodes;
		}
		if(i==nodes+1){
			xy[i]=xy[i-1]+L/(nodes*2);
		}
		if(info6){
			printf("%f,  ", xy[i] );
		}
	}
	if(info6)printf("\n");

}

void printfullmatrix(int n, int m, double **A, int info){
	
	/*---- Prints matrix of dimensions [n+2][m+2]	-------*/

	//Works for double values
	
	FILE * fPointer;
	fPointer = fopen("fullMATRIX.txt", "w" ); //use a for no overwriting files
	
	int i,j;
	for(j=0; j<m+2; j++){
		for(i=0; i<n+2; i++){


			if(info){
				if(i==0) printf("row %d \t:  ",j);
				printf("%f ", A[i][j]);				
			}
		
			fprintf(fPointer, "%f", A[i][j] );
			fprintf(fPointer, " ");

			//printf("Material[%d][%d] \t %p \t %d \n", i,j , &mat[i][j], mat[i][j]);
		}
		if(info)printf("\n");
		fprintf(fPointer, "\n" );		
	}	
	fclose(fPointer);	
}

void printinteriormatrix(int n, int m, double **A, int info){

	/*---- Prints matrix of dimensions [n][m]	-------*/
	//Works for int values
	
	FILE * fPointer;
	fPointer = fopen("nodeMATRIX.txt", "w" );
	
	int i,j;
	for(j=1; j<m+1; j++){
		for(i=1; i<n+1; i++){


			if(info){
				if(i==1) printf("row %d \t:  ",j);
				printf("%f ", A[i][j]);				
			}
		
			fprintf(fPointer, "%f", A[i][j] );
			fprintf(fPointer, " ");

			//printf("Material[%d][%d] \t %p \t %d \n", i,j , &mat[i][j], mat[i][j]);
		}
		if(info)printf("\n");
		fprintf(fPointer, "\n" );		
	}
	fclose(fPointer);
	
}

void intmatrixsum(int n, int m, int **A,int *divergence){

	/*-----  Sums int type matix elements   ----*/

	int i,j;
	int sum=0;
	for(j=1; j<m+1; j++){
		for(i=1; i<n+1; i++){			
		
			sum = sum + A[i][j];
		}	
	}

	if(sum == n*m) *divergence=0;
}

void temperaturemap(int n, int m, double **T, double Tini, double Tw1, double Tw2,double Tw3,double Tw4, int info2){

	/*-------- Fills initial temperature matrix  -----------*/
	/*------ with corresponding boundary conditions  -------*/

	int i,j;
	int posi;
	for (j=0; j<m+2; j++ ){
		for(i=0; i<n+2; i++){


			posi = position( j,  i,  n, m);
			if(posi ==1) T[i][j] = Tw1;			
			if(posi ==3) T[i][j] = Tw3; 
			if(posi ==0) T[i][j] = Tini;

			//Needed but not used
			if(posi ==4) T[i][j] = Tw4;
			if(posi ==2) T[i][j] = Tw2;

			/* Corner temperature
			if(posi ==11) T[i][j] = (Tw1+Tw2)/2;
			if(posi ==22) T[i][j] = (Tw2+Tw3)/2;
			if(posi ==33) T[i][j] = (Tw3+Tw4)/2; 
			if(posi ==44) T[i][j] = (Tw4+Tw1)/2;
			*/

			if (info2){
				if(i==0) printf("row %d \t:  ",j);
				printf("%f ",T[i][j]);
			}				
		} 
		if(info2) printf("\n");
	}

}

void temperaturebound(int n, int m, double **T, double Tw1, double Tw2,double Tw3,double Tw4){

	/*-------- Fills initial temperature matrix  -----------*/
	/*------ with corresponding boundary conditions  -------*/

	//For postprocessing could be a good tool

	int i,j;
	int posi;
	for (j=0; j<m+2; j++ ){
		for(i=0; i<n+2; i++){


			posi = position( j,  i,  n, m);
			if(posi ==1) T[i][j] = T[i][j+1];			
			if(posi ==3) T[i][j] = Tw3; 
			if(posi ==4) T[i][j] = Tw4;
			if(posi ==2) T[i][j] = Tw2;

		} 
	}

}


int position(int j, int i, int n, int m){

	/*--------- Finds walls position and corners -----------*/
	/*------ for a matrix of dimensions [n+2][m+2] ---------*/

	int position=0;			//In node

	if(i>0 && i<n+1 && j==0 ) position = 1;		//Wall 1 
	if(i==0 && j==0) position=11; 				//Corner 1

	if(i==0 && j>0 && j<m+1 ) position = 2;		//Wall 2
	if(i==0 && j==m+1) position=22; 			//Corner 2

	if(i>0 && i<n+1 && j==m+1 ) position = 3;	//Wall 3
	if(i==n+1 && j==m+1) position=33; 			//Corner 1

	if(j<m+1 && i==n+1 && j>0 ) position = 4;	//Wall 4
	if(i==n+1 && j==0) position=44; 			//Corner 1	

	return(position);	
}

int inposition(int j, int i, int n, int m){

	/*--------- Finds walls position and corners -----------*/
	/*------ for a matrix of dimensions [n][m] ---------*/

	int position=0;			//In node

	if(i>0 && i<n+1 && j==1 ) position = 1;		//Wall 1 
	if(i==1 && j==1) position=11; 				//Corner 1

	if(i==1 && j>1 && j<m ) position = 2;		//Wall 2
	if(i==1 && j==m) position=22; 			//Corner 2

	if(i>1 && i<n && j==m ) position = 3;	//Wall 3
	if(i==n && j==m) position=33; 			//Corner 1

	if(j<m && i==n && j>1 ) position = 4;	//Wall 4
	if(i==n && j==1) position=44; 			//Corner 1	

	return(position);	
}

void stabcriterion(double ae, double an, double aw, double as, double ap0){

	/*---- Checks if the stability criterion for ------*/
	/*----     the explicit method is followed   ------*/

	if(ap0 < ae+an+aw+as ) {
		printf("Stability criterion violated, reduce time step " );
		exit(5);
	}
}


void heatconduction(double *x, double *y, double **T, double **prop, int *divergence, double theta, int n, int m,double dx, double dy,double dt,int **mat,int **div, double time, int explicit, int implicit){

	/*-------  Solves 2D heat conduction equatio for conduction boundaries -----------*/

	/*--- Parameters needed ------*/

	//Node temperatures
	double Tp,Tp0,Tn,Ts,Te,Tw;

	//Equation coefficients
	double ap,ap0,an,as,ae,aw, b; 

	//Node geometry
	double dxe,dxw,dyn,dys;

	//Thermal conductivity
	double Kn,Ks,Ke,Kw,Kp;
	double kn,ks,ke,kw;

	//Material properties
	double density, cp;

	//Intern heat sources
	double Gc=0, Gp=0;

	//For implicit method convergence
	int totaldivergence=0;

	//Boundary nodes
	int posi;		//Auxiliar for position
	double Q1 = 60;	//Top boundary q flow
	double h2 = 9;



	

	int i,j;
	for (j=1; j<m+1; j++ ){
		for(i=1; i<n+1; i++){


			if(implicit==1){
				Tp0 =	T[i][j];
				Te 	=	T[i+1][j];
				Tw 	=	T[i-1][j];
				Tn 	=  	T[i][j-1];
				Ts 	=	T[i][j+1];				
			}
			
			dxe = x[i+1]-x[i];
			dxw = x[i]-x[i-1];
			dyn = y[j]-y[j-1];
			dys = y[j+1]-y[j];

			Kn 	= 	prop[2][mat[i][j+1]];
			Ks  =	prop[2][mat[i][j-1]];
			Ke  =   prop[2][mat[i+1][j]];
			Kw  =   prop[2][mat[i-1][j]];
			Kp  =   prop[2][mat[i][j]];

			kn 	=	(2*Kp*Kn)/(Kp+Kn);
			ks 	=	(2*Kp*Ks)/(Kp+Ks);
			kw 	=	(2*Kp*Kw)/(Kp+Kw);
			ke 	=	(2*Kp*Ke)/(Kp+Ke);

			ae 	=	(ke*dy)/dxe;
			an	=	(kn*dx)/dyn;
			aw 	=	(kw*dy)/dxw;
			as 	= 	(ks*dx)/dys;

			density = prop[0][mat[i][j]];
			cp 		= prop[1][mat[i][j]];

			//BOUNDARY CONDITIONS
			posi = inposition( j,  i,  n, m);
			if(posi == 1){
				an = 0;
				ap0 = (density*cp*dx*dy)/dt;
				b 	= Gc*dx*dy + ap0*Tp0 + Q1;	
			}
			if(posi == 2){
				aw = h2*dy;
			}	
			if(posi == 3){
				//no need
			}	
			if(posi == 4){
				//no need
			}
			if(posi == 11){
				an = 0;
				ap0 = (density*cp*dx*dy)/dt;
				b 	= Gc*dx*dy + ap0*Tp0 + Q1;

				aw = h2*dy;
			}
			if(posi == 22){
				aw = h2*dy;
				
			}	
			if(posi == 33){
				//no need
			}	
			if(posi == 44){

				an = 0;
				ap0 = (density*cp*dx*dy)/dt;
				b 	= Gc*dx*dy + ap0*Tp0 + Q1;
			}		


			if(implicit==1){

				if(posi != 1 || posi != 11 || posi != 44){
					ap0 = (density*cp*dx*dy)/dt;
					b 	= Gc*dx*dy + ap0*Tp0;
				}
							
				ap 	= ae+aw+an+as+ap0-Gp*dx*dy;

				T[i][j]	=	(ae*Te + aw*Tw + as*Ts + an*Tn  +b)/ap;


				if( abs(T[i][j]-Tp0) <= theta ){
					totaldivergence = totaldivergence +1;					
				}

				if(totaldivergence==(n*m)) *divergence=0; 

			}				
		}
	}
}


