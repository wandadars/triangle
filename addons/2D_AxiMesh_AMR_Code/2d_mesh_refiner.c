/**********************************************************************/
/*	Purpose: Code to read an unstructured 2D triangle mesh and    */
/*	         and create a .area file detailing specific areas for */
/*	         each triangle element for mesh refinement.           */
/*	         	                        		      */
/*	         	                                              */
/*	 Notes: Code reads 2 files: .node and .ele to produce .area   */
/*	 	file. This code is specifically designed to work only */
/*	 	with the 2D meshes that the program, Triangle, is     */
/*	 	used to create.  	                              */
/*	                                                              */
/*	Author 	: Christopher Neal				      */
/*	Date 	: August 19th, 2014   				      */
/*	Edited	: April 14th, 2015	        		      */
/**********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct //Structure to hold vertices of mesh
{
	double  	x, y; // Coordinates of vertex
} vertex ;


typedef struct //Structure to hold all mesh information
{
	int     nNodes, nElems, nNodeMarkers;
	vertex 	*nodes ;
	int	*nodeNums;
	int     *nodeMarkers;
	int     *nodeBoundaryMarkers;
   	int     **elems ;
	double  *ElemArea;
} trimesh ;


/* GLOBAL VARIABLE DECLARATIONS */
trimesh 	Trimesh ;
char *FLOAT_OUTPUT_FORMAT_SPEC = "%23.16E"; // field width is 26, 16 decimal digits, exponential notation with big E
char *INT_OUPTUT_FORMAT_SPEC = "%8d"; // Field width is 8, number is integer
char    casename[25];

/* DECLARATIONS FOR FUNCTIONS */
void 	output_triangle_area_file(void) ;
void    read_trimesh(void) ;
void	cramer_rule(double [3][4], double [3]);

void main()
{
	
   printf("\n\n----------WELCOME TO TRIANGLE MESH REFINEMENT PROGRAM----------\n\n");
   printf("This Code reads data from a 2D Triangle mesh and makes a .area file for producing a refined mesh\n") ;
   printf("Run this code in the directory where the Triangle mesh files are located.\n") ;

  // Read in trimesh format
  read_trimesh();
  
  // Print .area file
  output_triangle_area_file();
}




void read_trimesh()
{

        // Local Variable Declarations
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT,nAttributes ;
	int	debug = 0; // Debugging flag
	int	v1,v2,v3 ;
	double	x,y ;
	char	filename[30] ;
	FILE	*fp ;

	printf("----------STARTING TRIANGLE MESH READER----------\n\n");
	printf("The Code reads triangle meshes created by the program Triangle. \n\n") ;

        /* Store casename provided by user */
        printf("What is file casename?  i.e. <casename>.node:");
        scanf("%s",casename);

        /* Echo user input */
        printf("Selected casename is: %s \n",casename);

	/* Read node data */
  	sprintf(filename,"%s.node",casename) ;
	printf("Reading File: %s \n\n",filename) ; //Echo filename

      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d %d %d",&nNodes,&nDims,&nAttrs,&MarkerFlag) ; // Read first line of .node file

	/* Echo data read from .node file */
	printf("Node File Data\n");
	printf("# of Nodes: %d \n# of Dimensions: %d\nAttribute Flag: %d \nMarker Flag: %d \n\n",nNodes,nDims,nAttrs,MarkerFlag) ;

	Trimesh.nNodes = nNodes ; // Store # of nodes into Trimesh variable

	/* Allocate for nodes array(1D array) */
	Trimesh.nodes = (vertex *) malloc(nNodes*sizeof(vertex)) ;

	/* Allocate for node numbering array(1D array) */
	Trimesh.nodeNums = malloc(nNodes*sizeof(int));

	/* Allocate for nodeMarkers array(1D array) */
	Trimesh.nodeMarkers = malloc(nNodes*sizeof(int));
        
	/* Loop through file and read all node data */
	for(i=0; i<nNodes; i++)
	{
	
	  if(MarkerFlag == 1 ) // additional marker line is present
	  {
	  	fscanf(fp,"%d %lf %lf %d",&dummy,&x,&y,&dummy) ;
	  	//printf("%lf %lf \n",x,y) ; //Echo data from file
          	Trimesh.nodes[i].x = x ;
          	Trimesh.nodes[i].y = y ;
		Trimesh.nodeMarkers[i] = dummy;
		Trimesh.nodeNums[i] = i+1;
	  }
	  else
	  {
		fscanf(fp,"%d %lf %lf",&dummy,&x,&y) ;
                //printf("%lf %lf \n",x,y) ; //Echo data from file
                Trimesh.nodes[i].x = x ;
                Trimesh.nodes[i].y = y ;
		Trimesh.nodeNums[i] = i+1;

	  }


	}

	fclose(fp) ;

	/* Read elements */
  	sprintf(filename,"%s.ele",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename to be read

      	fp = fopen(filename,"r") ; // Open .ele file for reading elements

	fscanf(fp,"%d %d %d",&nElems,&nNpT,&nAttributes) ; // Read first line
	printf("Element File Data\n");
	printf("# of Elements: %d\n# of Verts Per Element:  %d\n# of Attributes: %d \n\n",nElems,nNpT,nAttributes) ; //Echo first line to user

	Trimesh.nElems = nElems ; // Store number of elements into Tetmesh variable

	/* Allocate for elements array ( 2D array) */
	Trimesh.elems = (int **) malloc(nElems*sizeof(int *)) ; // Allocate rows of array
        
        /* Loop through and store element vertices */
	for(i=0; i<nElems; i++)
	{
	  fscanf(fp,"%d %d %d %d",&dummy,&v1,&v2,&v3) ;
	  //printf("%d %d %d %d \n",dummy,v1,v2,v3) ;
	  Trimesh.elems[i] = (int *) malloc(nNpT*sizeof(int)) ; // Allocate columns of array

          /* Check this section if grid ouput is not correct(adjust clockwise/counter-clockwise read) */
	  Trimesh.elems[i][0] = v2 ; 
          Trimesh.elems[i][1] = v1 ;
          Trimesh.elems[i][2] = v3 ;
	}

	fclose(fp) ;

	printf("Finished Reading Triangle Mesh .node and .ele Data Files\n\n");
	
	/* DEBUGGING */
	if(debug == 1)//print entire contents of quadmesh data structure
        {

                printf("\n\nPrinting Contents of Trimesh Data Structure");
                printf("\nNumber of Nodes: %d",Trimesh.nNodes);
                printf("\nNumber of Elements: %d",Trimesh.nElems);

                printf("\n\n------------Printing Node Data------------");
                printf("\nNode # \t X Coord \t Y Coord");
                for(i=0;i<Trimesh.nNodes;i++)
                {
                        printf("\n %d %4.2e %4.2e",Trimesh.nodeNums[i],Trimesh.nodes[i].x,Trimesh.nodes[i].y);

                }

        }


}





void output_triangle_area_file(void)
{
	int	i, j,k;
	int 	AMR_Flag;
	double  xc, yc, xP, yP, r, R0, A0, Area_Percent,Amin_Factor,Amax_inner;
	double  Amin, Amax, rmin,rmax;
	double  *Area, *radial_pos;//area and radial location of elements
	char	filename[30];
	FILE	*fp ;

/*************************************************************/
/*               	Open file                            */
/*************************************************************/

	/* Read Initial Parameters from Input File */
	printf("\n%s\n","Reading Data from input files: amr_inputs.dat");
	printf("%s\n","Reading X,Y coordinates of circle Center and the AMR Flag");
        fp = fopen("amr_inputs.dat","r") ;
	
	fscanf(fp,"%lf %lf %d",&xc,&yc,&AMR_Flag);
	
	/* Read additional Parameters */
	if(AMR_Flag == 1)//linear element scaling, fixed areas
	{
		printf("\n%s\n","Linear Element Scaling with Mesh Dependent Minimum Area Selected");
	}
	else if(AMR_Flag == 2)//quadratic element scaling, fixed areas 
	{
		printf("\n%s\n","Quadratic Element Scaling with Mesh Dependent Minimum Area Selected");

		printf("\n%s\n%s\n%s","Reading:"," Refinement Radius","Area Percent Threshold");
		fscanf(fp,"%lf %lf",&R0,&Area_Percent);
	}
	else if(AMR_Flag == 3)// Prescribed Minimum Area With Linear Variation
	{
		printf("\n%s\n","Linear Element Scaling with Mesh Dependent Minimum Area, but Prescribed Minimum Area Scaling Selected");

		printf("\n%s\n%s\n%s","Reading:","Minimum Area Factor(0 to 1)");
		fscanf(fp,"%lf",&Amin_Factor);
	}
	else if(AMR_Flag == 4)//Qaudratic element scaling, prescribed minimum area
	{
		printf("\n%s\n","Quadratic Element Scaling with Mesh Dependent Minimum Area, but Prescribed Minimum Area Scaling Selected");

		printf("\n%s\n%s\n%s","Reading:"," Refinement Radius"," Area Percent Threshold");
                fscanf(fp,"%lf %lf",&R0,&Area_Percent);

		printf("\n%s\n"," Minimum Area Factor(0 to 1)");
                fscanf(fp,"%lf",&Amin_Factor);
	}
	else if(AMR_Flag == 5)//Quadratic element scaling, constant inner region element size
	{
	
		printf("\n%s\n","Quadratic Element Scaling with Prescribed Uniform Mesh Size inside A Prescribed Refinement Radius Selected");

                printf("\n%s\n%s","Reading:"," Refinement Radius");
                fscanf(fp,"%lf ",&R0);

                printf("\n%s\n"," Maximum Area Constraint");
                fscanf(fp,"%lf",&Amax_inner);

	}
	else if(AMR_Flag == 6)//Quadratic element scaling, linear inner region element scaling
        {

                printf("\n%s\n","Quadratic Element Scaling with Linear Mesh Scaling inside a Prescribed Refinement Radius Selected");

                printf("\n%s\n%s\n%s","Reading:"," Refinement Radius"," Area Percent Threshold");
                fscanf(fp,"%lf %lf",&R0,&Area_Percent);

        }


	//Echo User Inputs Back
	printf("\n%s\n","User Input Is:");
	printf("%s %lf\n","X-Center Coord: ",xc);
	printf("%s %lf\n","Y-Center Coord: ",yc);
	printf("%s %d\n","Element Scaling Flag: ",AMR_Flag);

	if(AMR_Flag == 2)
	{
	   printf("%s %lf\n","Refinement Radius: ",R0);
 	   printf("%s %lf\n","Percent Area Inside Refinement Radius: ",Area_Percent);
	}
	else if(AMR_Flag == 3)
	{
	   printf("%s %lf\n","Minimum Area Factor(0 to 1): ",Amin_Factor);
	}
	else if(AMR_Flag == 4)
	{
	   printf("%s %lf\n","Refinement Radius: ",R0);
           printf("%s %lf\n","Percent Area Inside Refinement Radius: ",Area_Percent);
	   printf("%s %lf\n","Minimum Area Factor(0 to 1): ",Amin_Factor);
	}
	else if(AMR_Flag == 5)
	{
	   printf("%s %lf\n","Refinement Radius: ",R0);
           printf("%s %lf\n","Maximum Area Constraint Inside Refinement Radius: ",Amax_inner);

	}
	else if(AMR_Flag == 6)
        {
           printf("%s %lf\n","Refinement Radius: ",R0);
           printf("%s %lf\n","Percent Area Inside Refinement Radius: ",Area_Percent);

        }


	fclose(fp); // Close input file
	
	printf("\n\nComputing Triangle area constraints for .area file...\n") ;

	sprintf(filename,"%s.area",casename) ;

	printf("Data will be written to: \"%s\" \n",filename) ;
	fp = fopen(filename,"w") ;


//===================================================================
//	Compute Element Areas and Print to File
//===================================================================

	/* Allocate array to hold element areas */
	Trimesh.ElemArea = malloc(Trimesh.nElems*sizeof(double));

	Area = malloc(Trimesh.nElems*sizeof(double));
	radial_pos = malloc(Trimesh.nElems*sizeof(double));

	printf("%s\n","Allocation of Arrays Complete");

	for(i=0;i<Trimesh.nElems;i++) // Loop over all elements
	{
		/* Initialize centroid location */
		xP = 0;
		yP = 0;
		double x[3],y[3]; // Array to hold node coordinates

		for(j=0;j<3;j++)
		{
		  /* Compute X and Y coordinates of element centroid */
		  xP = xP +(1.0/3.0)*( Trimesh.nodes[ Trimesh.elems[i][j]-1 ].x );
		  yP = yP +(1.0/3.0)*( Trimesh.nodes[ Trimesh.elems[i][j]-1 ].y );	

		  x[j] =  Trimesh.nodes[ Trimesh.elems[i][j]-1 ].x;
		  y[j] =  Trimesh.nodes[ Trimesh.elems[i][j]-1 ].y;	

		}

		
		/* Compute Element area using Heron's Formula */
		double a,b,c,s;
		a = sqrt( (x[0]-x[1])*(x[0]-x[1]) + (y[0] - y[1])*(y[0] - y[1]) );
		b = sqrt( (x[0]-x[2])*(x[0]-x[2]) + (y[0] - y[2])*(y[0] - y[2]) );
		c = sqrt( (x[1]-x[2])*(x[1]-x[2]) + (y[1] - y[2])*(y[1] - y[2]) );
		s = 0.5*(a+b+c);
		
		Area[i] = sqrt( s*(s-a)*(s-b)*(s-c) );	//Heron's Formula	
		
		Trimesh.ElemArea[i] = Area[i]; // Store element area;

		/* Compute element distance from center provided by user */
		radial_pos[i] = sqrt( (xc-xP)*(xc-xP)  + (yc-yP)*(yc-yP) );
	
		/* Compare computed values to find max and min of r and Area */
		if( i == 1 ) // first iteration assume values of max and min are same
		{
			Amin = Area[i];
			Amax = Area[i];
			rmin = radial_pos[i];
			rmax = radial_pos[i];
		}
		else // update max and min values
		{
			if(Area[i] < Amin)
			{
				Amin = Area[i];
			}

			if(Area[i] > Amax)
			{
				Amax = Area[i];
			}
			
			if(radial_pos[i] < rmin)
			{
				rmin = radial_pos[i];
			}
			
			if(radial_pos[i] > rmax)
			{
				rmax = radial_pos[i];
			}

		}

	}


	//Define fitting constants
	double C1,C2,C3,C4,C5,C6;

	if(AMR_Flag == 1) //Linear Variation
	{
	   
	   //Echo curve fit parameters
           printf("\n%s %4.2E \n%s %4.2E \n","Mininum Area: ",Amin,"Maximum Area: ",Amax);
           printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);
           
           //Linear Variation
           C1 = (Amax-Amin)/(rmax-rmin);
           C2 = Amin - C1*rmin;
           
           printf("\n%s\n","Coefficients to linear area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n","A(r) = C1r +C2","C1 = ",C1,"C2 = ",C2);
                                  
	}
	else if(AMR_Flag == 2) // Quadratic variation
	{

	   /*Compute intermediate area*/
	   A0 = Amin + (Amax-Amin)*(Area_Percent/100);
	          		
	   /*Echo curve fit parameters*/
	   printf("\n%s %4.2E \n%s %4.2E \n%s %4.2E \n","Mininum Area: ",Amin,"Critical Area: ",A0,"Maximum Area: ",Amax);
	   printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);
	   

	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n\n%s\n","Solving the following Equations to Scale Element Areas");
	   printf("%s\n%s\n%s\n%s\n","OUTER SOLUTION","2*C1*rmax + C2 = 0","C1*rmax^2 + C2*rmax + C3 = Amax","C1*R0^2 + C2*R0 + C3 = A0");

	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n%s\n%s\n%s\n%s\n","INNER SOLUTION","2*C1*rmin + C2 = 0","C1*rmin^2 + C2*rmin + C3 = Amin","C1*R0^2 + C2*R0 + C3 = A0");
	  

           double A[3][4];
	   double x[3];
		
	   /*Construct the augmented matrix for outer solution coefficients*/
	   A[0][0] = 2*rmax;
	   A[0][1] = 1;
	   A[0][2] = 0;
	   A[0][3] = 0;

	   A[1][0] = rmax*rmax;
           A[1][1] = rmax;
           A[1][2] = 1;
           A[1][3] = Amax;

	   A[2][0] = R0*R0;
	   A[2][1] = R0;
	   A[2][2] = 1;
	   A[2][3] = A0;

	   /*Call CRAMER'S RULE solver to get coefficients*/
	   cramer_rule(A,x); // x has the coefficients
	   C1 = x[0];
	   C2 = x[1];
	   C3 = x[2];
	   

	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for OUTER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C1 + [%4.2E]C2 + [%4.2E]C3 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



	   /*Construct Augmented Matrix for inner solution coefficients*/

	   A[0][0] = 2*rmin;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmin*rmin;
           A[1][1] = rmin;
           A[1][2] = 1;
           A[1][3] = Amin;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = A0;	   


	   /*Call CRAMER'S RULE solver to get coefficients*/
	   cramer_rule(A,x); // x has the coefficients
 	   
	   C4 = x[0];
	   C5 = x[1];
	   C6 = x[2];

	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
	   printf("\n\n%s\n","Coefficients for INNER SOLUTION...");
	   for(i=0;i<3;i++) // Print Matrix
	   { 
		printf("[%4.2E]C4 + [%4.2E]C5 + [%4.2E]C6 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
	   }
		

	   printf("\n%s\n","Coefficients to the OUTER quadratic area sizing function ");
 	   printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","A(r) = C1r^2 +C2r+C3","C1 = ",C1,"C2 = ",C2,"C3 = ",C3);

           printf("\n\n%s\n","Coefficients to the INNER quadratic area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","A(r) = C4r^2 +C5r+C6","C4 = ",C4,"C5 = ",C5,"C6 = ",C6);

	
	   printf("\n\n%s","Predicted Areas Using Sizing Function(Check to see if it matches input)");
	   printf("\n%s","OUTER area sizing function predictions");
	   printf("\n%s %4.2E","Max Area: ",C1*rmax*rmax+C2*rmax+C3);
	   printf("\n%s %4.2E\n","Critical Area: ",C1*R0*R0+C2*R0+C3);

           printf("\n\n%s","INNER area sizing function predictions");
           printf("\n%s %4.2E","Min Area: ",C4*rmin*rmin+C5*rmin+C6);
           printf("\n%s %4.2E\n","Critical Area: ",C4*R0*R0+C5*R0+C6);


	}
	else if(AMR_Flag == 3) // Prescribed Minimum Area With Linear Variation
	{
	
	   /*Update minimum area using user-defined factor*/
	   Amin = Amin*(1-Amin_Factor/100);
	
	   /*Echo curve fit parameters*/
	   printf("\n%s %4.2E \n%s %4.2E \n","Mininum Area: ",Amin,"Maximum Area: ",Amax);
           printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);

	   /*Linear Variation*/
	   C1 = (Amax-Amin)/(rmax-rmin);
           C2 = Amin - C1*rmin;

           printf("\n%s\n","Coefficients to linear area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n","A(r) = C1r +C2","C1 = ",C1,"C2 = ",C2);

		

	}
	else if(AMR_Flag == 4)
	{
	

	   Amin = Amin*Amin_Factor; // Update Minimum Area

	   /*Compute intermediate area*/
	   A0 = Amin + (Amax-Amin)*(Area_Percent/100);
	   
	   /*Echo curve fit parameters*/
	   printf("\n%s %4.2E \n%s %4.2E \n%s %4.2E \n","Mininum Area: ",Amin,"Critical Area: ",A0,"Maximum Area: ",Amax);
	   printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);
	   
	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n\n%s\n","Solving the following Equations to Scale Element Areas");
	   printf("%s\n%s\n%s\n%s\n","OUTER SOLUTION","2*C1*rmax + C2 = 0","C1*rmax^2 + C2*rmax + C3 = Amax","C1*R0^2 + C2*R0 + C3 = A0");
	   
	   /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
	   printf("\n%s\n%s\n%s\n%s\n","INNER SOLUTION","2*C1*rmax + C2 = 0","C1*rmin^2 + C2*rmin + C3 = Amin","C1*R0^2 + C2*R0 + C3 = A0");
	   
	   
	   double A[3][4];
	   double x[3];
	   
	   /*Construct the augmented matrix for outer solution coefficients*/
	   A[0][0] = 2*rmax;
	   A[0][1] = 1;
	   A[0][2] = 0;
	   A[0][3] = 0;
	   
	   A[1][0] = rmax*rmax;
	   A[1][1] = rmax;
	   A[1][2] = 1;
	   A[1][3] = Amax;
	   
	   A[2][0] = R0*R0;
	   A[2][1] = R0;
	   A[2][2] = 1;
	   A[2][3] = A0;
	   
	   /*Call CRAMER'S RULE solver to get coefficients*/
	   cramer_rule(A,x); // x has the coefficients
	   C1 = x[0];
	   C2 = x[1];
	   C3 = x[2];                                                                                                                                                                                    

  	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for OUTER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C1 + [%4.2E]C2 + [%4.2E]C3 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



           /*Construct Augmented Matrix for inner solution coefficients*/

           A[0][0] = 2*rmin;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmin*rmin;
           A[1][1] = rmin;
           A[1][2] = 1;
           A[1][3] = Amin;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = A0;


           /*Call CRAMER'S RULE solver to get coefficients*/
           cramer_rule(A,x); // x has the coefficients

           C4 = x[0];
           C5 = x[1];
           C6 = x[2];


	   
	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for INNER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C4 + [%4.2E]C5 + [%4.2E]C6 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



           printf("\n%s\n","Coefficients to the OUTER quadratic area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","A(r) = C1r^2 +C2r+C3","C1 = ",C1,"C2 = ",C2,"C3 = ",C3);

           printf("\n\n%s\n","Coefficients to the INNER quadratic area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","A(r) = C4r^2 +C5r+C6","C4 = ",C4,"C5 = ",C5,"C6 = ",C6);


           printf("\n\n%s","Predicted Areas Using Sizing Function(Check to see if it matches input)");
           printf("\n%s","OUTER area sizing function predictions");
           printf("\n%s %4.2E","Max Area: ",C1*rmax*rmax+C2*rmax+C3);
           printf("\n%s %4.2E\n","Critical Area: ",C1*R0*R0+C2*R0+C3);

           printf("\n\n%s","INNER area sizing function predictions");
           printf("\n%s %4.2E","Min Area: ",C4*rmin*rmin+C5*rmin+C6);
           printf("\n%s %4.2E\n","Critical Area: ",C4*R0*R0+C5*R0+C6);

 
	}
	else if(AMR_Flag ==5)
	{
	
	   Amin = Amax_inner; // Update Minimum Area

           /* Intermediate area is dictacted by minimum area now */
           A0 = Amin;

           /* Echo curve fit parameters */
           printf("\n%s %4.2E \n%s %4.2E \n%s %4.2E \n","Mininum Area: ",Amin,"Critical Area: ",A0,"Maximum Area: ",Amax);
           printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);

           /* Quadratic variation with coefficients solved by CRAMER'S RULE */
           printf("\n\n%s\n","Solving the following Equations to Scale Element Areas");
           printf("%s\n%s\n%s\n%s\n","OUTER SOLUTION","2*C1*rmax + C2 = 0","C1*rmax^2 + C2*rmax + C3 = Amax","C1*R0^2 + C2*R0 + C3 = A0");

           /* Constant Inner Mesh Size Variation Governed By User Input */
           printf("\n%s\n%s","INNER SOLUTION","A=Amin (constant)");


           double A[3][4];
           double x[3];

           /*Construct the augmented matrix for outer solution coefficients*/
           A[0][0] = 2*rmax;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmax*rmax;
           A[1][1] = rmax;
           A[1][2] = 1;
           A[1][3] = Amax;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = A0;

           /*Call CRAMER'S RULE solver to get coefficients*/
           cramer_rule(A,x); // x has the coefficients
           C1 = x[0];
           C2 = x[1];
           C3 = x[2];   
	 
	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for OUTER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C1 + [%4.2E]C2 + [%4.2E]C3 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }


	   /*Construct Augmented Matrix for inner solution coefficients*/

           A[0][0] = 2*rmin;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmin*rmin;
           A[1][1] = rmin;
           A[1][2] = 1;
           A[1][3] = Amin;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = Amin;


           /*Call CRAMER'S RULE solver to get coefficients*/
           cramer_rule(A,x); // x has the coefficients

           C4 = x[0];
           C5 = x[1];
           C6 = x[2];

	   
	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for INNER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C4 + [%4.2E]C5 + [%4.2E]C6 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



           printf("\n%s\n","Coefficients to the OUTER quadratic area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","A(r) = C1r^2 +C2r+C3","C1 = ",C1,"C2 = ",C2,"C3 = ",C3);

           printf("\n\n%s\n","Coefficients to the INNER quadratic area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","A(r) = C4r^2 +C5r+C6","C4 = ",C4,"C5 = ",C5,"C6 = ",C6);


           printf("\n\n%s","Predicted Areas Using Sizing Function(Check to see if it matches input)");
           printf("\n%s","OUTER area sizing function predictions");
           printf("\n%s %4.2E","Max Area: ",C1*rmax*rmax+C2*rmax+C3);
           printf("\n%s %4.2E\n","Critical Area: ",C1*R0*R0+C2*R0+C3);

           printf("\n\n%s","INNER area sizing function predictions");
           printf("\n%s %4.2E","Min Area: ",C4*rmin*rmin+C5*rmin+C6);
           printf("\n%s %4.2E\n","Critical Area: ",C4*R0*R0+C5*R0+C6);



	}
	else if(AMR_Flag ==6)
	{
	   

           /*Compute intermediate area*/
           A0 = Amin + (Amax-Amin)*(Area_Percent/100);

           /*Echo curve fit parameters*/
           printf("\n%s %4.2E \n%s %4.2E \n%s %4.2E \n","Mininum Area: ",Amin,"Critical Area: ",A0,"Maximum Area: ",Amax);
           printf("%s %f \n%s %f\n","Minimum Radius: ",rmin,"Maximum Radius: ", rmax);

           /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
           printf("\n\n%s\n","Solving the following Equations to Scale Element Areas");
           printf("%s\n%s\n%s\n%s\n","OUTER SOLUTION","2*C1*rmax + C2 = 0","C1*rmax^2 + C2*rmax + C3 = Amax","C1*R0^2 + C2*R0 + C3 = A0");

           /*Quadratic variation with coefficients solved by CRAMER'S RULE*/
           printf("\n%s\n%s\n%s\n%s\n","INNER SOLUTION","2*C1*rmin + C2 = 0","C1*rmin^2 + C2*rmin + C3 = Amin","C1*R0^2 + C2*R0 + C3 = A0");


           double A[3][4];
           double x[3];

           /*Construct the augmented matrix for outer solution coefficients*/
           A[0][0] = 2*rmax;
           A[0][1] = 1;
           A[0][2] = 0;
           A[0][3] = 0;

           A[1][0] = rmax*rmax;
           A[1][1] = rmax;
           A[1][2] = 1;
           A[1][3] = Amax;

           A[2][0] = R0*R0;
           A[2][1] = R0;
           A[2][2] = 1;
           A[2][3] = A0;

           /*Call CRAMER'S RULE solver to get coefficients*/
           cramer_rule(A,x); // x has the coefficients
           C1 = x[0];
           C2 = x[1];
           C3 = x[2];                                             

	
	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for OUTER SOLUTION...");
           for(i=0;i<3;i++) // Print Matrix
           {
                printf("[%4.2E]C1 + [%4.2E]C2 + [%4.2E]C3 = [%4.2E]\n",A[i][0],A[i][1],A[i][2],A[i][3]);
           }



           /*Solve for Linear  inner solution coefficients(Easy to do, no linear algebrea needed)*/
           C4 = (A0-Amin)/(R0-rmin);
           C5 = Amin-C4*rmin;


	   printf("\n%s","Numerical Coeffients for system of equations are the following:");
           printf("\n\n%s\n","Coefficients for INNER SOLUTION...");
          
	   printf("[%4.2E]C4 + [%4.2E]C5  = [%4.2E]\n",R0,1.0,A0);
	   printf("[%4.2E]C4 + [%4.2E]C5  = [%4.2E]\n",rmin,1.0,Amin);


           printf("\n%s\n","Coefficients to the OUTER quadratic area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E \n%s %10.3E","A_outer(r) = C1r^2 +C2r+C3","C1 = ",C1,"C2 = ",C2,"C3 = ",C3);

           printf("\n\n%s\n","Coefficients to the INNER quadratic area sizing function ");
           printf("%s \n%s %10.3E \n%s %10.3E ","A_inner(r) = C4r + C5","C4 = ",C4,"C5 = ",C5);


           printf("\n\n%s","Predicted Areas Using Sizing Function(Check to see if it matches input)");
           printf("\n%s","OUTER area sizing function predictions");
           printf("\n%s %4.2E","Max Area: ",C1*rmax*rmax+C2*rmax+C3);
           printf("\n%s %4.2E\n","Critical Area: ",C1*R0*R0+C2*R0+C3);

           printf("\n\n%s","INNER area sizing function predictions");
           printf("\n%s %4.2E","Min Area: ",C4*rmin + C5);
           printf("\n%s %4.2E\n","Critical Area: ",C4*R0 + C5);


	}


	/*************UPDATE ELEMENT AREAS***************************************/
	printf("\n%s","Updating Element Areas...");
	for(i = 0; i<Trimesh.nElems;i++) // Scale all areas between min and max according to r
	{

		if(AMR_Flag == 1) // Linear+Unaltered Min Area
		{
		  Area[i] = C1*radial_pos[i] + C2;
		}
		else if (AMR_Flag == 2)//Quadratic+Unaltered Min Area
		{
		
		  if(radial_pos[i] <= R0) //Scale Inner
		  {
		    Area[i] = C4*radial_pos[i]*radial_pos[i] + C5*radial_pos[i] + C6;
		  }
		  else // Scale Outer
		  {
		    Area[i] = C1*radial_pos[i]*radial_pos[i] + C2*radial_pos[i] + C3;
		  }

		}
		else if (AMR_Flag == 3)//Linear+Altered Min Area
		{
		   Area[i] = C1*radial_pos[i] + C2;
		
		}
		else if(AMR_Flag == 4) //Quadratic+Altered Min Area
		{

		  
             	   if(radial_pos[i] <= R0) // Scale Inner
                   {
                     Area[i] = C4*radial_pos[i]*radial_pos[i] + C5*radial_pos[i] + C6;
                   }
                   else // Scale Outer
                   {
                     Area[i] = C1*radial_pos[i]*radial_pos[i] + C2*radial_pos[i] + C3;
                   }

		}
		else if(AMR_Flag == 5) // Quadratic+Uniform inner
		{
		
		  if(radial_pos[i] <= R0) // Scale Inner
                   {
                     Area[i] = C4*radial_pos[i]*radial_pos[i] + C5*radial_pos[i] + C6;
                   }
                   else // Scale Outer
                   {
                     Area[i] = C1*radial_pos[i]*radial_pos[i] + C2*radial_pos[i] + C3;
                   }


		}
		else if(AMR_Flag == 6) // Quadratic+Uniform inner
                {

                  if(radial_pos[i] <= R0) // Scale Inner
                   {
                     Area[i] = C4*radial_pos[i]+ C5;
                   }
                   else // Scale Outer
                   {
                     Area[i] = C1*radial_pos[i]*radial_pos[i] + C2*radial_pos[i] + C3;
                   }


                }


	}

	fprintf(fp,"%8d\n",Trimesh.nElems);
	
	printf("%s\n","Printing to File Now");

	/* Print output to file */
	for(i=0;i<Trimesh.nElems;i++)
	{
	   fprintf(fp,"%8d %23.16E",i+1,Area[i]);
	   fprintf(fp,"\n");
	}

//**************************************************************************
// 	Close file
//**************************************************************************

	fclose(fp) ;

//**************************************************************************
//      Print trimesh object contents
//**************************************************************************

  	printf("Writing Adaptive Meshed Triangle .area file done.\n") ;

  	return ;
}


void    cramer_rule(double A[3][4], double x[3])
{

	/* Local variable declarations*/
	double D,D1,D2,D3;
	double first, second, third;;

	// The Ds are determinants, and first, second, and third are the
	// three parts of the 3x3 determinant

	first=    A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2]);
	second =  -1*A[0][1]*(A[1][0]*A[2][2]-A[2][0]*A[1][2]);
	third =   A[0][2]*(A[1][0]*A[2][1]-A[2][0]*A[1][1]);
	
	D = first + second + third;


	first=    A[0][3]*(A[1][1]*A[2][2]-A[2][1]*A[1][2]);
        second =  -1*A[0][1]*(A[1][3]*A[2][2]-A[2][3]*A[1][2]);
        third =   A[0][2]*(A[1][3]*A[2][1]-A[2][3]*A[1][1]);

        D1 = first + second + third;


	first=    A[0][0]*(A[1][3]*A[2][2]-A[2][3]*A[1][2]);
        second =  -1*A[0][3]*(A[1][0]*A[2][2]-A[2][0]*A[1][2]);
        third =   A[0][2]*(A[1][0]*A[2][3]-A[2][0]*A[1][3]);

        D2 = first + second + third;


	first=    A[0][0]*(A[1][1]*A[2][3]-A[2][1]*A[1][3]);
        second =  -1*A[0][1]*(A[1][0]*A[2][3]-A[2][0]*A[1][3]);
        third =   A[0][3]*(A[1][0]*A[2][1]-A[2][0]*A[1][1]);

        D3 = first + second + third;

	
	//Compute solutions using Cramer's Rule
	x[0] = D1/D;
	x[1] = D2/D;
	x[2] = D3/D;

}
