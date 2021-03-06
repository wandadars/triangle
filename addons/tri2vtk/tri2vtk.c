/***********************************************************************/
/*	Purpose: Code to convert an 2D unstructured triangle mesh that */
/*	         was generated by the mesh generation program Triangle */
/*	         to the VTK legacy format.                             */
/*	         	                                               */
/*	Author 	: Christopher Neal				       */
/*	Date 	: October 23rd, 2014 	        		       */
/***********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct
{
	double  	x, y; // Coordinates of vertex
} vertex ;


typedef struct
{
	int     nNodes, nElems, nEdges, nEdgeMarkers, nNodeMarkers;
	vertex 	*nodes ;
	int	nNpts;
	int	*nodeNums;
	int     *nodeMarkers;
	int     *nodeBoundaryMarkers;
	int	**nodeBoundaryCounts;
   	int     **elems ;
} trimesh ;


/* GLOBAL VARIABLE DECLARATIONS */
trimesh 	Trimesh ;
char *FLOAT_OUTPUT_FORMAT_SPEC = "%23.16E"; // field width is 26, 16 decimal digits, exponential notation with big E
char *INT_OUPTUT_FORMAT_SPEC = "%8d"; // Field width is 8, number is integer
double EPSILON = 1.0e-9;

char    casename[25]; //Filename used in both functions

/* DECLARATIONS FOR FUNCTIONS */
void 	output_vtk_format(void) ;
void    read_trimesh(void) ;


void main()
{
	

   printf("----------WELCOME TO TRI2VTK MESH FILTER PROGRAM----------\n\n");
   printf("The Code reads in data from a 2D Triangle  mesh \n") ;
   printf("and outputs the result to the legacy .vtk format for use in Paraview \n\n") ;

  /* Read in trimesh format */
  read_trimesh();

  /* Print hybrid mesh format in GRD format */
  output_vtk_format();
}




void read_trimesh()
{

        /* Local Variable Declarations */
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT,nAttributes ;
	int     nEdges ;
	int	debug = 0; // Debugging flag
	int	v1,v2,v3 ;
	double	x,y ;
	char	filename[30] ;
	char    line[128] ;
	FILE	*fp ;

	printf("----------STARTING TRIMESH READER----------\n\n");
	printf("The Code reads triangle meshes created by the program Triangle. \n\n") ;

        /*Store casename provided by user*/
        printf("What is file casename?  i.e. <casename>.node:");
        scanf("%s",casename);

        /*Echo user input*/
        printf("Selected casename is: %s \n",casename);

	/* Read nodes ****************************************************/
  	sprintf(filename,"%s.node",casename) ;
	printf("Reading File: %s \n\n",filename) ; //Echo filename

      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d %d %d",&nNodes,&nDims,&nAttrs,&MarkerFlag) ; // Read first line of .node file

	/*Echo data read from .node file*/
	printf("Node File Data\n");
	printf("# of Nodes: %d \n# of Dimensions: %d\nAttribute Flag: %d \nMarker Flag: %d \n\n",nNodes,nDims,nAttrs,MarkerFlag) ;

	Trimesh.nNodes = nNodes ; // Store # of nodes into Trimesh variable
	
	/*Store the number of nodes per element*/
	Trimesh.nNpts = nDims;

	/*Allocate for nodes array(1D array)*/
	Trimesh.nodes = (vertex *) malloc(nNodes*sizeof(vertex)) ;

	/*Allocate for node numbering array(1D array)*/
	Trimesh.nodeNums = malloc(nNodes*sizeof(int));

	/*Allocate for nodeMarkers array(1D array)*/
	Trimesh.nodeMarkers = malloc(nNodes*sizeof(int));
        
	/*Loop through file and read all node data*/
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

	/* Read elements *************************************************/
  	sprintf(filename,"%s.ele",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename to be read

      	fp = fopen(filename,"r") ; // Open .ele file for reading elements

	fscanf(fp,"%d %d %d",&nElems,&nNpT,&nAttributes) ; // Read first line
	printf("Element File Data\n");
	printf("# of Elements: %d\n# of Verts Per Element:  %d\n# of Attributes: %d \n\n",nElems,nNpT,nAttributes) ; //Echo first line to user

	Trimesh.nElems = nElems ; // Store number of elements into Tetmesh variable

	/*Allocate for elements array ( 2D array)*/
	Trimesh.elems = (int **) malloc(nElems*sizeof(int *)) ; // Allocate rows of array
        
        /*Loop through and store element vertices*/
	for(i=0; i<nElems; i++)
	{
	  fscanf(fp,"%d %d %d %d",&dummy,&v1,&v2,&v3) ;
	  //printf("%d %d %d %d \n",dummy,v1,v2,v3) ;
	  Trimesh.elems[i] = (int *) malloc(3*sizeof(int)) ; // Allocate columns of array
          /*Check this section if grid ouput is not correct(adjust clockwise/counter-clockwise read)*/
	  Trimesh.elems[i][0] = v3 ; 
          Trimesh.elems[i][1] = v2 ;
          Trimesh.elems[i][2] = v1 ;
	}

	fclose(fp) ;
	
	/*DEBUGGING PURPOSES*/
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

	/*Count Number of Elements with Boundary Markers*/
	if(MarkerFlag==1)
	{

	}

	/*Print Out Number of Elements with Boundary Markers*/


}




void output_vtk_format(void)
{
	int	i, j, k ;
	int	triangle_type=5; /*Consult VTK manual for cell type definitions*/
	int	count,sum ;
        int 	loc = 0;    //row index of repeated node(used only if node is repeated)
	char	filename[20];
	FILE	*fp ;

//*************************************************************
// 	Open file and read title
//*************************************************************

	printf("\n\nWriting Legacy VTK grid file...\n") ;

	sprintf(filename,"%s.vtk",casename) ;
	printf("Writing Legacy VTK format to: \"%s\" \n",filename) ;
	fp = fopen(filename,"w") ;


//===================================================================
//      Print VTK Version and Identifier
//===================================================================

	fprintf(fp,"%s \n","# vtk DataFile Version 3.1") ;


/*===================================================================*/
/*      Print VTK Version and Identifier                             */
/*===================================================================*/

        fprintf(fp,"%s \n","Triangle Mesh File From: ") ;

/*===================================================================*/
/*      Print File Type                                              */
/*===================================================================*/

        fprintf(fp,"%s \n","ASCII") ;


/*===================================================================*/
/*      Print File Type                                              */
/*===================================================================*/

        fprintf(fp,"%s \n\n","DATASET UNSTRUCTURED_GRID");


/*===================================================================*/
/*      Print Number of Nodes and Node Type                          */
/*===================================================================*/

        fprintf(fp,"%s %d %s \n","POINTS ",Trimesh.nNodes," FLOAT") ;


/*===================================================================*/
/*      Print Node Coordinates                                       */
/*===================================================================*/

	
	for(j=0;j<Trimesh.nNodes;j++)
	{
		fprintf(fp,"%lf %lf %lf \n",Trimesh.nodes[j].x, Trimesh.nodes[j].y, 0.0);
	}

	fprintf(fp,"\n");

/*===================================================================*/
/*      Print Cell Information                                       */
/*===================================================================*/

        fprintf(fp,"%s %d %d\n","CELLS", Trimesh.nElems, Trimesh.nElems*(Trimesh.nNpts+2)) ;

        for(j=0;j<Trimesh.nElems;j++)
        {
		fprintf(fp,"%d ",Trimesh.nNpts+1);

		for(k = 0; k <Trimesh.nNpts+1 ;k++)
		{
                  fprintf(fp,"%d ",Trimesh.elems[j][k]-1);
		}
		
		fprintf(fp,"\n");
        }

	fprintf(fp,"\n");

/*===================================================================*/
/*      Print Cell Types                                             */
/*===================================================================*/


        fprintf(fp,"%s %d \n","CELL_TYPES",Trimesh.nElems) ;

        for(j=0;j<Trimesh.nElems;j++)
        {
                fprintf(fp,"%d \n",triangle_type);
        }


	fprintf(fp,"\n");

//**************************************************************************
// 	Close file
//**************************************************************************

	fclose(fp) ;

//**************************************************************************
//      Print trimesh object contents
//**************************************************************************

  	printf("Writing Legacy VTK grid file done.\n") ;

  	return ;
}

