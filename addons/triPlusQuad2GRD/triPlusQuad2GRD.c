/**********************************************************************/
/*	Purpose: Code to combine a structured quadrilateral mesh with */
/*	         an unstructured triangle mesh, and convert the output*/
/*	         to the .grd mesh format for input to Rocflu's conflu */
/*	         program.                         		      */
/*	         	                                              */
/*	 Notes: The meshes are 2D meshes. The structured mesh is a    */
/*	 	concentric circle mesh, and the unstructured mesh     */
/*	 	fills the inner circle of the unstrucuted concentric  */
/*	 	circle mesh.                                          */
/*	                                                              */
/*	Author 	: Christopher Neal				      */
/*	Date 	: July 24th, 2014 	        		      */
/**********************************************************************/

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
	int	*nodeNums;
	int     *nodeMarkers;
	int     *nodeBoundaryMarkers;
   	int     **elems ;
    	int     **edges;
	int     *edgeMarkers ;	//all edge boundary markers
	int	*edgeCounts;
	int	*edgeBoundaryMarkers;	//unique edge boundary markers
} trimesh ;

typedef struct
{
	
	int     nNodes, nElems, nEdges, nEdgeMarkers;
        vertex  *nodes ;
	int	*nodeNums;
        int     **elems ;
        int     **edges;
        int     *edgeMarkers ;
        int     *edgeCounts;
        int     *edgeBoundaryMarkers;

} quadmesh ;

/* GLOBAL VARIABLE DECLARATIONS */
trimesh 	Trimesh ;
quadmesh	Quadmesh;
char *FLOAT_OUTPUT_FORMAT_SPEC = "%23.16E"; // field width is 26, 16 decimal digits, exponential notation with big E
char *INT_OUPTUT_FORMAT_SPEC = "%8d"; // Field width is 8, number is integer
double EPSILON = 1.0e-9;
int **matchArray;    //Array to hold respective nodes that match and their numbers in each mesh

/* DECLARATIONS FOR FUNCTIONS */
void 	output_grd_format(int) ;
void    read_trimesh(void) ;
void    read_quadmesh(void) ;
int 	stitch_meshes(void);


void main()
{
	
   int nMatches;

   printf("----------WELCOME TO TRI-QUAD MESH STITCHING AND OUTPUT PROGRAM----------\n\n");
   printf("The Code reads in data from a triangle and quadrilateral mesh and combines\n") ;
   printf("the two meshes and outputs the result to a .grd format for use in conflu \n\n") ;

  // Read in trimesh format
  read_trimesh();

  // Read in quadmesh format
   read_quadmesh();

  // Stitch meshes together
  nMatches=stitch_meshes();

  // Print hybrid mesh format in GRD format
  output_grd_format(nMatches);
}




void read_trimesh()
{

        // Local Variable Declarations
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT,nAttributes ;
	int     nEdges ;
	int	debug = 0; // Debugging flag
	int	v1,v2,v3,v4 ;
	double	x,y,z ;
	char	casename[25],filename[30] ;
	char    line[128] ;
	FILE	*fp ;

	printf("----------STARTING TRIMESH READER----------\n\n");
	printf("The Code reads triangle meshes created by the program Triangle. \n\n") ;

        //Store casename provided by user
        printf("What is file casename?  i.e. <casename>.node:");
        scanf("%s",casename);

        // Echo user input
        printf("Selected casename is: %s \n",casename);

	/* Read nodes ****************************************************/
  	sprintf(filename,"%s.node",casename) ;
	printf("Reading File: %s \n\n",filename) ; //Echo filename

      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d %d %d",&nNodes,&nDims,&nAttrs,&MarkerFlag) ; // Read first line of .node file

	// Echo data read from .node file
	printf("Node File Data\n");
	printf("# of Nodes: %d \n# of Dimensions: %d\nAttribute Flag: %d \nMarker Flag: %d \n\n",nNodes,nDims,nAttrs,MarkerFlag) ;

	Trimesh.nNodes = nNodes ; // Store # of nodes into Trimesh variable

	// Allocate for nodes array(1D array)
	Trimesh.nodes = (vertex *) malloc(nNodes*sizeof(vertex)) ;

	//Allocate for node numbering array(1D array)
	Trimesh.nodeNums = malloc(nNodes*sizeof(int));

	//Allocate for nodeMarkers array(1D array)
	Trimesh.nodeMarkers = malloc(nNodes*sizeof(int));
        
	// Loop through file and read all node data
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

	// Allocate for elements array ( 2D array)
	Trimesh.elems = (int **) malloc(nElems*sizeof(int *)) ; // Allocate rows of array
        
        //Loop through and store element vertices
	for(i=0; i<nElems; i++)
	{
	  fscanf(fp,"%d %d %d %d",&dummy,&v1,&v2,&v3) ;
	  //printf("%d %d %d %d \n",dummy,v1,v2,v3) ;
	  Trimesh.elems[i] = (int *) malloc(3*sizeof(int)) ; // Allocate columns of array
          //Check this section if grid ouput is not correct(adjust clockwise/counter-clockwise read)
	  Trimesh.elems[i][0] = v2 ; 
          Trimesh.elems[i][1] = v1 ;
          Trimesh.elems[i][2] = v3 ;
	}

	fclose(fp) ;

	/* Read Edges ****************************************************/
  	sprintf(filename,"%s.edge",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename being read

      	fp = fopen(filename,"r") ;

	fscanf(fp,"%d %d",&nEdges,&MarkerFlag) ;//MarkerFlag-->boundary markers flag
	printf("Edge File Data\n");
	printf("# of Edges: %d\nBoundary Marker Flag: %d \n",nEdges,MarkerFlag) ; // Echo first line back to screen

	Trimesh.nEdges = nEdges ; // Store total number of edges in mesh
	

	// Allocate for edges array( 2D array)
	Trimesh.edges = (int **) malloc(nEdges*sizeof(int *)) ; // Allocate rows of edge array
	Trimesh.edgeMarkers = (int *) malloc(nEdges*sizeof(int)) ;// Allocate 1D array for all edge markers
	
	Trimesh.nEdgeMarkers = 1;	// Initialize unique edge boundary marker counter
	Trimesh.edgeBoundaryMarkers = malloc(1*sizeof(int)); // Allocate first element of unique edge boundary marker array

	for(i=0; i<nEdges; i++)
	{
	  fscanf(fp,"%d %d %d %d",&dummy,&v1,&v2,&marker) ;
	  //printf("%d %d %d %d \n",dummy,v1,v2,marker) ;
	  Trimesh.edges[i] = (int *) malloc(2*sizeof(int)) ; // Allocate columns for row i

          //Order of read matters for Clockwise/Anti-Clockwise ordering
          Trimesh.edges[i][0] = v2 ;
          Trimesh.edges[i][1] = v1 ;

          Trimesh.edgeMarkers[i] = marker ;
	  
	  if(i == 0) // set currentMarker to the first marker that is read from file
	  {
		Trimesh.edgeBoundaryMarkers[i] = marker;
       	  }


	  int newMarker = 1;	//Boolean for tracking if new marker is found. 0=false, 1=true
	  for(j=0;j<Trimesh.nEdgeMarkers;j++) // Compare currently read marker with markers that have been stored
	  {
		if(marker == Trimesh.edgeBoundaryMarkers[j])
		{
			newMarker = 0;
		}

	  }

	  if( newMarker == 1  ) // A new unique boundary marker is found. 
	  {
		Trimesh.nEdgeMarkers=Trimesh.nEdgeMarkers+1;

 		// Dynamically resize edgeBoundaryMarkers array to hold newly found boundary marker
		int *temp = realloc(Trimesh.edgeBoundaryMarkers,Trimesh.nEdgeMarkers*sizeof(int));  // Store data in a larger temp array
   		if(temp != NULL) // realloc was successful
		{
			Trimesh.edgeBoundaryMarkers = temp; //Point original array to newly sized array
			Trimesh.edgeBoundaryMarkers[Trimesh.nEdgeMarkers-1] = marker; // Store newly read marker into array
		}
		else
		{
			free(Trimesh.edgeBoundaryMarkers);
			printf("Error allocating memory for Edge Boundary Marker Array!\n");
			return ;
		}

	  }
		

}

	// Display the number of unique Edge markers in .edge file
	printf("\nNumber of Unique Boundary Edge Markers: %d \n",Trimesh.nEdgeMarkers);

	// Allocate edgeCounts array
	Trimesh.edgeCounts = malloc(Trimesh.nEdgeMarkers*sizeof(int));

	// Initialize edgeCounts array to zero. This array will hold the number
	// of edges that have a specific boundary marker.
	for(i=0; i<Trimesh.nEdgeMarkers; i++)
	{
		Trimesh.edgeCounts[i] = 0;
	}
	


	//Loop through edgeMarkers array and count the number of edges that have all of the markers
	//that are stored in edgeBoundaryMarkers
	for(i=0;i<nEdges; i++) // Loop through every edge
	{
		for(j=0;j<Trimesh.nEdgeMarkers; j++) //Loop through all edge boundary markers
		{

			if(Trimesh.edgeMarkers[i] == Trimesh.edgeBoundaryMarkers[j] ) // see which edgeBoundaryMarker the edgeMarker matches
			{				
				Trimesh.edgeCounts[j] = Trimesh.edgeCounts[j] + 1; // Increment corresponding edgeCount
			}			

		}

	}


	fclose(fp) ;


	// Print the number of edges on each edge boundary that were counted in file
	printf("\nNumber of Edges On Each Boundary\n");	
	int sum = 0;
	for( i=0;i<Trimesh.nEdgeMarkers;i++)
	{
	  	printf("Boundary Marker %d(%d) = %d \n",Trimesh.edgeBoundaryMarkers[i],i+1,Trimesh.edgeCounts[i]);
		sum = sum + Trimesh.edgeCounts[i];  
	}

	//Error Check, Print face counts to screen
	printf("\nTotal Edges Counted:%d\n# Edges In File:%d\n\n",sum,Trimesh.nEdges);

	printf("Finished Reading Triangle Mesh Data File\n\n");
	
	//DEBUGGING PURPOSES
	if(debug == 1)//print entire contents of quadmesh data structure
        {

                printf("\n\nPrinting Contents of Trimesh Data Structure");
                printf("\nNumber of Nodes: %d",Trimesh.nNodes);
                printf("\nNumber of Elements: %d",Trimesh.nElems);
                printf("\nNumber of Edges: %d",Trimesh.nEdges);
                printf("\nNumber of Unique Edge Boundary Markers: %d",Trimesh.nEdgeMarkers);


                printf("\n\n------------Printing Node Data------------");
                printf("\nNode # \t X Coord \t Y Coord");
                for(i=0;i<Trimesh.nNodes;i++)
                {
                        printf("\n %d %4.2e %4.2e",Trimesh.nodeNums[i],Trimesh.nodes[i].x,Trimesh.nodes[i].y);

                }


        }


}






void read_quadmesh()
{

        // Local Variable Declarations
	int	i,j,marker,dummy ;
	int	nNodes,nDims,nAttrs,MarkerFlag ;
	int 	nElems,nNpT ;
	int     nEdges ;
	int	v1,v2,v3,v4 ;
	int 	debug = 0; 	//Debug flag for printing data
	unsigned long line_count;
	double	x,y ;
	char	casename[25],filename[30] ;
	char 	dummychar[15];
	char    line[128] ;
	FILE	*fp ;

	printf("\n\n----------STARTING QUADMESH READER----------\n\n");
	printf("The Code reads quadrilateral meshes in StarCD format. \n\n") ;

        //Store casename provided by user
        printf("What is file casename?  i.e. <casename>.node: ");
        scanf("%s",casename);

        // Echo user input
        printf("Selected casename is: %s \n",casename);

	/* Read nodes ****************************************************/
  	sprintf(filename,"%s.vrt",casename) ;
	printf("Reading File: %s \n\n",filename) ; //Echo filename

      	fp = fopen(filename,"r") ; // open file to count lines

	/* Line Counting For Files ***************************************/
    	line_count = 0;
	int ch;
    	if(fp == NULL)
	{
        	fclose(fp);
		printf("FAILED TO OPEN FILE: %s",filename);
		return;
    	}
	else
	{	
		do
    		{
  		 	ch = fgetc(fp);
   			if( ch == '\n')
			{
				 line_count++;  
			}
 
   		}while( ch != EOF );    		

    		fclose(fp);
	}
	/****************************************************************/
	
	fp = fopen(filename,"r");	//Open file for reading
	nNodes = line_count;
	nDims = 3;	//StarCD file has 3 dims for vertices, but this function only reads 2

	// Echo data read from .node file
	printf("Node File Data\n");
	printf("# of Nodes: %d \n# of Dimensions: %d \n\n",nNodes,nDims) ;

	Quadmesh.nNodes = nNodes ; // Store # of nodes into Trimesh variable

	// Allocate for nodes array(1D array)
	Quadmesh.nodes = (vertex *) malloc(nNodes*sizeof(vertex)) ;
	
	//Allocate for numNodes array(1D array)
	Quadmesh.nodeNums = malloc(nNodes*sizeof(int));
        
	// Loop through file and read all node data
	for(i=0; i<nNodes; i++)
	{
	
		fscanf(fp,"%d %lf %lf %lf",&dummy,&x,&y,&dummy) ;
               	//printf("%lf %lf \n",x,y) ; //Echo data from file
               	Quadmesh.nodes[i].x = x ;
               	Quadmesh.nodes[i].y = y ;
		Quadmesh.nodeNums[i] = i+1;

	}

	fclose(fp) ;

	/* Read elements *************************************************/
  	sprintf(filename,"%s.cel",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename to be read

      	fp = fopen(filename,"r") ; // Open .ele file for reading elements

	 /* Line Counting For Files ***************************************/
        line_count = 0;

        if(fp == NULL)
        {
                fclose(fp);
                printf("FAILED TO OPEN FILE: %s",filename);
                return;
        }
        else
        {
		do
                {
                        ch = fgetc(fp);
                        if( ch == '\n')
                        {
                                 line_count++;
                        }

                }while( ch != EOF );

                fclose(fp);
        }
        /****************************************************************/

        fp = fopen(filename,"r");       //Open file for reading

	nElems = line_count;
	nNpT = 4; //Number of nodes per element

	printf("Element File Data\n");
	printf("# of Elements: %d\n# of Verts Per Element:  %d \n\n",nElems,nNpT) ; //Echo data back to user

	Quadmesh.nElems = nElems ; // Store number of elements into Quadmesh variable

	// Allocate for elements array ( 2D array)
	Quadmesh.elems = (int **) malloc(nElems*sizeof(int *)) ; // Allocate rows of array
        
        //Loop through and store element vertices
	for(i=0; i<nElems; i++)
	{
	  fscanf(fp,"%d %d %d %d %d %d %d %d %d %d %d",&dummy,&v1,&v2,&v3,&v4,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy) ;
	  //printf("%d %d %d %d %d \n",i,v1,v2,v3,v4) ;
	  Quadmesh.elems[i] = (int *) malloc(4*sizeof(int)) ; // Allocate columns of array
          
	  //Check this section if grid ouput is not correct(adjust clockwise/counter-clockwise read)
	  Quadmesh.elems[i][0] = v1 ; 
          Quadmesh.elems[i][1] = v2 ;
          Quadmesh.elems[i][2] = v3 ;
	  Quadmesh.elems[i][3] = v4;
	}

	fclose(fp) ;

	/* Read Edges ****************************************************/
  	sprintf(filename,"%s.bnd",casename) ;
	printf("Reading File: %s \n\n",filename) ; // Echo filename being read

      	fp = fopen(filename,"r") ;

	 /* Line Counting For Files ***************************************/
        line_count = 0;

        if(fp == NULL)
        {
                fclose(fp);
                printf("FAILED TO OPEN FILE: %s",filename);
                return;
        }
        else
        {
		do
                {
                        ch = fgetc(fp);
                        if( ch == '\n')
                        {
                                 line_count++;
                        }

                }while( ch != EOF );

                fclose(fp);
        }
        /****************************************************************/

        fp = fopen(filename,"r");       //Open file for reading

	nEdges = line_count;
	MarkerFlag = 1;	// the StarCD file should always have markers

	printf("Edge File Data\n");
	printf("# of Edges: %d\nBoundary Marker Flag: %d \n",nEdges,MarkerFlag) ; // Echo first line back to screen

	Quadmesh.nEdges = nEdges ; // Store total number of edges in mesh
	

	// Allocate for edges array( 2D array)
	Quadmesh.edges = (int **) malloc(nEdges*sizeof(int *)) ; // Allocate rows of edge array
	Quadmesh.edgeMarkers = (int *) malloc(nEdges*sizeof(int)) ;// Allocate 1D array for all edge markers
	
	Quadmesh.nEdgeMarkers = 1;	// Initialize unique edge boundary marker counter
	Quadmesh.edgeBoundaryMarkers = malloc(1*sizeof(int)); // Allocate first element of unique edge boundary marker array

	for(i=0; i<nEdges; i++)
	{

	  fscanf(fp,"%d %d %d %d %d %d %d %s",&dummy,&v1,&v2,&dummy,&dummy,&marker,&dummy,&dummychar) ;
	  //printf("%d %d %d %s\n",i+1,v1,v2,dummychar) ;
	  Quadmesh.edges[i] = (int *) malloc(2*sizeof(int)) ; // Allocate columns for row i

          //Order of read matters for Clockwise/Anti-Clockwise ordering
          Quadmesh.edges[i][0] = v2 ;
          Quadmesh.edges[i][1] = v1 ;

          Quadmesh.edgeMarkers[i] = marker ;
	  
	  if(i == 0) // set currentMarker to the first marker that is read from file
	  {
		Quadmesh.edgeBoundaryMarkers[i] = marker;
       	  }


	  int newMarker = 1;	//Boolean for tracking if new marker is found. 0=false, 1=true
	  for(j=0;j<Quadmesh.nEdgeMarkers;j++) // Compare currently read marker with markers that have been stored
	  {
		if(marker == Quadmesh.edgeBoundaryMarkers[j])
		{
			newMarker = 0;
		}

	  }

	  if( newMarker == 1  ) // A new unique boundary marker is found. 
	  {
		Quadmesh.nEdgeMarkers=Quadmesh.nEdgeMarkers+1;

 		// Dynamically resize edgeBoundaryMarkers array to hold newly found boundary marker
		int *temp = realloc(Quadmesh.edgeBoundaryMarkers,Quadmesh.nEdgeMarkers*sizeof(int));  // Store data in a larger temp array
   		if(temp != NULL) // realloc was successful
		{
			Quadmesh.edgeBoundaryMarkers = temp; //Point original array to newly sized array
			Quadmesh.edgeBoundaryMarkers[Quadmesh.nEdgeMarkers-1] = marker; // Store newly read marker into array
		}
		else
		{
			free(Quadmesh.edgeBoundaryMarkers);
			printf("Error allocating memory for Edge Boundary Marker Array!\n");
			return ;
		}

	  }
		

}

	// Display the number of unique Edge markers in .edge file
	printf("\nNumber of Unique Boundary Edge Markers: %d \n",Quadmesh.nEdgeMarkers);

	// Allocate edgeCounts array
	Quadmesh.edgeCounts = malloc(Quadmesh.nEdgeMarkers*sizeof(int));

	// Initialize edgeCounts array to zero. This array will hold the number
	// of edges that have a specific boundary marker.
	for(i=0; i<Quadmesh.nEdgeMarkers; i++)
	{
		Quadmesh.edgeCounts[i] = 0;
	}
	


	//Loop through edgeMarkers array and count the number of edges that have all of the markers
	//that are stored in edgeBoundaryMarkers
	for(i=0;i<nEdges; i++) // Loop through every edge
	{
		for(j=0;j<Quadmesh.nEdgeMarkers; j++) //Loop through all edge boundary markers
		{

			if(Quadmesh.edgeMarkers[i] == Quadmesh.edgeBoundaryMarkers[j] ) // see which edgeBoundaryMarker the edgeMarker matches
			{				
				Quadmesh.edgeCounts[j] = Quadmesh.edgeCounts[j] + 1; // Increment corresponding edgeCount
			}			

		}

	}


	fclose(fp) ;


	// Print the number of edges on each edge boundary that were counted in file
	printf("\nNumber of Edges On Each Boundary\n");	
	int sum = 0;
	for( i=0;i<Quadmesh.nEdgeMarkers;i++)
	{
	  	printf("Boundary Marker %d(%d) = %d \n",Quadmesh.edgeBoundaryMarkers[i],i+1,Quadmesh.edgeCounts[i]);
		sum = sum + Quadmesh.edgeCounts[i];  
	}

	//Error Check, Print face counts to screen
	printf("\nTotal Edges Counted:%d\n# Edges In File:%d\n\n",sum,Quadmesh.nEdges);

	printf("Finished Reading Quadrilateral Mesh Data File\n\n");

	//DEBUGGING PURPOSES
	if(debug == 1)//print entire contents of quadmesh data structure
	{
	
		printf("\n\nPrinting Contents of quadmesh Data Structure");
		printf("\nNumber of Nodes: %d",Quadmesh.nNodes);
		printf("\nNumber of Elements: %d",Quadmesh.nElems);
		printf("\nNumber of Edges: %d",Quadmesh.nEdges);
		printf("\nNumber of Unique Edge Boundary Markers: %d",Quadmesh.nEdgeMarkers);


		printf("\n\n------------Printing Node Data------------");
		printf("\nNode # \t X Coord \t Y Coord");		
		for(i=0;i<Quadmesh.nNodes;i++)
		{
			printf("\n %d %4.2e %4.2e",Quadmesh.nodeNums[i],Quadmesh.nodes[i].x,Quadmesh.nodes[i].y);

		}


	}
}







int stitch_meshes()
{

   int i, j;
   int debug = 0; // set to 1 for verbose output, 2 for less output
   int mCount;	//Number of matching nodes between 2 meshes
	
   //Combine node lists. The first nodes are the quad mesh, then the tri mesh. There will
   //be overlap because of the shared boundary nodes


   printf("\n\n----------STITCHING MESHES TOGETHER----------\n\n");
   printf("The Code stitches the quad and tri  meshes together. \n\n") ;

   /*-----------Shift all Triangle Nodes to Start Numbering at end of Quadrilateral Node List--------*/
   printf("Shifting trimesh node data\n");
  //First shift all of the node data for the triangle mesh 
   for(i=0;i<Trimesh.nNodes;i++)
   {
  	Trimesh.nodeNums[i] = Trimesh.nodeNums[i] + Quadmesh.nNodes;
	if(debug ==1)
	{
		printf("\n %d %lf %lf",Trimesh.nodeNums[i],Trimesh.nodes[i].x,Trimesh.nodes[i].y);
	}
   }

   /*-------------Adjust all Triangle Mesh Cell Connectivies to reflect Node Numbering Translation-----*/

	for(i = 0;i<Trimesh.nElems;i++)
	{

	  Trimesh.elems[i][0] = Trimesh.elems[i][0] + Quadmesh.nNodes;
	  Trimesh.elems[i][1] = Trimesh.elems[i][1] + Quadmesh.nNodes;
   	  Trimesh.elems[i][2] = Trimesh.elems[i][2] + Quadmesh.nNodes;

	}




   /*---------Find Common Nodes Between the 2 Meshes and Store their Relation-----------------*/
   printf("\n\nMerging common nodes closer than %5.2E between the two meshes\n",EPSILON);
   mCount = 0;	//Initialize match count to 0

   //Find all matching nodes between the 2 meshes
   for(i=0; i<Trimesh.nNodes ;i++) // Loop over all triangle nodes
   {
	if(debug==2)
	{
	  printf("\n\nComparing Trimesh Node: %d to All Quadmesh Nodes",i+1);
	}

	//Compare each node in the trimesh array to the nodes in the quadmesh array
	for(j=0; j<Quadmesh.nNodes; j++) // Loop over all quadrilateral nodes
	{
	  if(debug==1)
	  {
	    printf("\n\nQuadmesh Node %d:  X[%4.2e]   Y[%4.2e]",j+1,Quadmesh.nodes[j].x,Quadmesh.nodes[j].y);
	    printf("\nTrimesh Node %d:  X[%4.2e]  Y[%4.2e]",i+1,Trimesh.nodes[i].x,Trimesh.nodes[i].y);
			
	    printf("\nDeltaX: %4.2e  DeltaY: %4.2e ",fabs(Quadmesh.nodes[j].x-Trimesh.nodes[i].x) ,fabs(Quadmesh.nodes[j].y-Trimesh.nodes[i].y) );
	  }
	
 	  if( fabs(Quadmesh.nodes[j].x - Trimesh.nodes[i].x) <=EPSILON && fabs( Quadmesh.nodes[j].y - Trimesh.nodes[i].y) <= EPSILON )
	  {
	    
	    if(debug==2)
	    {
	      printf("\n\nA Matching Node Pair Found Between Quad:%d  and Tri: %d",Quadmesh.nodeNums[j],Trimesh.nodeNums[i]);
	    }
	
	    //A match is found.
	    if(mCount == 0)//First match
	    {
	      mCount = mCount+1;
	      matchArray = malloc(mCount*sizeof(int *)) ; // Allocate rows of match array
	      matchArray[mCount-1] =  malloc(2*sizeof(int)) ; // Allocate columns of match array

              matchArray[mCount-1][0] = Quadmesh.nodeNums[j];
              matchArray[mCount-1][1] = Trimesh.nodeNums[i];
              
	      if(debug==2)
	      {
		printf("\nMatch #:%d  QuadNode:%d  TriNode:%d",mCount,matchArray[mCount-1][0],matchArray[mCount-1][1]);
	      }
	      
	      if(debug==1)
  	      {
	        printf("\nQuadmesh: XCOORD[%4.2e]    YCOORD[%4.2e]",Quadmesh.nodes[j].x,Quadmesh.nodes[j].y);
	        printf("\nTrimesh: XCOORD[%4.2e]    YCOORD[%4.2e]",Trimesh.nodes[i].x,Trimesh.nodes[i].y);
	      }

            }
	    else
	    {

  	      mCount = mCount+1;
	
	      //Resize matchArray to hold additional node
	      int **temp = (int **)realloc(matchArray,mCount*sizeof(int *));  // Store data in a larger temp array
	      if(temp != NULL) // realloc was successful
              {
  
                matchArray = temp; //Point original array to newly sized array

		matchArray[mCount-1] = malloc(2*sizeof(int));
		matchArray[mCount-1][0] = Quadmesh.nodeNums[j];
    	    	matchArray[mCount-1][1] = Trimesh.nodeNums[i];
	    	if(debug==2)
		{
		  printf("\nMatch #:%d  QuadNode:%d  TriNode:%d\n",mCount,matchArray[mCount-1][0],matchArray[mCount-1][1]);
		}
 		if(debug==1)
                {
                  printf("\nQuadmesh: XCOORD[%4.2e]    YCOORD[%4.2e]",Quadmesh.nodes[j].x,Quadmesh.nodes[j].y);
                  printf("\nTrimesh: XCOORD[%4.2e]    YCOORD[%4.2e]",Trimesh.nodes[i].x,Trimesh.nodes[i].y);
                }

              }
              else
              {
                free(matchArray);
                printf("Error allocating memory for Edge Boundary Marker Array!\n");
                return ;
              }		

	    }
	  } 

	}


   }


  //At this point the matching nodes have been located.
  //Print matching node data
  if(debug==2)
  {

	printf("\n\n--------Printing Node Matches Between The 2 Meshes--------");
	for(i=0;i<mCount;i++)
	{
	  printf("\nMatch # %d: Quad Node: %d   Tri Node: %d",i+1,matchArray[i][0],matchArray[i][1]);

	}
  }

  //Now that the matches have been located and stored. The code must go through the trimesh and update its elements
  //to reference the nodes from the quadmesh node list.

  return mCount;	//return # of node matches to main program

}





void output_grd_format(int mCount)
{
	int	i, j, k ;
	int	boundaryType=11;// Arbitrary boundary type(keep as is)
	int	count,sum ;
	int 	outerBndryMrkr;	//Outer boundary marker
	int 	repeated = 0; // 0 --> False, 1--> True
        int 	loc = 0;    //row index of repeated node(used only if node is repeated)
	char	casename[20], filename[20];
	FILE	*fp ;

//*************************************************************
// 	Open file and read title
//*************************************************************

	printf("\n\nWriting hybrid GRD grid file...\n") ;
	sprintf(casename,"cylds") ;

	sprintf(filename,"%s.grd",casename) ;
	printf("Writing hybrid GRD format to: \"%s\" \n",filename) ;
	fp = fopen(filename,"w") ;

	//Outer Boundary Marker  provided by user
	printf("What is boundary marker for the outer boundary?\nChoices are:");
	for(i = 0;i<Quadmesh.nEdgeMarkers;i++)
	{
		printf("\nMarker #: %d",Quadmesh.edgeBoundaryMarkers[i]);

	}
	
	printf("\n\nSelection: ");
	scanf("%d",&outerBndryMrkr);
	
	//Echo user input back to user
	printf("\nBoundary Selection: %d ",outerBndryMrkr);

//===================================================================
//	Print Node Coordinates
//===================================================================

	fprintf(fp,"%8d %8d %8d",Quadmesh.nElems+Trimesh.nElems,Quadmesh.nNodes+Trimesh.nNodes-mCount,1) ;

	printf("\nPrinting Node Coordinates...") ;
	//Print quad mesh nodes
	for(i=0;i<Quadmesh.nNodes;i++)
	{
		fprintf(fp,"\n%23.16e \t %23.16e",Quadmesh.nodes[i].x,Quadmesh.nodes[i].y);

	}
	

	//Print Triangle Nodes(only non-repeated ones though)
	for(i=0;i<Trimesh.nNodes;i++)
	{
	
		//Check if each node about to be printed is on the repeated node list
		repeated = 0; // 0 --> False, 1--> True
               	loc = 0;    //row index of repeated node(used only if node is repeated)
		for(j=0;j<mCount;j++)
		{
			if(Trimesh.nodeNums[i] == matchArray[j][1])
			{
				repeated =1;	//Node is repeated from quad mesh
			}

		}

		if(repeated == 0)
		{
			fprintf(fp,"\n%23.16e \t %23.16e",Trimesh.nodes[i].x,Trimesh.nodes[i].y);

		}

	}
	

//=====================================================================
//	Cell connectivity
//=====================================================================


	/*----------First Print All Quad Mesh Cell connectivities---------------*/
	for(i=0;i<Quadmesh.nElems;i++)
	{
		fprintf(fp,"\n%8d %8d %8d %8d",Quadmesh.elems[i][0],Quadmesh.elems[i][1],Quadmesh.elems[i][2],Quadmesh.elems[i][3]);

	}

	/*-------Now Print Triangle Mesh Cell Connectivities---------------------*/
	fprintf(fp,"\n");
	for(i=0;i<Trimesh.nElems;i++)
	{

		for(j=0;j<3;j++)//Loop through all 3 nodes per element
		{

			//If a node on the list of repeated nodes is to be printed, print the
			//corresponding quad mesh's node in place
			repeated = 0; // 0 --> False, 1--> True
			loc = 0;	//row index of repeated node(used only if node is repeated)
			for( k=0;k<mCount;k++)//Loop through repeated node list
			{
				if(Trimesh.elems[i][j] == matchArray[k][1])
				{
					repeated = 1;
					loc = k;
				}

			}

			if(repeated == 1)// node is repeated therefore use quad mesh's node instead
			{
	
				fprintf(fp,"%8d ",matchArray[loc][0]);
	
			}
			else
			{

				//Just print trimesh element
				//Note: For triangles, print 4 nodes for connectivity, just make the last one equal to the first
				fprintf(fp,"%8d ",Trimesh.elems[i][j]);

			}

		}
		
		fprintf(fp,"%8d",Trimesh.elems[i][0]);
		fprintf(fp,"\n");


	}


//==============================================================================
//	Boundary Edge Connectivities & Type
//==============================================================================
	
	for(i=0;i<Quadmesh.nEdgeMarkers;i++)
	{
		if(Quadmesh.edgeBoundaryMarkers[i] == outerBndryMrkr)
		{
		  //printf("\nEdge Marker: %d    Selected Marker: %d",Quadmesh.edgeBoundaryMarkers[i],outerBndryMrkr);
		  fprintf(fp,"%8d %8d",boundaryType,Quadmesh.edgeCounts[i]);

		}
	}

	for(i=0;i<Quadmesh.nEdges;i++)
	{
		if(Quadmesh.edgeMarkers[i] == outerBndryMrkr)
		{
			fprintf(fp,"\n%8d  %8d",Quadmesh.edges[i][0],Quadmesh.edges[i][1]);
		}		
	
	}


//**************************************************************************
// 	Close file
//**************************************************************************

	fclose(fp) ;

//**************************************************************************
//      Print trimesh object contents
//**************************************************************************

  	printf("Writing hybrid GRD grid file done.\n") ;

  	return ;
}

