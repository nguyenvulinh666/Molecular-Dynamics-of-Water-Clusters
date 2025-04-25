#include <stdio.h>
#include <stdlib.h>
#include "aceplug.h"
#include "metadyn.h"
//IDX can be any private data of the plugin
#define IDX 1

int aceplug_init( 
	struct aceplug_sim_t *s,
	int argc,
	char **argkey,
	char **argval)
{   
  int i;
  int nrepl,repl;
  real timeunit, timestep, dt;
  real rte0,rteio;
  real *mass;
  real *charge;
  char *metainp;
  real box[3];

  charge    = (real *)calloc(s->natoms,sizeof(real));
  mass      = (real *)calloc(s->natoms,sizeof(real));
       
  if (argc!=4) { printf("Wrong argument number. \n You must provide the name of inputfile and boxsize %i \n", argc ); }
  else if(argc==4) {  
	printf("I have %d arguments: \n", argc );
	for(i=0; i < argc; i++ ) {
		printf("key %s value %s\n", argkey[i], argval[i]);
		metainp=argval[0];
		box[0]=strtod(argval[1],NULL);
		box[1]=strtod(argval[2],NULL);
		box[2]=strtod(argval[3],NULL);
	}
	dt = 1.0; 
	nrepl = 1;
	repl=1;

 //     Target_temperature controllare che non sia in kcal/mol (*0.59227/298)
	rte0  = 298.0;
	rteio = rte0; 

	for(i=0;i<s->natoms;i++) {
	  mass[i]   = (double) s->mass[i].x;
	  charge[i] = (double) s->charge[i].x;
	}

        init_metadyn(s->natoms, charge, mass, 
		     dt, repl, nrepl, rte0, rteio, metainp, box );
	  }
 

	return 0;
}

int aceplug_update( struct aceplug_sim_t *s) {
    meta_force_calculation(s);
	return 0;
}

int aceplug_terminate( void *privdata ) {
	printf("Hello! I am plugin-example %d\n", IDX );
	return 0;
}
