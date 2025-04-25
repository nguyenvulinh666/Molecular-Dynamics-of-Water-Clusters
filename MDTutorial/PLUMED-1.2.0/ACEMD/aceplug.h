#ifndef __aceplug_PLUGIN
#define __aceplug_PLUGIN 1

// We need the CUDA vector types
#include "vector_types.h"

#define aceplug_PLUGIN_VER 1

  /* Time unit in femtoseconds */
#define   TIMEFACTOR   48.88821 
  /* Coulomb's constant for electrostatics, units kcal*A/mol/e^2 */
#define MD_COULOMB  332.0636
  /* Boltzman's constant for temperature, units kcal/A/K */
#define MD_BOLTZMAN  0.001987191

#define MD_KCAL_TO_KJ  4.1868


struct aceplug_sim_t {
	float4 *pos;
	float4 *frc;
	float1  *mass;
	float1  *charge;
	int    natoms;
	int    step;
	void *privdata;
} aceplug_sim_t;


#ifdef __cplusplus
extern "C" {
#endif

int aceplug_init( struct aceplug_sim_t*, int argc, char **argkey, char **argval );

int aceplug_update( struct aceplug_sim_t* );

int aceplug_terminate( void * privdata );

#ifdef __cplusplus
}
#endif



#endif
