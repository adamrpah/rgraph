#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>
#include "tools.h"
#include "graph.h"
#include "modules.h"
#include "bipartite.h"

main(int argc, char **argv)
{
  FILE *outF, *inF;
  int seed = 1111;
  struct binet *binet = NULL;
  struct group *part = NULL;
  gsl_rng *randGen;
  double Ti, Tf, Ts, fac;
  char file_name[100];
  int invert;

  /*
    ------------------------------------------------------------
    Prompt for user-defined parameters
    ------------------------------------------------------------
  */
  if (argc < 6) {
    printf("\nUse: bipartmod_cl seed net_file_name iteration_factor cooling_factor column");
    return -1; 
  }

  seed = atoi(argv[1]);
  file_name = argv[2];
  fac = atof(argv[3]);
  Ts = atof(argv[4]);
  invert = atoi(argv[5]);

  /*
    ------------------------------------------------------------------
    Initialize the random number generator
    ------------------------------------------------------------------
  */
  randGen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(randGen, seed);

  /*
    ------------------------------------------------------------------
    Build the network
    ------------------------------------------------------------------
  */
  inF = fopen(file_name, "r");
  binet = FBuildNetworkBipart(inF, 0, 0);
  fclose(inF);
  if (invert == 1)
    InvertBipart(binet);

  /*
    ------------------------------------------------------------------
    Find the modules using the bipartite network
    ------------------------------------------------------------------
  */
  Ti = 1. / (double)CountNodes(binet->net1);
  Tf = 0.;

  part = SACommunityIdentBipart(binet,
				Ti, Tf, Ts, fac,
				0, 'o', 1, 'm',
				randGen);
  outF = fopen("modules_bipart.dat", "w");
  FPrintPartition(outF, part, 0);
  fclose(outF);
  RemovePartition(part);

  // Free memory
  // ------------------------------------------------------------
  gsl_rng_free(randGen);
  RemoveBipart(binet);
}
