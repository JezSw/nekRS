/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "ogstypes.h"

#include "gslib.h"

void *ogsHostSetup(MPI_Comm meshComm,
                   dlong NuniqueBases,
                   hlong *gatherGlobalNodes,
                   int nonsymm,
                   int verbose){

  /* gslib stuff */
  comm_ext world;
  struct comm com;
  
  /*  MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*) &world); */
  world = (comm_ext)meshComm; // MPI_COMM_WORLD;
  
  comm_init(&com, world);

  /* for the moment borrow gslib array */
  slong *id = tmalloc(slong, NuniqueBases);
  
  dlong n;
  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    id[n] = (slong) gatherGlobalNodes[n];
  }

  struct gs_data *gsh = gs_setup(id, NuniqueBases, &com, nonsymm, gs_pairwise, 0); // gs_auto, gs_crystal_router, gs_pw

  free(id);

  return gsh;
}


void ogsHostFree(void *gsh){

  gs_free(static_cast<gs_data*>(gsh));

}


void ogsGsUnique(hlong *gatherGlobalNodes,
                 dlong NuniqueBases,
                 MPI_Comm meshComm){

  /* gslib stuff */
  comm_ext world;
  struct comm com;
  
  /*  MPI_Comm_dup(MPI_COMM_WORLD, (MPI_Comm*) &world); */
  world = (comm_ext)meshComm; // MPI_COMM_WORLD;
  
  comm_init(&com, world);

  /* for the moment borrow gslib array */
  slong *id = tmalloc(slong, NuniqueBases);
  
  dlong n;
  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    id[n] = (slong) gatherGlobalNodes[n];
  }

  gs_unique(id, NuniqueBases, &com);

  for(n=0;n<NuniqueBases;++n){ /* at some point need to choose int */
    gatherGlobalNodes[n] = (hlong) id[n];
  }

  free(id);
}
