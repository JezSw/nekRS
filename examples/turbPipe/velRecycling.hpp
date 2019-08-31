/*
   copy velocity data of a given slab (slabIdSrc) to another slab
   (slideIdDst) also known as recycling

   Note: This implementation relies on a special global element 
         numbering which is only true for extruded meshes in z from nek!
*/

#include <nekInterfaceAdapter.hpp>

namespace velRecycling {

static ogs_t *ogs;
static ins_t *ins;

static occa::memory o_wrk;

static dfloat *flux, *area;
static occa::memory o_flux, o_area;

static dfloat *tmp1, *tmp2;
static occa::memory o_tmp1, o_tmp2;

static occa::kernel setValueBCKernel;
static occa::kernel getBCFluxKernel;
static occa::kernel sumReductionKernel; 
static occa::kernel scalarMultiplyKernel; 

bool buildKernelCalled = 0;
bool setupCalled = 0;

void buildKernel(ins_t *ins)
{
  mesh_t *mesh = ins->mesh; 

  string fileName;
  int rank = mesh->rank;
  fileName.assign(getenv("NEKRS_INSTALL_DIR"));
  fileName += "/okl/recycling.okl";
  occa::properties& kernelInfo = *ins->kernelInfo;
  for (int r=0;r<2;r++){
    if ((r==0 && rank==0) || (r==1 && rank>0)) {
       setValueBCKernel     =  mesh->device.buildKernel(fileName.c_str(), "setBCVectorValue", kernelInfo);
       getBCFluxKernel      =  mesh->device.buildKernel(fileName.c_str(), "getBCFlux", kernelInfo);
       sumReductionKernel   =  mesh->device.buildKernel(fileName.c_str(), "sumReduction", kernelInfo);
       scalarMultiplyKernel =  mesh->device.buildKernel(fileName.c_str(), "scalarMultiply", kernelInfo);
    }
    MPI_Barrier(mesh->comm);
  }
} 

void copy()
{
  mesh_t *mesh = ins->mesh; 
  const int bc = 2;
  const dfloat wbar = 1.0; 

  // copy recycling plane in interior to inlet
  o_wrk.copyFrom(ins->o_U, ins->NVfields*ins->Ntotal*sizeof(dfloat));
  setValueBCKernel(mesh->Nelements, 0.0, bc, ins->fieldOffset,
                   o_wrk, mesh->o_vmapM, mesh->o_EToB);

  //ogsGatherScatterMany(o_wrk, ins->NVfields, ins->fieldOffset,
  //                     ogsDfloat, ogsAdd, ogs);
  for(int k=0;k<ins->dim;++k)
    ogsGatherScatter(o_wrk+k*ins->fieldOffset*sizeof(dfloat),
                     ogsDfloat, ogsAdd, ogs);
  
  // rescale
  getBCFluxKernel(mesh->Nelements, bc, ins->fieldOffset, o_wrk,
                  mesh->o_vmapM, mesh->o_EToB, mesh->o_sgeo, o_area, o_flux);
 
  const int NfpTotal = mesh->Nelements*mesh->Nfaces*mesh->Nfp; 
  sumReductionKernel(NfpTotal, o_area, o_flux, o_tmp1, o_tmp2); 
  o_tmp1.copyTo(tmp1);
  o_tmp2.copyTo(tmp2);
  dfloat sbuf[2] = {0,0};
  for(int n=0; n<ins->Nblock; n++){
    sbuf[0] += tmp1[n]; 
    sbuf[1] += tmp2[n]; 
  }
  MPI_Allreduce(MPI_IN_PLACE, sbuf, 2, MPI_DFLOAT, MPI_SUM, mesh->comm); 

  const dfloat scale = -wbar*sbuf[0] / sbuf[1]; 
  scalarMultiplyKernel(ins->NVfields*ins->Ntotal, scale, o_wrk);
}


void setup(ins_t *ins_, occa::memory o_wrk_, const int nelSlab, const int slabIdSrc)
{
  ins = ins_;
  o_wrk = o_wrk_; 
  mesh_t *mesh = ins->mesh;

  const int slabIdDst = 0;
  const int offset = (slabIdSrc - slabIdDst) * nelSlab;
  const int egStartSrc = slabIdSrc * nelSlab;   
  const int egStartDst = slabIdDst * nelSlab;   
  const int Ntotal = mesh->Np * mesh->Nelements;

  hlong *ids = (hlong *) calloc(Ntotal, sizeof(hlong));

  for (int e=0; e < mesh->Nelements; e++){
    for (int n=0; n < mesh->Np; n++) 
      ids[e*mesh->Np + n] = 0; // gs noop

    const int eg = nek_lglel(e); // 0-based

    if (eg >= egStartSrc && eg < egStartSrc + nelSlab)
      for (int n=0; n < mesh->Np; n++) {
        ids[e*mesh->Np + n] = (eg-offset)*mesh->Np+n+1;
      }
    if (eg >= egStartDst && eg < egStartDst + nelSlab)
      for (int n=0; n < mesh->Np; n++)  {
        ids[e*mesh->Np + n] = (eg*mesh->Np+n+1); 
      }
  }

  ogs = ogsSetup(Ntotal, ids, mesh->comm, 1, mesh->device);
  free(ids);

  tmp1   = (dfloat *) calloc(ins->Nblock, sizeof(dfloat));
  o_tmp1 = mesh->device.malloc(ins->Nblock*sizeof(dfloat), tmp1);
  tmp2   = (dfloat *) calloc(ins->Nblock, sizeof(dfloat));
  o_tmp2 = mesh->device.malloc(ins->Nblock*sizeof(dfloat), tmp2);

  const int NfpTotal = mesh->Nelements*mesh->Nfaces*mesh->Nfp;
  flux   = (dfloat *)calloc(NfpTotal, sizeof(dfloat));
  o_flux = mesh->device.malloc(NfpTotal*sizeof(dfloat), flux);
  area   = (dfloat *)calloc(NfpTotal, sizeof(dfloat));
  o_area = mesh->device.malloc(NfpTotal*sizeof(dfloat), area);
}


} // namespace 