/*
 * Small C shim for NVTX ranges used by the Fortran timing module.
 *
 * NVTX v3 is header-only and dynamically connects to the profiler when one
 * is attached. If the CUDA toolkit header is unavailable, keep these entry
 * points as no-ops so a non-NVTX build remains linkable.
 */

#if defined(__has_include)
#  if __has_include(<nvtx3/nvToolsExt.h>)
#    include <nvtx3/nvToolsExt.h>
#    define SPECFEM_HAVE_NVTX 1
#  elif __has_include("/usr/local/cuda/include/nvtx3/nvToolsExt.h")
#    include "/usr/local/cuda/include/nvtx3/nvToolsExt.h"
#    define SPECFEM_HAVE_NVTX 1
#  endif
#endif

#ifndef SPECFEM_HAVE_NVTX
#  define SPECFEM_HAVE_NVTX 0
#endif

void specfem_nvtx_acoustic_push(void)
{
#if SPECFEM_HAVE_NVTX
  (void)nvtxRangePushA("SPECFEM::acoustic");
#endif
}

void specfem_nvtx_acoustic_pop(void)
{
#if SPECFEM_HAVE_NVTX
  (void)nvtxRangePop();
#endif
}

void specfem_nvtx_viscoelastic_push(void)
{
#if SPECFEM_HAVE_NVTX
  (void)nvtxRangePushA("SPECFEM::viscoelastic");
#endif
}

void specfem_nvtx_viscoelastic_pop(void)
{
#if SPECFEM_HAVE_NVTX
  (void)nvtxRangePop();
#endif
}
