#ifndef __TIO_PPC_METHODS_H__
#define __TIO_PPC_METHODS_H__ 1

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

enum
{
   _pTIO_PPC_METHOD_ALT  = 0,
   _pTIO_PPC_METHOD_ZERO = 1,
   _pTIO_PPC_METHOD_SET  = 2
};

extern int _pTIO_ppc_f32_bitmask (int nsd, int ppc_method, int sz, float *val,
                                  const float *pmiss_val);

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
