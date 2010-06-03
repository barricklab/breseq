#define PARI_TUNE

#ifdef PARI_TUNE
extern long KARATSUBA_SQRI_LIMIT;
extern long KARATSUBA_MULI_LIMIT;
extern long KARATSUBA_MULR_LIMIT;
extern long MONTGOMERY_LIMIT;
extern long REMIIMUL_LIMIT;
extern long INVMOD_GMP_LIMIT;
extern long DIVRR_GMP_LIMIT;
extern long Flx_INVMONTGOMERY_LIMIT;
extern long Flx_POW_MONTGOMERY_LIMIT;
extern long EXPNEWTON_LIMIT;
extern long LOGAGM_LIMIT;
extern long LOGAGMCX_LIMIT;
extern long AGM_ATAN_LIMIT;
extern long Flx_SQR_LIMIT;
extern long Flx_MUL_LIMIT;
extern long RgX_SQR_LIMIT;
extern long RgX_MUL_LIMIT;
#else
#  define KARATSUBA_SQRI_LIMIT     __KARATSUBA_SQRI_LIMIT
#  define KARATSUBA_MULI_LIMIT     __KARATSUBA_MULI_LIMIT
#  define KARATSUBA_MULR_LIMIT     __KARATSUBA_MULR_LIMIT
#  define MONTGOMERY_LIMIT         __MONTGOMERY_LIMIT
#  define REMIIMUL_LIMIT           __REMIIMUL_LIMIT
#  define INVMOD_GMP_LIMIT         __INVMOD_GMP_LIMIT
#  define DIVRR_GMP_LIMIT          __DIVRR_GMP_LIMIT
#  define EXPNEWTON_LIMIT          __EXPNEWTON_LIMIT
#  define LOGAGM_LIMIT             __LOGAGM_LIMIT
#  define LOGAGMCX_LIMIT           __LOGAGMCX_LIMIT
#  define AGM_ATAN_LIMIT           __AGM_ATAN_LIMIT
#  define Flx_INVMONTGOMERY_LIMIT  __Flx_INVMONTGOMERY_LIMIT
#  define Flx_POW_MONTGOMERY_LIMIT __Flx_POW_MONTGOMERY_LIMIT
#  define Flx_SQR_LIMIT            __Flx_SQR_LIMIT
#  define Flx_MUL_LIMIT            __Flx_MUL_LIMIT
#  define RgX_SQR_LIMIT            __RgX_SQR_LIMIT
#  define RgX_MUL_LIMIT            __RgX_MUL_LIMIT
#endif
