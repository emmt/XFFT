/*
 * xfft.c --
 *
 * Implement eXtended FFT for Yorick, using Swarztrauber FFT routines
 * compiled with Yorick, or FFTW (version 2 or 3).
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2009-2015: Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *-----------------------------------------------------------------------------
 */

/* RATIONALE AND IMPLEMENTATION NOTES

   FFTW assumes row-major order while Yorick is column-major order.  Hence the
   dimension lists are reversed when creating FFTW plans.  For real-to-complex
   or complex-to-real transforms, this means that the "half" dimension of
   complex Hermitian arrays corresponds to the first Yorick dimension (the
   fastest varying one).

   For real-to-complex or complex-to-real transforms, in terms of number of
   real values (a complex counting as 2 reals), the first dimension of the
   complex array is 2*(N1/2 + 1) where N1 is the first dimension of the real
   array.  Whatever is N1, 2*(N1/2 + 1) is always larger than N1, hence
   input and output arrays must be different.

   In the benchmarks, out-of-place transforms are slightly faster than
   in-place ones.  It is expected that, most of the time, the user of the
   FFTW-Yorick interface wants to perform out-of-place transforms (preserving
   the input array).  However, (1) using aligned memory imposes to allocate
   internals workspaces and to copy arrays from/to this workspaces; (2) not
   all out-of-place transforms (notably the complex-to-real ones) are able to
   preserve the input array.  Hence, to minimize the number of copies (which
   have a time cost for large arrays) and achieve quasi-optimal performances
   for most cases.  All plans have their private workspace of size suitable to
   store any real/complex array with dimensions corresponding to the plan.  If
   memory alignement is required (option ALIGN=1 in fftw_new), the workspace
   is allocated with fftw_malloc and destroyed with fftw_free; otherwise,
   p_malloc and p_free are used.  If memory alignement is required, in-place
   transforms are always used, on entry the input array is copied (and perhaps
   converted) into the workspace, then the in-place transform is computed,
   then the contents of the workspace is copied into the output array (taking
   care of zero padding and data cconversion).  If memory alignement is not
   required, out-of-place transforms are always used and the internal
   workspace may be used to avoid destroying the contents of the input array.

   In FFTW, a complex-to-real transform has the following constraints:
    - input and output arrays cannot have the same sizes;
    - preserving input for out-of-place transforms is not possible;
    - out-of-place transform is slightly faster;

   Hence the strategy for complex-to-real transforms is:
   - FFTW2: use out-of-place transforms;
   - FFTW3: if memory alignement is imposed, use in-place transforms (to
     limit the number of copies); otherwise, use out-of-place transforms;
   - for out-of-place transforms, if input cannot be destroyed, the internal
     buffer is used as a scratch array.
*/

#include <string.h>
#include <stdio.h>
#include <yapi.h>
#include <pstdlib.h>
#include <play.h>
#include <ydata.h>

#ifdef YORICK_CHAR_IS_SIGNED
# define BYTE   signed char
#else
# define BYTE unsigned char
#endif

/* Define implementation name and perform minimal checks. */
#if defined(CFFT)
# if defined(FFTW3) || defined(FFTW2) || defined(MKL)
#  define _XFFT_MULTIPLE
# endif
#elif defined(FFTW3)
# if defined(CFFT) || defined(FFTW2) || defined(MKL)
#  define _XFFT_MULTIPLE
# endif
#elif defined(FFTW2)
# if defined(CFFT) || defined(FFTW3) || defined(MKL)
#  define _XFFT_MULTIPLE
# endif
#elif defined(MKL)
# if defined(CFFT) || defined(FFTW2) || defined(FFTW3)
#  define _XFFT_MULTIPLE
# endif
#endif
#ifdef _XFFT_MULTIPLE
# error multiple implementations
#endif
#ifndef XFFT_impl
# error macro XFFT_impl must be defined
#endif
#ifndef XFFT_IMPL
# error macro XFFT_IMPL must be defined
#endif

#define _XFFT_STRINGIFY(a)              #a
#define _XFFT_EXPAND(a)                 a
#define _XFFT_JOIN2(a1,a2)              a1##a2
#define _XFFT_JOIN3(a1,a2,a3)           a1##a2##a3
#define _XFFT_JOIN4(a1,a2,a3,a4)        a1##a2##a3##a4
#define _XFFT_JOIN5(a1,a2,a3,a4,a5)     a1##a2##a3##a4##a5

#define XFFT_STRINGIFY(a)               _XFFT_STRINGIFY(a)
#define XFFT_JOIN(a,b)                  _XFFT_JOIN2(a,b)
#define XFFT_JOIN2(a1,a2)               _XFFT_JOIN2(a1,a2)
#define XFFT_JOIN3(a1,a2,a3)            _XFFT_JOIN3(a1,a2,a3)
#define XFFT_JOIN4(a1,a2,a3,a4)         _XFFT_JOIN4(a1,a2,a3,a4)
#define XFFT_JOIN5(a1,a2,a3,a4,a5)      _XFFT_JOIN5(a1,a2,a3,a4,a5)

#define XFFT_IMPL_NAME  XFFT_STRINGIFY(XFFT_impl)
#define XFFT_PKG_NAME   "xfft_" XFFT_IMPL_NAME

/* Load definitions for FFTW. */
#ifdef FFTW3
# include <fftw3.h>
# define CHOICE(a,b) b
#else
# if defined(FFTW_PREFIX) && (FFTW_PREFIX != 0)
#  ifdef USE_THREADS
#   include <drfftw_threads.h>
#  else
#   include <drfftw.h>
#  endif
# else
#  ifdef USE_THREADS
#   include <rfftw_threads.h>
#  else
#   include <rfftw.h>
#  endif
# endif
# ifdef FFTW_ENABLE_FLOAT
#  error only double precision real supported
# endif
# define CHOICE(a,b) a
#endif

#define FALSE 0
#define TRUE  1


typedef struct _xform xform_t;

static void setup_dimlist(xform_t *xform, long dims[], int forward);
static void create_plan(xform_t *xform, int forward);

/* Indices of keywords. */
static long align_index = -1L;
static long dims_index = -1L;
static long impl_index = -1L;
static long nevals_index = -1L;
static long nthreads_index = -1L;
static long planning_index = -1L;
static long rdims_index = -1L;
static long real_index = -1L;
static long zdims_index = -1L;

/* Methods. */
static void    free_xform(void *);
static void   print_xform(void *);
static void    eval_xform(void *, int);
static void extract_xform(void *, char *);

static y_userobj_t xform_class = {
  /* type_name:  */ "xfft_" XFFT_IMPL_NAME,
  /* on_free:    */    free_xform,
  /* on_print:   */   print_xform,
  /* on_eval:    */    eval_xform,
  /* on_extract: */ extract_xform,
  /* uo_ops:     */ (void *)0
};

#define XFFT_ESTIMATE    0
#define XFFT_MEASURE     1
#define XFFT_PATIENT     2
#define XFFT_EXHAUSTIVE  3

/* These constants are defined (without the leading "Y") in "fftw.i". */
#define XFFT_DIRECT                       0
#define XFFT_CONJUGATE_TRANSPOSE          1
#define XFFT_INVERSE                      2
#define XFFT_INVERSE_CONJUGATE_TRANSPOSE  3
#define XFFT_FORWARD                      XFFT_DIRECT
#define XFFT_BACKWARD                     XFFT_CONJUGATE_TRANSPOSE

struct _xform {
  double scale; /* scaling factor to normalize the FFT */
  long dims[Y_DIMSIZE];
  long r_size; /* number of real's in the real array */
  long z_size; /* number of real's in the complex array */
  long nevals;
  void *forward; /* FFTW plan for forward transform */
  void *backward; /* FFTW plan for backward transform */
  void *ws; /* workspace */
  int nthreads;
  int planning;
  int real;
  int align;
};

static void free_xform(void *ptr)
{
  xform_t *xform = (xform_t *)ptr;
#ifdef FFTW3
  if (xform->forward != NULL) {
    fftw_destroy_plan(xform->forward);
  }
  if (xform->backward != NULL) {
    fftw_destroy_plan(xform->backward);
  }
  if (xform->ws != NULL) {
    if (xform->align) {
      fftw_free(xform->ws);
    } else {
      p_free(xform->ws);
    }
  }
#else /* FFTW2 */
  if (xform->real) {
    if (xform->forward != NULL) {
      rfftwnd_destroy_plan(xform->forward);
    }
    if (xform->backward != NULL) {
      rfftwnd_destroy_plan(xform->backward);
    }
  } else {
    if (xform->forward != NULL) {
      fftwnd_destroy_plan(xform->forward);
    }
    if (xform->backward != NULL) {
      fftwnd_destroy_plan(xform->backward);
    }
  }
  if (xform->ws != NULL) {
    p_free(xform->ws);
  }
#endif /* FFTW2/FFTW3 */
}

static void print_xform(void *ptr)
{
  y_print("xfft operator with "XFFT_IMPL_NAME" implementation", 1);
}

/* Implement the on_extract method to query a member of the object. */
static void extract_xform(void *ptr, char *member)
{
  xform_t *xform = (xform_t *)ptr;
  long index = yget_global(member, 0);

  if (index == dims_index ||
      index == rdims_index ||
      index == zdims_index) {
    long dims_of_dims[2], *xform_dims, *dims, j, ndims;
    xform_dims = xform->dims;
    ndims = xform_dims[0];
    if (ndims < 0) {
      ypush_nil();
    } else {
      dims_of_dims[0] = 1;
      dims_of_dims[1] = ndims + 1;
      dims = ypush_l(dims_of_dims);
      for (j = 0; j <= ndims; ++j) {
	dims[j] = xform_dims[j];
      }
      if (ndims >= 1 && index == zdims_index && xform->real) {
	dims[1] = dims[1]/2L + 1L;
      }
    }
  } else if (index == real_index) {
    ypush_int(xform->real ? TRUE : FALSE);
  } else if (index == align_index) {
    ypush_int(xform->align ? TRUE : FALSE);
  } else if (index == planning_index) {
    ypush_long(xform->planning);
  } else if (index == nthreads_index) {
    ypush_long(xform->nthreads);
  } else if (index == nevals_index) {
    ypush_long(xform->nevals);
  } else if (index == impl_index) {
    ypush_q(NULL)[0] = p_strcpy(XFFT_IMPL_NAME);
  } else {
    ypush_nil();
  }
}

/*---------------------------------------------------------------------------*/
/* CONVERSION FUNCTIONS */

/* For optimal speed, these function perform, at the same time, copy,
   conversion, zero padding and rescaling of the values. */

#if (Y_CHAR != 0)
# error code assumes that Y_CHAR = 0
#endif
#if (Y_SHORT != 1)
# error code assumes that Y_SHORT = 1
#endif
#if (Y_INT != 2)
# error code assumes that Y_INT = 2
#endif
#if (Y_LONG != 3)
# error code assumes that Y_LONG = 3
#endif
#if (Y_FLOAT != 4)
# error code assumes that Y_FLOAT = 4
#endif
#if (Y_DOUBLE != 5)
# error code assumes that Y_DOUBLE = 5
#endif
#if (Y_COMPLEX != 6)
# error code assumes that Y_COMPLEX = 6
#endif

#define IS_INTEGER(id) ((id) >= Y_CHAR && (id) <= Y_LONG)

/* ZPAD_R2H - Zero-pad and convert a real array to form an half-Hermitian
              array.

   Arguments are:
   ID is the data type of the input array;
   INP is the address of the input array;
   OUT is the address of the output array;
   NTOT is the number of elements of the real array (not the half Hermitian one);
   N1 is the leading dimension of the real array (not the half Hermitian one);
   A is the scale factor (1.0 for no scaling).

   The real array is N1-by-N2 (with N2 the product of the trailing dimensions)
   while the half-Hermitian array is H1-by-N2 with H1 = (N1/2 + 1).  Hence,
   the half-Hermitian array has 2*H1*N2 = K1*N2 reals while the real array has
   NTOT = N1*N2 reals.  If N1 is odd, K1 = 2*H1 = N1+1; otherwise
   K1 = 2*H1 = N1+2.
*/
#ifdef FFTW3
#define ZPAD_R2H(id, out, inp, n1, ntot, ptr)  \
  zpad_r2h[id](out, inp, n1, ntot, ptr)

#define FUNCTION(X, TYPE)                                       \
  static void zpad_r2h_##X(double *dst,                         \
                           const void *inp,                     \
                           const long n1,                       \
                           const long ntot,                     \
                           const double a)                      \
  {                                                             \
    const double zero = 0.0;                                    \
    const TYPE *src = (const TYPE *)inp;                        \
    long i, j, k1, n2 = ntot/n1;                                \
    if ((n1 & 1L) != 0L) {                                      \
      k1 = n1 + 1;                                              \
      if (a != 1.0) {                                           \
        for (j = 0; j < n2; ++j, src += n1, dst += k1) {        \
          for (i = 0; i < n1; ++i) {                            \
            dst[i] = ((double)src[i])*a;                        \
          }                                                     \
          dst[n1] = zero;                                       \
        }                                                       \
      } else {                                                  \
        for (j = 0; j < n2; ++j, src += n1, dst += k1) {        \
          for (i = 0; i < n1; ++i) {                            \
            dst[i] = src[i];                                    \
          }                                                     \
          dst[n1] = zero;                                       \
        }                                                       \
      }                                                         \
    } else {                                                    \
      k1 = n1 + 2;                                              \
      if (a != 1.0) {                                           \
        for (j = 0; j < n2; ++j, src += n1, dst += k1) {        \
          for (i = 0; i < n1; ++i) {                            \
            dst[i] = ((double)src[i])*a;                        \
          }                                                     \
          dst[n1] = zero;                                       \
          dst[n1 + 1] = zero;                                   \
        }                                                       \
      } else {                                                  \
        for (j = 0; j < n2; ++j, src += n1, dst += k1) {        \
          for (i = 0; i < n1; ++i) {                            \
            dst[i] = src[i];                                    \
          }                                                     \
          dst[n1] = zero;                                       \
          dst[n1 + 1] = zero;                                   \
        }                                                       \
      }                                                         \
    }                                                           \
  }

FUNCTION(c, BYTE);
FUNCTION(s, short);
FUNCTION(i, int);
FUNCTION(l, long);
FUNCTION(f, float);
FUNCTION(d, double);
#undef FUNCTION
static void (*zpad_r2h[])(double *, const void *, const long,
                          const long, const double) = {
  zpad_r2h_c, zpad_r2h_s, zpad_r2h_i, zpad_r2h_l, zpad_r2h_f, zpad_r2h_d
};
#endif /* FFTW3 */

/* In-place scaling of an array (assuming SCALE != 1). */

static void scale_d(double *arr, const long ntot, const double scale)
{
  long i;
  for (i = 0; i < ntot; ++i) {
    arr[i] *= scale;
  }
}

static void scale_z(double *arr, const long ntot, const double scale)
{
  long i;
  for (i = 0; i < ntot; ++i) {
    arr[2*i] *= scale;
    arr[2*i + 1] *= scale;
  }
}


/* Out-of-place copy, convert and scale an array of any type to a double
   array. */

#define FUNCTION(X,TYPE)                                                \
  static void scale_##X##_to_d(double *dst, const void *inp,            \
                               const long ntot, const double scale)     \
  {                                                                     \
    const TYPE *src = (const TYPE *)inp;                                \
    long i;                                                             \
    if (scale != 1.0) {                                                 \
      for (i = 0; i < ntot; ++i) {                                      \
        dst[i] = src[i]*scale;                                          \
      }                                                                 \
    } else {                                                            \
      for (i = 0; i < ntot; ++i) {                                      \
        dst[i] = src[i];                                                \
      }                                                                 \
    }                                                                   \
  }
FUNCTION(c, BYTE)
FUNCTION(s, short)
FUNCTION(i, int)
FUNCTION(l, long)
FUNCTION(f, float)
FUNCTION(d, double)
#undef FUNCTION

static void (*scale_x_to_d[])(double *, const void *,
                              const long, const double) = {
  scale_c_to_d, scale_s_to_d, scale_i_to_d, scale_l_to_d,
  scale_f_to_d, scale_d_to_d };

#define SCALE_X_TO_D(id, inp, out, ntot, scale) \
  scale_x_to_d[id](inp, out, ntot, scale)


/* Out-of-place copy, convert and scale an array of any type to a complex
   array. */

#define FUNCTION(X,TYPE)                                                \
static void scale_##X##_to_z(double *dst, const void *inp,              \
                             const long ntot, const double scale)       \
{                                                                       \
  const double zero = 0.0;                                              \
  const TYPE *src = (const TYPE *)inp;                                  \
  long i;                                                               \
  if (scale != 1.0) {                                                   \
    for (i = 0; i < ntot; ++i) {                                        \
      dst[2*i] = src[i]*scale;                                          \
      dst[2*i + 1] = zero;                                              \
    }                                                                   \
  } else {                                                              \
    for (i = 0; i < ntot; ++i) {                                        \
      dst[2*i] = src[i];                                                \
      dst[2*i + 1] = zero;                                              \
    }                                                                   \
  }                                                                     \
}
FUNCTION(c, BYTE)
FUNCTION(s, short)
FUNCTION(i, int)
FUNCTION(l, long)
FUNCTION(f, float)
FUNCTION(d, double)
#undef FUNCTION
static void scale_z_to_z(double *dst, const void *inp,
                         const long ntot, const double scale)
{
  const double *src = (const double *)inp;
  long i;
  if (scale != 1.0) {
    for (i = 0; i < ntot; ++i) {
      dst[2*i] = src[2*i]*scale;
      dst[2*i + 1] = src[2*i + 1]*scale;
    }
  } else {
    memcpy(dst, src, (2*sizeof(double))*ntot);
  }
}

static void (*scale_x_to_z[])(double *, const void *,
                              const long, const double) = {
  scale_c_to_z, scale_s_to_z, scale_i_to_z, scale_l_to_z,
  scale_f_to_z, scale_d_to_z, scale_z_to_z };

#define SCALE_X_TO_Z(id, inp, out, ntot, scale) \
  scale_x_to_z[id](inp, out, ntot, scale)

/*---------------------------------------------------------------------------*/

/* Macros for calling out-of-place FFT routines. */

#ifdef FFTW3

# define OUT_OF_PLACE_FFT(INP, OUT) \
  fftw_execute_dft(plan, INP, OUT)

# define OUT_OF_PLACE_FFT_R2C(INP, OUT) \
  fftw_execute_dft_r2c(plan, INP, OUT)

# define OUT_OF_PLACE_FFT_C2R(INP, OUT) \
  fftw_execute_dft_c2r(plan, INP, OUT)

#else /* FFTW2 */

# if USE_THREADS
#  define OUT_OF_PLACE_CALL(PFX,SFX,INP,OUT)                       \
     if (xform->nthreads > 1)                                      \
       PFX##_threads##SFX(xform->nthreads, plan, INP, OUT);        \
     else                                                          \
       PFX##SFX(plan, INP, OUT)
# else
#  define OUT_OF_PLACE_CALL(PFX,SFX,INP,OUT) \
       PFX##SFX(plan, INP, OUT)
# endif

# define OUT_OF_PLACE_FFT(INP, OUT)                             \
  OUT_OF_PLACE_CALL(fftwnd,_one,INP,OUT)

# define OUT_OF_PLACE_FFT_R2C(INP, OUT)                         \
  OUT_OF_PLACE_CALL(rfftwnd,_one_real_to_complex,INP,OUT)

# define OUT_OF_PLACE_FFT_C2R(INP, OUT)                         \
  OUT_OF_PLACE_CALL(rfftwnd,_one_complex_to_real,INP,OUT)

#endif /* FFTW3 or FFTW2 */

/* iarg = argc - k is stack index of k-th argument of the builtin function */

/*
 *   rank = 0, scalar array, FFT is a no-op
 *   rank = 1
 */
#define COMPLEX_TO_COMPLEX  0
#define COMPLEX_TO_REAL     1
#define REAL_TO_COMPLEX     2

static void eval_xform(void *ptr, int argc)
{
  double scale;
  xform_t *xform = (xform_t *)ptr;
  void *inp, *out, *tmp, *plan;
  long ntot, index;
#ifdef FFTW3
  long inp_dim1, out_dim1;
  int rank;
#endif
  int arg_type, xform_type, forward, overwrite, rescale, job, scratch;
  long dims[Y_DIMSIZE];

  /* Get the direction of the transform and check number of arguments. */
  if (argc == 2) {
    arg_type = yarg_typeid(0);
    if (IS_INTEGER(arg_type) && yarg_rank(0) == 0) {
      job = ygets_i(0);
    } else if (arg_type == Y_VOID) {
      job = 0;
    } else {
      job = -1;
    }
    switch (job) {
    case XFFT_DIRECT:
      forward = TRUE;
      rescale = FALSE;
      break;
    case XFFT_CONJUGATE_TRANSPOSE:
      forward = FALSE;
      rescale = FALSE;
      break;
    case XFFT_INVERSE:
      forward = FALSE;
      rescale = TRUE;
      break;
    case XFFT_INVERSE_CONJUGATE_TRANSPOSE:
      forward = TRUE;
      rescale = TRUE;
      break;
    default:
      y_error("bad job");
      return; /* to avoid compiler warnings */
    }
    yarg_drop(1);
  } else {
    /* Default is same as DIRECT. */
    if (argc != 1) {
      y_error("syntax: op(a) or op(a, job) with op the XFFT operator");
    }
    forward = TRUE;
    rescale = FALSE;
  }

  /* Figure out whether or not to perform in-place operation and, if needed,
     set INDEX to save result in a variable.  FIXME: The rationale of
     yarg_scratch is strange: it returns 1 (hence true) when the argument is a
     variable reference; this is why, I combine the result of yarg_scratch and
     yget_ref to properly figure out whether an argument can be used
     in-place. */
  index = yget_ref(0);
  scratch = (index < 0 ? yarg_scratch(0) : FALSE);
  if (yarg_subroutine()) {
    /* FIXME: it should be possible to perform in-place operation for other
       kind of arguments such as hash table members. */
    if (index < 0) {
      y_error("when called as a subroutine, argument must be " \
              "a simple variable, not a temporary expression");
    }
    overwrite = TRUE;
  } else {
    overwrite = (scratch != 0);
    index = -1L; /* avoids setting global symbol at the end */
  }

  /* Check data type of input array. */
  inp = ygeta_any(0, &ntot, dims, &arg_type);
  if (arg_type < 0 || arg_type > Y_COMPLEX) {
    y_error("invalid non-numerical data type");
  }
  if (xform->real) {
    if (forward) {
      if (arg_type == Y_COMPLEX) {
	y_error("invalid complex input for forward real-complex transform");
      }
      xform_type = REAL_TO_COMPLEX;
    } else {
      xform_type = COMPLEX_TO_REAL;
    }
  } else {
    xform_type = COMPLEX_TO_COMPLEX;
  }

  /* Check dimension list of the input array and get dimension list of the
     result.  Must be done prior to workspace allocation and plan creation
     because this may set the dimension list of the transform. */
#ifdef FFTW3
  rank = dims[0];
  inp_dim1 = (rank > 0 ? dims[1] : 0); /* save 1st dimension of input */
#endif
  setup_dimlist(xform, dims, forward);
#ifdef FFTW3
  out_dim1 = (rank > 0 ? dims[1] : 0); /* save 1st dimension of output */
#endif
  if (ntot == 1) {
    /* Transform of a scalar is the identity so, at most, it is just a matter
       of converting the data type of the input. */
    if (xform_type == COMPLEX_TO_REAL) {
      out = ygeta_d(0, NULL, NULL);
    } else {
      out = ygeta_z(0, NULL, NULL);
    }
    goto done;
  }

  /* Choose the plan, allocate it if necessary. */
  if (forward) {
    if (xform->forward == NULL) create_plan(xform, TRUE);
    plan = xform->forward;
  } else {
    if (xform->backward == NULL) create_plan(xform, FALSE);
    plan = xform->backward;
  }
  scale = (rescale ? xform->scale : 1.0);

  if (xform_type == COMPLEX_TO_COMPLEX) {

    /**************************/
    /* COMPLEX-TO-COMPLEX FFT */
    /**************************/

#ifdef FFTW3
    if (xform->align) {
      /* Perform in-place transform with internal workspace. */
      tmp = xform->ws;
      SCALE_X_TO_Z(arg_type, tmp, inp, ntot, scale);
      fftw_execute(plan);
      if (overwrite && arg_type == Y_COMPLEX) {
        out = inp;
        index = -1L; /* no needs to set global symbol at the end */
      } else {
        yarg_drop(1);
        out = ypush_z(dims);
      }
      memcpy(out, tmp, (2*sizeof(double))*ntot);
      goto done;
    }
#endif

    /* Perform out-of-place transform with destroyable input. */
    if (overwrite && arg_type == Y_COMPLEX) {
      tmp = inp;
      if (scale != 1.0) {
        scale_z(tmp, ntot, scale);
      }
    } else {
      tmp = xform->ws;
      SCALE_X_TO_Z(arg_type, tmp, inp, ntot, scale);
      yarg_drop(1);
    }
    out = ypush_z(dims);
    OUT_OF_PLACE_FFT(tmp, out);

  } else if (xform_type == REAL_TO_COMPLEX) {

    /***********************/
    /* REAL-TO-COMPLEX FFT */
    /***********************/

#ifdef FFTW3
    if (xform->align) {
      tmp = xform->ws;
      ZPAD_R2H(arg_type, tmp, inp, inp_dim1, ntot, scale);
      yarg_drop(1);
      fftw_execute(plan); /* FIXME: fftw_execute_dft_r2c(plan, tmp, tmp); */
      out = ypush_z(dims);
      memcpy(out, tmp, xform->z_size*sizeof(double));
      goto done;
    }
#endif

    if (overwrite && arg_type == Y_DOUBLE) {
      tmp = inp;
      if (scale != 1.0) {
        scale_d(tmp, ntot, scale);
      }
    } else {
      tmp = xform->ws;
      SCALE_X_TO_D(arg_type, tmp, inp, ntot, scale);
      yarg_drop(1);
    }
    out = ypush_z(dims);
    OUT_OF_PLACE_FFT_R2C(tmp, out);

  } else {

    /***********************/
    /* COMPLEX-TO-REAL FFT */
    /***********************/

#ifdef FFTW3
    if (xform->align) {
      /* Perform in-place transform, then perform a packed copy of the
         half-Hermitian result into a real array. */
      double *dst;
      const double *src;
      long i, j, k1, n1, n2;
      tmp = xform->ws;
      SCALE_X_TO_Z(arg_type, tmp, inp, ntot, scale);
      fftw_execute(plan); /* FIXME: fftw_execute_dft_c2r(plan, tmp, tmp); */
      src = (const double *)tmp;
      dst = (double *)ypush_d(dims);
      n1 = out_dim1;
      n2 = ntot/inp_dim1;
      k1 = 2L*inp_dim1;
      for (j = 0L; j < n2; ++j, src += k1, dst += n1) {
        for (i = 0L; i < n1; ++i) {
          dst[i] = src[i];
        }
      }
      goto done;
    }
#endif

    /* Perform out-of-place transform. */
    if (overwrite && arg_type == Y_COMPLEX) {
      tmp = inp;
      if (scale != 1.0) {
        scale_z(tmp, ntot, scale);
      }
    } else {
      tmp = xform->ws;
      SCALE_X_TO_Z(arg_type, tmp, inp, ntot, scale);
      yarg_drop(1);
    }
    out = ypush_d(dims);
    OUT_OF_PLACE_FFT_C2R(tmp, out);

  }

 done:
  /* Increment the number of evaluations and, possibly, save the result. */
  ++xform->nevals;
  if (index >= 0L) {
    yput_global(index, 0);
  }
}

static void create_plan(xform_t *xform, int forward)
{
#ifdef FFTW3
  void *inp, *out;
  double dummy[5];
#endif
  int j, rank, n[Y_DIMSIZE - 1];
  unsigned int flags;

  /* Shortcut. */
  if ((forward ? xform->forward : xform->backward) != NULL) {
    return;
  }

  /* Copy the list of dimensions.  FFTW dimensions are in row-major format
     (the first dimension's index varies most slowly and the last dimension's
     index varies most quickly) whereas Yorick's dimensions are in
     column-major format. Also remember that the first element of a Yorick's
     dimension list is the "rank"; that is, the number of dimensions. */
  rank = xform->dims[0];
  for (j = 0; j < rank; ++j) {
    n[j] = (int)xform->dims[rank - j];
  }

  /* Set the flags for the transform. */
#ifdef FFTW3
  if (xform->align) {
    /* Will always use the internal buffer for in-place transforms. */
    flags =  FFTW_DESTROY_INPUT;
  } else {
    /* Will always use a temporary input for out-of-place transforms. */
    flags =  FFTW_UNALIGNED | FFTW_DESTROY_INPUT;
  }
  if (xform->planning == XFFT_ESTIMATE) {
    flags |= FFTW_ESTIMATE;
  } else if (xform->planning == XFFT_MEASURE) {
    flags |= FFTW_MEASURE;
  } else if (xform->planning == XFFT_PATIENT) {
    flags |= FFTW_PATIENT;
  } else if (xform->planning == XFFT_EXHAUSTIVE) {
    flags |= FFTW_EXHAUSTIVE;
  }
#else /* FFTW2 */
  if (xform->planning == XFFT_ESTIMATE) {
    flags = FFTW_OUT_OF_PLACE | FFTW_ESTIMATE;
  } else {
    flags = FFTW_OUT_OF_PLACE | FFTW_MEASURE;
  }
#endif /* FFTW3 or FFTW2 */

  /* Allocate workspace for the transforms (see top comments for the rationale
     about the in-place / out-of-place transforms and alignement of
     memory). */
  if (xform->ws == NULL) {
#ifdef FFTW3
    if (xform->align) {
      xform->ws = fftw_malloc(xform->z_size*sizeof(double));
    } else {
      xform->ws = p_malloc(xform->z_size*sizeof(double));
    }
#else /* FFTW2 */
    xform->ws = p_malloc(xform->z_size*sizeof(double));
#endif /* FFTW2 or FFTW3 */
    if (xform->ws == NULL) {
      y_error("failed to allocate FFTW workspace");
    }
  }

  /* Create the plan. */
#ifdef FFTW3
  if (xform->align) {
    inp = out = xform->ws;
  } else if (xform->planning == XFFT_ESTIMATE) {
    /* Use different addresses with enough spacing for a complex and different
       alignments. */
    inp = &dummy[0];
    out = &dummy[3];
  } else {
    /* Other planning methods will require to perform some FFT's so input and
       output arrays must be given.  Here we allocated different arrays since
       we will use out-of-place transforms. */
    inp = ypush_scratch(xform->z_size*2*sizeof(double), NULL);
    out = ((double *)inp) + xform->z_size;
  }
# if USE_THREADS
  fftw_plan_with_nthreads(xform->nthreads);
# endif
#endif
  if (forward) {
    if (xform->real) {
      xform->forward = CHOICE(rfftwnd_create_plan(rank, n, FFTW_FORWARD, flags),
                              fftw_plan_dft_r2c(rank, n, inp, out, flags));
    } else {
      xform->forward = CHOICE(fftwnd_create_plan(rank, n, FFTW_FORWARD, flags),
                              fftw_plan_dft(rank, n, inp, out, FFTW_FORWARD, flags));
    }
    if (xform->forward == NULL) {
      y_error("failed to allocate FFTW plan for forward transform");
    }
  } else {
    if (xform->real) {
      xform->backward = CHOICE(rfftwnd_create_plan(rank, n, FFTW_BACKWARD, flags),
                               fftw_plan_dft_c2r(rank, n, inp, out, flags));
    } else {
      xform->backward = CHOICE(fftwnd_create_plan(rank, n, FFTW_BACKWARD, flags),
                               fftw_plan_dft(rank, n, inp, out, FFTW_BACKWARD, flags));
    }
    if (xform->backward == NULL) {
      y_error("failed to allocate FFTW plan for backward transform");
    }
  }
#ifdef FFTW3
  /* Drop temporary arrays used to benchmark the plans. */
  if (! xform->align && xform->planning != XFFT_ESTIMATE) {
    yarg_drop(1);
  }
#endif
}

/* Usage:
 *
 *    setup_dimlist(xform, dims, forward);
 *
 * with XFORM = FFT object, DIMS = on entry (on return), dimension list of input
 * (output) array, FORWARD = true for forward transform.
 *
 * This private function perform several tasks:
 *
 *   1. Initialize dimension list built-in the FFT object, if not yet done;
 *      otherwise, check that dimension list of input array matches the
 *      built-in one.
 *
 *   2. Store the dimension list of the result of the transform in DIMS.
 *
 * In case of error, y_error is called so that there is no return.
 */
static void setup_dimlist(xform_t *xform, long dims[], int forward)
{
  long j, ndims, r_dim1, z_dim1, *xform_dims, r_size, z_size;

  xform_dims = xform->dims;

  if (xform_dims[0] < 0L) {
    /* Dimension list never specified. */
    if (! forward && xform->real) {
      y_error("cannot guess real-complex dimension list " \
              "from backward transform");
    }
    ndims = dims[0];
    if (ndims >= Y_DIMSIZE) {
      y_error("too many dimensions");
    }
    xform_dims[0] = ndims;
    r_size = 1L;
    for (j = 1; j <= ndims; ++j) {
      xform_dims[j] = dims[j];
      r_size *= dims[j];
    }
    if (ndims >= 1L && xform->real) {
      r_dim1 = dims[1];
      z_dim1 = r_dim1/2L + 1L;
      z_size = 2L*(r_size/r_dim1)*z_dim1;
      dims[1] = z_dim1;
    } else {
      z_size = 2L*r_size;
    }
    xform->scale = 1.0/r_size;
    xform->z_size = z_size;
    xform->r_size = r_size;
    return;
  }

  /* Check that dimension lists match. */
  ndims = xform_dims[0];
  if (dims[0] != ndims) {
    goto mismatch;
  }
  if (ndims >= 1L) {
    if (xform->real) {
      /* Real-complex transform: do something special with the first
	 dimension. */
      z_dim1 = xform_dims[1]/2L + 1L;
      if (forward) {
	if (dims[1] != xform_dims[1]) {
	  goto mismatch;
	}
	dims[1] = z_dim1;
      } else {
	if (dims[1] != z_dim1) {
	  goto mismatch;
	}
	dims[1] = xform_dims[1];
      }
      j = 1L;
    } else {
      /* Complex-complex transform: nothing special with the first
	 dimension. */
      j = 0L;
    }
    while (++j <= ndims) {
      if (dims[j] != xform_dims[j]) {
	goto mismatch;
      }
    }
  }
  return;

 mismatch:
  y_error("incompatible dimension list");
}

static int first_time = TRUE;
static void initialize(void)
{
  if (first_time) {
    align_index = yget_global("align", 0);
    dims_index = yget_global("dims", 0);
    impl_index = yget_global("impl", 0);
    nevals_index = yget_global("nevals", 0);
    nthreads_index = yget_global("nthreads", 0);
    planning_index = yget_global("planning", 0);
    rdims_index = yget_global("rdims", 0);
    real_index = yget_global("real", 0);
    zdims_index = yget_global("zdims", 0);
#ifdef USE_THREADS
    if (CHOICE(fftw_threads_init() != 0, fftw_init_threads() == 0)) {
      y_error("initialization of FFTW threads failed");
    }
#endif
    first_time = FALSE;
  }
}

/*---------------------------------------------------------------------------*/
/* BUILTIN FUNCTIONS */

/* Macro for name of builtin functions. */
#define XFFT_BUILTIN(name) XFFT_JOIN4(xfft_,XFFT_impl,_,name)

static void XFFT_BUILTIN(new)(int argc)
{
  xform_t *xform;
  long index, planning, *dims, dimlist[Y_DIMSIZE];
  int iarg, id, flag, j, real, align, nthreads;

  /* Setup internals. */
  if (first_time) {
    initialize();
  }

  /* Parse arguments. */
  if (yarg_subroutine()) {
    y_error("xfft_new must be called as a function");
  }
  planning = XFFT_ESTIMATE;
  real = FALSE;
  nthreads = 1;
  align = FALSE;
  dims = NULL;
  if (argc % 2 == 1) {
    if (argc != 1 || ! yarg_nil(0)) {
      goto expectingKeywords;
    }
  } else {
    for (iarg = argc - 1; iarg >= 0; --iarg) {
      index = yarg_key(iarg);
      if (index < 0L) {
      expectingKeywords:
        y_error("all arguments must be keywords");
        return;
      }
      if (index == align_index) {
	id = yarg_typeid(--iarg);
	if (IS_INTEGER(id)) {
	  align = (ygets_l(iarg) != 0L ? TRUE : FALSE);
	} else if (id != Y_VOID) {
	  y_error("bad value for ALIGN keyword");
	}
      } else if (index == real_index) {
	id = yarg_typeid(--iarg);
	if (IS_INTEGER(id)) {
	  real = (ygets_l(iarg) != 0L ? TRUE : FALSE);
	} else if (id != Y_VOID) {
	  y_error("bad value for REAL keyword");
	}
      } else if (index == nthreads_index) {
	id = yarg_typeid(--iarg);
	if (IS_INTEGER(id)) {
	  nthreads = ygets_i(iarg);
        } else if (id != Y_VOID) {
          nthreads = -1; /* make sure to trigger an error */
        }
        if (nthreads <= 0) {
	  y_error("bad value for NTHREADS keyword");
        }
      } else if (index == planning_index) {
	id = yarg_typeid(--iarg);
	flag = TRUE; /* used to trigger an error */
	if (IS_INTEGER(id)) {
	  planning = ygets_l(iarg);
	  if (planning == XFFT_ESTIMATE || planning == XFFT_MEASURE ||
	      planning == XFFT_PATIENT || planning == XFFT_EXHAUSTIVE) {
#ifndef FFTW3
	    /* Only ESTIMATE or MEASURE implemented by FFTW2. */
	    if (planning != XFFT_ESTIMATE) {
	      planning = XFFT_MEASURE;
	    }
#endif
	    flag = FALSE;
	  }
	} else if (id == Y_VOID) {
	  flag = FALSE;
	}
	if (flag) {
	  y_error("bad value for PLANNING keyword");
	}
      } else if (dims == NULL && index == dims_index) {
	id = yarg_typeid(--iarg);
        if (id != Y_VOID) {
          if (IS_INTEGER(id)) {
            dims = ygeta_l(iarg, NULL, dimlist);
            if (dimlist[0] != 1 || dimlist[1] != dims[0] + 1L) {
              dims = NULL;
            } else {
              /* Check and copy dimension list (a copy is needed because
                 setup_dimlist may overwrite the first dimension). */
              dimlist[0] = dims[0];
              for (j = 1; j <= dims[0]; ++j) {
                if (dims[j] < 1L) {
                  dims = NULL;
                  break;
                }
                dimlist[j] = dims[j];
              }
              dims = dimlist;
            }
          }
          if (dims == NULL) {
            y_error("bad value for DIMS keyword");
          }
        }
      } else {
	y_error("unknown keyword");
      }
    }
  }
#ifndef FFTW3
  align = FALSE;
#endif
#ifndef USE_THREADS
  nthreads = 1;
#endif

  /* Create FFTW object and setup dimension list. */
  xform = (xform_t *)ypush_obj(yfunc_obj(&xform_class), sizeof(xform_t));
  xform->scale = 1.0;
  xform->dims[0] = -1L; /* not yet initialized */
  xform->r_size = 0L;
  xform->z_size = 0L;
  xform->nevals = 0L;
  xform->forward = NULL;
  xform->backward = NULL;
  xform->ws = NULL;
  xform->nthreads = nthreads;
  xform->planning = planning;
  xform->real = real;
  xform->align = align;
  if (dims != NULL) {
    setup_dimlist(xform, dims, 1);
  }
}

static void XFFT_BUILTIN(version)(int argc)
{
  ypush_q(NULL)[0] = p_strcpy(XFFT_VERSION);
}

static void XFFT_BUILTIN(threads)(int argc)
{
#ifdef USE_THREADS
  ypush_int(TRUE);
#else
  ypush_int(FALSE);
#endif
}

static void XFFT_BUILTIN(implementation)(int argc)
{
  ypush_q(NULL)[0] = p_strcpy(XFFT_IMPL_NAME);
}

static void XFFT_BUILTIN(best_dim)(int argc)
{
  /* FFTW is best at handling sizes of the form 2^a 3^b 5^c 7^d 11^e 13^f,where e+f
     is either 0 or 1, and the other exponents are arbitrary. Other sizes are
     computed by means of a slow, general-purpose algorithm (which
     nevertheless retains O(n log n) performance even for prime sizes). It is
     possible to customize FFTW for different array sizes; see Installation
     and Customization. Transforms whose sizes are powers of 2 are especially
     fast. */
  long len, best, k, i2, i3, i5, i7;

  if (argc != 1) y_error("xfft_best_dim takes exactly one argument");
  len = ygets_l(0);
  if (len <= 0) y_error("invalid dimension length");
  best = 2*len;
  k = 1;
  while (1) {
    for (i7 = k; i7 < best; i7 *= 7) {
      for (i5 = i7; i5 < best; i5 *= 5) {
	for (i3 = i5; i3 < best; i3 *= 3) {
	  /* last loop (power of 2) is exited as soon as N >= LEN */
	  for (i2 = i3; i2 < len; i2 *= 2)
	    ; /* empty loop body */
	  if (i2 == len) {
            ypush_long(len);
            return;
          }
	  if (i2 < best) {
            best = i2;
          }
	}
      }
    }
    if (k == 1) {
      k = 11;
    } else if (k == 11) {
      k = 13;
    } else {
      ypush_long(best);
      return;
    }
  }
}

static void XFFT_BUILTIN(indgen)(int argc)
{
  long *p, dims[2];
  long i, n, k;
  int half;

  if (argc != 1 && argc != 2) y_error("xfft_indgen takes one or two argument");
  if (argc == 2) {
    half = yarg_true(0);
  } else {
    half = FALSE;
  }
  n = ygets_l(argc - 1);
  if (n < (half ? 0L : 1L)) y_error("invalid dimension");
  k = n/2L + 1L;
  dims[0] = 1;
  dims[1] = (half ? k : n);
  p = ypush_l(dims);
  for (i = 0; i < k; ++i) {
    p[i] = i;
  }
  if (! half) {
    for (i = k; i < n; ++i) {
      p[i] = i - n;
    }
  }
}

#if 0
static void XFFT_BUILTIN(dummy)(int argc)
{
  fprintf(stderr, "argc = %d\n", argc);
  ypush_nil();
}

#endif

static void XFFT_BUILTIN(export_wisdom)(int argc)
{
  char *filename, *native;
  int err;
  if (argc != 1 || yarg_rank(0) != 0 || yarg_typeid(0) != Y_STRING) {
    y_error("xfft_export_wisdom takes a single scalar string argument");
  }
  filename = ygets_q(0);
  if (filename == NULL || filename[0] == 0) {
    y_error("bad filename");
  }
  native = p_native(filename);
#ifdef FFTW3
  err = (fftw_export_wisdom_to_filename(native) == 0);
#else
  {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
      err = TRUE;
    } else {
      fftw_export_wisdom_to_file(file);
      fclose(file);
      err = FALSE;
    }
  }
#endif
  p_free(native);
  if (err) {
    y_error("failed to export wisdom to file");
  }
  ypush_nil();
}

static void XFFT_BUILTIN(import_wisdom)(int argc)
{
  char *filename, *native;
  int err;
  if (argc != 1 || yarg_rank(0) != 0 || yarg_typeid(0) != Y_STRING) {
    y_error("xfft_import_wisdom takes a single scalar string argument");
  }
  filename = ygets_q(0);
  if (filename == NULL || filename[0] == 0) {
    y_error("bad filename");
  }
  native = p_native(filename);
#ifdef FFTW3
  err = (fftw_import_wisdom_from_filename(native) == 0);
#else
  {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
      err = TRUE;
    } else {
      err = (fftw_import_wisdom_from_file(file) != FFTW_SUCCESS);
      fclose(file);
    }
  }
#endif
  p_free(native);
  if (err) {
    y_error("failed to import wisdom from file");
  }
  ypush_nil();
}

static void XFFT_BUILTIN(forget_wisdom)(int argc)
{
  if (! yarg_subroutine() || argc > 1 || (argc == 1 && ! yarg_nil(0))) {
    y_error("xfft_forget_wisdom" \
            " must be called as a subroutine with no arguments");
  }
  fftw_forget_wisdom();
}

/*---------------------------------------------------------------------------*/
/* PACKAGE INITIALIZATION */

/* Package name and list of include files */

static char *y0_pkgname = XFFT_PKG_NAME;

static char *y0_includes[] = {
  XFFT_PKG_NAME ".i",
  0,
  0
};

/* Collect pointers and names for built-in functions and data. */

#define XFFT_FUNC_LIST                     \
  XFFT_FUNC_ITEM(new),                     \
  XFFT_FUNC_ITEM(version),                 \
  XFFT_FUNC_ITEM(threads),                 \
  XFFT_FUNC_ITEM(implementation),          \
  XFFT_FUNC_ITEM(best_dim),                \
  XFFT_FUNC_ITEM(indgen),                  \
  XFFT_FUNC_ITEM(export_wisdom),           \
  XFFT_FUNC_ITEM(import_wisdom),           \
  XFFT_FUNC_ITEM(forget_wisdom),

#define XFFT_DATA_LIST /* no exported data */

#define XFFT_FUNC_ITEM(n) &XFFT_BUILTIN(n)
static BuiltIn *y0_routines[] = {XFFT_FUNC_LIST 0};
#undef XFFT_FUNC_ITEM

#define XFFT_DATA_ITEM(n) &n
static void *y0_values[] = {XFFT_DATA_LIST 0};
#undef XFFT_DATA_ITEM

#define XFFT_FUNC_ITEM(n) "xfft_" XFFT_STRINGIFY(n)
#define XFFT_DATA_ITEM(n) "xfft_" XFFT_STRINGIFY(n)
static char *y0_names[] = {XFFT_FUNC_LIST XFFT_DATA_LIST 0};
#undef XFFT_FUNC_ITEM
#undef XFFT_DATA_ITEM

/* Define package initialization function. */

#define YK_XFFT XFFT_JOIN(yk_xfft_,XFFT_impl)

PLUG_EXPORT char *YK_XFFT(char ***, BuiltIn ***,
                          void ***, char ***);

char *YK_XFFT(char ***ifiles, BuiltIn ***code,
                          void ***data, char ***varname)
{
  *ifiles = y0_includes;
  *code = y0_routines;
  *data = y0_values;
  *varname = y0_names;
  return y0_pkgname;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * End:
 */
