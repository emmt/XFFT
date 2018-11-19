/*
 * xfft.i -
 *
 * Interface for eXtended FFT in Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2009-2018: Éric Thiébaut <https://github.com/emmt/XFFT>
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

local XFFT_DIRECT, XFFT_FORWARD, XFFT_CONJUGATE_TRANSPOSE, XFFT_BACKWARD;
local XFFT_INVERSE, XFFT_INVERSE_CONJUGATE_TRANSPOSE;
local XFFT_ESTIMATE, XFFT_MEASURE, XFFT_PATIENT, XFFT_EXHAUSTIVE;
local xfft_new;
/* DOCUMENT op = xfft_new(...);

     This function creates a new FFT operator, all arguments are given as
     keywords.  Once created, the operator OP can be used as a function:

       z = op(a, job); // compute the FFT of A
       op, a, job;     // compute in-place FFT

     where optional argument JOB is:

       0 or XFFT_DIRECT or XFFT_FORWARD, for the direct transform;
       1 or XFFT_CONJUGATE_TRANSPOSE or XFFT_BACKWARD, for the backward
            transform;
       2 or XFFT_INVERSE, for the inverse transform;
       3 or XFFT_INVERSE_CONJUGATE_TRANSPOSE, for the conjugate
            transpose of the inverse transform;

     By default JOB=0 (so the forward or direct transform is performed).  The
     inverse transform is just the backward transform normalized (i.e. divided
     by the number of elements of the array).  The conjugate transpose of the
     inverse transform is just the backward transform normalized.

     The operator OP can be queried as a structured object:

       op.align       // 1 if OP effectively uses aligned memory, 0 else
       op.real        // 1 if OP is for real-complex transform, 0 else
       op.planning    // value of the planning strategy (see below)
       op.dims        // dimension list of the transform
       op.rdims       // dimension list of real arrays
       op.zdims       // dimension list of complex arrays
       op.nthreads    // effective number of threads used by OP
       op.nevals      // number of evaluation of the transform
       op.impl        // name of the FFT implementation for OP

     For a complex transform the members RDIMS, ZDIMS and DIMS are identical;
     for a real-complex transform, RDIMS and DIMS are identical but the first
     dimension of ZDIMS may be different: ZDIMS(2) = DIMS(2)/2 + 1, other
     dimensions are identical.

     Keyword DIMS can be used to preset the dimension list of the
     transform. By default, the dimension list is obtained the first time the
     operator is used.  However for a real-complex transform, if the first
     operation is a backward transform it is not possible to correctly guess
     the dimension list and this keyword muts be used.

     Keyword REAL must be set true to create an operator for real-complex
     transforms.

     Keyword NTHREADS can be set to specify the number of threads to use.  If
     multi-threading is not supported --- you can check this with
     xfft_threads() --- the number of threads will always be 1 whatever the
     value of the NTHREADS keyword.

     Keyword ALIGN can be set true to apply FFT on aligned memory.  This
     involves using workspaces specially allocated to be properly aligned and,
     hence, some overheads to copy the values between Yorick arrays and these
     workspaces.  If using aligned memory is not supported by the FFT
     implementation, the value of ALIGN keyword is ignored.

     Keyword PLANNING can be set to specify the strategy to build the plan
     when FFTW is used to compute the transform.  The possible planning values
     are: XFFT_ESTIMATE (default), XFFT_MEASURE, XFFT_PATIENT, or
     XFFT_EXHAUSTIVE.  With FFTW-2, XFFT_PATIENT and XFFT_EXHAUSTIVE are the
     same as XFFT_MEASURE.


   SEE ALSO: fft, xfft_threads, xfft_export_wisdom.
 */

/* These constants should match those in "linop.i" and in "xfft.c". */
XFFT_ESTIMATE   = 0;
XFFT_MEASURE    = 1;
XFFT_PATIENT    = 2;
XFFT_EXHAUSTIVE = 3;
XFFT_FORWARD = XFFT_DIRECT = 0;
XFFT_BACKWARD = XFFT_CONJUGATE_TRANSPOSE = 1;
XFFT_INVERSE = 2;
XFFT_INVERSE_CONJUGATE_TRANSPOSE = 3;

func xfft_info(..)
/* DOCUMENT xfft_info;
         or xfft_info, op1, op2, ...;
      Without any arguments, this subroutine prints some informations about
      the current XFTT plugin.  Otherwise, this subroutine prints informations
      about the XFFT operators OP1, OP2, etc.

   SEE ALSO xfft_new.
 */
{
  local op;
  n = 0;
  while (more_args()) {
    eq_nocopy, op, next_arg();
    if (! is_void(op)) {
      p = op.planning;
      if (p == XFFT_ESTIMATE) {
        p = "XFFT_ESTIMATE";
      } else if (p == XFFT_MEASURE) {
        p = "XFFT_MEASURE";
      } else if (p == XFFT_PATIENT) {
        p = "XFFT_PATIENT";
      } else if (p == XFFT_EXHAUSTIVE) {
        p = "XFFT_EXHAUSTIVE";
      } else {
        p = "*UNKNOWN*";
      }
      write, format="implentation: %s, nevals: %d, nthreads: %d, dims: %s, real: %d, align: %d, planning: %s\n",
        op.impl, op.nevals, op.nthreads, sum(print(op.dims)), op.real, op.align, p;
      ++n;
    }
  }
  if (n < 1) {
    write, format="XFFT version: %s\n", xfft_version();
    write, format="XFFT implementation: %s\n", xfft_implementation();
    write, format="XFFT multi-threading: %s\n", (xfft_threads() ? "yes" : "no");
    return;
  }
}

local xfft_best_dim;
/* DOCUMENT xfft_best_dim(len);
     Returns the smallest integer which is greater or equal LEN and which is a
     multiple of suitable powers of 2, 3, 5, 7, 11 and 13.
   SEE ALSO xfft_indgen, xfft_new. */

local xfft_indgen;
/* DOCUMENT xfft_indgen(len);
         or xfft_indgen(len, half);
     Return FFT frequencies along a dimension of length LEN.  If optional
     argument HALF is true, only returns the first half range corresponding to
     the positive frequencies (same as indgen(0:LEN/2)).  This is to deal with
     the first dimension of real to complex transforms.

   SEE ALSO xfft_best_dim, xfft_new. */

local xfft_version;
local xfft_implementation;
local xfft_threads;
/* DOCUMENT xfft_version();
         or xfft_implementation();
	 or xfft_threads();

     The function xfft_version() returns a string with the version of XFFT
     plug-in.

     The function xfft_implementation() returns a string with the name of the
     FFT implementation currently used by XFFT.

     The function xfft_threads() returns a boolean value indicating whether
     multi-threading is supported by the FFT implementation currently used by
     XFFT.


   SEE ALSO: xfft_new. */

local xfft_export_wisdom;
local xfft_import_wisdom;
local xfft_forget_wisdom;
/* DOCUMENT xfft_export_wisdom, FILENAME;
         or xfft_import_wisdom, FILENAME;
         or xfft_forget_wisdom;

     These routines deal with the wisdom memorized by FFTW implementation
     about the best strategy found so far for how to compute FFT of different
     sizes.  If FFT is not implemented by FFTW, these functions does nothing.

     The subroutine xfft_export_wisdom saves FFTW wisdom currently accumulated
     into file FILENAME.

     The subroutine xfft_import_wisdom loads FFTW wisdom saved into file
     FILENAME.

     The subroutine xfft_forget_wisdom discards any FFTW wisdom accumulated so
     far.

   SEE ALSO: wfft_new.
 */

func xfft_try_load(name)
/* DOCUMENT xfft_try_load(name);

     The function xfft_try_load() attempts to load a Yorick plug-in catching
     errors.  Return 0 on success, -1 on failure.  Argument NAME is the name
     of the plug-in, it may have directory separators '/' to force setting
     (temporarily) the plug-in directory.

   SEE ALSO: plug_in, plug_dir, catch. */
{
  index = strfind("/", name, back=1n)(2);
  if (index >= 0) {
    old_dirs = plug_dir(strpart(name, 1:index));
    name = strpart(name, index+1:0);
  }
  if (catch(0x08)) {
    if (index >= 0) {
      plug_dir, old_dirs;
    }
    return -1;
  }
  plug_in, name;
  if (index >= 0) {
    plug_dir, old_dirs;
  }
  return 0;
}

func xfft_try_include(impl)
{
  extern xfft_new;
  old_xfft_new = xfft_new;
  xfft_new = [];
  include, "xfft_"+impl+".i", 3;
  if (is_func(xfft_new) == 2) {
    return 0;
  }
  xfft_new = old_xfft_new;
  return -1;
}

if (is_func(xfft_new) != 2 &&
    xfft_try_include("fftw3_threads") != 0 &&
    xfft_try_include("fftw2_threads") != 0 &&
    xfft_try_include("fftw3") != 0 &&
    xfft_try_include("fftw2") != 0 &&
    xfft_try_include("cfft") != 0) {
  error, "no implementation for XFFT plug-in found";
}
