/*
 * xfft_generic.i --
 *
 * Yorick eXtended FFT using implementation @IMPL@.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2009 Éric Thiébaut <thiebaut@obs.univ-lyon1.fr>
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can use, modify
 * and/or redistribute the software under the terms of the CeCILL-C license
 * as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated with
 * loading, using, modifying and/or developing or reproducing the software by
 * the user in light of its specific status of free software, that may mean
 * that it is complicated to manipulate, and that also therefore means that it
 * is reserved for developers and experienced professionals having in-depth
 * computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling
 * the security of their systems and/or data to be ensured and, more
 * generally, to use and operate it in the same conditions as regards
 * security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 *
 *-----------------------------------------------------------------------------
 *
 * $Id$
 * $Log$
 */

if (is_func(plug_in) && is_func(_xfft_@impl@_new) != 2) {
  plug_in, "xfft_@impl@";
  _xfft_@impl@_new = xfft_new;
  _xfft_@impl@_version = xfft_version;
  _xfft_@impl@_threads = xfft_threads;
  _xfft_@impl@_implementation = xfft_implementation;
  _xfft_@impl@_best_dim = xfft_best_dim;
  _xfft_@impl@_indgen = xfft_indgen;
  _xfft_@impl@_export_wisdom = xfft_export_wisdom;
  _xfft_@impl@_import_wisdom = xfft_import_wisdom;
  _xfft_@impl@_forget_wisdom = xfft_forget_wisdom;
} else {
  xfft_new = _xfft_@impl@_new;
  xfft_version = _xfft_@impl@_version;
  xfft_threads = _xfft_@impl@_threads;
  xfft_implementation = _xfft_@impl@_implementation;
  xfft_best_dim = _xfft_@impl@_best_dim;
  xfft_indgen = _xfft_@impl@_indgen;
  xfft_export_wisdom = _xfft_@impl@_export_wisdom;
  xfft_import_wisdom = _xfft_@impl@_import_wisdom;
  xfft_forget_wisdom = _xfft_@impl@_forget_wisdom;
}  

///* The following declarations are hooks for Yorick's codger to specifically
//   link with built-in functions provided by the plug-in. */
//extern xfft_new;
//extern xfft_version;
//extern xfft_threads;
//extern xfft_implementation;
//extern xfft_export_wisdom;
//extern xfft_import_wisdom;
//extern xfft_forget_wisdom;

/* Make sure common variables, interpreted functions and documentation are
   defined. */
if (is_void(XFFT_DIRECT)) {
  include, "xfft.i", 1;
}

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 78
 * coding: utf-8
 * End:
 */
