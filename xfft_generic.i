/*
 * xfft_generic.i --
 *
 * Yorick eXtended FFT using implementation @IMPL@.
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
