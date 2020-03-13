# History of changes in XFFT

## Version 0.2.0 (released on 2020/03/13)
- real-real transforms of FFTW3 are implemented (Hartley transform, discrete
  cosine/sine transforms, ...);
- fix sed command for MacOS.

## Version 0.1.0 (released on 2015/06/19)
- add pseudo-configure script;
- plug-in can now be compiled from a remote build directory;
- add autoload file;
- change to GNU-GPL license;
- XFFT is now hosted by GitHub (https://github.com/emmt/XFFT).

## Version 0.0.4 (released on 2013/05/14)
- make "unsigned char" the default C type for Yorick "char" (unless
  macro YORICK_CHAR_IS_SIGNED is defined).

## Version 0.0.3 (released on 2012/09/18)
- fixed a bug due to improperly detecting scratch arrays;
- fixed trailing space in version number in Makefile;

## Version 0.0.2 (released on 2012/06/12)
- fixed installation of xfft.i;
- added preprocessor directives in xfft.c to avoid compilation warnings;
- fixed text in error messages.

## Version 0.0.1 (released on 2011/11/25)
