XFFT: eXtended support for the Fast Fourier Transform in Yorick
===============================================================

Description
-----------

**XFFT** (for *eXtended Fast Fourier Transform*) is a
[Yorick](http://yorick.github.com/) extension to uniformely interface
between various implementations of FFT (notably
[FFTW](http://www.fftw.org), the Fastest Fourier Transform in the West,
and, in the future, CFFT based on Swarztrauber routines) and benefit from
multi-threading, complex-to-real and real-to-complex transforms, etc.

The idea is to simplify the use of FFT transforms (XFFT provides direct,
inverse, complex-to-real and real-to-complex transforms) and write the
same code whatever the particular FFT implementation used at runtime.

Note that, for the moment, the ability of Yorick's FFT to compute the
forward or backward transform along some (not all) directions of the input
array are not supported by XFFT.


Installation
------------

In short, building and installing the plug-in can be as quick as:
````{.sh}
cd $BUILD_DIR
$SRC_DIR/configure
make
make install
````
where `$BUILD_DIR` is the build directory (at your convenience) and
`$SRC_DIR` is the source directory of the plug-in code.  The build and
source directories can be the same in which case, call `./configure` to
configure for building.

If the plug-in has been properly installed, it is sufficient to use any
of its functions (like `xfft_new`) to automatically load the plug-in.
You may force the loading of the plug-in by something like:
````{.cpp}
#include "xfft.i"
````
or
````{.cpp}
require, "xfft.i";
````
in your code.

More detailled installation explanations are given below.


0. You must have [Yorick](http://yorick.github.com/) installed on your
   machine and the various FFT libraries.

1. Unpack the plug-in code somewhere.

2. Configure for compilation.  There are two possibilities:

   * For an **in-place build**, go to the source directory of the plug-in 
     code and run the configuration script:
     ````{.sh}
     cd $SRC_DIR
     ./configure
     ````
     To see the configuration options, call:
     ````{.sh}
     ./configure --help
     ````
	 See below for more informations about the options.

   * To compile in a **different build directory**, say `$BUILD_DIR`, create
     the build directory, go to the build directory and run the
     configuration script:
     ````{.sh}
     mkdir -p $BUILD_DIR
     cd $BUILD_DIR
     $SRC_DIR/configure
     ````
     where `$SRC_DIR` is the path to the source directory of the plug-in code.
     To see the configuration options, call:
     ````{.sh}
     $SRC_DIR/configure --help
     ````
	 See below for more informations about the options.

3. Compile the code:
   ````{.sh}
   make clean
   make
   ````

4. Install the plug-in in Yorick directories:
   ````{.sh}
   make install
   ````


Configuration Options
---------------------

The configuration is the opportunity to specify which implementation(s) of
FFT you want to use and the corresponding compilation flags and libraries.
You may specify as many implementations as you want (this will result in
several different plugins) providing you have the corresponding libraries
(and header files).  The possibilities are:
```
    fftw2           - FFTW version 2
    fftw2-threads   - FFTW version 2 with multi-threading
    fftw3           - FFTW version 3
    fftw3-threads   - FFTW version 3 with multi-threading
```
In the future, the following front-ends will be implemented:
```
    cfft            - Swarztrauber FFT routines
    mkl             - FFTW from the Math kernel Library with multi-threading
    mkl-threads     - FFTW from the Math kernel Library with multi-threading
```
CFFT implementation uses Swarztrauber FFT routines already compiled in Yorick
so this one should be always available.

For each specific implementation, say `IMPL`, the configuration can be done
via options of the `configure` script:

* `--with-IMPL` or `--with-IMPL=yes` to compile support for `IMPL`.

* `--without-IMPL` or `--with-IMPL=no` to not compile support for `IMPL`.

* `--IMPL-cflags=...` for additional compiler (or preprocessor) flags.

* `--IMPL-ldflags=...` for additional linker flags.

* `--IMPL-deplibs=...` for additional libraries.

For instance, the default configuration settings correspond to:
```sh
    configure \
      --with-fftw2=yes \
      --fftw2-cflags="" \
      --fftw2-ldflags="" \
      --fftw2-deplibs="-lrfftw -lfftw" \
      --with-fftw2-threads=yes \
      --fftw2-threads-cflags="" \
	  --fftw2-threads-ldflags="" \
      --fftw2-threads-deplibs="-lrfftw_threads -lfftw_threads -lrfftw -lfftw -lpthread" \
      --with-fftw3=yes \
      --fftw3-cflags="" \
      --fftw3-ldflags="" \
      --fftw3-deplibs="-lfftw3" \
      --with-fftw3-threads=yes \
      --fftw3-threads-cflags="" \
      --fftw3-threads-ldflags="" \
      --fftw3-threads-deplibs="-lfftw3_threads -lfftw3 -lpthread"
```


Implementation
--------------

FFT transforms are implemented by XFFT in Yorick as special objects which
behave as functions to compute the forward or backward transforms and which
take care themself of all the housekeeping with (de-)allocation of required
FFTW plans and workspaces.  These ressources are only allocated once and only
when required.  That is, assuming you have compiled XFFT to use FFTW to
compute the transforms, if you create an XFFT object and only compute forward
FFT with it, it will never create a FFTW plan for backward transform.  The
drawback of this hiding is that an XFFT object can only be used for a given
dimension list of the arrays to transform, similarly XFFT objects for
real-complex and complex-complex transforms must be different (even if they
apply to arrays of same dimensions).  However the same XFFT object can be used
for forward and backward transform.  For instance:
```
    f = xfft_new();          // create a new XFFT object with default options
    z = f(a);                // compute the forward transform of A
    z = f(a, XFFT_FORWARD);  // idem
    b = f(z, XFFT_BACKWARD); // compute the backward transform of Z
    f, a, job;               // in-place transform, optional job is
                             // XFFT_FORWARD, or XFFT_BACKWARD, etc.
    f = [];                  // destroys the XFFT object (and release
                             // ressources)
```
The main difference between FFTW2 and FFTW3 is that FFTW3 makes use of
workspaces that are properly aligned to boost the performances (e.g., by
using
[SIMD]'https://fr.wikipedia.org/wiki/Single_instruction_multiple_data)
instructions).  Each FFTW3 plan is aware of its input and output workspaces
(which can be the same for an in-place transform).  This is why each FFTW
object owns its workspaces (though they are the same for the forward and
backward plans of complex transforms).  In the interfacing with Yorick,
input array must be copied into the input workspace, then the FFT is
computed, then the output workspace is copied into the output Yorick array.
These copy operations have a slight impact on the performances but
destroying the contents of the input workspace in FFTW transform is not an
issue.  For 1-D transforms, the out-of-place transform is slighlty faster
than the in-place one.  For N-D transforms (N > 1), the performances are
very similar or better for the in-place transforms.  To save memory without
sacrificing performances, FFTW objects use out-of-place transform for 1-D
data and in-place transform otherwise.

C2R destroys its input even for out-of-place transform (unless in 1D and if
flag `FFTW_PRESERVE_INPUT` is set, but this yields to some sacrifice of
performance, therefore I always consider that input is not preserved).

The same interface is implemented for FFTW2 and FFTW3 (and will perhaps be
implemented for Swarztrauber's FFT which is built into Yorick).  However,
because some symbols have the same names in the two libraries, the two plugins
cannot be loaded in the same Yorick session.  The recommended way to load FFTW
plugin is:
```C
   #include "xfft.i"
```
which will load FFTW3 if possible and FFTW2 otherwise, *etc.*  Thanks to the
`autoload` comman of Yorick, this should be automatically done as soon as
you use one of the XFFT functions.  To explicitly use FFTW2:
```C
   #include "xfft_fftw2.i"
```
and to explicitly use FFTW3 with threads:
```C
   #include "xfft_fftw3_threads.i"
```

Example of use:
```C
   f = xfft_new();
   input = array(double, f.dims);
   if (f.real) {
      // Real-complex transform.
      dims = f.dims;
      dims(2) = dims(2)/2 + 1;
   } else {
      // Complex-complex transform.
      dims = f.dims;
   }
   output = array(complex, dims);
```


Dimension List
--------------

> Following Yorick conventions, dimensions are given here in column major
> format, the first array index varies most rapidly; however beware that
> FFTW documentation is for row major format.

For complex-complex transforms, input and output arrays have the same
dimensions whether the transform is done in-place or not.  Input and
output workspaces can be the same that is array(s) of `N1 x N2 x ... x Nd`
complex values.

For real-complex transforms, the real array is `N1 x N2 x ... x Nd` whereas
the complex array is Hermitian and of size `(N1/2 + 1) x N2 x ... x Nd`.  For
out-of-place transforms, the workspace for the real array must be at least
`N1 x N2 x ... x Nd` real data, the workspace for the complex array must be
at least `(N1/2 + 1) x N2 x ... x Nd` complex data.  For in-place transforms,
the single workspace must be large enough to store the real or the complex
array, that is: `2*(N1/2 + 1) x N2 x ... x Nd` real data.


Choosing a Specific Implementation
----------------------------------

If you want to use a specific FFT implementation, directly include the
corresponding library file, for instance:
```C
    #include "xfft_fftw2_threads.i"
```
to use FFTW (version 2) with multi-threading.  By default, when you include
"xfft.i" the first time, it tries to load a multi-threaded version and by
prefrence FFTW (version 3).  Thus the order of preference is:
```C
    mkl_threads     - FFTW from the Math kernel Library with multi-threading
    fftw3_threads   - FFTW version 3 with multi-threading
    fftw2_threads   - FFTW version 2 with multi-threading
    mkl             - FFTW from the Math kernel Library with multi-threading
    fftw3           - FFTW version 3
    fftw2           - FFTW version 2
    cfft            - Swarztrauber FFT routines
```
Some provisions have been made to allow for different implementations of the
FFT to co-exist in a Yorick session but this is not always possible.  For
instance, FFTW 2 and 3 cannot be safely loaded at the same time.

If multiple implementations are loaded, `fftw_new` will create an FFT transform
object which uses the last implementation that have nee loaded.  To choose a
specific implementation, then use:
```C
    op = _xfft_IMPL_new(...)
```
with `IMPL` the name of the implementation (e.g., `fftw2_threads`) to create an
object which uses FFTW (version 2) with multi-threading.
