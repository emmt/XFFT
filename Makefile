#
# Makefile --
#
# Make rules for XFFT plug-in.
#
# -----------------------------------------------------------------------------

# Current version (a final 'x' indicates a development version)
VERSION = 0.2.0

#---- BEGIN OF CUSTOMIZABLE PART ----------------------------------------------
# Edit this part for customization.

# Choose which package to build (fftw2 and/or fftw3 with/without threads):
TARGETS = fftw2 fftw2_threads fftw3 fftw3_threads

# Set compiler and (rarely) loader flags and libraries for FFTW2 (you can add
# -DFFTW_PREFIX=1 to FFTW2_CFLAGS if you have compiled FFTW with option
# --enable-prefix):
FFTW2_CFLAGS  =
FFTW2_LDFLAGS =
FFTW2_DEPLIBS = -lrfftw -lfftw

# Similarly for FFTW2 with threads.
FFTW2_THREADS_CFLAGS =
FFTW2_THREADS_LDFLAGS =
FFTW2_THREADS_DEPLIBS = -lrfftw_threads -lfftw_threads -lrfftw -lfftw -lpthread

# Similarly for FFTW3.
FFTW3_CFLAGS  =
FFTW3_LDFLAGS =
FFTW3_DEPLIBS = -lfftw3

# Similarly for FFTW3 with threads.
FFTW3_THREADS_CFLAGS  =
FFTW3_THREADS_LDFLAGS =
FFTW3_THREADS_DEPLIBS = -lfftw3_threads -lfftw3 -lpthread

# Similarly for CFFT:
CFFT_CFLAGS  =
CFFT_LDFLAGS =
CFFT_DEPLIBS =

#---- END OF CUSTOMIZABLE PART ------------------------------------------------

# where are the sources? (automatically filled in by configure script)
srcdir=.

# these values filled in by yorick -batch make.i
Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

all: $(TARGETS)

help:
	@echo "The targets are:"
	@echo "    make all                    # to build the software"
	@echo "    make clean                  # to remove some files"
#	@echo "    make distclean              # to remove non-distribution files"
	@echo "    make distrib [VERSION=###]  # to pack a distribution"
	@echo "    make install                # to install the software"

# These common sed rules must be the last:
SED_COMMON = \
  -e 's,@VERSION@,$(VERSION),g' \
  -e 's,@srcdir@,$(srcdir),g' \
  -e 's,@Y_MAKEDIR@,$(Y_MAKEDIR),g' \
  -e 's,@Y_EXE@,$(Y_EXE),g' \
  -e 's,@Y_EXE_PKGS@,$(Y_EXE_PKGS),g' \
  -e 's,@Y_EXE_HOME@,$(Y_EXE_HOME),g' \
  -e 's,@Y_EXE_SITE@,$(Y_EXE_SITE),g' \
  -e 's,@Y_HOME_PKG@,$(Y_HOME_PKG),g' \
  -e 's,%\(impl\|IMPL\)%,@\1@,g'

SED_CFFT = \
  -e 's,@IMPL@,CFFT,g' \
  -e 's,@impl@,cfft,g' \
  -e 's,xfft_generic\.mkf,xfft_cfft.mkf,g' \
  -e 's,@PKG_DEPLIBS@,$(CFFT_DEPLIBS),g' \
  -e 's,@PKG_CFLAGS@,$(CFFT_CFLAGS) -DCFFT,g' \
  -e 's,@PKG_LDFLAGS@,$(CFFT_LDFLAGS),g' \

SED_FFTW2 = \
  -e 's,@IMPL@,FFTW2,g' \
  -e 's,@impl@,fftw2,g' \
  -e 's,xfft_generic\.mkf,xfft_fftw2.mkf,g' \
  -e 's,@PKG_DEPLIBS@,$(FFTW2_DEPLIBS),g' \
  -e 's,@PKG_CFLAGS@,$(FFTW2_CFLAGS) -DFFTW2,g' \
  -e 's,@PKG_LDFLAGS@,$(FFTW2_LDFLAGS),g'

SED_FFTW2_THREADS = \
  -e 's,@IMPL@,FFTW2_THREADS,g' \
  -e 's,@impl@,fftw2_threads,g' \
  -e 's,xfft_generic\.mkf,xfft_fftw2_threads.mkf,g' \
  -e 's,@PKG_DEPLIBS@,$(FFTW2_THREADS_DEPLIBS),g' \
  -e 's,@PKG_CFLAGS@,$(FFTW2_THREADS_CFLAGS) -DFFTW2 -DUSE_THREADS,g' \
  -e 's,@PKG_LDFLAGS@,$(FFTW2_THREADS_LDFLAGS),g'

SED_FFTW3 = \
  -e 's,@IMPL@,FFTW3,g' \
  -e 's,@impl@,fftw3,g' \
  -e 's,xfft_generic\.mkf,xfft_fftw3.mkf,g' \
  -e 's,@PKG_DEPLIBS@,$(FFTW3_DEPLIBS),g' \
  -e 's,@PKG_CFLAGS@,$(FFTW3_CFLAGS) -DFFTW3,g' \
  -e 's,@PKG_LDFLAGS@,$(FFTW3_LDFLAGS),g'

SED_FFTW3_THREADS = \
  -e 's,@IMPL@,FFTW3_THREADS,g' \
  -e 's,@impl@,fftw3_threads,g' \
  -e 's,xfft_generic\.mkf,xfft_fftw3_threads.mkf,g' \
  -e 's,@PKG_DEPLIBS@,$(FFTW3_THREADS_DEPLIBS),g' \
  -e 's,@PKG_CFLAGS@,$(FFTW3_THREADS_CFLAGS) -DFFTW3 -DUSE_THREADS,g' \
  -e 's,@PKG_LDFLAGS@,$(FFTW3_THREADS_LDFLAGS),g'


# Rules for CFFT:

cfft: xfft_cfft

xfft_cfft: xfft_cfft.mkf
	make -f xfft_cfft.mkf

xfft_cfft.mkf: ${srcdir}/xfft_generic.mkf Makefile
	sed <"$<" >"$@" $(SED_CFFT) $(SED_COMMON)


# Rules for FFTW2:

fftw2: xfft_fftw2

xfft_fftw2: xfft_fftw2.mkf
	make -f xfft_fftw2.mkf

xfft_fftw2.mkf: ${srcdir}/xfft_generic.mkf Makefile
	sed <"$<" >"$@" $(SED_FFTW2) $(SED_COMMON)


# Rules for FFTW2 with threads:

fftw2_threads: xfft_fftw2_threads

xfft_fftw2_threads: xfft_fftw2_threads.mkf
	make -f xfft_fftw2_threads.mkf

xfft_fftw2_threads.mkf: ${srcdir}/xfft_generic.mkf Makefile
	sed <"$<" >"$@" $(SED_FFTW2_THREADS) $(SED_COMMON)


# Rules for FFTW3:

fftw3: xfft_fftw3

xfft_fftw3: xfft_fftw3.mkf
	make -f xfft_fftw3.mkf

xfft_fftw3.mkf: ${srcdir}/xfft_generic.mkf Makefile
	sed <"$<" >"$@" $(SED_FFTW3) $(SED_COMMON)


# Rules for FFTW3 with threads:

fftw3_threads: xfft_fftw3_threads

xfft_fftw3_threads: xfft_fftw3_threads.mkf
	make -f xfft_fftw3_threads.mkf

xfft_fftw3_threads.mkf: ${srcdir}/xfft_generic.mkf Makefile
	sed <"$<" >"$@" $(SED_FFTW3_THREADS) $(SED_COMMON)


# General rules:

clean:
	rm -f core *~
	for impl in $(TARGETS); do \
	  mkf=xfft_$$impl.mkf; \
	  if test -f $$mkf; then \
	    make -f $$mkf clean; \
	  fi; \
	done
	rm -f xfft_cfft.i xfft_cfft.mkf
	rm -f xfft_fftw*.i xfft_fftw*.mkf

install:
	for impl in $(TARGETS); do \
	  mkf=xfft_$$impl.mkf; \
	  if test -f $$mkf; then \
	    make -f $$mkf install; \
	  fi; \
	done

#distclean: clean-here
#	rm -rf fftw2 fftw3

DISTRIB_FILES = Makefile AUTHORS LICENSE.md README.md NEWS.md FAQ.txt configure \
                xfft.c xfft.i xfft-start.i xfft_generic.i xfft_generic.mkf

bump-version:
	@oldversion=`sed < Makefile -e '/^ *VERSION *=/!d;s/^ *VERSION *= *\(.*[^ ]\) *$$/\1/'`; \
	newversion=`echo $(VERSION) | sed -e 's/^ *//;s/ *$$//'`; \
	if test "$$oldversion" = "$$newversion"; then \
	  echo "No version change (VERSION=$$oldversion)."; \
	  echo "To change the version, try:"; \
	  echo "   make $@ VERSION=..."; \
	else \
	  sed <Makefile >Makefile.tmp -e "s/^ *VERSION *=.*$$/VERSION = $$newversion/"; \
	  mv -f Makefile.tmp Makefile; \
	  echo "New version is $$newversion (previous version was $$oldversion)."; \
	  echo "You may make a new release with:"; \
	  echo "    git commit -m 'New version v$$newversion' Makefile"; \
	  echo "    git tag 'v$$newversion'"; \
	  echo "    git push all"; \
	  echo "    git push --tags"; \
	  echo "    make distrib"; \
	fi; \
	return 0

distrib:
	@if test "x$(VERSION)" = "x"; then \
	  echo >&2 "error: bad VERSION"; \
	  return 1; \
	else \
	  version=$(VERSION); \
	fi; \
	pkgdir=xfft-$$version; \
	archive=$$pkgdir.tar.bz2; \
	if test -e "$$pkgdir"; then \
	  echo >&2 "error: $$pkgdir already exists, bump or override version number"; \
	  echo >&2 "error: for instance: make VERSION=XXXX $@"; \
	  return 1; \
	fi; \
	if test -e "$$archive"; then \
	  echo >&2 "error: $$archive already exists, bump or override version number"; \
	  echo >&2 "error: for instance: make VERSION=XXXX $@"; \
	  return 1; \
	fi; \
	mkdir "$$pkgdir"; \
	for file in $(DISTRIB_FILES); do \
	  cp -a "$$file" "$$pkgdir/."; \
	done; \
	rm -f "$$pkgdir/Makefile"; \
	sed <Makefile >"$$pkgdir/Makefile" -e 's/^ *VERSION *=.*$$/VERSION = $(VERSION)/'; \
	tar cf - "$$pkgdir" | bzip2 -9 > "$$archive"; \
	rm -rf "$$pkgdir"; \
	echo "archive $$archive created"; \
	return 0

.PHONY: clean distrib install \
        xfft_cfft \
        xfft_fftw2 xfft_fftw2_threads \
        xfft_fftw3 xfft_fftw3_threads

# ------------------------------------------------------------- end of Makefile
# Local Variables:
# mode: Makefile
# tab-width: 8
# fill-column: 78
# coding: utf-8
# End:
# -----------------------------------------------------------------------------
