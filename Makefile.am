#
# Char0Gauss: Linear algebra stuff
#
# This file is part of the build system of a GAP kernel extension.
# Requires automake.
#
ACLOCAL_AMFLAGS = -I m4

BINARCHDIR = bin/$(GAPARCH)
GAPINSTALLLIB = $(BINARCHDIR)/Char0Gauss.so

lib_LTLIBRARIES = Char0Gauss.la

Char0Gauss_la_SOURCES = src/Char0Gauss.cc
Char0Gauss_la_CPPFLAGS = $(GAP_CPPFLAGS) -DCONFIG_H
Char0Gauss_la_LDFLAGS = -module -avoid-version
if SYS_IS_CYGWIN
Char0Gauss_la_LDFLAGS += -no-undefined -version-info 0:0:0 
Char0Gauss_la_LDFLAGS += -Wl,$(GAPROOT)/bin/$(GAPARCH)/gap.dll
endif

all-local: $(GAPINSTALLLIB)

$(GAPINSTALLLIB): Char0Gauss.la
	$(mkdir_p) $(BINARCHDIR)
if SYS_IS_CYGWIN
	cp .libs/Char0Gauss.dll $(GAPINSTALLLIB)
else
	cp .libs/Char0Gauss.so $(GAPINSTALLLIB)
endif

clean-local:
	rm -f $(GAPINSTALLLIB)

distclean-local:
	rm -rf bin/*

doc:
	($(GAPROOT)/bin/gap.sh -A makedoc.g)
