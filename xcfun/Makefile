# Compile time options:
# -DNDEBUG Turn off run time checks in production calculations to gain speed
#NOTE: Using -fast gives wrong output from pbex, disappears at -O2
#level (without ipo) for icpc 10.1.017

# C++ compiler and flags
CXX = icx
CC = icx
# Turn of messages that hide real errors
#CFLAGS = -Wall -g -O3 -wd981 -wd279 -wd383 -wd1572 -wd177 -DNDEBUG #-vec-report0
CFLAGS = -Wall -g -O3  -DNDEBUG #-vec-report0
# icx takes the C++ headers of the newest GCC in /usr/lib/gcc, which may lack
# libstdc++-dev; use the GCC of g++ instead (empty: icx default)
GCC_INSTALL_DIR?=$(shell g++ -print-libgcc-file-name 2>/dev/null | xargs -r dirname)
ifneq ($(GCC_INSTALL_DIR),)
CFLAGS+=--gcc-install-dir=$(GCC_INSTALL_DIR)
endif
CXXFLAGS = $(CFLAGS) -fno-rtti -fno-exceptions

# Fortran compiler and flags
FC = ifx

intel : lib \
 		fortran/xcfun_autogen.o \
		fortran/xcfun_module.o

include Makefile.common
