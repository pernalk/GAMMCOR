# Compile time options:
# -DNDEBUG Turn off run time checks in production calculations to gain speed
#NOTE: Using -fast gives wrong output from pbex, disappears at -O2
#level (without ipo) for icpc 10.1.017

# C++ compiler and flags
CXX = icx
CC = icx
# Turn of messages that hide real errors
CFLAGS = -Wall -g -O3 -wd981 -wd279 -wd383 -wd1572 -wd177 -DNDEBUG #-vec-report0 
CXXFLAGS = $(CFLAGS) -fno-rtti -fno-exceptions

# Fortran compiler and flags
FC = ifort

intel : lib \
 		fortran/xcfun_autogen.o \
		fortran/xcfun_module.o

include Makefile.common
