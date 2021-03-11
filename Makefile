S = SOURCE/
O = OBJ/

PROG = gammcor

OBJ = $(O)mainp.o $(O)initia.o $(O)dmscf.o $(O)misc.o $(O)optocc.o \
      $(O)eigenno.o $(O)optnorb.o\
      $(O)matvec.o $(O)hf.o $(O)matrix.o $(O)cpdmft.o \
      $(O)nonadia.o $(O)dftgrid.o $(O)lsd_sr.o $(O)dftfun_exerfpbe.o \
      $(O)dftfun_exerf.o $(O)dftfun_ecerfpbe.o $(O)dftfun_ecerf.o \
      $(O)dftacg_pw92c.o $(O)projector.o $(O)ekt.o \
      $(O)xcfun.o \
      $(O)gridmolpro.o \
      $(O)sorter.o $(O)tran.o $(O)systemdef.o \
      $(O)types.o $(O)inputfill.o $(O)abmats.o $(O)abfofo.o \
      $(O)sapt_main.o $(O)sapt_exch.o $(O)sapt_pol.o $(O)sapt_utils.o \
      $(O)exmisc.o $(O)exdpino.o $(O)exi.o $(O)exappr.o \
      $(O)srefex.o $(O)diis.o \
      $(O)timing.o \
      $(O)srlrdynamic.o $(O)erpa.o $(O)interpa.o  $(O)exact2el.o $(O)optapsg.o $(O)newton.o $(O)acfd.o $(O)accas.o \
      $(O)caspidft.o $(O)ac_exact_2el.o $(O)vv10.o \
      xcfun_intel/fortran/xcfun_module.o xcfun_intel/fortran/xcfun_autogen.o

FCC = ifort -assume byterecl
FFLAGS = -mkl -heap-arrays  -O3 -I xcfun_intel/fortran
LIBS = - -L./xcfun_intel/lib -lxcfun 

$(PROG) :  $(OBJ) 
	$(FCC) $(FFLAGS) -o $(PROG) $(OBJ) $(LIBS)
$(O)mainp.o : $(S)mainp.f $(S)commons.inc $(O)types.o $(O)inputfill.o $(O)systemdef.o $(O)sapt_main.o
	$(FCC) $(FFLAGS)  -c $(S)mainp.f -o $(O)mainp.o
$(O)initia.o : $(S)initia.f $(S)commons.inc $(O)types.o $(O)sorter.o $(O)tran.o $(O)abmats.o 
	$(FCC) $(FFLAGS)  -c $(S)initia.f -o $(O)initia.o 
$(O)misc.o : $(S)misc.f $(S)commons.inc
	$(FCC) $(FFLAGS)  -c $(S)misc.f -o $(O)misc.o
$(O)dmscf.o : $(S)dmscf.f $(S)commons.inc $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)dmscf.f -o $(O)dmscf.o
$(O)optocc.o : $(S)optocc.f $(S)commons.inc
	$(FCC) $(FFLAGS)  -c $(S)optocc.f -o $(O)optocc.o
$(O)eigenno.o : $(S)eigenno.f $(S)commons.inc
	$(FCC) $(FFLAGS)  -c $(S)eigenno.f -o $(O)eigenno.o
$(O)optnorb.o : $(S)optnorb.f $(S)commons.inc
	$(FCC) $(FFLAGS)  -c $(S)optnorb.f -o $(O)optnorb.o
$(O)matvec.o : $(S)matvec.f
	$(FCC) $(FFLAGS)  -c $(S)matvec.f -o $(O)matvec.o
$(O)hf.o : $(S)hf.f
	$(FCC) $(FFLAGS)  -c $(S)hf.f -o $(O)hf.o
$(O)matrix.o : $(S)matrix.f
	$(FCC) $(FFLAGS)  -c $(S)matrix.f -o $(O)matrix.o
$(O)cpdmft.o : $(S)cpdmft.f
	$(FCC) $(FFLAGS)  -c $(S)cpdmft.f -o $(O)cpdmft.o
$(O)nonadia.o : $(S)nonadia.f
	$(FCC) $(FFLAGS)  -c $(S)nonadia.f -o $(O)nonadia.o
$(O)dftgrid.o : $(S)dftgrid.f $(S)commons.inc $(O)abmats.o
	$(FCC) $(FFLAGS)  -c $(S)dftgrid.f -o $(O)dftgrid.o
$(O)lsd_sr.o : $(S)lsd_sr.f
	$(FCC) $(FFLAGS)  -c $(S)lsd_sr.f  -o $(O)lsd_sr.o
$(O)dftfun_exerfpbe.o : $(S)dftfun_exerfpbe.f
	$(FCC) $(FFLAGS)  -c $(S)dftfun_exerfpbe.f  -o $(O)dftfun_exerfpbe.o
$(O)dftfun_exerf.o : $(S)dftfun_exerf.f
	$(FCC) $(FFLAGS)  -c $(S)dftfun_exerf.f  -o $(O)dftfun_exerf.o
$(O)dftfun_ecerfpbe.o : $(S)dftfun_ecerfpbe.f
	$(FCC) $(FFLAGS)  -c $(S)dftfun_ecerfpbe.f  -o $(O)dftfun_ecerfpbe.o
$(O)dftfun_ecerf.o : $(S)dftfun_ecerf.f
	$(FCC) $(FFLAGS)  -c $(S)dftfun_ecerf.f  -o $(O)dftfun_ecerf.o
$(O)dftacg_pw92c.o : $(S)dftacg_pw92c.f
	$(FCC) $(FFLAGS)  -c $(S)dftacg_pw92c.f  -o $(O)dftacg_pw92c.o
$(O)projector.o : $(S)projector.f
	$(FCC) $(FFLAGS)  -c $(S)projector.f  -o $(O)projector.o
$(O)ekt.o : $(S)ekt.f
	$(FCC) $(FFLAGS)  -c $(S)ekt.f  -o $(O)ekt.o
$(O)srlrdynamic.o : $(S)srlrdynamic.f
	$(FCC) $(FFLAGS)  -c $(S)srlrdynamic.f -o $(O)srlrdynamic.o
$(O)erpa.o : $(S)erpa.f
	$(FCC) $(FFLAGS)  -c $(S)erpa.f -o $(O)erpa.o
$(O)interpa.o : $(S)interpa.f $(O)abmats.o $(O)abfofo.o
	$(FCC) $(FFLAGS)  -c $(S)interpa.f -o $(O)interpa.o
$(O)exact2el.o : $(S)exact2el.f
	$(FCC) $(FFLAGS)  -c $(S)exact2el.f -o $(O)exact2el.o
$(O)optapsg.o : $(S)optapsg.f
	$(FCC) $(FFLAGS)  -c $(S)optapsg.f -o $(O)optapsg.o
$(O)newton.o : $(S)newton.f
	$(FCC) $(FFLAGS)  -c $(S)newton.f -o $(O)newton.o
$(O)acfd.o : $(S)acfd.f $(O)abmats.o $(O)abfofo.o
	$(FCC) $(FFLAGS)  -c $(S)acfd.f -o $(O)acfd.o
$(O)accas.o : $(S)accas.f
	$(FCC) $(FFLAGS)  -c $(S)accas.f -o $(O)accas.o
$(O)xcfun.o : $(S)xcfun.f90
	$(FCC) $(FFLAGS)  -c $(S)xcfun.f90  -o $(O)xcfun.o
$(O)gridmolpro.o : $(S)gridmolpro.f90
	$(FCC) $(FFLAGS)  -c $(S)gridmolpro.f90 -o $(O)gridmolpro.o
$(O)timing.o : $(S)timing.f90
	$(FCC) $(FFLAGS)  -c $(S)timing.f90 -o $(O)timing.o
$(O)types.o : $(S)types.f90
	$(FCC) $(FFLAGS)  -c $(S)types.f90 -o $(O)types.o
$(O)inputfill.o : $(S)inputfill.f90 $(O)types.o  
	$(FCC) $(FFLAGS)  -c $(S)inputfill.f90 -o $(O)inputfill.o
$(O)systemdef.o : $(S)systemdef.f90 $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)systemdef.f90 -o $(O)systemdef.o
$(O)sorter.o : $(S)sorter.f90 $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)sorter.f90 -o $(O)sorter.o
$(O)tran.o : $(S)tran.f90 $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)tran.f90 -o $(O)tran.o
$(O)diis.o : $(S)diis.f90 $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)diis.f90 -o $(O)diis.o
$(O)abmats.o : $(S)abmats.f90 $(O)types.o $(O)tran.o
	$(FCC) $(FFLAGS)  -c $(S)abmats.f90 -o $(O)abmats.o
$(O)abfofo.o : $(S)abfofo.f90 $(O)types.o $(O)tran.o
	$(FCC) $(FFLAGS)  -c $(S)abfofo.f90 -o $(O)abfofo.o
$(O)srefex.o : $(S)srefex.f90 $(O)types.o $(O)exmisc.o
	$(FCC) $(FFLAGS)  -c $(S)srefex.f90 -o $(O)srefex.o
$(O)exi.o : $(S)exi.f90
	$(FCC) $(FFLAGS)  -c $(S)exi.f90 -o $(O)exi.o
$(O)exmisc.o : $(S)exmisc.f90 $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)exmisc.f90 -o $(O)exmisc.o
$(O)exappr.o : $(S)exappr.f90 $(O)exmisc.o
	$(FCC) $(FFLAGS)  -c $(S)exappr.f90 -o $(O)exappr.o
$(O)exdpino.o : $(S)exdpino.f90 $(O)types.o $(O)tran.o $(O)timing.o $(O)exmisc.o
	$(FCC) $(FFLAGS)  -c $(S)exdpino.f90 -o $(O)exdpino.o
$(O)sapt_main.o : $(S)sapt_main.f90 $(O)types.o $(O)systemdef.o $(O)tran.o $(O)sorter.o $(O)sapt_exch.o $(O)sapt_pol.o $(O)sapt_utils.o $(O)abmats.o $(O)abfofo.o $(O)exdpino.o
	$(FCC) $(FFLAGS)  -c $(S)sapt_main.f90 -o $(O)sapt_main.o
$(O)sapt_utils.o : $(S)sapt_utils.f90 $(O)types.o $(O)tran.o $(O)diis.o
	$(FCC) $(FFLAGS)  -c $(S)sapt_utils.f90 -o $(O)sapt_utils.o
$(O)sapt_pol.o : $(S)sapt_pol.f90 $(O)types.o $(O)tran.o $(O)sapt_utils.o
	$(FCC) $(FFLAGS)  -c $(S)sapt_pol.f90 -o $(O)sapt_pol.o
$(O)sapt_exch.o : $(S)sapt_exch.f90 $(O)types.o $(O)tran.o $(O)timing.o $(O)sapt_utils.o $(O)exmisc.o $(O)exi.o $(O)exappr.o $(O)srefex.o
	$(FCC) $(FFLAGS)  -c $(S)sapt_exch.f90 -o $(O)sapt_exch.o
$(O)caspidft.o : $(S)caspidft.f $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)caspidft.f -o $(O)caspidft.o
$(O)vv10.o : $(S)vv10.f $(O)types.o
	$(FCC) $(FFLAGS)  -c $(S)vv10.f -o $(O)vv10.o
$(O)ac_exact_2el.o : $(S)ac_exact_2el.f
	$(FCC) $(FFLAGS)  -c $(S)ac_exact_2el.f -o $(O)ac_exact_2el.o


.PHONY : clean
clean :
	rm $(O)*.o *.mod $(PROG)

