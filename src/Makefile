LIBNAME=funfits
#LD=ld
#FFLAGS=-g

OBJS=  addreg.o    css.o       TPS.o       EXP.o     FEXP.o    GASP.o   \
       nkreg.o     nkden.o     spreg.o     qspreg.o

NNREG= bfgsfmin.o  dblas_le.o  linpakle.o  netev.o   nnutil.o  objfun.o

LLE=   dblas_le.o  dr_lle.o    liapexp.o   linpakle.o  llest.o

GCVPACK=  DNreps.o   ddcom.o    dmakek.o   dreps.o    dsuy.o     prmut.o \
DNtpss.o   ddevf.o    dmaket.o   drsap.o    dtpss.o    dzdc.o    dvmin.o\
dcfcr.o    ddiag.o    dpdcr.o    dset.o     dtrace.o   fact.o \
dcfr1.o    dftkf.o    dpmse.o    dsetup.o   dtsvdc.o   linpack.o \
dcpmut.o   dgcv.o     dpred.o    dsgdc.o    duni.o     mkkxys.o \
dcrtz.o    dgcv1.o    dprmut.o   dsgdc1.o   dvl.o      mkpoly.o \
dctsx.o    dgemv.o    dptpss.o   dsnsm.o    dvlop.o    mktpar.o




$(LIBNAME):	$(OBJS) lle nnreg tpsreg 
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS) # -lf2c -lm

clean:
	@rm -f *.o *.so nnreg.x nnregci.x

realclean:
	@rm -f Makefile *.o *.so

lle: $(LLE) 
	@$(F77) $(FFLAGS) $(LLE) -o ../exec/lle.x 

nnreg: $(NNREG) nnreg.f nnregci.f
	@$(F77) $(FFLAGS) nnreg.f -o ../exec/nnreg.x $(NNREG)
	@$(F77) $(FFLAGS) nnregci.f -o ../exec/nnregci.x $(NNREG)

tpsreg: $(GCVPACK)
	@$(F77) $(FFLAGS) ev.f rest.f $(GCVPACK) -o ../exec/tpsregev.x
	@$(F77) $(FFLAGS) tps.f rest.f $(GCVPACK) -o ../exec/tpsreg.x

