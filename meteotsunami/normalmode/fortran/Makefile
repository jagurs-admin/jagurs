FC=ifort -DJAGURS #-DNONAN_NOFILT
FFLAGS=-cpp -check all -g -traceback -qmkl -I$(MKLROOT)/include/fftw -fp-model precise
FFLAGS=-cpp -O3                      -qmkl -I$(MKLROOT)/include/fftw -fp-model precise
LIBS=-L./amos -lamos

all: PHcalc_wave6 printp

PHcalc_wave6: mod_PHcalc_wave6.f90
	$(FC) -o $@ $(FFLAGS) $^ $(LIBS)

printp: printp.f90
	$(FC) -o $@ $(FFLAGS) $^

clean:
	rm -f PHcalc_wave6 printp *.mod
