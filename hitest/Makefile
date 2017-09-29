

all: class like_fourier.o
	gcc -Wno-missing-braces -Wno-missing-field-initializers -shared \
	    -o like_fourier.so -fPIC -fopenmp ../theory/like_fourier.c \
	    -L../hi_class -lfftw3 -lgsl -lgslcblas -lclass -lm
# Make sure that -lm is last, otherwise weird run-time linking errors occur!

test: class like_fourier.o
	gcc -o test -fPIC -fopenmp like_fourier.o \
	    -L../hi_class -lfftw3 -lgsl -lgslcblas -lclass -lm

like_fourier.o: ../theory/like_fourier.c ../theory/cosmo3D.c
	gcc -Wno-missing-braces -Wno-missing-field-initializers -g \
	    -fPIC -fopenmp -c ../theory/like_fourier.c

class:
	cd ../hi_class; $(MAKE)

clean:
	rm -f *.o *.so *.exe test

tim:
	gcc -Wno-missing-braces -Wno-missing-field-initializers -g \
	    -o test.exe -fPIC like_fourier.c \
	    -L../../hi_class -lfftw3 -lgsl -lgslcblas -lclass -lm	    
	gcc -Wno-missing-braces -Wno-missing-field-initializers -shared \
	    -o like_fourier.so -fPIC like_fourier.c \
	    -L../../hi_class -lfftw3 -lgsl -lgslcblas -lclass -lm
