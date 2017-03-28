FFLAGS = -O3 -Wall -fPIC
FC = gfortran 


objects = double.o pdf_calculator.o

pdf_calculator.o : pdf_calculator.f90 double.o
		$(FC) -c $(FFLAGS) $<

double.o : double.f90
	$(FC) -c $(FFLAGS) $<

python: pdf_calculator.o double.o pdf_calculator.f90
		f90wrap -vvvvv -m pdf_calculator_funcs pdf_calculator.f90 -k kind_map_ajm -P
		f2py-f90wrap -c -m _pdf_calculator_funcs *.o f90wrap*.f90

clean:
	rm -f *.o *.mod *.MOD

pyclean:
	rm -f *.o *.mod *.MOD *.so *.pyc */*.pyc
	rm -f f90wrap* _pdf_calculator_funcs
	rm -rf pdf_calculator_funcs
