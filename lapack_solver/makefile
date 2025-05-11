LAPACK_LIBDIR = -L/usr/local/opt/lapack/lib
LAPACK_INCDIR = -I/usr/local/opt/lapack/include
#was -I$(brew --prefix lapack)/include and -L$(brew --prefix lapack)/lib
#	I replaced (brew --prefix lapack) with the output but $(brew --prefix lapack) would have done it automatically.
#	in my solution, I should have erased the $ as well



TUTORIAL: lapack_solver.o
	gcc -Wall -Wextra -Wpedantic -Ofast -flto aufgabe35.o $(LAPACK_LIBDIR) -l lapacke -lm -o TUTORIAL

lapack_solver.o: lapack_solver.c ../../my_numerics.h
	gcc -Wall -Wextra -Wpedantic -I.. $(LAPACK_INCDIR) -c -flto -Ofast lapack_solver.c

run: TUTORIAL
	./TUTORIAL

clean:
	@rm -rf *.o
	@rm -f TUTORIAL
