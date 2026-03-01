cppobj = main.cpp mpilib.cpp mpilib.h  mytime.cpp mytime.h
obj = main.o mpilib.o
comp = mpicxx
exec = a.out ga.out gp.out gpO3.out tsan
mpicxx = -isystem /usr/lib/x86_64-linux-gnu/mpich/include/ -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

wo3 = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
testfold =matrixtests/
testexe = ./a.out



a.out: $(cppobj)

	$(comp) $(mpicxx) $(cppobj) -o a.out


mainzip:
	zip Frolov_PS.zip $(cppobj) makefile 

zip:
	zip mpig.zip $(cppobj) makefile  


clean:
	rm -fr $(exec)

git:
	../sequentialGauss/commit.sh
	
	

push:
	git push

test:
	./test.sh  $(testexe) $(testfold)

all: $(exec)
