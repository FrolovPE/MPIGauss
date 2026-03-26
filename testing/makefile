cppobj = main.cpp mpilib.cpp mpilib.h  mytime.cpp mytime.h
obj = main.o mpilib.o
comp = mpicxx
exec = a.out ga.out san.out
flags = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
mpicxx = -isystem /usr/lib/x86_64-linux-gnu/mpich/include/ -O3 $(flags)
gmpicxx = -isystem /usr/lib/x86_64-linux-gnu/mpich/include/  $(flags)
sanmpicxx = -isystem /usr/lib/x86_64-linux-gnu/mpich/include/ -O0 -g -fsanitize=address -fno-omit-frame-pointer  $(flags)
wo3 = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
testfold =matrixtests/



a.out: $(cppobj)

	$(comp) $(mpicxx) $(cppobj) -o $@

ga.out:$(cppobj)

	$(comp) $(gmpicxx) $(cppobj) -o $@

san.out:$(cppobj)

	$(comp) $(sanmpicxx) $(cppobj) -o $@

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
