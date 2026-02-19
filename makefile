obj = main.cpp mpilib.cpp mpilib.h  mytime.cpp mytime.h
comp = mpicxx
exec = a.out ga.out gp.out gpO3.out tsan
mpicxx = -isystem /opt/impi-5.1.3.223/intel64/include -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
parallel = -pthread  $(cxx)
tsan =  $(parallel) -fsanitize=thread
gprO3 = -pg $(cxx)
gpr = -pg $(wo3)
wo3 = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
testfold =matrixtests/
testexe = ./a.out



a.out: $(obj)

	$(comp) $(cxx) $(obj) -o a.out

tsan: $(obj)

	$(comp) $(tsan) $(obj) -o $@

mainzip:
	zip Frolov_PS.zip $(obj) makefile 

zip:
	zip mpig.zip $(obj) makefile  


gp.out: $(obj)

	$(comp) $(gpr) $(obj) -o gp.out

gpO3.out: $(obj)

	$(comp) $(gprO3) $(obj) -o gpO3.out

ga.out: $(obj)

	$(comp) $(wo3) $(obj) -o ga.out

clean:
	rm -fr $(exec)

git:
	../sequentialGauss/commit.sh
	
	

push:
	git push

test:
	./test.sh  $(testexe) $(testfold)

all: $(exec)
