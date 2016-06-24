
#ARCH    = i386-apple-darwin8.10.0
ARCH     = powerpc-apple-darwin8.10.0
SYSLIB   = -L/usr/X11R6/lib -L/sw/lib -lreadline -lncurses -lm

#ARCH    = x86_64-unknown-linux-gnu
#ARCH    = i686-pc-linux-gnu
#SYSLIB  = -L/usr/lib64 -lreadline -lncurses -lm

#HOME     = /usr/people/zeyun
FETK_INC = -I${HOME}/include
FETK_LIB = -L${HOME}/lib/${ARCH} -ltriangle -lmaloc

default : go

tricall.o : tricall.c
	g++ -c -g -O2 -I. ${FETK_INC} tricall.c

go : tricall.o
	g++ -g -O2 -o go tricall.o ${FETK_LIB} ${SYSLIB}

clean :
	-rm -f *.o go

