
include config/make.def

OBJS = dataCompression.o himenoBMTxps_comp.o 
#OBJS = himenoBMTxps.o

CFLAGS = -lm 
#CFLAGS = -lm -g 

RM = /bin/rm
PROG = himeno
PROGDIR=bin

${PROGDIR}/${PROG}: all

${PROG} : all

all: ${OBJS} param.h
	${CLINK} ${CLINKFLAGS} -o ${PROGDIR}/${PROG} ${OBJS} ${CMPI_LIB}${CFLAGS}

#  ${CLINK} ${CLINKFLAGS} -o ${PROGDIR}/${PROG} ${OBJS} ${CMPI_LIB}

.c.o : param.h
	 ${MPICC}  ${CFLAGS} ${CMPI_INC} -c ${CFLAGS} $<

clean :
	${RM} -f ${PROGDIR}/${PROG} ${OBJS} core log*.txt param.h
	
