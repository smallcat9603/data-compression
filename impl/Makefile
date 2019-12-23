OBJS = himenoBMTxps.o dataCompression.o

CC = mpicc
RM = /bin/rm
PROG = bmt

CFLAGS = -O3
LDFLAGS = ${CFLAGS}

all : ${PROG}

${PROG} : ${OBJS} param.h dataCompression.h
	${CC} -o $@ ${OBJS} ${LDFLAGS}

.c.o : param.h dataCompression.h
	${CC} -c ${CFLAGS} $<

clean :
	${RM} -f ${PROG} ${OBJS} core






