OBJS = pingpong.o himenoBMTxps.o k-means.o dataCompression.o

CC = mpicc
RM = /bin/rm
PROG = pingpong bmt k-means

CFLAGS = -O3 -lm
LDFLAGS = ${CFLAGS}

all : ${PROG}

pingpong : pingpong.o dataCompression.o 
	${CC} -o $@ pingpong.o dataCompression.o ${LDFLAGS}

bmt : himenoBMTxps.o dataCompression.o 
	${CC} -o $@ himenoBMTxps.o dataCompression.o ${LDFLAGS}	

k-means : k-means.o dataCompression.o 
	${CC} -o $@ k-means.o dataCompression.o ${LDFLAGS}

.c.o : param.h dataCompression.h
	${CC} -c ${CFLAGS} $<

clean :
	${RM} -f ${PROG} ${OBJS}






