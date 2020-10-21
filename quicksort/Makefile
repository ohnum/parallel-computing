CC         =  mpicc
CCFLAGS    =  -O3
LIBS       =  -lmpi -lm

p_qsort: p_qsort.c
	$(CC) $(CCFLAGS) -o p_qsort p_qsort.c $(LIBS)

s_qsort: Squicksort.c
	$(CC) $(CCFLAGS) -o s_qsort Squicksort.c $(LIBS)

clean:
	$(RM) p_qsort
