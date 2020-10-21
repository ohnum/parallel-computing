CC         =  mpicc
CCFLAGS    =  -O3
LIBS       =  -lmpi -lm

mm_mul: mm_mul.c
	$(CC) $(CCFLAGS) -o mm_mul mm_mul.c $(LIBS)

clean:
	$(RM) mm_mul
