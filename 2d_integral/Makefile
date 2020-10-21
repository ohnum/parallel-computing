CC         =  mpicc
CCFLAGS    =  -O3 -Wall
LIBS       =  -lmpi -lm

integral2d:             integral2d.c
	$(CC) $(CCFLAGS) -o integral2d integral2d.c $(LIBS)

clean:
	$(RM) integral2d
