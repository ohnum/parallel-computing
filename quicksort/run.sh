DIR="./";
CORE="2";
LEN=300000000;

make clean
make
(mpirun -np ${CORE} ${DIR}p_qsort ${LEN})
