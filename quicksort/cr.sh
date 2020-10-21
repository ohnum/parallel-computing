DIR="./";
CORE="32";
ITER=4;  #no of times to run simulation and time it
OUT=./time_32P.txt
LEN=1000000000;

printf "_________________________________________________________\n" >> ${OUT}
make -C ${DIR} >> ${OUT}
printf "CORE=${CORE}, LEN=${LEN}\n" >> ${OUT}
for (( num=1; num <= ${ITER}; num++ ))
do
  (mpirun -np ${CORE} ${DIR}p_qsort ${LEN}) | tee -a ${OUT}
done

make clean -C ${DIR}
less ${OUT}
