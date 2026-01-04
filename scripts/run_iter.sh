#!/bin/bash

BIN=./bin/tpPoisson1D_iter
N=100          # problem size (change to 1000 if needed)
OUTDIR=results_iter

mkdir -p ${OUTDIR}

echo "Running iterative methods for N = ${N}"
echo "-------------------------------------"

# Richardson (optimal alpha)
echo ">> Richardson (ALPHA)"
${BIN} 0 ${N} > ${OUTDIR}/richardson.log
mv RESVEC.dat ${OUTDIR}/RESVEC_richardson.dat

# Jacobi
echo ">> Jacobi"
${BIN} 1 ${N} > ${OUTDIR}/jacobi.log
mv RESVEC.dat ${OUTDIR}/RESVEC_jacobi.dat

# Gauss-Seidel
echo ">> Gauss-Seidel"
${BIN} 2 ${N} > ${OUTDIR}/gs.log
mv RESVEC.dat ${OUTDIR}/RESVEC_gs.dat

echo "-------------------------------------"
echo "All runs completed."
echo "Results saved in ${OUTDIR}/"
