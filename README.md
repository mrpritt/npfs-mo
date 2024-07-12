# Non-permutation flow shop scheduling with missing operations

This repository contains detailed results and code for the paper [Non-permutation flow shop scheduling with missing operations](https://doi.org/10.1016/j.cor.2024.106742).

## Detailed results

The detailed results can be found in folder [data](data). In the folder you also can find a [R notebook](data/tables.Rmd) to reproduce the tables, including a [rendered version in HTML](https://raw.githack.com/mrpritt/npfs-mo/main/data/tables.html).

## Code

The code is contained in the folder `src`. To compile, clone the repo, and use
```bash
cd npfs-mo/build
git submodule update --init 
cmake -DCMAKE_BUILD_TYPE=Release ../src
make
```
You will need [Boost](https://www.boost.org), and if you want to build the exact solver, also a [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) installation at `$CPLEX_ROOT_DIR`.

To run the experiments from the paper, for example on instance 0.2_10_05_02.txt, do the following.
```bash
./npfsmo --flowtime --timelimit -1 --iterfactor 0.1 --npfs 0.2_10_05_02.txt
```
This will produce a couple of output lines. The next-to-last value in the line tagged INFO is the flowtime found by the IGA. By default all parameters are fixed to the settings of the paper, and the random seed is fixed to 1. Therefore, since the stopping criterion is the number of iterations, and not time, you should be able to exactly reproduce the values from the tables.

## How to cite
```bibtex
@Article{Ritt.Rossit/2024,
  author =  {Marcus Ritt and Daniel Alejandro Rossit},
  title =   {Effective heuristics for permutation and non-permutation ﬂow shop scheduling with missing operations},
  journal = {Comp. Oper. Res.},
  year =    {2024},
  volume =  {106742}
}
```
