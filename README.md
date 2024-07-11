# Non-permutation flow shop scheduling with missing operations

This repository contains detailed results and code for the paper [Non-permutation flow shop scheduling with missing operations](https://doi.org/10.1016/j.cor.2024.106742).

## Detailed results

The detailed results can be found in folder [data](data). In the folder you also can find a [R notebook](data/tables.Rmd) to reproduce the tables, including a [rendered version in HTML](data/tables.html).

## Code

The code is contained in the folder `src`. To compile use
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../src
make
```

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
