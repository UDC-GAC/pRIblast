# pRIblast
![gnu workflow](https://img.shields.io/github/actions/workflow/status/UDC-GAC/pRIblast/compile-gnu.yml?label=gnu)
[![issues](https://img.shields.io/github/issues/UDC-GAC/pato)](https://github.com/UDC-GAC/pRIblast/issues)
[![doi](https://img.shields.io/badge/doi-j.future.2022.08.014-blue)](https://doi.org/10.1016/j.future.2022.08.014)
![version](https://img.shields.io/badge/release-v0.0.3-blue)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

pRIblast is a high efficient, parallel application for extensive lncRNA-RNA interaction prediction. pRIblast is based on the work of T. Fukunaga and M. Hamada, [RIblast](https://github.com/fukunagatsu/RIblast/), and it has been fully optimized to reduce I/O latencies and memory usage to the bare minimum.

## Version
Version 0.0.3.

## Requirements
To compile and execute pRIblast, the following software is required:
* GNU Make.
* C++ compiler (with support for OpenMP and the `c++17` standard).
* MPI implementation (MPI-3 compliant).

For instance, a valid combination of these tools may be: GNU Make v3.82, GCC v9.3.0 and OpenMPI v3.1.4.

## Compilation
Download the source code from this repository, either use Git or download a copy from GitHub, and let GNU Make automatically compile pRIblast for you. As a result, there will be a newly created binary file named pRIblast in the `target` folder of your current working directory.

## Execution
To execute pRIblast, fetch the MPI runtime interface as follows
```
mpirun -np <p> -x OMP_NUM_THREADS=<t> pRIblast <options>
```
where `<p>` is the numnber of processes that will exist in the MPI group and `<t>` the number of threads spawned per MPI process.

As for the program options, [RIblast](https://github.com/fukunagatsu/RIblast/)'s official repository provides a detailed list of the available execution modes (i.e. database construction and RNA interaction search) and per mode parameters. However, pRIblast implements new options to have fine grained control over the execution of the parallel algorithm. Those options are:
```
 (db) -a  <std>, sets the parallel algorithm used to distribute data among processes (block | heap | dynamic)
 (db) -p  <str>, sets a per process local path for fast writing of temporary output files
 (db) -c  <int>, sets the database page size (smaller page implies less memory usage)
(ris) -a  <str>, sets the parallel algorithm used to distribute data among processes (block | area | dynamic)
(ris) -p  <str>, sets a per process local path for fast writing of temporary output files
```

### Execution example
Suppose you want to execute pRIblast (both the `db` and `ris` steps using the `dynamic` algorithm) on a 16-node multicore cluster using a FASTA file `db.fa`, which contains RNA sequences to construct a database, a FASTA file `ris.fa`, which contains the RNA sequences you want to predict interactions against the database, a page size of 500 sequences, and 1 process per node with 16 threads each. Furthermore, there exist a local, temporary disk attached to every node located in `/tmp/scratch` that allows fast writing of temporary output files.

First, create the target RNA database running the pRIblast database construction step as follows
```bash
mpirun -np 16 -x OMP_NUM_THREADS=16 \
       pRIblast db -i db.fa -o rna-db -a dynamic -p /tmp/scratch -c 500
```
And then, predict interactions against the database running the pRIblast RNA interaction search step as follows
```bash
mpirun -np 16 -x OMP_NUM_THREADS=16 \
       pRIblast ris -i ris.fa -o predictions.txt -d /path/to/rna-db -a dynamic -p /tmp/scratch
```

Note that the `-p` option is not mandatory, but it is highly recommended to use it if there exists a local, temporary disk attached to every node, as this will reduce I/O latencies. And also, note that the `-c` option is only available for the database construction step. It sets the page size of the database, i.e. the number of RNA sequences that will be loaded into memory at once. The smaller the page size, the less memory will be used.

### Configuration of threads, processes and algorithms
To achieve maximum performance, avoid running the `pure-block` algorithm. Its only purpose is to benchmark. Instead, use the `heap` (database construction step) and the `area-sum` (RNA interaction search step) algorithms if computing nodes have a high number of CPU cores available to take advantage of the multithreading performance optimization heuristics developed within the tool. Spawn one process per socket and run as many threads as cores it has. Otherwise, use the `dynamic` algorithm if the number of available nodes is low and/or the number of CPU cores per node is low. Spawn one process per core.

## Cite us
If you use pRIblast in your research, please cite our work using the following reference:
```
Iñaki Amatria-Barral, Jorge González-Domínguez and Juan Touriño. "pRIblast: A highly efficient parallel application for comprehensive lncRNA-RNA interaction prediction". Future Generation Computer Systems, 138, pages 270-279 (January 2023)
```

## License
pRIblast is free software and as such it is distributed under the [MIT License](LICENSE). However, pRIblast makes use of several modules which are not original pieces of work. Therefore, its usage is subject to their correspoding [THIRDPARTYLICENSE](THIRDPARTYLICENSES) and all rights are reserved to their authors.
