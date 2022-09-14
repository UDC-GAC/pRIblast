# pRIblast
[![doi badge](https://badgen.net/badge/DOI/j.future.2022.08.014/blue)](https://doi.org/10.1016/j.future.2022.08.014)
[![compile pRIblast (release)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-release.yml/badge.svg)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-release.yml)
[![compile pRIblast (debug)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-debug.yml/badge.svg)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-debug.yml)

pRIblast is a high efficient, parallel application for extensive lncRNA-RNA interaction prediction. It is based on the work of T. Fukunaga and M. Hamada, [RIblast](https://github.com/fukunagatsu/RIblast/), and it has been fully optimized to reduce I/O latencies and memory usage to the bare minimum.

## Version
Version 0.0.2.

## Requirements
To compile and execute pRIblast, the following software is required:
* GNU Make.
* C++ compiler (with support for OpenMP directives).
* MPI implementation (MPI-3 compliant).

For instance, a valid combination of these tools may be: GNU Make v3.82, GCC v8.3.0 and OpenMPI v3.1.4.

## Compilation
Download the source code from this repository, either use Git or download a copy from GitHub, and let GNU Make automatically compile pRIblast for you. As a result, there will be a newly created binary file named pRIblast in the `target` folder of your current working directory.

If you are interested in benchmarking the application, compile with `make BUILD=debug`.

## Execution
To execute pRIblast, fetch the MPI runtime interface as follows
```
foo@bar:target$ mpirun -np <p> -x OMP_NUM_THREADS=<t> pRIblast <options>
```
where `<p>` is the numnber of processes that will exist in the MPI group and `<t>` the number of threads spawned per MPI process.

As for the program options, [RIblast](https://github.com/fukunagatsu/RIblast/)'s official repository provides a detailed list of the available execution modes (i.e. database construction and RNA interaction search) and per mode parameters. However, pRIblast implements new options to have fine grained control over the execution of the parallel algorithm. Those options are:
```
(db)  -a  <std>, sets the parallel algorithm used to distribute data among processes (block|heap|dynamic)
(db)  -p  <str>, sets a per process local path for fast writing of temporary output files
(db)  -c <uint>, sets a database chunk size (smaller chunk implies less memory usage)
(ris) -a  <str>, sets the parallel algorithm used to distribute data among processes (block|area|dynamic)
(ris) -p  <str>, sets a per process local path for fast writing of temporary output files
```

### Execution example
Suppose you want to execute pRIblast (both the `db` and `ris` steps using the `heap` algorithm and the `area-sum` algorithm, respectively) in a 16-node multicore cluster using the `drosophila` dataset, a database chunk of 500 sequences and 1 process per node with 16 threads each. Furthermore, there exist a local, temporary disk attached to every node located in `/tmp/scratch`. First, download the `drosophila` dataset from the [Ensembl genome browser](ftp://ftp.ensembl.org/pub/release-97/fasta/). And then, use the script `misc/preprocess.sh`
```
foo@bar:misc$ ./preprocess.sh /path/to/drosophila.fa
```
to split the FASTA file into two different datasets: one with all RNAs (used to build the target database) and a lncRNA-only file (used to predict interactions against the database). Secondly, create the fragmented database running the pRIblast database construction step
```
foo@bar:target$ mpirun -np 16 -x OMP_NUM_THREADS=16 \
                       pRIblast db -i /path/to/db-drosophila.fa \
                                   -o /path/to/db-drosophila \
                                   -a heap \
                                   -p /tmp/scratch \
                                   -c 500
```
And finally, execute the interaction search step on 16 nodes issuing the following command
```
foo@bar:target$ mpirun -np 16 -x OMP_NUM_THREADS=16 \
                       pRIblast ris -i /path/to/ris-drosophila.fa \
                                    -o /path/to/out-drosophila.txt \
                                    -d /path/to/db-drosophila \
                                    -a area \
                                    -p /tmp/scratch
```

### Configuration of threads, processes and algorithms
To achieve maximum performance, use the pRIblast parallel algorithms as follows:
* Do not use the `pure-block` algorithm (i.e. `-a block`). Its only purpose is to benchmark.
* Use the `heap` (database construction step) and the `area-sum` (RNA interaction search step) algorithms (i.e. `-a heap` and `-a area`) if computing nodes have a high number of CPU cores available to take advantage of the multithreading performance optimization heuristics developed within the tool. Don't use hyperthreading. Spawn one process per socket and run as many threads as cores it has.
* Use the `dynamic` algorithm (i.e. -a `dynamic`) if computing nodes are single processors or the tool is not run exclusively in the nodes. Spawn one process per core. Hyperthreading may help.

## Cite us
```
Iñaki Amatria-Barral, Jorge González-Domínguez and Juan Touriño. "pRIblast: A highly efficient parallel application for comprehensive lncRNA-RNA interaction prediction". Future Generation Computer Systems (2022)
```

## License
pRIblast is free software and as such it is distributed under the [MIT License](LICENSE). However, pRIblast makes use of several modules which are not original pieces of work. Therefore, its usage is subject to their correspoding [THIRDPARTYLICENSES](THIRDPARTYLICENSES) License and all rights are reserved to their authors.
