# pRIblast
[![doi badge](https://badgen.net/badge/DOI/j.future.2022.08.014/blue)](https://doi.org/10.1016/j.future.2022.08.014)
[![compile pRIblast (release)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-release.yml/badge.svg)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-release.yml)
[![compile pRIblast (debug)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-debug.yml/badge.svg)](https://github.com/UDC-GAC/pRIblast/actions/workflows/compile-priblast-debug.yml)

pRIblast is a high efficient, parallel application for extensive lncRNA-RNA interaction prediction. It is based on the work of T. Fukunaga and M. Hamada, [RIblast](https://github.com/fukunagatsu/RIblast/), and it has been fully optimized to reduce I/O latencies and memory usage to the bare minimum.

## Version
Version 0.0.1.

## Requirements
To compile and execute pRIblast, the following packages are required:
* GNU Make (>=v3.82).
* GCC (>=v7.3.0).
* OpenMPI (>=v3.1.0).
* OpenMP runtime (>=v3.1).

Note that there exist different MPI implementations, such as OpenMPI or MPICH. Please, make sure your library has support for RDMA communications (i.e. one-sided operations) before compiling pRIblast.
## Compilation
Download the source code from this repository, either use Git or download a copy from GitHub, and let GNU Make automatically compile pRIblast for you. As a result, there will be a newly created binary file named pRIblast in the `target` folder of your current working directory.

If you are interested in benchmarking the application, compile with `make BUILD=debug`.

## Execution
To execute pRIblast, fetch the MPI runtime interface as follows
```
foo@bar:target$ mpirun -np <p> -x OMP_NUM_THREADS=<t> pRIblast <options>
```
where `<p>` is the numnber of processes that will exist in the MPI group and `<t>` the number of threads spawned per MPI process.

As for the program options, [RIblast](https://github.com/fukunagatsu/RIblast/)'s official repository provides a fairly detailed list of available execution modes (i.e. database construction and RNA interaction search) and per mode flags. However, pRIblast implements new options to have fine grained control over the execution of the parallel algorithm. Those options are:
```
(db)  -c <uint>, sets a database chunk size (smaller chunk implies less memory usage)
(ris) -p  <str>, sets a per process local path for fast writing of temporary output files
(ris) -a  <str>, sets the parallel algorithm used to distribute data among processes (block|area|dynamic)
```

### Execution example
Suppose you want to execute the `area-sum` algorithm in a 16-node multicore cluster using the `drosophila` dataset, a database chunk of 500 sequences and 1 process per node with 16 threads each. Furthermore, there exist a local, temporary disk attached to every node located in `/tmp/scratch`. First, download the `drosophila` dataset from the [Ensembl genome browser](ftp://ftp.ensembl.org/pub/release-97/fasta/). Secondly, use the script `misc/preprocess.sh`
```
foo@bar:misc$ ./preprocess.sh /path/to/drosophila.fa
```
to split the FASTA file into two different datasets: one with all RNAs (used to build the target database) and a lncRNA-only file (used to predict interactions against the database). Third, create the fragmented database running the pRIblast database construction step
```
foo@bar:target$ ./pRIblast db -i /path/to/db-drosophila.fa \
                              -o /path/to/db-drosophila \
                              -c 500
```
And finally, execute the interaction search step on 16 nodes issuing the following command
```
foo@bar:target$ mpirun -np 16 -x OMP_NUM_THREADS=16 \
                       pRIblast ris -i /path/to/ris-drosophila.fa \
                                    -o /path/to/out-drosophila.txt \
                                    -d /path/to/db-drosophila \
                                    -p /tmp/scratch \
                                    -a area
```

### Configuration of threads, processes and algorithms
To achieve maximum performance, use the pRIblast parallel algorithms as follows:
* Do not use the `pure-block` algorithm (i.e. `-a block`). Its only purpose is to benchmark.
* Use the `area-sum` algorithm (i.e. `-a area`) if there exist plenty of compute resources with respect to the dataset size. No hyperthreading needed. Use as many threads as processes.
* Use the `dynamic-decomp` (i.e. `-a dynamic`) algorithm if the number of computing resources is small with respect to the dataset size. Hyperthreading will not hurt. More processes than threads works best here.

## License
pRIblast is free software and as such it is distributed under the [MIT License](licenses/MIT.txt). However, pRIblast makes use of a header-only module (see [minmaxheap.h](src/minmaxheap.h)) which is not an original work. Therefore, its usage is subject to the [GPLv3](licenses/GPL3.txt) License and all rights are reserved to [Kilian Gebhardt](https://github.com/kilian-gebhardt).
