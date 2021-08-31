# pRIblast
pRIblast is a high efficient, parallel application for extensive lncRNA-RNA interaction prediction. It is based on the work of T. Fukunaga and M. Hamada, [RIblast](https://github.com/fukunagatsu/RIblast/), and it has been fully optimized to reduce I/O latencies and memory usage to the bare minimum.

## Version
Version 0.0.0.

## Requirements
To compile and execute pRIblast, the following packages are required:
* GNU Make (>=v3.82).
* GCC (>=v7.3.0).
* OpenMPI (>=v3.1.0).
* OpenMP runtime (>=3.1).

Note that there exist different MPI implementations, such as OpenMPI or MPICH. Please, make sure your library has support for RDMA (i.e. one-sided operations) before compiling pRIblast. Also, it is possible to compile the program with older GCC and OpenMP versions, but they will not be officially supported.

## Compilation
First of all, download the source code from this repository. Either use Git or download a copy from GitHub. Secondly, navigate to the `src` directory and choose which version of the parallel algorithm you would like to use (see below for further advice):
* `pure-block`: naive pure block decomposition algorithm
* `area-sum`: greedy algorithm which decomposes sequences according to their size
* `dynamic-decomp`: share all algorithm, all processes see the same batch of sequences

And lastly, navigate to the folder of the algorithm you wish to use and let GNU Make automatically compile it for you. As a result, there will be a newly created binary file named pRIblast in your current working directory.

## Execution
To execute pRIblast, fetch the MPI runtime interface as follows
```
foo@bar:src/alg$ mpirun -np <p> -x OMP_NUM_THREADS=<t> pRIblast <options>,
```
where `<p>` is the numnber of processes that will exist in the MPI group and `<t>` the number of threads spawned per MPI process.

As for the program options, [RIblast](https://github.com/fukunagatsu/RIblast/)'s official repository provides a fairly detailed list of available execution modes (i.e. database construction and RNA interaction search) and per mode flags. However, pRIblast implements new options to have fine grained control over the execution of the parallel algorithm. Those options are:
```
(db)  -c <uint>, sets a database chunk size (smaller chunk implies less memory usage)
(ris) -p <path>, sets a per process local path for fast writing of temporary output files
(ris) -t <0|1>, debug execution
```

### Execution example
Suppose you want to execute the `area-sum` algorithm in a 16-node multicore cluster using the `drosophila` dataset, a database chunk of 500 sequences and 1 process per node with 16 threads each. Furthermore, there exist a local, temporary disk attached to every node located in `/tmp/scratch`. First, download the `drosophila` dataset from the [Ensembl genome browser](ftp://ftp.ensembl.org/pub/release-97/fasta/). Secondly, use the script `misc/preprocess.sh`
```
foo@bar:misc$ ./preprocess.sh /path/to/drosophila.fa,
```
to split the FASTA file into two different datasets: one with all RNAs (used to build the target database), and a lncRNAs only (used to predict interactions against the database). Third, create the fragmented database running the pRIblast database construction step
```
foo@bar:src/area-sum$ ./pRIblast db -i /path/to/db-drosophila.fa \
                                    -o /path/to/db-drosophila \
                                    -c 500.
```
And finally, execute the interaction search step on 16 nodes issuing the following command
```
foo@bar:src/area-sum$ mpirun -np 16 -x OMP_NUM_THREADS=16 \
                             pRIblast ris -i /path/to/ris-drosophila.fa \
                                          -d /path/to/db-drosophila \
                                          -o /path/to/out-drosophila.txt \
                                          -p /tmp/scratch.
```

### Configuration of threads, processes and algorithms
To achieve maximum performance, use the pRIblast algorithm as follows:
* Do not use the `pure-block` algorithm. Its only purpose is to benchmark.
* Use the `area-sum` algorithm if there exist plenty of compute resources with respect to the dataset size. No hyperthreading needed. Use as many threads as processes.
* Use the `dynamic-decomp` algorithm if the number of computing resources is small with respect to the dataset size. Hyperthreading will not hurt. More processes than threads works best here.
