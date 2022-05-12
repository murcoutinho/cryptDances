# cryptDances
Tool to perform cryptanalysis of ChaCha, Salsa and Forro in high performance environments

## Dependencies
This project was developed with OpenMPI 2.1.1 and CUDA 11.2, there shouldn't be any major issue with more recent versions, but it **could**.

## Build

The first lines of the Makefile specify the full path of the compilers we're using since they're not in out PATH and are not installed in their default locations. In order to build those have to be changed to match your environment, including nvcc's arch flag. That being done it's a simple:

```
make
```

## Running

In order to run in a computer with a single gpu:

```
mpirun -np 2 [executable]
```

The number of process is the number of GPUs you want to use +1. In order to run in a cluster with MPI you'll need to specify a hostfile. Supposing that the setup has 3 computers: *front*, *compute-0-1* and *compute-0-2*. Where front does not have a GPU, and *compute-0-x* has 2 (each).

hosts:

---
front slots=1
compute-0-1 slots=2
compute-0-2 slots=2
---

By running the following command, *front* will run a process gathering data from the other cluster nodes, while *compute-0-x* will run a process for each GPU it has doing the actual calculations and sending back the results to *front*.

```
mpirun -np 5 -hostfile hosts [executable]
```
