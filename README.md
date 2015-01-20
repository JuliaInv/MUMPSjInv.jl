#MUMPS

Efficiently solving large and sparse linear systems is essential in many areas of scientific computing. This project provides an interface to [MUMPS ( a MUltifrontal Massively Parallel sparse direct Solver)](http://mumps.enseeiht.fr/) in Julia. To improve the memory-efficiency of the method the software package  [METIS (Serial Graph Partitioning and Fill-reducing Matrix Ordering)](http://glaros.dtc.umn.edu/gkhome/views/metis) is used. Please refer to the MUMPS and metis websites for more information and updated versions.

## Copyright

Please make sure you have read and understood the licenses of [MUMPS  4.10.0](http://graal.ens-lyon.fr/MUMPS/index.php?page=dwnld)  and [Metis (4.0.3)](http://www.filewatcher.com/m/metis-4.0.3.tar.gz.522624-0.html). Due to the use of Metis, this code can be freely used for educational and research purposes by non-profit institutions and US government agencies only. The wrapper itself is provided as free software under MIT License in the hope it is useful.

Note that this software is provided "as is" with absolutely no warranty. Use at your own risk.


## Benchmark Example

In our experience MUMPS is considerably faster than Julia's backslash both in  solving real or complex systems. It also seems to consume less memory. Another benefit is that once the factorization is generated it can be applied to multiple right-hand sides. Here are results of two benchmarks ran on a MacBookPro with 2.6 GHz Intel Core i7 Processor and 16 GB RAM.

### Solving Poisson System
We are solving linear systems arising from discretization of Poisson's equation for different discretization sizes for a random right hand side. We compare the runtime for MUMPS solver and the backslash. 

| Grid size | MUMPS (sec) | Julia (sec) | Speedup|
| --------- | ------------|-------------|--------|
| 8^3  		|  0.002      |  0.003      | 1.3 	 |
| 16^3 		|  0.025      | 0.043       | 1.2 	 |
| 24^3 		|  0.161      | 0.286       | 1.7 	 |
| 32^3 		|  0.545      | 1.747       | 2.9 	 |
| 48^3 		|  3.515      | 20.738      | 5.5 	 |
| 64^3 		|  13.286     | 125.902     | 7.2 	 |

### Solving a complex system
As above, we obtain  discretizations of Poisson's equation for different grid sizes,  add an imaginary identity matrix and solve for one complex random right hand side. We compare the runtime for MUMPS solver and the backslash. 

| Grid size | MUMPS (sec) | Julia (sec) | Speedup|
|----------:|------------:|------------:|-------:|
| 8^3  		|0.004        | 0.065       | 15.8   |
| 16^3 		|0.044        | 0.073       | 1.7    |
| 24^3 		|0.354        | 0.629       | 1.8    |
| 32^3 		|1.418        | 4.106       | 2.9    |
| 48^3 		|12.150       |  53.757     | 4.4    |
| 64^3 		|59.465       |  MEM_EX     |   	 |

MEM_EX: On the 64^3 grid computation was interrupted, memory was insufficient and Julia began swapping.

Some more tests are provided in MUMPS/tests. 

## Installing MUMPS

This package was tested on Mac and Ubuntu only. Here, open a terminal navigate to MUMPS/src and run make. Please do a pull-request if you have improvements. 


**ToDo: Make MUMPS available from Julia's package manager! **





