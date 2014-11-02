#A combined MPI and OpenMP implementation of the simulation of a chain of particles

This project is a continuation from the first project. The aim of the project is
to implement the program from the first project into a combined MPI-OpenMP
framework. This is quite a common approach for exploiting both coarse and fine
grain parallelism in programs. A problem is partitioned coarsely at the top
level and finely within each individual part. The coarse level partitioning is
done using MPI and the finer level partitioning is usually done on a multi-core
machine or on a Graphics Processing Unit (GPU).

The aim of the project is to partition the linear chain into a few smaller parts
and distribute these parts to different computers by using MPI. The computation
on each part would now occur within the individual machines using the cores
available on those machines.

The coding in the project will be minimal, it is the case quite often that
parallelizing a piece of sequential code requires only small but well thought
out modifications. Also, though in principle your code should be executable
across a cluster of machines, it will be sufficient to test it across two
machines.

The main challenge in the project would be the communication and data
distribution. Each particle in the chain gets a new position, acceleration,
kinetic energy and potential energy. However, the latter two are dependent on
the first two. Assume that you have divided the chain of particles into two
parts, P1 and P2, from left to right. In other words, the first 50 particles are
in P1 and the last 50 particles are in P2. Assume that the last particle in P1
(the 50th particle) is called p and the first particle in P2 is called q. We
need the position of p for calculating the aceleration of q. This calculation is
done in the accel(.,.) function. Hence there is a need to send the position of p
to the computer that has been allocated P2 for a correct computation of the
acceleration, position and velocities of the particles in P2 (the latter two are
updated elsewhere in the code).

Hence the computation should proceed in rounds. The position and accelerations
from a chain need to be sent to the chain to its right after each round. The
computation within each part will be done using OpenMP.

**The deliverables:**

- The first deliverable is of course your modified C code with appropriate MPI
and OpenMP directives. You should also comment your code suitably, so that the
code itself can be read and your modifications can be understood.
- You have to also submit a (reasonably) small document where you should explain
how you have implemented the parallelism in the code and why. Include also all
the decisions related to your implementation in the first project, so that the
document can be read without checking your first project again.
- You should develop your code on a single machine with multiple MPI processes.
Performance improvement is not the primary focus of this project. Hence, you can
experiment with finer partitions of the chain. Though the finer the partition,
the communication overhead will be higher. However, in reality the amount of
data that need to be sent is small. You should eventually run your code on at
least two machines. OpenMPI is now installed in the machines in Lab 2.03 and you
will need to ssh to multiple machines without password. I have explained in the
lab sheet how this can be done.
- You should conduct a performance analysis of your program, irrespective of a
speed-up or slow down. First, analyse the performance on a single machine with
different partitions of the chain. It will be better to use much longer chains
(the original chain is of length 100) and relatively smaller values of dt for
faster turnaround. The second performance analysis should be done on multiple
machines, you should use at least two machines. You can also use a varying
number of threads in OpenMP for your performance analysis.
- Of course the output from your parallelized program should match the output of
the sequential program.

**Note:** The project can be done either individually or in a group consisting
of a maximum of two students.

**Amitava Datta  
September 2014**

Taken From:  
http://undergraduate.csse.uwa.edu.au/units/CITS3402/labs/project2.html
