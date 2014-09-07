#Programming Project or Programming Challenge

The aim of the programming project is to parallelize a C code on multi-core machines to gain as much speed-up as possible. The project will be evaluated on the degree of speed-up that your program achieves. I will provide you with a machine where you can analyze the performance of your program. As such, I expect the coding to be minimal, as understanding is the key to the success in this project. Also, your code will be very useful in my research as I may use the code in future in the Amazon cloud, deploying many cores.

**Background:** Though understanding of the problem is not strictly necessary, here is a background on the problem. This is a modified version of a famous simulation known as the Fermi-Pasta-Ulam (FPU) problem in non-linear dynamics. These three scientists used this simulation to implement one of the first compute intensive tasks in a digital computer for understanding nonlinear phenomena. I use it for studying two special kinds of non-linear waves, called solitons and breathers with some modifications in the settings of the FPU problem. I will not go into the details of the physics behind it, but just describe the simulation setting.

It is a one-dimensional string of particles connected by bonds. When a particle is given a displacement, its neighbours get displaced as well, and a valid intuition is that the disturbance will propagate through the chain. The disturbance of a particle can be measured by noting its displacement, or kinetic energy. However, depending on the two parameters alpha and beta, it is possible to localize the energy in this chain and that is the focus of this study and the code.

The dynamics of the system, i.e., the displacements of the particles, is computed using Newton's second law. However, we need to numerically integrate Newton's second law to compute the displacement and energy of each particle at every step. You have already seen an example of numerical integration in the lab when we approximated the value of Pi. You have also seen that the the accuracy of the integration depends on the step size used. The accuracy in this project is the most important aspect of the simulation, as we want to understand the physical phenomenon and any numerical inaccuracy can give spurious results that are not related to the physical phenomenon. Hence, we have to use very small step sizes, typically 10^(-7) or 0.00000001 or even lower. And this is the main problem. The program takes a long time to execute when we use such small step sizes. But we have to simulate many scenarios within a reasonable amount of time and hence we need speed-up.

**The Program:** Here is the code for the simulation.

The program has about 225 lines of C code and it is your responsibility to understand the flow of the program. The integration method used is called the "Velocity Verlet" integration. It is not necessary to understand how it works, but you can lookup the Wikipedia article if you are interested. Here are some details of the program.

 - The program is written in a bit of old style C and there are many global variables. The most important parts for us are the main() and the accel() functions. The latter function calculates the acceleration at every step and is called from the main.
 - The simulation is done in discrete steps and the variable nprntstps indicates the number of simulation steps. It is currently set to 10001, but you can change (reduce) it to make the running time faster.
 - The key variable is dt which determines the step size for the Verlet integration. It is currently set to 0.001, but ideally it should be much smaller to avoid numerical errors. As I have mentioned, it should be 0.0000001 or even 0.00000001 or 0.000000001. I usually get acceptable numerical errors with the last value of dt. You should use only larger values for dt while developing your parallel code and testing performance, but eventually should use small values for dt to check your results. The other important parameter is chainlngth, which is set to 100 mainly because the simulation time becomes prohibitive if I use more than 100 particles. I hope to increase this parameter to higher values when I get a very fast implementation from you.
 - There two other very important parameters alpha and beta that set the ratio of linearity and non-linearity in the system. However, these two parameters have no effect on the run time of the program. Still you should run your final program with a few different values of these two parameters to make sure that your program is correct (more on correctness below).
 - The real action is inside the big while loop (while (n < nprntstps)) in the main() function. I will not explain the logic as it is your task to understand the control flow in the code, however, the running time of the program is almost entirely due to this while loop. Note that the overall running time of the program depends on both nprntstps and dt.
 - I don't have any suggestion as to how you should parallelize the code and you can use whatever method you like using of course OpenMP, or directly using pthreads. However, you have to make sure that the resulting program gives correct output. You will notice the program writes to several files. These are mainly displacement, velocity, kinetic and potential energies of the particles in the chain. You have to make sure the corresponding files from the original sequential program and your modified program match exactly. You can write a small program or use some linux utility to compare the files. I will compare the files while marking the project and the files from your code must exactly match the files from the sequential code for the same values of all the parameters that I have mentioned above. Hence you must test your code in the machine provided.

**The machine:**

We have set up a machine hpc.csse.uwa.edu.au for this project. Hyperthreading is already disabled in this machine. It has four cores and hence, it is possible that optimal performance can be gained with four threads. However, you are welcome to use any number of threads as long as you get an improvement in performance. We are also implementing a 'batch system' so that only one of your programs will be executed at any time in this computer. This is important as otherwise your performance figures will not be correct. I will soon inform you about how to submit your program for execution. You can access the machine this way:

 - ssh -l xxxxxxxx hpc.csse.uwa.edu.au  
where xxxxxxxx is your 8-digit student number. Use your pheme password to login.
 - It will create a temporary home area for you where you can edit, compile and run your program. You can also access your linux home area from there.

**The deliverables:**

 - The first deliverable is of course your modified C code with appropriate OpenMP directives. You should also comment your code suitably, so that the code itself can be read and your modifications can be understood.
 - You have to also submit a (reasonably) small document where you should explain how you have implemented the parallelism in the code and why. It should also contain extensive performance evaluation of your code with different parameters. I am not looking for anything specific here, rather well designed performance evaluation and speed-up results that highlight your work.
Note: The project can be done either individually or in a group consisting of a maximum of two students.

**Amitava Datta  
September 2014**

Taken From:  
http://undergraduate.csse.uwa.edu.au/units/CITS3402/labs/project1.html
