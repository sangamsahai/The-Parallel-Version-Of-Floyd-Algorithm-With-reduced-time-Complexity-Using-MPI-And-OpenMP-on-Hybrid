# The-Parallel-Version-Of-Floyd-Algorithm-With-reduced-time-Complexity-Using-MPI-And-OpenMP-on-Hybrid
This Algorithm reduces the time complexity of Floyd's Algorithm by a huge factor. 
The complexity comes down from ( to O(N Log_2 P / Sqrt(P)) )N^3) where P is the number of nodes used in the cluster.
This is a hybrid approach in which we use MPI for inter process communication and OpenMP for Multithreading.
 MPI is used for communication between the nodes(distributed memory architecture).
 And OpenMP facilitates the multithreading on different cores of the same node, which have shared memory architecture.

The Document 'Parallel Programming- Final Project Report.docx' discusses the problem in detail.
It then elaborates the technique and the Algorithm used.
And finally , it briefs the steps to compile and run the C code on a cluster.
