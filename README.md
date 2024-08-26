# ToolLoadingProblem
Comparison of Max Pipe Construction Algorithm (MPCA) and Keep Tool Needed Soonest (KTNS) algorithms for the Tool Loading Problem (TLP)

Source code, benchmark files and README file also available at https://github.com/cherniavskii-mikhail/ToolLoadingProblem.

A draft of the paper available at https://arxiv.org/abs/2205.06042.

Abbreviations:
1) TLP - the Tool Loading Problem
2) SSP - the job Sequencing and tool Switching Problem
3) KTNS - Keep Tool Needed Soonest algorithm (proposed by Christopher Tang, Eric Denardo, see https://doi.org/10.1287/opre.36.5.767)
4) MPCA - Max Pipe Construction Algorithm (proposed by the authors of the manuscript)

Structure:

1) MPCA_vs_KTNS_comparison.cpp contains C++ code for comparison of MPCA and KTNS algorithms.
2) Folder Catanzaro contains classical benchmark instances for the SSP and TLP introduced in https://doi.org/10.1016/j.ejor.2015.02.018
3) README.md (this file)

To reproduce the computational results:

1) Compile MPCA_vs_KTNS_comparison.cpp. To compile run: g++ MPCA_vs_KTNS_comparison.cpp -o MPCA_vs_KTNS_comparison.exe -O3
2) Run compiled MPCA_vs_KTNS_comparison.exe.
3) Open results.txt to see results. Compare results with the Table in Section Experimental results of the manuscript.


Format of benchmark Tool Loading Problem instances:
1) The first line contains a single integer N - the number of jobs.
2) The second line contains a single integer M - the number of tools.
3) The third line contains a single integer C - the capacity of the magazine.
4) After the third line there is a 0-1 matrix of size M x N, where element in pisition (i, j) equals 1, if tool i is needed for job j and 0, otherwise. 

Other well known benchmark instances (Catanzaro, Crama, Yanasse, Mecler) for the SSP and TLP are available at https://github.com/jordanamecler/HGS-SSP/tree/master.

MPCA_vs_KTNS_comparison.cpp description:

1) KTNS function is an implementation of the KTNS algortihm by Jordana Mecler, Anand Subramanian, Thibaut Vidal, see https://doi.org/10.1016/j.cor.2020.105153. We use KTNS as a "base-line" algorithm to compare it with our algorithms MPCA, MPCA_bitwise and our bitwise implementation of KTNS called KTNS_bitwise.
2) MPCA function is our new MPCA proposed in the manuscript.
3) MPCA_bitwise function is our bitwise implementation of MPCA.
4) KTNS_bitwise function is our bitwise implementation of KTNS.
5) algorithms_verefication function runs 10^5 random tests to verify that KTNS, MPCA, KTNS_bitwise, MPCA_bitwise give the same output.
6) Experiments structure is designed to run computational experiments, i.e. enumerate problem instances, read data from files, execute algorithms, measure the time spent on calculations, write cumputational results to file results.txt.
7) TLP_Data structure is a container for a TLP instance data.
8) remove_ones is an auxiliary function for KTNS_bitwise.

