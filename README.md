# ToolLoadingProblem
Comparison of MPCA and KTNS algorithms for the Tool Loading Problem

Code files, benchmark files and README file also available at https://github.com/cherniavskii-mikhail/ToolLoadingProblem

A draft of the paper available at https://arxiv.org/abs/2205.06042

Structure:

1) IGA_vs_KTNS_comparison.cpp contains C++ code for comparison of MPCA and KTNS algorithms.
2) Folder Catanzaro contains classical benchmark instances for the Job Sequencing and Tool Switching Problem introduced in https://doi.org/10.1016/j.ejor.2015.02.018
3) README.md file

To reproduce the computational results:

1) Compile IGA_vs_KTNS_comparison.cpp. To compile run: g++ MPCA_vs_KTNS_comparison.cpp -o MPCA_vs_KTNS_comparison.exe -O3
2) Run compiled IGA_vs_KTNS_comparison.exe.
3) Open results.txt to see results. Compare results with Table 6 in Section EXPERIMENTAL RESULTS of the manuscript.


Format of benchmark Tool Loading Problem instances:
1) The first line contains a single integer N - the number of jobs.
2) The second line contains a single integer M - the number of tools.
3) The third line contains a single integer C - the capacity of the magaz
4) After the third line there is a 0-1 matrix of size M x N, where element in pisition (i,j) equals 1 if tool j is needed for job i, and 0 otherwise. 

Other well known benchmark datasets (Catanzaro, Crama, Yanasse, Mecler) for the SSP and TLP are available at https://github.com/jordanamecler/HGS-SSP/tree/master.

