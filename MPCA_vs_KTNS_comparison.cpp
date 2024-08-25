#include <iostream>
#include <fstream>
#include <vector>
#include <string> 
#include <stdio.h> 
#include <time.h> 
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <iomanip>

#define vec vector

using namespace std;

const string OUTPUT_FILE = "results.txt"; 

struct Data {
    vec<vec<int> > jobsToolsMatrix;
    vec<vec<int> > L;
    vec<vec<int> > loadedMatrix;
    int numJobs, numTools, magCapacity;
    vec<int> used;
    vec<int> W_n;

    int start_dist;
    vec<int> filled_slots;
    vec<int> job_capacity;
    vec<int> cromosome_free;
    vec<unsigned long long> cromosome_bytes;
    vec<int> jobnum_to_free;
    vec<unsigned long long> jobnum_to_bytes;
    vec<vec<int> > jobs_pairs_intersect;

    vec<vec<int> > jobnum_to_tools_list;
    vec<vec<int> > cromosome_jobnum_to_tools_list;
    vec<int> last_seen_tool;
    vec<int> usedToFullMag;
    vec<vec<int> > MAGAZINE;
    vec<vec<int> > H_1_MAG;
    vec<vec<int> > tool_in_i_not_in_prev;
    vec<unordered_set<int> > MAG_SETS;
    vec<unordered_set<int> > T_SETS;
    vec<int> empty_vec;
    vec<vec<int> > tool_in_i_not_in_prev_start;
    vec<vec<int> > MN;
    Data(int numJobs_arg, int numTools_arg, int magCapacity_arg, vec<vec<int> >& jobsToolsMatrix_arg)
        :numJobs(numJobs_arg),
        numTools(numTools_arg),
        magCapacity(magCapacity_arg),
        jobsToolsMatrix(jobsToolsMatrix_arg) {

        used.assign(numTools, 0);
        L.assign(numTools, vec<int>(numJobs, -1));
        W_n.assign(numTools, -1);
        loadedMatrix.assign(numJobs, vec<int>(numTools, -1));


        jobnum_to_bytes.assign(numJobs, 0LL);
        for (int i = 0; i < numJobs; i++) {
            string str = "";
            for (auto v : jobsToolsMatrix[i]) { str += to_string(v); }
            jobnum_to_bytes[i] = bitset<64>(str).to_ullong();
        }
        cromosome_bytes.assign(numJobs, 0LL);

        start_dist = magCapacity * (numJobs - 1);
        jobnum_to_free.assign(numJobs, 0);
	    job_capacity.assign(numJobs, 0);
        for (int i = 0; i < numJobs; i++) {
	    job_capacity[i] = (int)count(jobsToolsMatrix[i].begin(), jobsToolsMatrix[i].end(), 1);
            jobnum_to_free[i] = magCapacity - (int)count(jobsToolsMatrix[i].begin(), jobsToolsMatrix[i].end(), 1);
            start_dist -= jobnum_to_free[i];
        }
        cromosome_free.assign(numJobs, 0);
	    filled_slots.assign(numJobs, 0);
	
        jobs_pairs_intersect.assign(numJobs, vec<int>(numJobs, 0));
        for (int i = 0; i < numJobs; i++) {
            for (int j = 0; j < numJobs; j++) {
                jobs_pairs_intersect[i][j] = (int)__builtin_popcountll(jobnum_to_bytes[i] & jobnum_to_bytes[j]);
            }
        }

        jobnum_to_tools_list.assign(numJobs, vec<int>(0));
        cromosome_jobnum_to_tools_list.assign(numJobs, vec<int>(0));
        for (int i = 0; i < numJobs; i++) {
            for (int j = 0; j < numTools; j++) {
                if (jobsToolsMatrix[i][j] == 1) {
                    jobnum_to_tools_list[i].push_back(j);
                    cromosome_jobnum_to_tools_list[i].push_back(j);
                }

            }
        }
        last_seen_tool.assign(numTools, -1);
        usedToFullMag.assign(numTools, 0);
        MAG_SETS.assign(numTools, unordered_set<int>());
        T_SETS.assign(numTools, unordered_set<int>());
        for (int i = 0; i < numJobs; i++) {
            for (auto& tool : jobnum_to_tools_list[i]) {
                T_SETS[i].insert(tool);
            }
        }
        empty_vec.assign(0,0);
        tool_in_i_not_in_prev_start.assign(numJobs, vec<int>(0));
        tool_in_i_not_in_prev = tool_in_i_not_in_prev_start;
        MAGAZINE = tool_in_i_not_in_prev_start;
        H_1_MAG = tool_in_i_not_in_prev_start;
        usedToFullMag.assign(numTools, 0);
        MN.assign(numJobs, vec<int>(numTools, 0));
    }
    
};


unsigned int(*func_to_test)(Data&, vec<int>&);
string func_to_test_name;

unsigned int KTNS(Data& data, vec<int>& chromosome) {
    int jump = -1;

    for (int i = (data.numJobs - 1); i >= 0; i--) {
        for (unsigned int j = 0; j < data.numTools; j++) {
            if (i != jump && data.jobsToolsMatrix[chromosome[i]][j] == 1) {
                data.L[j][i] = (unsigned int)i;
            }
            else if (i < (int)data.numJobs - 1) {
                data.L[j][i] = data.L[j][i + 1];
            }
            else {
                data.L[j][i] = data.numJobs;
            }
            data.used[j] = false;
        }
    }



    unsigned int switches = 0;
    unsigned int capacity = 0;
    unsigned int tool = 0;
    double minVal;

    for (unsigned int i = 0; i < data.numTools; i++) {
        if (data.L[i][0] == 0) {
            data.W_n[i] = 1;
            data.used[i] = true;
            capacity++;
        }
        else {
            data.W_n[i] = 0;
        }
    }

    while (capacity < data.magCapacity) {
        minVal = numeric_limits<double>::infinity();
        for (unsigned int i = 0; i < data.numTools; i++) {
            if (!data.used[i] && (data.L[i][0] < minVal)) {
                tool = i;
                minVal = data.L[i][0];
            }
        }
        data.used[tool] = true;
        data.W_n[tool] = 1;
        capacity++;
    }

    data.loadedMatrix[0] = data.W_n;

    unsigned int maxVal;
    for (unsigned int n = 1; n < data.numJobs; n++) {
        for (unsigned int i = 0; i < data.numTools; i++) {
            if (data.W_n[i] != 1 && data.L[i][n] == n) {
                data.W_n[i] = 1;
                capacity++;
            }
        }
        while (capacity > data.magCapacity) {
            maxVal = n;
            for (unsigned int i = 0; i < data.numTools; i++) {
                if (data.W_n[i] == 1 && data.L[i][n] > maxVal) {
                    tool = i;
                    maxVal = data.L[i][n];
                }
            }
            data.W_n[tool] = 0;
            capacity--;
            switches++;
        }
        data.loadedMatrix[n] = data.W_n;
    }
    return switches;
}



unsigned int MPCA_bitwise(Data& data, vec<int>& chromosome) {
    unsigned long long end_tools, intersect;

    for (int i = 0; i < data.numJobs; i++) {
        data.cromosome_bytes[i] = data.jobnum_to_bytes[chromosome[i]];
        data.cromosome_free[i] = data.jobnum_to_free[chromosome[i]];
    }

    int dist = data.start_dist;
    for (int i = 0; i < data.numJobs - 1; i++) {
        dist -= data.jobs_pairs_intersect[chromosome[i]][chromosome[i + 1]];
    }
    int fullmag = -1;


    for (int end = 1; end < data.numJobs; end++) {
        int end_minus_1 = end - 1;
        if (data.cromosome_free[end_minus_1] == 0) { fullmag = end_minus_1; }
        end_tools = data.cromosome_bytes[end] & (~data.cromosome_bytes[end_minus_1]);
        int min_free = 99999;
        for (int start = end - 2; start >= 0; start--) {
            if (fullmag > start || end_tools == 0LL) { break; }
            if (min_free > data.cromosome_free[start + 1]) { min_free = data.cromosome_free[start + 1]; }
            intersect = data.cromosome_bytes[start] & end_tools;
            if (intersect > 0) {
                end_tools = end_tools & (~intersect);
                int tubes_count = (int)__builtin_popcountll(intersect);
                if (tubes_count > min_free) { tubes_count = min_free; }
                if (tubes_count > 0) {
                    for (int j = start + 1; j < end; j++) {
                        data.cromosome_free[j] -= tubes_count;
                        if (data.cromosome_free[j] == 0) { fullmag = j; }
                    }
                    dist -= tubes_count;
                    min_free -= tubes_count;
                }
            }
        }
    }

    return dist;
}


unsigned int MPCA(Data& data, vec<int>& chromosome) {
    int i,end, pipes_count = 0, last_full = -1;
    data.last_seen_tool.assign(data.numTools, -2);

    for (end = 0; end < data.numJobs; end++) {
        for (int tool : data.jobnum_to_tools_list[chromosome[end]]) {
            if (last_full <= data.last_seen_tool[tool]) {
                for (i = end - 1; i > data.last_seen_tool[tool]; i--) {
                    if ((++data.filled_slots[i]) == data.magCapacity) {
                        last_full = i;
                        break;
                    }
                }
                pipes_count++;
            }
            data.last_seen_tool[tool] = end;
        }
        data.filled_slots[end] = data.job_capacity[chromosome[end]];
        if (data.filled_slots[end] == data.magCapacity) { last_full = end; }
    }
    return data.start_dist - pipes_count;
}

unsigned long long remove_ones(unsigned long long x, int ones_to_remove) {
    for (int i = 0; i < ones_to_remove; i++) {
        x &= x - 1;
    }
    return x;
}


unsigned int KTNS_bitwise(Data& data, vec<int>& chromosome) {
    int remove_count, to_keep_count, intersect_count, to_add_count;
    unsigned long long now, candidates, intersect;
    vec<unsigned long long> M(data.numJobs, 0);
    for (int i = 1; i < data.numJobs; i++) {
        data.cromosome_bytes[i] = data.jobnum_to_bytes[chromosome[i]];
        M[i] = data.jobnum_to_bytes[chromosome[i]];
    }
    candidates = 0;
    int candidates_count = 0;
    for (int i = 0; i < data.numJobs; i++) {
        to_add_count = (int)__builtin_popcountll(data.cromosome_bytes[i] & ~candidates);
        if (candidates_count + to_add_count < data.magCapacity) {
            candidates |= data.cromosome_bytes[i];
            candidates_count += to_add_count;
        }
        else {
            candidates |= remove_ones(data.cromosome_bytes[i] & ~candidates, candidates_count + to_add_count - data.magCapacity);
            M[0] = candidates;
            break;
        }
    }

    for (int s = 1; s < data.numJobs; s++) {
        now = M[s - 1] | data.cromosome_bytes[s];
        remove_count = (int)__builtin_popcountll(now) - data.magCapacity;
        if (remove_count == 0) {
            M[s] = now;
            continue;
        }
        else if (remove_count < 0) { cout << "error remove_count < 0 \n"; throw 1; }
        candidates = now & ~data.cromosome_bytes[s];
        to_keep_count = (int)__builtin_popcountll(candidates) - remove_count;
        for (int e = s + 1; e <= data.numJobs; e++) {
            if (e == data.numJobs) {
                M[s] |= remove_ones(candidates, (int)__builtin_popcountll(candidates) - to_keep_count);
                break;
            }
            intersect = candidates & data.cromosome_bytes[e];
            intersect_count = (int)__builtin_popcountll(intersect);
            if (intersect_count < to_keep_count) {
                to_keep_count -= intersect_count;
                M[s] |= intersect;
                candidates = candidates & ~intersect;
            }
            else {
                M[s] |= remove_ones(intersect, intersect_count - to_keep_count);
                break;
            }
        }
    }

    int sw = 0;
    for (int i = 0; i < data.numJobs - 1; i++) {
        sw += (int)__builtin_popcountll(M[i + 1] & ~M[i]);
    }
    return sw;
}


double ktns_time(string dataset_path, string prem_file_path) {
    int n = -1, m = -1, C = -1;
    std::cout << dataset_path << "<<<<<\n";
    ifstream fin(dataset_path);
    assert(fin && "dataset file not found");

    fin >> n >> m >> C;
    cout << "alg=" << func_to_test_name << ", n=" << n << ", m=" << m << ", C=" << C << ", file=" << dataset_path << ".\n";
    vec<vec<int> > jobsToolsMatrix(n, vec<int>(m, -1));
    int a;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fin >> jobsToolsMatrix[j][i];
        }
    }
    fin.close();

    

    Data data(n, m, C, jobsToolsMatrix);

    ifstream fin2(prem_file_path);
    assert(fin2 && "permutation file not found");
    int n_premuts;
    fin2 >> n_premuts;
    vec<vec<int> > quest(n_premuts, vec<int>(n, -1));

    for (int i = 0; i < n_premuts; i++) {
        for (int j = 0; j < n; j++) {
            fin2 >> quest[i][j];
        }
    }
    fin2.close();

    int result;
    clock_t start = clock();
    for (int i = 0; i < n_premuts; i++) {
        result = func_to_test(data, quest[i]);
    }
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    return seconds;
}


double dataset_time(string folde_name, string dataset) {
    string dir = ".\\";
    dir += folde_name + "\\";
    string ABCD = string() + dataset[0];
    string capacity_banchmark_index = string() + dataset[1];


    double summary_time = 0;
    for (int instance_num = 1; instance_num <= 10; instance_num++) {
        string dataset_path = dir + "Tabela" + capacity_banchmark_index + "\\dat" + ABCD + to_string(instance_num);
        string prem_file_path = dir + "Tabela" + capacity_banchmark_index + "\\permuts" + ABCD;
        summary_time += ktns_time(dataset_path, prem_file_path);
    }
    return summary_time;
}



void ktns_tca_time() {
    vec<string> datasets;
    for(auto v : "ABCD")
        for(auto u : "1234")
            datasets.push_back(string() + v + u);
    unordered_map<string, double> time_dict;

    
    ofstream output_file;
    output_file.open("results.txt", ios::app);

    output_file << "algorithm: " << func_to_test_name << "\n";
    for (string dataset : datasets) {
        double summary_time = dataset_time("Catanzaro", dataset);
        time_dict[dataset] = summary_time;
        cout << "'" << dataset << "': " << time_dict[dataset] << ",\n";
        output_file << "'" << dataset << "': " << time_dict[dataset] << ",\n";
    }
    output_file.close();
}



int main() {
    
    
    std::cout << std::fixed << std::setprecision(5);
    std::ofstream output_file; output_file.open(OUTPUT_FILE, std::ofstream::out | std::ofstream::trunc); output_file.close();

    func_to_test_name = "MPCA_bitwise";
    func_to_test = MPCA_bitwise;
    ktns_tca_time();

    func_to_test_name = "KTNS_bitwise";
    func_to_test = KTNS_bitwise;
    ktns_tca_time();

    func_to_test_name = "MPCA";
    func_to_test = MPCA;
    ktns_tca_time();

    func_to_test_name = "KTNS";
    func_to_test = KTNS;
    ktns_tca_time();


    return 0;
}




