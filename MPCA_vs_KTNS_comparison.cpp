#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <bitset>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include <iomanip>

const std::string OUTPUT_FILE = "results.txt";
const int MAX_INT = (int)1e9;

struct TLP_Data
{
    /*
        Container for a TLP instance tlp_data.
    */

    // TLP input data
    int numJobs, numTools, magCapacity;
    std::vector<std::vector<int>> jobsToolsMatrix;

    // auxiliary structures
    std::vector<std::vector<int>> L;
    std::vector<std::vector<int>> loadedMatrix;
    std::vector<int> used;
    std::vector<int> W_n;
    int number_of_tools_with_repetitions;
    std::vector<int> filled_slots;
    std::vector<int> job_capacity;
    std::vector<int> job_seqence;
    std::vector<int> job_seqence_free;
    std::vector<unsigned long long> job_seqence_bytes;
    std::vector<int> jobnum_to_free;
    std::vector<unsigned long long> jobnum_to_bytes;
    std::vector<std::vector<int>> jobnum_to_tools_list;
    std::vector<int> last_seen_tool;

    TLP_Data(int numJobs_arg, int numTools_arg, int magCapacity_arg, std::vector<std::vector<int>> &jobsToolsMatrix_arg)
        : numJobs(numJobs_arg),
          numTools(numTools_arg),
          magCapacity(magCapacity_arg),
          jobsToolsMatrix(jobsToolsMatrix_arg)
    {

        used.assign(numTools, 0);
        L.assign(numTools, std::vector<int>(numJobs, -1));
        W_n.assign(numTools, -1);
        loadedMatrix.assign(numJobs, std::vector<int>(numTools, -1));

        jobnum_to_bytes.assign(numJobs, 0LL);
        for (int i = 0; i < numJobs; i++)
        {
            std::string str = "";
            for (int v : jobsToolsMatrix[i])
            {
                str += std::to_string(v);
            }
            jobnum_to_bytes[i] = std::bitset<64>(str).to_ullong();
        }
        job_seqence_bytes.assign(numJobs, 0LL);

        number_of_tools_with_repetitions = magCapacity * (numJobs - 1);
        jobnum_to_free.assign(numJobs, 0);
        job_capacity.assign(numJobs, 0);
        for (int i = 0; i < numJobs; i++)
        {
            job_capacity[i] = (int)count(jobsToolsMatrix[i].begin(), jobsToolsMatrix[i].end(), 1);
            jobnum_to_free[i] = magCapacity - (int)count(jobsToolsMatrix[i].begin(), jobsToolsMatrix[i].end(), 1);
            number_of_tools_with_repetitions -= jobnum_to_free[i];
        }
        job_seqence_free.assign(numJobs, 0);
        filled_slots.assign(numJobs, 0);

        jobnum_to_tools_list.assign(numJobs, std::vector<int>(0));
        for (int i = 0; i < numJobs; i++)
        {
            for (int j = 0; j < numTools; j++)
            {
                if (jobsToolsMatrix[i][j] == 1)
                {
                    jobnum_to_tools_list[i].push_back(j);
                }
            }
        }
        last_seen_tool.assign(numTools, -1);
    }
};

unsigned int KTNS(TLP_Data &tlp_data)
{
    /*
        Main implementation of Keep Tool Needed Soonest algortihm.
            -proposed by Christopher Tang, Eric Denardo, see https://doi.org/10.1287/opre.36.5.767
            -implemented by Jordana Mecler, Anand Subramanian, Thibaut Vidal, see https://doi.org/10.1016/j.cor.2020.105153
            -available at https://github.com/jordanamecler/HGS-SSP
        Returns the minimal number of switches for a given instance of the Tool Loading Problem.
    */
    int jump = -1;
    for (int i = (tlp_data.numJobs - 1); i >= 0; i--)
    {
        for (unsigned int j = 0; j < tlp_data.numTools; j++)
        {
            if (i != jump && tlp_data.jobsToolsMatrix[tlp_data.job_seqence[i]][j] == 1)
            {
                tlp_data.L[j][i] = (unsigned int)i;
            }
            else if (i < (int)tlp_data.numJobs - 1)
            {
                tlp_data.L[j][i] = tlp_data.L[j][i + 1];
            }
            else
            {
                tlp_data.L[j][i] = tlp_data.numJobs;
            }
            tlp_data.used[j] = false;
        }
    }

    unsigned int switches = 0;
    unsigned int capacity = 0;
    unsigned int tool = 0;
    double minVal;

    for (unsigned int i = 0; i < tlp_data.numTools; i++)
    {
        if (tlp_data.L[i][0] == 0)
        {
            tlp_data.W_n[i] = 1;
            tlp_data.used[i] = true;
            capacity++;
        }
        else
        {
            tlp_data.W_n[i] = 0;
        }
    }

    while (capacity < tlp_data.magCapacity)
    {
        minVal = std::numeric_limits<double>::infinity();
        for (unsigned int i = 0; i < tlp_data.numTools; i++)
        {
            if (!tlp_data.used[i] && (tlp_data.L[i][0] < minVal))
            {
                tool = i;
                minVal = tlp_data.L[i][0];
            }
        }
        tlp_data.used[tool] = true;
        tlp_data.W_n[tool] = 1;
        capacity++;
    }

    tlp_data.loadedMatrix[0] = tlp_data.W_n;

    unsigned int maxVal;
    for (unsigned int n = 1; n < tlp_data.numJobs; n++)
    {
        for (unsigned int i = 0; i < tlp_data.numTools; i++)
        {
            if (tlp_data.W_n[i] != 1 && tlp_data.L[i][n] == n)
            {
                tlp_data.W_n[i] = 1;
                capacity++;
            }
        }
        while (capacity > tlp_data.magCapacity)
        {
            maxVal = n;
            for (unsigned int i = 0; i < tlp_data.numTools; i++)
            {
                if (tlp_data.W_n[i] == 1 && tlp_data.L[i][n] > maxVal)
                {
                    tool = i;
                    maxVal = tlp_data.L[i][n];
                }
            }
            tlp_data.W_n[tool] = 0;
            capacity--;
            switches++;
        }
        tlp_data.loadedMatrix[n] = tlp_data.W_n;
    }
    return switches;
}

unsigned long long remove_ones(unsigned long long x, int ones_to_remove)
{
    /*
        Auxiliary function for KTNS_bitwise.
    */
    for (int i = 0; i < ones_to_remove; i++)
    {
        x &= x - 1;
    }
    return x;
}

unsigned int KTNS_bitwise(TLP_Data &tlp_data)
{
    /*
        Bitwise implementation of Keep Tool Needed Soonest algorithm (proposed by the authors).
        Returns the minimal number of switches for a given instance of the Tool Loading Problem.
        Implementation limitations: number of tools is at most 64 (long long is used to represent magazine states).
    */
    assert(tlp_data.numTools <= 64);
    int remove_count, to_keep_count, intersect_count, to_add_count;
    unsigned long long now, candidates, intersect;
    std::vector<unsigned long long> magazine_bit(tlp_data.numJobs, 0);
    for (int i = 1; i < tlp_data.numJobs; i++)
    {
        tlp_data.job_seqence_bytes[i] = tlp_data.jobnum_to_bytes[tlp_data.job_seqence[i]];
        magazine_bit[i] = tlp_data.jobnum_to_bytes[tlp_data.job_seqence[i]];
    }
    candidates = 0;
    int candidates_count = 0;
    for (int i = 0; i < tlp_data.numJobs; i++)
    {
        to_add_count = (int)__builtin_popcountll(tlp_data.job_seqence_bytes[i] & ~candidates);
        if (candidates_count + to_add_count < tlp_data.magCapacity)
        {
            candidates |= tlp_data.job_seqence_bytes[i];
            candidates_count += to_add_count;
        }
        else
        {
            candidates |= remove_ones(tlp_data.job_seqence_bytes[i] & ~candidates, candidates_count + to_add_count - tlp_data.magCapacity);
            magazine_bit[0] = candidates;
            break;
        }
    }

    for (int s = 1; s < tlp_data.numJobs; s++)
    {
        now = magazine_bit[s - 1] | tlp_data.job_seqence_bytes[s];
        remove_count = (int)__builtin_popcountll(now) - tlp_data.magCapacity;
        if (remove_count == 0)
        {
            magazine_bit[s] = now;
            continue;
        }

        candidates = now & ~tlp_data.job_seqence_bytes[s];
        to_keep_count = (int)__builtin_popcountll(candidates) - remove_count;
        for (int e = s + 1; e <= tlp_data.numJobs; e++)
        {
            if (e == tlp_data.numJobs)
            {
                magazine_bit[s] |= remove_ones(candidates, (int)__builtin_popcountll(candidates) - to_keep_count);
                break;
            }
            intersect = candidates & tlp_data.job_seqence_bytes[e];
            intersect_count = (int)__builtin_popcountll(intersect);
            if (intersect_count < to_keep_count)
            {
                to_keep_count -= intersect_count;
                magazine_bit[s] |= intersect;
                candidates = candidates & ~intersect;
            }
            else
            {
                magazine_bit[s] |= remove_ones(intersect, intersect_count - to_keep_count);
                break;
            }
        }
    }

    int switches = 0;
    for (int i = 0; i < tlp_data.numJobs - 1; i++)
    {
        switches += (int)__builtin_popcountll(magazine_bit[i + 1] & ~magazine_bit[i]);
    }
    return switches;
}

unsigned int MPCA(TLP_Data &tlp_data)
{
    /*
        Main implementation of Max Pipe Construction Algorithm (proposed by the authors).
        Returns the minimal number of switches for a given instance of the Tool Loading Problem.
    */
    int i, end, pipes_count = 0, last_full = -1;
    tlp_data.last_seen_tool.assign(tlp_data.numTools, -2);

    for (end = 0; end < tlp_data.numJobs; end++)
    {
        for (int tool : tlp_data.jobnum_to_tools_list[tlp_data.job_seqence[end]])
        {
            if (last_full <= tlp_data.last_seen_tool[tool])
            {
                for (i = end - 1; i > tlp_data.last_seen_tool[tool]; i--)
                {
                    if ((++tlp_data.filled_slots[i]) == tlp_data.magCapacity)
                    {
                        last_full = i;
                        break;
                    }
                }
                pipes_count++;
            }
            tlp_data.last_seen_tool[tool] = end;
        }
        tlp_data.filled_slots[end] = tlp_data.job_capacity[tlp_data.job_seqence[end]];
        if (tlp_data.filled_slots[end] == tlp_data.magCapacity)
        {
            last_full = end;
        }
    }
    int switches = tlp_data.number_of_tools_with_repetitions - pipes_count;
    return switches;
}

unsigned int MPCA_bitwise(TLP_Data &tlp_data)
{
    /*
        Bitwise implementation of Max Pipe Construction Algorithm (proposed by the authors).
        Returns the minimal number of switches for a given instance of the Tool Loading Problem.
        Implementation limitations: number of tools is at most 64 (long long is used to represent magazine states).
    */
    assert(tlp_data.numTools <= 64);
    unsigned long long end_tools, intersect;

    for (int i = 0; i < tlp_data.numJobs; i++)
    {
        tlp_data.job_seqence_bytes[i] = tlp_data.jobnum_to_bytes[tlp_data.job_seqence[i]];
        tlp_data.job_seqence_free[i] = tlp_data.jobnum_to_free[tlp_data.job_seqence[i]];
    }

    int max_switches = tlp_data.number_of_tools_with_repetitions;
    for (int i = 0; i < tlp_data.numJobs - 1; i++)
    {
        max_switches -= (int)__builtin_popcountll(tlp_data.jobnum_to_bytes[tlp_data.job_seqence[i]] & tlp_data.jobnum_to_bytes[tlp_data.job_seqence[i + 1]]);
    }
    int total_pipes_count = 0;

    int fullmag = -1;
    for (int end = 1; end < tlp_data.numJobs; end++)
    {
        if (tlp_data.job_seqence_free[end - 1] == 0)
        {
            fullmag = end - 1;
        }
        end_tools = tlp_data.job_seqence_bytes[end] & (~tlp_data.job_seqence_bytes[end - 1]);
        int min_empty = MAX_INT;
        for (int start = end - 2; start >= 0; start--)
        {
            if (fullmag > start || end_tools == 0LL)
            {
                break;
            }
            if (min_empty > tlp_data.job_seqence_free[start + 1])
            {
                min_empty = tlp_data.job_seqence_free[start + 1];
            }
            intersect = tlp_data.job_seqence_bytes[start] & end_tools;
            if (intersect > 0)
            {
                end_tools = end_tools & (~intersect);
                int pipes_count = (int)__builtin_popcountll(intersect);
                if (pipes_count > min_empty)
                {
                    pipes_count = min_empty;
                }
                if (pipes_count > 0)
                {
                    for (int j = start + 1; j < end; j++)
                    {
                        tlp_data.job_seqence_free[j] -= pipes_count;
                        if (tlp_data.job_seqence_free[j] == 0)
                        {
                            fullmag = j;
                        }
                    }
                    total_pipes_count += pipes_count;
                    min_empty -= pipes_count;
                }
            }
        }
    }
    int switches = max_switches - total_pipes_count;
    return switches;
}

void algorithms_verefication()
{
    /*
        Verifies KTNS, MPCA, KTNS_bitwise, MPCA_bitwise by running
        the algorithms on random Tool Loading Problem insances and
        checking whether they give the same result.
        A random TLP instance is selected 100 times
        and 1000 job sequences are generated for each instance.
        In summary we have 10^5 tests.
    */
    std::cout << "Verifying KTNS, MPCA, KTNS_bitwise, MPCA_bitwise...\n";
    for (int rnd_paths = 0; rnd_paths < 100; rnd_paths++)
    {
        // select a random Tool Loading Problem instance
        std::string capacity_banchmark_index = std::string() + (char)('1' + std::rand() % 4);
        std::string ABCD = std::string() + (char)('A' + std::rand() % 4);
        int instance_num = std::rand() % 10 + 1;
        std::string rand_TLP_instance_path = "./Catanzaro/Tabela" + capacity_banchmark_index + "/dat" + ABCD + std::to_string(instance_num);

        // read data from file that was selected randomly
        int n = -1, m = -1, C = -1;
        std::ifstream fin(rand_TLP_instance_path);
        fin >> n >> m >> C;
        assert(n > 0);
        assert(m > 0);
        assert(C > 0);
        std::vector<std::vector<int>> jobsToolsMatrix(n, std::vector<int>(m, -1));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                fin >> jobsToolsMatrix[j][i];
        fin.close();
        TLP_Data tlp_data(n, m, C, jobsToolsMatrix);

        std::vector<int> rand_job_seqence(n);
        for (int i = 0; i < n; i++)
            rand_job_seqence[i] = i;

        for (int i = 0; i < 1000; i++)
        {
            // generate a random sequence of jobs
            std::random_shuffle(rand_job_seqence.begin(), rand_job_seqence.end());
            tlp_data.job_seqence = rand_job_seqence;
            // run algorithms
            int MPCA_bitwise_result = MPCA_bitwise(tlp_data);
            int KTNS_bitwise_result = KTNS_bitwise(tlp_data);
            int MPCA_result = MPCA(tlp_data);
            int KTNS_result = KTNS(tlp_data);

            // check whether algorithms give the same result
            assert(MPCA_bitwise_result == KTNS_result);
            assert(KTNS_bitwise_result == KTNS_result);
            assert(MPCA_result == KTNS_result);
        }
    }
    std::cout << "KTNS, MPCA, KTNS_bitwise, MPCA_bitwise passed verification and return same outputs for 10^5 TLP instances.\n"
              << "\n";
}

struct Experiments
{
    /*
        Run computational experiments for the given algorithm, i.e. enumerate problem instances,
        read data from files, execute algorithm, measure the time spent on calculations,
        write cumputational results to file results.txt.
    */
    unsigned int (*alg_for_test)(TLP_Data &);
    std::string alg_for_test_name;

    Experiments(unsigned int (*alg_for_test_arg)(TLP_Data &), std::string alg_for_test_name_arg)
    {
        this->alg_for_test = alg_for_test_arg;
        this->alg_for_test_name = alg_for_test_name_arg;
    }
    double instance_time(std::string TLP_instance_path)
    {
        /*
            For the given TLP instance runs an algorithm for 10^5
            generated random job sequences and returns total time spent.
        */

        // read TLP data from file
        int n = -1, m = -1, C = -1;
        std::ifstream fin(TLP_instance_path);
        fin >> n >> m >> C;
        assert(n > 0);
        assert(m > 0);
        assert(C > 0);
        std::vector<std::vector<int>> jobsToolsMatrix(n, std::vector<int>(m, -1));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                fin >> jobsToolsMatrix[j][i];
        fin.close();
        TLP_Data tlp_data(n, m, C, jobsToolsMatrix);

        // generate 10^5 random job sequences
        int NUMBER_OF_JOB_SEQUENCES = 100000;
        std::vector<std::vector<int>> quest(NUMBER_OF_JOB_SEQUENCES);
        std::vector<int> rand_job_seqence(n);
        for (int i = 0; i < n; i++)
            rand_job_seqence[i] = i;
        for (int i = 0; i < NUMBER_OF_JOB_SEQUENCES; i++)
        {
            quest[i] = rand_job_seqence;
            std::random_shuffle(rand_job_seqence.begin(), rand_job_seqence.end());
        }

        // execute an algorithm for each job sequence
        std::cout << "alg=" << alg_for_test_name << ", n=" << n << ", m=" << m << ", C=" << C << ", file=" << TLP_instance_path << ".\n";
        clock_t start = clock();
        for (int i = 0; i < NUMBER_OF_JOB_SEQUENCES; i++)
        {
            tlp_data.job_seqence = quest[i];
            alg_for_test(tlp_data);
        }
        clock_t end = clock();
        double seconds = (double)(end - start) / CLOCKS_PER_SEC;
        return seconds;
    }

    double dataset_time(std::string dataset)
    {
        /*
            Triggers tests for 10 instances of the given dataset and returns the total time spent.
        */
        std::string ABCD = std::string() + dataset[0];
        std::string capacity_banchmark_index = std::string() + dataset[1];

        double total_time = 0;
        for (int instance_num = 1; instance_num <= 10; instance_num++)
        {
            std::string TLP_instance_path = "./Catanzaro/Tabela" + capacity_banchmark_index + "/dat" + ABCD + std::to_string(instance_num);
            total_time += instance_time(TLP_instance_path);
        }
        return total_time;
    }

    void run_tests()
    {
        /*
            Triggers running tests and writes test results to a file results.txt.
        */

        // list of datasets
        std::vector<std::string> datasets;
        datasets.push_back("A1");
        datasets.push_back("A2");
        datasets.push_back("A3");
        datasets.push_back("A4");
        datasets.push_back("B1");
        datasets.push_back("B2");
        datasets.push_back("B3");
        datasets.push_back("B4");
        datasets.push_back("C1");
        datasets.push_back("C2");
        datasets.push_back("C3");
        datasets.push_back("C4");
        datasets.push_back("D1");
        datasets.push_back("D2");
        datasets.push_back("D3");
        datasets.push_back("D4");

        std::ofstream output_file;
        output_file.open("results.txt", std::ios::app);
        output_file << "algorithm: " << alg_for_test_name << "\n";
        output_file.close();

        for (std::string dataset_name : datasets)
        {
            // trigger tests
            double summary_time = dataset_time(dataset_name);
            std::cout << "'" << dataset_name << "': " << summary_time << ",\n";

            // write results
            output_file.open("results.txt", std::ios::app);
            output_file << "'" << dataset_name << "': " << summary_time << ",\n";
            output_file.close();
        }
    }
};

int main()
{
    srand(time(NULL));
    std::cout << std::fixed << std::setprecision(3);
    std::ofstream output_file;
    output_file.open(OUTPUT_FILE, std::ofstream::out | std::ofstream::trunc);
    output_file.close();

    algorithms_verefication();

    Experiments test_MPCA_bitwise(MPCA_bitwise, "MPCA_bitwise");
    test_MPCA_bitwise.run_tests();

    Experiments test_KTNS_bitwise(KTNS_bitwise, "KTNS_bitwise");
    test_KTNS_bitwise.run_tests();

    Experiments test_MPCA(MPCA, "MPCA");
    test_MPCA.run_tests();

    Experiments test_KTNS(KTNS, "KTNS");
    test_KTNS.run_tests();

    return 0;
}