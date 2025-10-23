#include <iostream>
#include "RunProtal.h"
#include "Test.h"

//#include <htslib/sam.h>
//#include "Test.h"
#include "MetagenomeSimulator.h"

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);



    std::cout << "Running protal tests..." << std::endl;

    MetagenomeSimulator simulator { 
        "/home/joachim/Data/data/protal/datasets/random_forest_test_data/index_sub.txt", 
        "/home/joachim/Data/data/protal/datasets/random_forest_test_data/test_acc_species.txt" 
    };


    auto powerlaw1 = power_law_distribution(1000);
    std::cout << "Powerlaw samples:" << std::endl;
    for (size_t i = 0; i < powerlaw1.size(); ++i) {
        std::cout << powerlaw1[i] << ", ";
    }
    std::cout << std::endl << std::endl;    

    auto powerlaw2 = power_law_distribution(1000);
    std::cout << "Powerlaw samples:" << std::endl;
    for (size_t i = 0; i < powerlaw2.size(); ++i) {
        std::cout << powerlaw2[i] << ", ";
    }
    std::cout << std::endl << std::endl;    

    auto powerlaw3  = power_law_distribution(100);
    std::cout << "Powerlaw samples:" << std::endl;
    for (size_t i = 0; i < powerlaw3.size(); ++i) {
        std::cout << powerlaw3[i] << ", ";
    }
    std::cout << std::endl << std::endl;    


    auto powerlaw4 = power_law_distribution(10);
    std::cout << "Powerlaw samples:" << std::endl;
    for (size_t i = 0; i < powerlaw4.size(); ++i) {
        std::cout << powerlaw4[i] << ", ";
    }
    std::cout << std::endl << std::endl;    

    auto [counts, rel_abundances, sequencing_depths] = simulator.simulate(
        SimulationOptions{ 10'000'000, 0.01, 300, 150 }, 
        { 2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL,
          2500000UL + static_cast<unsigned long>(rand()) % 2500001UL }
    );

    for (size_t i = 0; i < counts.size(); ++i) {
        std::cout << "Genome " << i << ": "
                  << "Counts = " << counts[i] << ", "
                  << "Relative Abundance = " << rel_abundances[i] << ", "
                  << "Sequencing Depth = " << sequencing_depths[i] << std::endl;
    }

    exit(0);

    protal::Run(argc, argv);

    return 0;
}
