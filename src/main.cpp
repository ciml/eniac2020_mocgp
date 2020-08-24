#include "templates/cgp.tpp"
#include <ctime>

int main(int argc, char const *argv[])
{
    std::clock_t c_start = std::clock();

    int seed = atoi(argv[1]);
    srand(seed);
    std::cout << "Seed: " << seed << std::endl;
    bdd_init(5000000, 50000);

    Cgp<int, bdd> cgp(argc, argv);

    #ifdef MO
    if(cgp.evolve()){
        cgp.optimize_MO();
    }
    #endif
    #ifndef MO
    if(cgp.evolve()){
        cgp.optimize();
    }
    #endif

    bdd_done();

    std::clock_t c_end = std::clock();
    long double time_elapsed_s = (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "CPU time used: " << time_elapsed_s << "s\n";

    return 0;
}