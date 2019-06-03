#include <random>
#include <omp.h>
#include <iostream>

static std::vector<std::mt19937_64> engines;
static std::vector<std::uniform_real_distribution<double>> distribs;

extern "C" {

void cxx11_init_random(int* seeds)
{
  for (int i=0; i<omp_get_max_threads(); ++i) {
    engines.emplace_back(seeds[i]);
    distribs.emplace_back(0.0,1.0);
  }
}

double cxx11_random_number()
{
  int tid = omp_get_thread_num();
  return (distribs[tid])(engines[tid]);
}

}
