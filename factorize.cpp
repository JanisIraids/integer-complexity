#include "factorize.h"

#include <iostream>
#include <fstream>
#include <msieve.h>


// Some copy-paste code from msieve 1.50
// All credit to the original authors
/*--------------------------------------------------------------------*/
int get_random_seeds(uint32 *seed1, uint32 *seed2) {

  uint32 tmp_seed1, tmp_seed2;

  /* In a multithreaded program, every msieve object
     should have two unique, non-correlated seeds
     chosen for it */

  std::ifstream rand_device("/dev/urandom", std::ios::binary);
  int rv = -1;

  if (rand_device.is_open()) {

    /* Yay! Cryptographic-quality nondeterministic randomness! */

    rand_device.read(reinterpret_cast<char*>(&tmp_seed1), sizeof(uint32));
    rand_device.read(reinterpret_cast<char*>(&tmp_seed2), sizeof(uint32));
    rand_device.close();
  }

  /* The final seeds are the result of a multiplicative
     hash of the initial seeds */

  (*seed1) = tmp_seed1 * ((uint32)40499 * 65543);
  (*seed2) = tmp_seed2 * ((uint32)40499 * 65543);
  return rv;
}


void factor(mpz_t number, mpz_t* factors)
{
  msieve_obj *g_curr_factorization = NULL;

  static uint32 seed1 = 0, seed2 = 0;
  char *savefile_name = NULL;
  char *logfile_name = NULL;
  char *nfs_fbfile_name = NULL;
  uint32 flags;
  uint32 max_relations = 0;
  enum cpu_type cpu;
  uint32 cache_size1; 
  uint32 cache_size2; 
  uint32 num_threads = 0;
  uint32 which_gpu = 0;
		
  get_cache_sizes(&cache_size1, &cache_size2);
  cpu = get_cpu_type();

  flags = MSIEVE_FLAG_USE_LOGFILE;
  flags |= MSIEVE_FLAG_DEEP_ECM;
  flags &= ~(MSIEVE_FLAG_USE_LOGFILE |
	     MSIEVE_FLAG_LOG_TO_STDOUT);

  get_random_seeds(&seed1, &seed2);

  char int_start[2048];
  gmp_sprintf(int_start, "%Zd", number);
  msieve_obj *obj;
  msieve_factor *mfactor;

  /* point to the start of the integer or expression;
     if the start point indicates no integer is present,
     don't try to factor it :) */

  g_curr_factorization = msieve_obj_new(int_start, flags,
					savefile_name, logfile_name,
					nfs_fbfile_name,
					seed1, seed2, max_relations,
					cpu,
					cache_size1, cache_size2,
					num_threads, which_gpu,
				        NULL);
  if (g_curr_factorization == NULL) {
    std::cout << "factoring initialization failed" << std::endl;
    return;
  }

  msieve_run(g_curr_factorization);

  if (!(g_curr_factorization->flags & MSIEVE_FLAG_FACTORIZATION_DONE)) {
    std::cout << "current factorization was interrupted" << std::endl;
    exit(0);
  }

  /* If no logging is specified, at least print out the
     factors that were found */

  mfactor = g_curr_factorization->factors;
  
  int factornum = 0;
  mpz_set_d(factors[factornum++], 1);
  while (mfactor != NULL) {
    gmp_sscanf(mfactor->number, "%Zd", factors[factornum++]);
    mfactor = mfactor->next;
  }
  mpz_set_d(factors[factornum++], 0);

  /* save the current value of the random seeds, so that
     the next factorization will pick up the pseudorandom
     sequence where this factorization left off */

  seed1 = g_curr_factorization->seed1;
  seed2 = g_curr_factorization->seed2;

  /* free the current factorization struct. The following
     avoids a race condition in the signal handler */

  obj = g_curr_factorization;
  g_curr_factorization = NULL;
  if (obj)
    msieve_obj_free(obj);
  
}

void factor(mpz_class number, std::vector<mpz_class>& factors)
{
  msieve_obj *g_curr_factorization = NULL;

  static uint32 seed1 = 0, seed2 = 0;
  char *savefile_name = NULL;
  char *logfile_name = NULL;
  char *nfs_fbfile_name = NULL;
  uint32 flags;
  uint32 max_relations = 0;
  enum cpu_type cpu;
  uint32 cache_size1; 
  uint32 cache_size2; 
  uint32 num_threads = 0;
  uint32 which_gpu = 0;
		
  get_cache_sizes(&cache_size1, &cache_size2);
  cpu = get_cpu_type();

  flags = MSIEVE_FLAG_USE_LOGFILE;
  flags |= MSIEVE_FLAG_DEEP_ECM;
  flags &= ~(MSIEVE_FLAG_USE_LOGFILE |
	     MSIEVE_FLAG_LOG_TO_STDOUT);

  get_random_seeds(&seed1, &seed2);

  char int_start[2048];
  gmp_sprintf(int_start, "%Zd", number.get_mpz_t());
  msieve_obj *obj;
  msieve_factor *mfactor;

  /* point to the start of the integer or expression;
     if the start point indicates no integer is present,
     don't try to factor it :) */

  g_curr_factorization = msieve_obj_new(int_start, flags,
					savefile_name, logfile_name,
					nfs_fbfile_name,
					seed1, seed2, max_relations,
					cpu,
					cache_size1, cache_size2,
					num_threads, which_gpu,
				        NULL);
  if (g_curr_factorization == NULL) {
    std::cout << "factoring initialization failed" << std::endl;
    return;
  }

  msieve_run(g_curr_factorization);

  if (!(g_curr_factorization->flags & MSIEVE_FLAG_FACTORIZATION_DONE)) {
    std::cout << "current factorization was interrupted" << std::endl;
    exit(0);
  }

  /* If no logging is specified, at least print out the
     factors that were found */

  mfactor = g_curr_factorization->factors;
  
  factors.push_back(mpz_class(1));
  while (mfactor != NULL) {
    mpz_class f;
    gmp_sscanf(mfactor->number, "%Zd", f.get_mpz_t());
    factors.push_back(f);
    mfactor = mfactor->next;
  }
  factors.push_back(mpz_class(0));

  /* save the current value of the random seeds, so that
     the next factorization will pick up the pseudorandom
     sequence where this factorization left off */

  seed1 = g_curr_factorization->seed1;
  seed2 = g_curr_factorization->seed2;

  /* free the current factorization struct. The following
     avoids a race condition in the signal handler */

  obj = g_curr_factorization;
  g_curr_factorization = NULL;
  if (obj)
    msieve_obj_free(obj);
  
}
