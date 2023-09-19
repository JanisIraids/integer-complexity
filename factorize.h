#pragma once

#include <gmp.h>
#include <gmpxx.h>
#include <vector>

void factor(mpz_t number, mpz_t* factors);
void factor(mpz_class number, std::vector<mpz_class>& factors);
