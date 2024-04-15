#pragma once

#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <limits>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <gmp.h>
#include <gmpxx.h>
#include <png++/png.hpp>
#include <ncurses.h>

#include "fretriever.h"
#include "heap.h"
#include "factorize.h"
#include "helpers.h"

class Experiment {
protected:
    std::ostream &out; // output stream
public:
    explicit Experiment(std::ostream &outs);

    virtual ~Experiment() = default;

    virtual void printUsage() = 0;

    virtual void perform(std::vector<std::string> args);

    virtual std::string id() = 0;
};

class ExperimentWRet : public Experiment {
protected:
    FRetriever *ret;

    virtual void _perform(std::vector<std::string> args) = 0;

public:
    explicit ExperimentWRet(std::ostream &outs);

    void perform(std::vector<std::string> args) override;
};

template<bool PRINT>
class GenHypo : public ExperimentWRet {
public:

    explicit GenHypo(std::ostream &outs) : ExperimentWRet(outs) {};

    void printUsage() override {
        std::cout << "Test the hypothesis, that f(2^a * 3^b * 5^c) = 2a + 3b + 5c for all a, b, c, with c < 6"
                  << std::endl;
    }

    std::string id() override {
        return "gh";
    }

protected:
    void _perform(const std::vector<std::string> args) override {
        out << "Counterexamples to the hypothesis that f(2^a * 3^b * 5^c) = 2a + 3b + 5c" << std::endl;
        if (PRINT)
            std::cout << "Checking..." << std::endl;
        for (ull i = 0, i2 = 1; i2 < ret->size; i++, i2 *= 2)
            for (ull j = 0, j3 = 1; i2 * j3 < ret->size; j++, j3 *= 3)
                for (ull k = 0, k5 = 1; i2 * j3 * k5 < ret->size && k < 6; k++, k5 *= 5)
                    if (ret->f(i2 * j3 * k5) < 2 * i + 3 * j + 5 * k)
                        out << "f(2^" << i << "*3^" << j << "*5^" << k << ") = " << ret->f(i2 * j3 * k5) << " < "
                            << 2 * i + 3 * j + 5 * k << std::endl;
        if (PRINT)
            std::cout << "Done!" << std::endl;
    }
};

template <bool PRINT>
class SecondSmallest : public ExperimentWRet {
public:

    explicit SecondSmallest(std::ostream &outs) : ExperimentWRet(outs) {};

    void printUsage() override {
        std::cout << "Print the second smallest integer of complexity =n or * otherwise" << std::endl;
    }

    std::string id() override {
        return "ss";
    }

protected:
    void _perform(const std::vector<std::string> args) override {
        out << "Second smallest integer of complexity n" << std::endl;
        const ui MAXCOMPL = 100;
        ull freq[MAXCOMPL];
        ull second[MAXCOMPL];
        ull E[MAXCOMPL];
        E[0] = 0;
        E[1] = 1;
        E[2] = 2;
        E[3] = 3;
        E[4] = 4;
        for (ui i = 5; i < MAXCOMPL; i++)
            E[i] = E[i - 3] * 3;
        for (auto & f : freq)
            f = 0;
        for (ull i = 1; i < ret->size; i++) {
            ui c = ret->bf(i);
            freq[c]++;
            if (freq[c] == 2)
                second[c] = i;
        }
        ui i = 1;
        while (freq[i] > 1 || ret->size > E[i]) {
            if (freq[i] > 1)
                out << i << " " << second[i] << std::endl;
            else
                out << i << " *" << std::endl;
            i++;
        }
    }
};

class ShortestExpr : public Experiment {

    void printShortest(long long n);

protected:
    FRetriever *ret;

public:
    explicit ShortestExpr(std::ostream &outs);

    void printUsage() override;

    std::string id() override;

    void perform(std::vector<std::string> args) override;
};

bool grdbl(const double &lhs, const double &rhs);

bool cmpll(const unsigned long long &lhs, const unsigned long long &rhs);

bool lll(const ui &lhs, const ui &rhs);

template<bool PRINT>
class GoodNums : public Experiment {
    FRetriever *ret = nullptr;

public:
    explicit GoodNums(std::ostream &outs) : Experiment(outs) {};

    void printUsage() override {
        std::cout << "2 parameters: numbers n, m" << std::endl;
        std::cout << "Prints n numbers with the smallest logarithmic complexity" << std::endl;
        std::cout << "All numbers are less than m" << std::endl;
        std::cout << "Numbers are ordered by increasing logarithmic complexity" << std::endl;
    }

    std::string id() override {
        return "good";
    }

    void perform(const std::vector<std::string> args) override {
        const ui CS = 1000000000;
        ret = new FRetriever(1, CS);
        ui n, m;
        const ui MAXFACTORS = 100;
        std::istringstream(args[0]) >> n;
        std::istringstream(args[1]) >> m;
        unsigned long long total = 0, curr = 0;
        ui maxN = m;
        ui maxNsqrt = 10 + static_cast<unsigned long long>(sqrt(static_cast<double>(maxN)));
        MinHeap<double, unsigned long long, grdbl> logcompl(2.0, 5.0, n + 1);
        MinHeap<ui, ui, lll> divs(0xFFFFFFFFU, 0, maxNsqrt);
        ui fs[MAXFACTORS];
        ui divord[MAXFACTORS];
        ui cumdiv[MAXFACTORS];
        ui currord[MAXFACTORS];
        for (ui i = 2; i <= maxN; i++) {
            int icompl = ret->f(i);
            double logc = lc(i, icompl);
            bool hasMult = false;
            if (divs.size() > 0) {
                ui dti = divs.top();
                if (divs.getK(dti) == i) {
                    ui q = i;
                    int factornumber = 0;
                    if (logcompl.size() < n || logc < logcompl.getK(logcompl.top())) {
                        while (divs.getK(dti) == i) {
                            ui factor = divs.getV(dti);
                            if (divs.getK(dti) + factor > factor)
                                divs.increase(dti, divs.getK(dti) + factor);
                            else
                                divs.removeTop();
                            fs[factornumber] = factor;
                            cumdiv[factornumber] = 1;
                            currord[factornumber] = 0;
                            divord[factornumber] = 1;
                            q /= factor;
                            while (q % factor == 0) {
                                divord[factornumber]++;
                                q /= factor;
                            }
                            dti = divs.top();
                            factornumber++;
                        }
                        if (q > 1) {
                            fs[factornumber] = q;
                            divord[factornumber] = 1;
                            cumdiv[factornumber] = 1;
                            currord[factornumber] = 0;
                            factornumber++;
                        }
                        ui divisor = 1;
                        ui o = (divord[0] & 1 ? divord[0] / 2 + 1 : divord[0] / 2);
                        while (currord[0] < o) {
                            divisor *= fs[0];
                            currord[0]++;
                        }
                        cumdiv[0] = divisor;

                        // try all combinations
                        int j = 0;
                        while (true) {
                            j++;
                            while (j < factornumber) {
                                cumdiv[j] = divisor;
                                j++;
                            }
                            // !!! first check if this is the last one and prepare for the next !!!
                            j--;
                            while ((j >= 0 ? currord[j] == divord[j] : 0)) {
                                currord[j] = 0;
                                j--;
                            }
                            if (j < 0)
                                break;
                            else {
                                int fd = divisor < CS ? ret->bf(divisor) : ret->f(divisor);
                                ui divided = i / divisor;
                                int fid = divided < CS ? ret->bf(divided) : ret->f(divided);
                                if (fid + fd <= icompl) {
                                    hasMult = true;
                                    break;
                                }
                                currord[j]++;
                                divisor = cumdiv[j] = cumdiv[j] * fs[j];
                            }
                        }
                    } else
                        while (divs.getK(dti) == i) {
                            ui factor = divs.getV(dti);
                            if (divs.getK(dti) + factor > factor)
                                divs.increase(dti, divs.getK(dti) + factor);
                            else
                                divs.removeTop();
                            dti = divs.top();
                        }
                } else if (i <= maxNsqrt)
                    divs.insert(i + i, i);
            } else
                divs.insert(i + i, i);
            if (!hasMult) {
                if (logcompl.size() == n) {
                    if (logc < logcompl.getK(logcompl.top())) {
                        logcompl.removeTop();
                        logcompl.insert(logc, i);
                    }
                } else
                    logcompl.insert(logc, i);
            }
            if (PRINT) {
                curr++;
                if (curr == 1000000) {
                    total++;
                    std::cout << total << "M done" << std::endl;
                    curr = 0;
                }
            }
        }
        out << "Numbers that cannot be written as products ordered by logarithmic complexity" << std::endl;
        out << logcompl.size() << std::endl;
        while (logcompl.size() > 0) {
            long long ti = logcompl.top();
            out << logcompl.getV(ti) << " " << logcompl.getK(ti) << std::endl;
            logcompl.removeTop();
        }
        delete ret;
    }
};

template<bool PRINT>
class Collapse : public ExperimentWRet {
protected:
    std::vector<std::pair<ui, ui> > good;
    std::vector<ui> coll;
    mpz_t m_fsize{};
    static const ui MAXBITS = 150; // seems enough

    void loadCollapsable(const std::string& fn) {
        std::ifstream in(fn.c_str());
        ui n;
        while (!in.eof()) {
            in >> n;
            coll.push_back(n);
        }
        in.close();
    }

    void loadGoodNumbers(const std::string& fn) {
        std::ifstream in(fn.c_str());
        char header[1000];
        in.getline(header, 1000);
        in.getline(header, 1000);
        ui n;
        in >> n;
        while (!in.eof()) {
            good.emplace_back(n, ret->f(n));
            double d;
            in >> d;
            in >> n;
        }
        in.close();
        for (ui i = 0; i < good.size() / 2; i++) {
            std::pair<ui, ui> temp = good[i];
            good[i] = good[good.size() - i - 1];
            good[good.size() - i - 1] = temp;
        }
    }

    void loadfsize() {
        unsigned long long fsize = ret->size;
        std::ostringstream oss;
        oss << fsize;
        mpz_init_set_str(m_fsize, oss.str().c_str(), 10);
    }

    virtual ui complexity(mpz_t m_n) {
        if (mpz_cmp(m_n, m_fsize) < 0)
            return ret->fbig(m_n);
        ui total = 0;
        for (auto & i : good)
            while (mpz_divisible_ui_p(m_n, i.first)) {
                mpz_divexact_ui(m_n, m_n, i.first);
                total += i.second;
                if (mpz_cmp(m_n, m_fsize) < 0)
                    return total + ret->fbig(m_n);
            }
        mpz_sub_ui(m_n, m_n, 1);
        return total + 1 + complexity(m_n);
    }

public:
    explicit Collapse(std::ostream &outs) : ExperimentWRet(outs) {};

    void printUsage() override {
        std::cout << "2 parameters: file from which to read good numbers and file that contains numbers to collapse"
                  << std::endl;
        std::cout
                << "Prints the power in which a number collapses by repeatedly attempting division with good numbers or subtracting 1"
                << std::endl;
    }

    std::string id() override {
        return "collapse";
    }

    void _perform(const std::vector<std::string> args) override {
        loadGoodNumbers(args[0]);
        if (PRINT)
            std::cout << "Good numbers loaded" << std::endl;
        loadCollapsable(args[1]);
        if (PRINT)
            std::cout << "Collapsable numbers loaded" << std::endl;
        loadfsize();
        for (unsigned int & i : coll) {
            ui icompl = ret->f(i);
            mpz_t m_prime;
            mpz_init_set_d(m_prime, i);
            mpz_t m_pow;
            mpz_init_set_d(m_pow, i);
            mpz_mul(m_pow, m_pow, m_prime);
            ui pow = 2;
            mpz_t m_minus1;
            mpz_init(m_minus1);
            ui collapsesAt = 0;
            ui powcompl;
            while (mpz_sizeinbase(m_pow, 2) < MAXBITS) {
                mpz_sub_ui(m_minus1, m_pow, 1);
                powcompl = 1 + complexity(m_minus1);
                if (PRINT)
                    std::cout << "checked power " << pow << " for " << i << std::endl;
                if (powcompl < pow * icompl) {
                    collapsesAt = pow;
                    break;
                }
                mpz_mul(m_pow, m_pow, m_prime);
                pow++;
            }
            if (collapsesAt == 0)
                out << i << " " << icompl << " ?(>" << pow - 1 << ")" << std::endl;
            else
                out << i << " " << icompl << " " << collapsesAt << " " << powcompl << " < " << pow * icompl
                    << std::endl;
            mpz_clear(m_prime);
            mpz_clear(m_pow);
            mpz_clear(m_minus1);
        }
        mpz_clear(m_fsize);
    }
};

template<bool PRINT>
class Collapse2 : public Collapse<PRINT> {
protected:
    static const ui FCACHESIZE = 1000000000;
    static const ui MAXDEPTH = 100;
    static const ui FACTORS = 4;
    static constexpr double EPS = 1.0e-6;
    mpz_t factors[MAXDEPTH][Collapse<PRINT>::MAXBITS + 1]{};

    virtual ui complexity(mpz_t m_n, int depth) {
        if (mpz_cmp(m_n, this->m_fsize) < 0)
            return this->ret->fbig(m_n);
        ui total = 0;
        for (ui i = 0; i < this->good.size(); i++)
            while (mpz_divisible_ui_p(m_n, this->good[i].first)) {
                mpz_divexact_ui(m_n, m_n, this->good[i].first);
                total += this->good[i].second;
                if (mpz_cmp(m_n, this->m_fsize) < 0)
                    return total + this->ret->fbig(m_n);
            }
        factor(m_n, factors[depth]);
        ui totalFactors = 1;
        while (mpz_cmp_d(factors[depth][totalFactors], 0) > 0)
            totalFactors++;
        totalFactors--;

        // first deal with the case when there are a lot of small factors
        // here the dummy factor is used
        while (totalFactors > FACTORS) {
            ui factori[Collapse<PRINT>::MAXBITS + 1]; // factor indices
            ui totali = 0; // total indices
            ui c = 0;
            double minlogc = 5;
            mpz_t currprod; // current product
            mpz_init_set_d(currprod, 1);
            ui curri = 0; // last factor index
            ui currfactori[Collapse<PRINT>::MAXBITS + 1]; // current factor indices
            currfactori[0] = 0;
            ui j = currfactori[curri] + 1;
            mpz_t temp;
            mpz_init(temp);
            do { // j == next possible factor index
                bool found = false;
                for (; j <= totalFactors; j++) // try to find another factor
                {
                    if (totalFactors - curri >= FACTORS + 1 &&
                        !(mpz_cmp(factors[depth][j - 1], factors[depth][j]) == 0 && currfactori[curri] != j - 1)) {
                        mpz_mul(temp, currprod, factors[depth][j]);
                        if (mpz_cmp(temp, this->m_fsize) < 0) {
                            mpz_set(currprod, temp);
                            curri++;
                            currfactori[curri] = j;
                            j++;
                            char bignumstr[Collapse<PRINT>::MAXBITS];
                            mpz_get_str(bignumstr, 10, currprod);
                            long long int n;
                            sscanf(bignumstr, "%lld", &n);
                            ui ncompl;
                            if (n < FCACHESIZE)
                                ncompl = this->ret->bf(n);
                            else
                                ncompl = this->ret->f(n);
                            double logcompl = ((double) ncompl) * log(3) / log((double) n);
                            if (logcompl <= minlogc + EPS && ncompl > c) {
                                minlogc = logcompl;
                                c = ncompl;
                                for (ui i = 1; i <= curri; i++)
                                    factori[i - 1] = currfactori[i];
                                totali = curri;
                            }
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) { // remove last factor and try to find another with greater index
                    mpz_divexact(currprod, currprod, factors[depth][currfactori[curri]]);
                    j = currfactori[curri] + 1;
                    curri--;
                }
            } while (curri + 1 != 0);
            mpz_clear(currprod);
            mpz_clear(temp);
            // divide by the factors
            total += c;
            for (int i = static_cast<int>(totali) - 1; i >= 0; i--) {
                for (ui k = factori[i]; k < totalFactors; k++)
                    mpz_set(factors[depth][k], factors[depth][k + 1]);
                totalFactors--;
            }
        }

        // now onto the brute force part
        // precalculate the complexity of all products
        ui all_subsets[1 << FACTORS];
        mpz_t product;
        mpz_init(product);
        ui doexceptall = total == 0 ? 1 : 0;
        for (ui i = 1; i < (1 << totalFactors) - doexceptall; i++) {
            mpz_set_d(product, 1);
            for (ui j = 0; j < totalFactors; j++)
                if (i & (1 << j))
                    mpz_mul(product, product, factors[depth][j + 1]);
            if (mpz_cmp(product, this->m_fsize) < 0)
                all_subsets[i] = this->ret->fbig(product);
            else {
                mpz_sub_ui(product, product, 1);
                all_subsets[i] = complexity(product, depth + 1) + 1;
            }
        }
        mpz_clear(product);

        // check all partitions of the product set
        ui set_index[FACTORS];
        ui max_set_index[FACTORS];
        for (ui i = 0; i < totalFactors; i++) {
            set_index[i] = 0;
            max_set_index[i] = 0;
        }
        ui minsofar = 100000;
        if (!doexceptall) {
            // calc total with all factors lumped together
            ui indices[FACTORS];
            for (ui i = 0; i <= max_set_index[totalFactors - 1]; i++)
                indices[i] = 0;
            for (ui i = 0; i < totalFactors; i++)
                indices[set_index[i]] ^= (1 << i);
            ui subtotal = 0;
            for (ui i = 0; i <= max_set_index[totalFactors - 1]; i++)
                subtotal += all_subsets[indices[i]];
            if (subtotal < minsofar)
                minsofar = subtotal;
        }
        if (totalFactors > 1) {
            while (max_set_index[totalFactors - 1] < totalFactors - 1) {
                // next partition
                ui j = totalFactors - 1;
                while (j > 0 ? (set_index[j] == max_set_index[j - 1] + 1) : false)
                    j--;
                set_index[j]++;
                ui msi = max(set_index[j], max_set_index[j - 1]);
                max_set_index[j] = msi;
                for (ui i = j + 1; i < totalFactors; i++) {
                    set_index[i] = 0;
                    max_set_index[i] = msi;
                }
                // calc total
                ui indices[FACTORS];
                for (ui i = 0; i <= max_set_index[totalFactors - 1]; i++)
                    indices[i] = 0;
                for (ui i = 0; i < totalFactors; i++)
                    indices[set_index[i]] |= (1 << i);
                ui subtotal = 0;
                for (ui i = 0; i <= max_set_index[totalFactors - 1]; i++)
                    subtotal += all_subsets[indices[i]];
                if (subtotal < minsofar)
                    minsofar = subtotal;
            }
        }
        total += minsofar;
        return total;
    }

    virtual ui complexity(mpz_t m_n) {
        return std::min(Collapse<PRINT>::complexity(m_n), complexity(m_n, 0));
    }

public:
    explicit Collapse2(std::ostream &outs) : Collapse<PRINT>(outs) {};

    virtual void printUsage() {
        std::cout << "2 parameters: file from which to read good numbers and file that contains numbers to collapse"
                  << std::endl;
        std::cout
                << "Prints the power in which a number collapses by repeatedly attempting division with good numbers and then factorizing and doing a full search by trying each subdivision of factors into factor groups where each group is subtracted 1"
                << std::endl;
    }

    virtual std::string id() {
        return "collapse2";
    }

    virtual void perform(const std::vector<std::string> args) {
        this->ret = new FRetriever(1, FCACHESIZE);
        for (auto & factor : factors)
            for (ui j = 0; j < Collapse<PRINT>::MAXBITS; j++)
                mpz_init(factor[j]);
        this->_perform(args);
        for (auto & factor : factors)
            for (ui j = 0; j < Collapse<PRINT>::MAXBITS; j++)
                mpz_clear(factor[j]);
        delete this->ret;
    }
};

template<bool PRINT>
class Defect : public ExperimentWRet {
public:

    explicit Defect(std::ostream &outs) : ExperimentWRet(outs) {};

    void printUsage() override {
        std::cout << "Print the statistics of defect |A_k(x)| for integer values of k and n values of x" << std::endl;
    }

    std::string id() override {
        return "defect";
    }

protected:
    void _perform(const std::vector<std::string> args) override {
        const ui MAXK = 14;
        const double EPS = 0.001;
        const long long N = ret->size - 1;
        ui xbins;
        std::istringstream(args[0]) >> xbins;
        long long xbinsize = N / xbins;
        out << "The next " << MAXK + 1 << " lines contain " << xbins << " numbers - the size of A_k(x)." << std::endl;
        out << "Each line corresponds to k from 0 to " << MAXK << "." << std::endl;
        out << "Each column corresponds to x from " << xbinsize << " to " << xbinsize * xbins << "." << std::endl;
        std::vector<std::vector<long long> > bins(MAXK + 1, std::vector<long long>(xbins, 0));
        long long nextl = 2;
        long long l = 0;

        long long xbin = 0;
        long long nextx = xbinsize + 1;

        if (PRINT)
            std::cout << "Starting..." << std::endl;
        for (long long i = 1; i <= N; i++) {
            if (i == nextx) {
                xbin++;
                nextx += xbinsize;
                if (PRINT)
                    std::cout << "At i = " << i << std::endl;
            }
            if (i == nextl) {
                l++;
                while (pow(3.0, static_cast<double>(l + 1) / 3.0) - EPS < static_cast<double>(i))
                    l++;
                double d = pow(3.0, static_cast<double>(l + 1) / 3.0);
                if (fabs(d - round(d)) < EPS)
                    nextl = (long long) round(d);
                else
                    nextl = ceil(d);
            }
            bins[ret->bf(i) - l][xbin]++;
        }
        for (long long i = 1; i <= MAXK; i++)
            bins[i][0] += bins[i - 1][0];
        for (long long j = 1; j < xbins; j++)
            bins[0][j] += bins[0][j - 1];
        for (long long i = 1; i <= MAXK; i++)
            for (long long j = 1; j < xbins; j++)
                bins[i][j] += bins[i - 1][j] + bins[i][j - 1] - bins[i - 1][j - 1];
        for (long long i = 0; i <= MAXK; i++) {
            out << bins[i][0];
            for (long long j = 1; j < xbins; j++)
                out << " " << bins[i][j];
            out << std::endl;
        }
        if (PRINT)
            std::cout << "Done!" << std::endl;
    }
};

template <bool PRINT>
class SophiePrimes : public Experiment {
    FRetriever *ret = nullptr;

public:
    explicit SophiePrimes(std::ostream &outs) : Experiment(outs) {};

    void printUsage() override {
        std::cout << "1 parameter: number n" << std::endl;
        std::cout << "Prints primes p <= n, such that (p-1+||p||)/log(p)<(p-2+||2p-1||)/log(2p-1)" << std::endl;
    }

    std::string id() override {
        return "primessop";
    }

    void perform(const std::vector<std::string> args) override {
        const ui CS = 1000000000; // 1e9
        ret = new FRetriever(1, CS);
        unsigned long long n;
        std::istringstream(args[0]) >> n;
        std::vector<unsigned long long> prs;
        std::vector<unsigned char> comps;
        for (unsigned long long i = 2; i <= n; i++)
            if (ret->f(i - 1) + 1 > ret->f(i))
                if (is_prime(i))
                    if (is_prime(i * 2 - 1)) {
                        prs.push_back(i);
                        comps.push_back(ret->f(i));
                    }
        for (unsigned long long i = 0; i < prs.size(); i++)
            if ((double)(prs[i] - 2 + ret->f(2 * prs[i] - 1)) / log((double)(prs[i] * 2 - 1)) > ((double)(prs[i] - 1 + comps[i])) * log((double)prs[i]))
                out << i << std::endl;
        delete ret;
    }
};

template<bool PRINT>
class AverageDigit : public ExperimentWRet {
public:

    explicit AverageDigit(std::ostream &outs) : ExperimentWRet(outs) {};

    void printUsage() override {
        std::cout
                << "Print the histogram of average digit in k-nary and logarithmic complexity up to l in a m x n bw image printed to file x"
                << std::endl;
    }

    std::string id() override {
        return "avehist";
    }

protected:
    void _perform(const std::vector<std::string> args) override {
        ui base;
        std::istringstream(args[0]) >> base;
        long long maxn;
        std::istringstream(args[1]) >> maxn;
        int imgdimx, imgdimy;
        std::istringstream(args[2]) >> imgdimy;
        std::istringstream(args[3]) >> imgdimx;

        std::vector<std::vector<long long> > bins(imgdimy, std::vector<long long>(imgdimx, 0));
        double minavedig = 0;
        double maxavedig = base - 1;
        double minlogc = 3;
        double maxlogc = 26 / log(1439) * log(3);

        double xbinsize = (maxavedig - minavedig) / imgdimx;
        double ybinsize = (maxlogc - minlogc) / imgdimy;

        // printing histogram
        std::vector<ui> krep(50, 0);
        if (PRINT)
            std::cout << "Setup done, starting..." << std::endl;
        ui digsum = 0;
        ui dignum = 0;
        long long maxfr = 0;
        for (long long i = 1; i <= maxn; i++) {
            ui j = 0;
            while (krep[j] + 1 == base) {
                krep[j] = 0;
                digsum -= (base - 1);
                j++;
            }
            krep[j] += 1;
            digsum += 1;
            dignum = max(dignum, j + 1);
            double avedig = (static_cast<double>(digsum)) / dignum;
            double logc = (static_cast<double>(ret->bf(i))) * log(3) / log((double)i);
            ui xbin = max(std::min(static_cast<int>(floor((avedig - minavedig) / xbinsize)), imgdimx - 1), 0);
            ui ybin = max(std::min(static_cast<int>(floor((logc - minlogc) / ybinsize)), imgdimy - 1), 0);
            bins[ybin][xbin] += 1;
            maxfr = max(maxfr, bins[ybin][xbin]);
        }

        if (PRINT)
            std::cout << "Now generating image..." << std::endl;

        png::image<png::gray_pixel> im(imgdimx, imgdimy);
        for (int row = 0; row < imgdimy; row++)
            for (int col = 0; col < imgdimx; col++) {
                im[row][col] = png::gray_pixel(max(min(static_cast<int>(floor(
                        static_cast<double>(bins[imgdimy - 1 - row][col]) / static_cast<double>(maxfr) * 256)), 255), 0));
            }
        im.write(args[4]);

        if (PRINT)
            std::cout << "Done!" << std::endl;
    }
};


template<bool PRINT>
class ComplexityBruteForce : public Experiment {

protected:
    FRetriever *ret = nullptr;
    ui depth{}; // for displaying the recursion tree using ncurses
    // magical constants
    // we wish this code to work for up to 100 bits
    static const ui MAX_ECOMPL = 300;
    static const ull CS = 1024 * 1024 * 1024; // hold the complexity of smallest CS numbers in memory
    static const ui MAX_LB = 10000; //
    const std::string BLANK_LINE = "                                                ";
    mpz_class E[MAX_ECOMPL];
    mpz_class m_fsize;
    ui LBS[MAX_LB]{}; // LBS[i] = minimum complexity over numbers >= i

    ull max_divisors = 0; // curious about the largest number of divisors encountered

    struct hash_mpz {
        static const ui PRIME = 2147480197;

        size_t operator()(const mpz_class &n) const {
            return mpz_fdiv_ui(n.get_mpz_t(), PRIME);
        }
    };

    struct cmp_mpz {
        bool operator()(const mpz_class &a, const mpz_class &b) const {
            return (cmp(a, b) == 0);
        }
    };

    std::unordered_map<mpz_class, std::pair<ui, bool>, hash_mpz, cmp_mpz> fmulcache;
    std::unordered_map<mpz_class, std::pair<ui, bool>, hash_mpz, cmp_mpz> faddcache;


    virtual void init() {
        ret = new FRetriever(1, CS);
        ret->bf(0); // fill cache
        E[0] = 0;
        E[1] = 1;
        E[2] = 2;
        E[3] = 3;
        E[4] = 4;
        for (ui i = 5; i < MAX_ECOMPL; i++) {
            E[i] = E[i - 3];
            E[i] *= 3;
        }

        mpz_class n;
        n = 0;
        for (unsigned int & lb : LBS) {
            lb = trivialLowerBound(n);
            n += 1;
        }
        if (PRINT) {
            depth = 0;
            initscr();
            noecho();
            curs_set(0);
        }
    }

    virtual ui trivialLowerBound(const mpz_class& m_n) {
        ui i = 1;
        while (cmp(E[i], m_n) < 0)
            i++;
        return i;
    }

    virtual void clear() {
        delete ret;
        if (PRINT) {
            endwin();
        }
    }

    void loadfsize() {
        ull fsize = ret->size;
        std::ostringstream oss;
        oss << fsize;
        m_fsize = oss.str();
    }

public:
    explicit ComplexityBruteForce(std::ostream &outs) : Experiment(outs) {};

    void printUsage() override {
        std::cout << "2 parameters: n, k. The experiment tests if the complexity of n does not exceed k" << std::endl;
        std::cout << "Outputs if n has complexity not exceeding k" << std::endl;
    }

    std::string id() override {
        return "complexitybf";
    }

    void perform(const std::vector<std::string> args) override {
        init();
        loadfsize();
        mpz_class number;
        number = args[0];
        ui maxComplexity = std::stoul(args[1], nullptr, 10);
        ui cm = complexityMul(number, maxComplexity);
        ui ca = complexityAdd(number, maxComplexity);
        clear();
        if (cm == 0 && ca == 0) {
            out << "Integer complexity of " << number << " > " << maxComplexity << std::endl;
        } else {
            if (cm > 0 && ca > 0)
                out << "Integer complexity of " << number << " = " << min(cm, ca) << std::endl;
            else
                out << "Integer complexity of " << number << " = " << max(cm, ca) << std::endl;
        }
    }


    // return 0 if complexity > complexityUpperBound, otherwise return complexity
    // NB: this function only works correctly on sufficiently large numbers!
    virtual ui complexityAdd(const mpz_class &number, const ui complexityUpperBound) {
        if (faddcache.count(number) == 1) {
            std::pair<ui, bool> v = faddcache[number];
            if (v.first > complexityUpperBound)
                return 0;
            if (v.second)
                return v.first;
        } else {
            faddcache[number] = std::make_pair(trivialLowerBound(number), false);
        }
        if (PRINT) {
            std::stringstream ss;
            ss << "+" << number << " " << complexityUpperBound;
            move(depth, 0);
            addstr(ss.str().c_str());
            refresh();
            depth++;
        }
        const ui MAXCOMPL = std::numeric_limits<ui>::max()/2;
        ui complexity = MAXCOMPL;
        mpz_class a;
        mpz_class b;
        a = number;
        a--;
        b = 1;
        ui addend = 1;

        while (trivialLowerBound(a) + LBS[addend] < complexityUpperBound + 2) {
            ui bcomplexity = ret->fbig(b);
            ui c = complexityMul(a, min(complexity - 1, complexityUpperBound) - bcomplexity);
            if (c > 0) {
                complexity = c + bcomplexity;
//                if (PRINT)
//                    std::cout << a << "+" << b << std::endl;
            }
            a--;
            b++;
            addend++;
        }
        if (complexity == MAXCOMPL) {
            faddcache[number] = std::make_pair(max(faddcache[number].first, complexityUpperBound + 1), false);
            complexity = 0;
        } else {
            faddcache[number] = std::make_pair(complexity, true);
        }
        if (PRINT) {
            depth--;
            move(depth, 0);
            addstr(BLANK_LINE.c_str());
        }
        return complexity;
    }

    // return 0 if complexity > maxCompl, otherwise return complexity
    virtual ui complexityMul(const mpz_class &number, const ui maxCompl) {
        const ui MAXCOMPL = std::numeric_limits<ui>::max()/2;
        ui complexity = MAXCOMPL;
        if (cmp(number, m_fsize) < 0) {
            complexity = ret->fbig(number);
            if (complexity <= maxCompl)
                return complexity;
            else
                return 0;
        }
        if (fmulcache.count(number) == 1) {
            std::pair<ui, bool> v = fmulcache[number];
            if (v.first > maxCompl)
                return 0;
            if (v.second)
                return v.first;
        } else {
            fmulcache[number] = std::make_pair(trivialLowerBound(number), false);
        }
        if (PRINT) {
            std::stringstream ss;
            ss << "*" << number << " " << maxCompl << "                                      ";
            move(depth, 0);
            addstr(ss.str().c_str());
            refresh();
            depth++;
        }

        std::vector<mpz_class> fs;
        std::vector<mpz_class> pfs;
        std::vector<ui> pfpowers;

        factor(number, fs);

        ui nfs = 1, distpf = 0;
        while (cmp(fs[nfs], 0) > 0) {
            pfs.push_back(fs[nfs]);
            pfpowers.push_back(1);
            while (cmp(fs[nfs], fs[nfs + 1]) == 0) {
                pfpowers[distpf]++;
                nfs++;
            }
            nfs++;
            distpf++;
        }

        // have to be careful about overflow here. For numbers up to 2^100 divisor count of 2^23 seems sufficient
        ull ndivisors = 1;
        for (ui i = 0; i < distpf; i++)
            ndivisors *= pfpowers[i] + 1;
        if (ndivisors > max_divisors) {
            max_divisors = ndivisors;
            if (PRINT)
            {
                move(0, 60);
                addstr("Max divisors: ");
                addstr(std::to_string(max_divisors).c_str());
                move(depth, 0);
                refresh();
            }
        }
        if (ndivisors == 2) // true iff the number is prime
            complexity = MAXCOMPL;
        else {
            std::vector<mpz_class> divisors(ndivisors); // all the divisors
            std::vector<ui> divlb(ndivisors); // divisor lower bound so far (we fill with trivial at the start)
            std::vector<bool> lbopt(ndivisors); // true if lower bound is optimal

            std::vector<ui> divindex(nfs, 0); // holds the index for the divisor in divisors and divlb
            std::vector<std::vector<ui> > pfindices(nfs, std::vector<ui>(distpf,
                                                                         0)); // holds the prime decomposition indices for each divisor (for easier recalculation of divindex)
            std::vector<ui> cumindex(distpf + 1); // holds the cumulative index for ith prime factor
            std::vector<ui> totalpfindex(distpf, 0); // holds the total number of pfs used

            // precalculate cumulative index and divisors and their lower bounds
            ui di = 0;
            std::vector<ui> pfindex(distpf, 0);
            cumindex[0] = 1;
            for (ui i = 0; i < distpf; i++)
                cumindex[i + 1] = cumindex[i] * (pfpowers[i] + 1);

            for (ui i = 0; i < ndivisors - 1; i++) {
                divisors[di] = 1;
                for (ui j = 0; j < distpf; j++)
                    for (ui k = 0; k < pfindex[j]; k++)
                        divisors[di] *= pfs[j];
                if (cmp(divisors[di], m_fsize) < 0) {
                    divlb[di] = ret->fbig(divisors[di]);
                    lbopt[di] = true;
                } else {
                    divlb[di] = trivialLowerBound(divisors[di]);
                    lbopt[di] = false;
                }

                ui l = distpf - 1;
                while (pfindex[l] == pfpowers[l]) {
                    di -= pfindex[l] * cumindex[l];
                    pfindex[l] = 0;
                    l--;
                }
                pfindex[l]++;
                di += cumindex[l];
            }
            divisors[ndivisors - 1] = number;
            divlb[ndivisors - 1] = trivialLowerBound(number);
            lbopt[ndivisors - 1] = false;


            // here the magic happens
            // the initial state is such that every prime divisor is a separate multiplier


            ui fi = 0;
            for (ui i = 0; i < distpf; i++) {
                while (totalpfindex[i] < pfpowers[i]) {
                    pfindices[fi][i]++;
                    divindex[fi] += cumindex[i];
                    fi++;
                    totalpfindex[i]++;
                }
            }

            while (fi > 1) {
                // do stuff for the multiplicative partition, e.g., print for debugging purposes :)
                ui totallb = 0;
                for (ui i = 0; i < fi; i++)
                    totallb += divlb[divindex[i]];
                ui quotient = ndivisors - 1;
                ui cumlb = 0;
                if (totallb <= std::min(complexity - 1, maxCompl)) {
                    for (ui i = 0; i < fi; i++) {
                        if (!lbopt[divindex[i]]) {
                            ui maxSC = divlb[divindex[i]] + std::min(maxCompl, complexity - 1) - totallb;
                            ui c = complexityAdd(divisors[divindex[i]], maxSC);
                            if (c > 0) {
                                ui delta = c - divlb[divindex[i]];
                                divlb[divindex[i]] = c;
                                lbopt[divindex[i]] = true;
                                for (ui j = i; j < fi; j++) {
                                    if (divindex[j] == divindex[i])
                                        totallb += delta;
                                }
                                if (totallb > std::min(maxCompl, complexity - 1))
                                    break;
                            } else {
                                divlb[divindex[i]] = maxSC + 1;
                                totallb = maxCompl + 1; // if c == 0, this partition can be discarded quickly
                                break;
                            }
                        }
                        cumlb += divlb[divindex[i]];
                        quotient -= divindex[i];
                        if (quotient > 0 && cmp(divisors[quotient], m_fsize) <
                                            0) // we're also done if we know the complexity of the remaining stuff
                        {
                            totallb = cumlb + ret->fbig(divisors[quotient]);
                            if (totallb <= maxCompl)
                                complexity = std::min(complexity, totallb);
                            break;
                        }
                    }
                    if (totallb <= maxCompl) {
                        complexity = std::min(complexity, totallb);
//                        if (PRINT) {
//                            for (ui i = 0; i < fi - 1; i++)
//                                std::cout << divisors[divindex[i]] << "*";
//                            std::cout << divisors[divindex[fi - 1]] << std::endl;
//                        }
                    }
                }

                // find next multiplicative partition
                // algorithm: repeatedly scan from right to left to find rightmost largest used pf

                for (ui i = distpf; i > 0; i--) {
                    ui pfi = i - 1;
                    ui j = fi - 1;
                    while (pfindices[j][pfi] == 0)
                        j--;
                    totalpfindex[pfi] -= pfindices[j][pfi];
                    divindex[j] -= pfindices[j][pfi] * cumindex[pfi];
                    pfindices[j][pfi] = 0;
                    if (j > 0) {
                        j--;
                        while (j > 0) {
                            bool isMatching = true;
                            for (ui k = 0; k < pfi; k++)
                                if (pfindices[j - 1][k] != pfindices[j][k]) {
                                    isMatching = false;
                                    break;
                                }
                            if (!isMatching || (pfindices[j - 1][pfi] > pfindices[j][pfi]))
                                break;
                            else {
                                totalpfindex[pfi] -= pfindices[j][pfi];
                                divindex[j] -= pfindices[j][pfi] * cumindex[pfi];
                                pfindices[j][pfi] = 0;
                            }
                            j--;
                        }
                        // dump all the stuff and break from the outer for loop
                        pfindices[j][pfi]++;
                        divindex[j] += cumindex[pfi];
                        totalpfindex[pfi]++;
                        bool isNonempty = true;
                        while (isNonempty) {
                            j++;
                            isNonempty = false;
                            for (ui k = 0; k < distpf; k++)
                                if (pfindices[j][k] > 0) {
                                    isNonempty = true;
                                    break;
                                }
                        }
                        fi = j;
                        if (totalpfindex[pfi] < pfpowers[pfi])
                            j = pfi;
                        else
                            j = pfi + 1;
                        while (j < distpf) {
                            for (ui k = totalpfindex[j]; k < pfpowers[j]; k++) {
                                pfindices[fi][j]++;
                                divindex[fi] += cumindex[j];
                                totalpfindex[j]++;
                                fi++;
                            }
                            j++;
                        }
                        break;
                    }
                }
            }


            // end of magic

        }

        if (complexity == MAXCOMPL) {
            fmulcache[number] = std::make_pair(max(fmulcache[number].first, maxCompl + 1), false);
            complexity = 0;
        } else {
            fmulcache[number] = std::make_pair(complexity, true);
        }
        if (PRINT) {
            depth--;
            move(depth, 0);
            addstr(BLANK_LINE.c_str());
        }
        return complexity;
    }

};

template<bool PRINT>
class FCalculator : public Experiment {
    const ui PROGPER = 1000000; // print a progress message every PROGPER numbers
    const ui MAXDIVNUM = 13; // maximum number of distinct prime factors, 13 is sufficient for up to ~ 7*10^17
    const ui DIVCACHE = 100; // the values f[n/2], f[n/3], f[n/4], ..., f[n/DIVCACHE] will be cached
    const ui CACHESIZE = 1000000; // cache CACHESIZE numbers
    const ull MN = 16534727299ULL; // smallest number for which +8 or greater is necessary (+9 for MN)

    // Globals
    std::vector<std::vector<unsigned char> > cache;
    std::vector<ull> cachestart;
    std::vector<unsigned char> finit;
    ull finitsize{}; // compute complexity for numbers up to finitsize-1 using the sieve method
    std::vector<unsigned char> ftail;
    ull ftailsize{};
    ull ftailoff{}; // ftail[x] == f[x+ftailoff]
    std::fstream ffile;
    ull n{};
    std::vector<ull> addends; // numbers not possible to write as sums
    std::vector<ull> E;
    MinHeap<ull, ull, cmpll> divisors; // we should have enough for pi(sqrt(n))

    unsigned char max_f(const ull& x) {
        if (x > 1)
            return __builtin_popcountll(x) - 1 + 2 * (63 - __builtin_clzll(x));
        else
            return x;
    }

//    ull max_add_check(const ull& x) {
//        unsigned char f = max_f(x);
//        return max_add_check_with_fbound(x, f);
//    }

    ull max_add_check_up_to(const ull& x) {
        unsigned char f = 1;
        ull y;
        for (y = x; y > 2; y /= 2)
            f += 3;
        f += y;
        return max_add_check_with_fbound(x, f);
    }

    ull max_add_check_with_fbound(const ull& x, const unsigned char& f) {
        unsigned char i = f / 2;
        while (E[i] + E[f - i] < x)
            i--;
        return E[i];
    }

    void init_divisors() {
        // sieve primes
        ull largestprime = std::min(isqrt(n), finitsize);
        for (ull i = 2; i <= largestprime; i++) {
            bool hasdiv = false;
            while (divisors.getK(divisors.top()) == i) {
                divisors.increase(divisors.top(), divisors.getK(divisors.top()) + divisors.getV(divisors.top()));
                hasdiv = true;
            }
            if (!hasdiv)
                divisors.insert(2 * i, i);
        }

        // now increase the keys until >= finitsize;
        for (ui i = divisors.size(); i > 0; i--)
            divisors.increase(i, pmultiple(finitsize, divisors.getV(i)) * divisors.getV(i));
    }

    void init_caches() {
        cachestart.resize(DIVCACHE + 1);
        cache.resize(DIVCACHE + 1);
        for (ull i = 2; i <= DIVCACHE; i++) {
            cachestart[i] = 0;
            cache[i].resize(CACHESIZE);
        }
    }

    void refresh_caches(ull start) {
        for (ull i = 2; i <= DIVCACHE; i++) {
            if (i * (cachestart[i] + CACHESIZE - 1) < start + 2 * CACHESIZE) {
                cachestart[i] = start / i;
                if (start % i)
                    cachestart[i]++;
                ffile.seekg(cachestart[i], std::ios_base::beg);
                ffile.read((char *) cache[i].data(), sizeof(unsigned char) * CACHESIZE);
            }
        }
    }


    void init_ftail() {
        ftailsize = max_add_check_up_to(n) + 1;
        ftail.resize(ftailsize);
        ffile.seekg(finitsize - ftailsize, std::ios_base::beg);
        ffile.read((char *) ftail.data(), sizeof(unsigned char) * ftailsize);
        ftailoff = finitsize;
    }

    void calculate_f(std::vector<unsigned char> &f) {
        bool changed = true;
        ull fsize = f.size();
        for (ull i = 0; i < fsize; i++)
            f[i] = max_f(i);
        while (changed) {
            changed = false;
            // multiplication
            for (ull i = 2; i < fsize; i++) {
                for (ull j = 2 * i, k = 2; j < fsize; j += i, k++)
                    if (f[j] > f[k] + f[i]) {
                        changed = true;
                        f[j] = f[k] + f[i];
                    }
            }
            // addition
            ull easy = std::min(fsize, MN);
            for (ull i = 9; i < easy; i++) {
                if (f[i - 1] + 1 < f[i]) {
                    changed = true;
                    f[i] = f[i - 1] + 1;
                }
                if (f[i - 6] + 5 < f[i]) {
                    changed = true;
                    f[i] = f[i - 6] + 5;
                }
            }
            for (ull i = easy; i < fsize; i++) {
                if (f[i - 1] + 1 < f[i]) {
                    changed = true;
                    f[i] = f[i - 1] + 1;
                }
                ull mc = max_add_check_with_fbound(i, f[i]);
                for (ull j = 2; j <= mc; j++)
                    if (f[i] > f[j] + f[i - j]) {
                        changed = true;
                        f[i] = f[j] + f[i - j];
                        mc = max_add_check_with_fbound(i, f[i]);
                    }
            }
        }
    }

    void calc_addends() {
        ull addendsize = max_add_check_up_to(n) + 2;
        if (PRINT && addendsize >= MN)
            std::cout << "Incorrect addends" << std::endl;
        std::vector<unsigned char> f(addendsize, 0);
        calculate_f(f);
        addends.push_back(1);
        addends.push_back(6);

        for (ull i = 7; i <= addendsize; i++) {
            if (f[i] < std::min(f[i - 1] + 1, f[i - 6] + 5))
                addends.push_back(i);
        }
        addends.push_back(std::numeric_limits<ull>::max() / 4);
    }

    void calc_E() {
        for (ui i = 0; i <= 4; i++)
            E.push_back(i);
        for (ui i = 5; i < 121; i++)
            E.push_back(3 * E[i - 3]);
    }


    void finit_thread(const ull threadnum, bool &changed, const ull start, const ull end) {
        const ull finitsqrt = isqrt(finitsize);
        // multiplication
        for (ull i = 2; i < finitsqrt; i++) {
            for (ull k = pmultiple(start, i), j = k * i; j < end; j += i, k++)
                if (finit[j] > finit[k] + finit[i]) {
                    changed = true;
                    finit[j] = finit[k] + finit[i];
                }
        }
        if (PRINT)
            std::cout << "Multiplication done [" << threadnum << "]" << std::endl;
        // addition
        for (ull i = std::max(9ULL, start); i < std::min(MN, end); i++) {
            if (finit[i - 1] + 1 < finit[i]) {
                changed = true;
                finit[i] = finit[i - 1] + 1;
            }
            if (finit[i - 6] + 5 < finit[i]) {
                changed = true;
                finit[i] = finit[i - 6] + 5;
            }
        }
        for (ull i = MN; i < end; i++) {
            if (finit[i - 1] + 1 < finit[i]) {
                changed = true;
                finit[i] = finit[i - 1] + 1;
            }
            ull mc = max_add_check_with_fbound(i, finit[i]);
            for (ull j = 1; addends[j] <= mc; j++)
                if (finit[i] > finit[addends[j]] + finit[i - addends[j]]) {
                    changed = true;
                    finit[i] = finit[addends[j]] + finit[i - addends[j]];
                    mc = max_add_check_with_fbound(i, finit[i]);
                }
        }
        if (PRINT)
            std::cout << "Addition done [" << threadnum << "]" << std::endl;
    }

    void calculate_finit_threaded() {
        const ui THREADS = 8;
        ull start = 9;
        ull bucketsize = (finitsize - start) / THREADS;
        bool changed = true;
        for (ull i = 0; i < finitsize; i++)
            finit[i] = max_f(i);
        if (PRINT)
            std::cout << "Setting base 2 upper bounds done." << std::endl;
        ull passes = 0;

        while (changed) {
            changed = false;
            std::vector<std::thread> threads;
            for (ui i = 0; i < THREADS - 1; i++)
                threads.push_back(std::thread(&FCalculator<PRINT>::finit_thread, this, i, std::ref(changed),
                                              start + bucketsize * i, start + bucketsize * (i + 1)));
            threads.push_back(std::thread(&FCalculator<PRINT>::finit_thread, this, THREADS - 1, std::ref(changed),
                                          start + bucketsize * (THREADS - 1), finitsize));
            for (ui i = 0; i < THREADS; i++)
                threads[i].join();
            passes++;
            if (PRINT)
                std::cout << "passes = " << passes << std::endl;
        }
    }

    void save_init() {
        ffile.seekp(0, std::ios_base::beg);
        ffile.write((char *) finit.data(), finitsize);
    }


    void init() {
        finitsize = std::min(get_total_system_memory() / 2, n + 1);
        if (PRINT) {
            std::cout << "Init size: " << finitsize << std::endl;
        }
        finit.resize(finitsize);

        calc_E();
        calc_addends();
        //calculate_finit();
        calculate_finit_threaded();
        save_init();
    }

    unsigned char& ftailf(const ull& i) {
        if (i < ftailoff)
            return (ftail[i + ftailsize - ftailoff]);
        else
            return (ftail[i - ftailoff]);
    }

    void calculate() {
        std::vector<ull> div(MAXDIVNUM, 0); // divisors
        std::vector<ull> divord(MAXDIVNUM, 0); // divisor ord
        std::vector<ull> cumdiv(MAXDIVNUM, 0);
        std::vector<ull> currdivord(MAXDIVNUM, 0);
        long divi; // used div
        ull cacheperiod = 0;
        ull prev = 0;
        ull sqrtn = isqrt(n);
        for (ull i = finitsize; i <= n; i++) {
            // flush/refresh caches
            if (cacheperiod == 2 * CACHESIZE) {
                cacheperiod = 0;
                refresh_caches(i);
                if (PRINT) {
                    if (i - prev > PROGPER) {
                        std::cout << "at " << i << std::endl;;
                        prev = i;
                    }
                }
            }
            if (i >= ftailoff + ftailsize) {
                ffile.seekp(ftailoff, std::ios_base::beg);
                ffile.write((char *) ftail.data(), ftailsize);
                ftailoff = i;
            }

            // set a reasonable bound
            ftailf(i) = ftailf(i - 1) + 1;
            // multiplications
            // first retrieve all prime divisors
            ull j = i;
            divi = 0;
            while (divisors.getK(divisors.top()) == i) {
                ull di = divisors.top();
                div[divi] = divisors.getV(di);
                cumdiv[divi] = 1;
                currdivord[divi] = 0;
                divisors.increase(di, divisors.getK(di) + divisors.getV(di));
                divord[divi] = 1;
                j /= div[divi];
                while (j % div[divi] == 0) {
                    j /= div[divi];
                    divord[divi]++;
                }
                divi++;
            }
            if (j != 1) {
                div[divi] = j;
                cumdiv[divi] = 1;
                currdivord[divi] = 0;
                divord[divi] = 1;
                divi++;
                if (j <= sqrtn)
                    divisors.insert(2 * j, j);
            }
            ull divisor = 1;
            ull o = (divord[0] & 1 ? divord[0] / 2 + 1 : divord[0] / 2);
            while (currdivord[0] < o) {
                divisor *= div[0];
                currdivord[0]++;
            }
            cumdiv[0] = divisor;

            // try all combinations
            long k = 0;
            while (true) {
                k++;
                while (k < divi) {
                    cumdiv[k] = divisor;
                    k++;
                }
                // !!! first check if this is the last one and prepare for the next !!!
                k--;
                while ((k >= 0 ? currdivord[k] == divord[k] : 0)) {
                    currdivord[k] = 0;
                    k--;
                }
                if (k < 0) {
                    break;
                } else {
                    ull q = i / divisor;
                    unsigned char fdivisor = (q <= DIVCACHE ? cache[q][divisor - cachestart[q]] : finit[divisor]);
                    unsigned char fq = (divisor <= DIVCACHE ? cache[divisor][q - cachestart[divisor]] : finit[q]);
                    unsigned char alt = fq + fdivisor;
                    if (alt < ftailf(i))
                        ftailf(i) = alt;
                    currdivord[k]++;
                    divisor = cumdiv[k] = cumdiv[k] * div[k];
                }
            }

            // sums
            ull mc = max_add_check_with_fbound(i, ftailf(i));
            for (ull l = 0; addends[l] <= mc; l++) {
                unsigned char alt = finit[addends[l]] + ftailf(i - addends[l]);
                if (alt < ftailf(i)) {
                    ftailf(i) = alt;
                    mc = max_add_check_with_fbound(i, alt);
                }
            }
            cacheperiod++;
        }
        ffile.seekp(ftailoff, std::ios_base::beg);
        ffile.write((char *) ftail.data(), sizeof(unsigned char) * (n + 1 - ftailoff));
        ftailoff = n + 1;
    }

public:
    explicit FCalculator(std::ostream &outs) : Experiment(outs), divisors(std::numeric_limits<ull>::max(), 0, 6000000) {};

    void printUsage() override {
        std::cout << "2 parameters: number n, and filename f" << std::endl;
        std::cout << "Computes integer complexity for all numbers 0 to n" << std::endl;
        std::cout << "Stores in i-th byte of file f the complexity of number i" << std::endl;
    }

    std::string id() override {
        return "complexity";
    }

    void perform(const std::vector<std::string> args) override {
        std::istringstream(args[0]) >> n;
        ffile.open(args[1], std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::trunc);
        if (PRINT && !ffile.is_open()) {
            std::cout << "Failed to open file: " << strerror(errno) << std::endl;
        }
        init();
        if (finitsize - 1 < n) {
            init_divisors();
            if (PRINT)
                std::cout << "Primes done." << std::endl;
            init_ftail();
            if (PRINT)
                std::cout << "Tail done." << std::endl;
            init_caches();
            refresh_caches(finitsize);
            if (PRINT)
                std::cout << "Caches done." << std::endl;
            calculate();
        }
        ffile.close();
    }
};
