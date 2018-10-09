#ifndef LAB2_MULTIPLICATIVEPRNG_H
#define LAB2_MULTIPLICATIVEPRNG_H

#include "PRNG.h"


class MultiplicativePRNG : public PRNG
{
private:
    long long module;
    long long seed;
    long long last;
    int multiplier;
public:
    MultiplicativePRNG(long long module, long long seed, int multiplier);
    explicit MultiplicativePRNG(const MultiplicativePRNG *prng);
    ~MultiplicativePRNG();

    double next() override;
    void reset() override;

    MultiplicativePRNG* clone() const override;
};


#endif //LAB2_MULTIPLICATIVEPRNG_H
