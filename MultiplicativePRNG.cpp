#include "MultiplicativePRNG.h"

MultiplicativePRNG::MultiplicativePRNG(long long module, long long seed, int multiplier)
{
    this->module = module;
    this->multiplier = multiplier;
    this->seed = seed;
    this->last = seed;
}

MultiplicativePRNG::MultiplicativePRNG(const MultiplicativePRNG* prng)
{
    this->module = prng->module;
    this->multiplier = prng->multiplier;
    this->seed = prng->seed;
    this->last = prng->last;
}

MultiplicativePRNG::~MultiplicativePRNG()
= default;

double MultiplicativePRNG::next()
{
    last = (last * multiplier) % module;
    return (double)last / module;
}

void MultiplicativePRNG::reset()
{
    last = seed;
}

MultiplicativePRNG* MultiplicativePRNG::clone() const
{
    return new MultiplicativePRNG(this);
}
