#include "NormalDistributionModel.h"
#include <cmath>

NormalDistributionModel::NormalDistributionModel(double mean, double variance, const PRNG* prng):
mean(mean),
variance(variance),
prng((PRNG*)(prng)->clone())
{
    cachedValue = 0.0;
    hasCachedValue = false;
}

NormalDistributionModel::NormalDistributionModel(const NormalDistributionModel* normalDistributionModel):
prng((PRNG*)(normalDistributionModel->prng->clone())),
mean(normalDistributionModel->mean),
variance(normalDistributionModel->variance)
{
    this->cachedValue = normalDistributionModel->cachedValue;
    this->hasCachedValue = normalDistributionModel->hasCachedValue;
}

NormalDistributionModel::~NormalDistributionModel()
{
    delete prng;
}

double NormalDistributionModel::next()
{
    double result;

    if(!hasCachedValue)
    {
        double squareRoot = sqrt(-2.0 * log(prng->next()));
        double angle = 2.0 * M_PI * prng->next();
        double deviation = sqrt(variance);

        result = mean + deviation * (squareRoot * cos(angle));
        cachedValue = mean + deviation * (squareRoot * sin(angle));
        hasCachedValue = true;
    }
    else
    {
        result = cachedValue;
        hasCachedValue = false;
    }

    return result;
}

void NormalDistributionModel::reset()
{
    prng->reset();
    hasCachedValue = false;
}

NormalDistributionModel* NormalDistributionModel::clone() const
{
    return new NormalDistributionModel(this);
}
