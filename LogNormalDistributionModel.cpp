#include "LogNormalDistributionModel.h"
#include <cmath>
#include "NormalDistributionModel.h"

LogNormalDistributionModel::LogNormalDistributionModel(double mu, double quadSigma, const PRNG* prng):
        mu(mu),
        quadSigma(quadSigma),
        prng((PRNG*)(prng)->clone())
{
    standartNormalDistributionModel = new NormalDistributionModel(0, 1, prng);
}

LogNormalDistributionModel::LogNormalDistributionModel(const LogNormalDistributionModel* logNormalDistributionModel):
        prng((PRNG*)(logNormalDistributionModel->prng->clone())),
        mu(logNormalDistributionModel->mu),
        quadSigma(logNormalDistributionModel->quadSigma)
{

}

LogNormalDistributionModel::~LogNormalDistributionModel()
{
    delete standartNormalDistributionModel;
    delete prng;
}

double LogNormalDistributionModel::next()
{
    return exp(mu + sqrt(quadSigma) * standartNormalDistributionModel->next());
}

void LogNormalDistributionModel::reset()
{
    prng->reset();
}

LogNormalDistributionModel* LogNormalDistributionModel::clone() const
{
    return new LogNormalDistributionModel(this);
}
