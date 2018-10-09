#include "CauchyDistributionModel.h"
#include <cmath>
#include "NormalDistributionModel.h"

CauchyDistributionModel::CauchyDistributionModel(double location, double scale, const PRNG* prng):
        location(location),
        scale(scale),
        prng((PRNG*)(prng)->clone())
{

}

CauchyDistributionModel::CauchyDistributionModel(const CauchyDistributionModel* cauchyDistributionModel):
        prng((PRNG*)(cauchyDistributionModel->prng->clone())),
        location(cauchyDistributionModel->location),
        scale(cauchyDistributionModel->scale)
{

}

CauchyDistributionModel::~CauchyDistributionModel()
{
    delete prng;
}

double CauchyDistributionModel::next()
{
    return location + scale * tan(2 * M_PI * prng->next());
}

void CauchyDistributionModel::reset()
{
    prng->reset();
}

CauchyDistributionModel* CauchyDistributionModel::clone() const
{
    return new CauchyDistributionModel(this);
}
