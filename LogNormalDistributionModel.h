#ifndef LAB3_LOGNORMALDISTRIBUTIONMODEL_H
#define LAB3_LOGNORMALDISTRIBUTIONMODEL_H

#include "PRNG.h"

class LogNormalDistributionModel : public PRNG {
private:
    PRNG *standartNormalDistributionModel;
    PRNG *prng;
    double mu;
    double quadSigma;
public:
    LogNormalDistributionModel(double mu, double quadSigma, const PRNG *prng);
    explicit LogNormalDistributionModel(const LogNormalDistributionModel *logNormalDistributionModel);
    ~LogNormalDistributionModel();

    double next() override;
    void reset() override;

    LogNormalDistributionModel* clone() const override;
};

#endif //LAB3_LOGNORMALDISTRIBUTIONMODEL_H
