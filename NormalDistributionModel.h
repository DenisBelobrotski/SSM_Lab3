#ifndef LAB3_NORMALDISTRIBUTIONMODEL_H
#define LAB3_NORMALDISTRIBUTIONMODEL_H

#include "PRNG.h"


class NormalDistributionModel : public PRNG {
private:
    PRNG *prng;
    double mean;
    double variance;
    double cachedValue;
    bool hasCachedValue;
public:
    NormalDistributionModel(double mean, double variance, const PRNG *prng);
    explicit NormalDistributionModel(const NormalDistributionModel *normalDistributionModel);
    ~NormalDistributionModel();

    double next() override;
    void reset() override;

    NormalDistributionModel* clone() const override;
};


#endif //LAB3_NORMALDISTRIBUTIONMODEL_H
