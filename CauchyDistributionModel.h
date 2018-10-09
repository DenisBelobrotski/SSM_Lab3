#ifndef LAB3_CAUCHYDISTRIBUTIONMODEL_H
#define LAB3_CAUCHYDISTRIBUTIONMODEL_H

#include "PRNG.h"

class CauchyDistributionModel : public PRNG {
private:
    PRNG *prng;
    double location;
    double scale;
public:
    CauchyDistributionModel(double location, double scale, const PRNG *prng);
    explicit CauchyDistributionModel(const CauchyDistributionModel *cauchyDistributionModel);
    ~CauchyDistributionModel();

    double next() override;
    void reset() override;

    CauchyDistributionModel* clone() const override;
};


#endif //LAB3_CAUCHYDISTRIBUTIONMODEL_H
