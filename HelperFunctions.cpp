#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "MultiplicativePRNG.h"
#include "NormalDistributionModel.h"
#include "LogNormalDistributionModel.h"
#include "CauchyDistributionModel.h"

using namespace std;

FILE* outputFile;

bool checkPearsonTestNormal(const double* quantiles, double* sequence, int num, int cellNum, int erfStepNum,
        double mean, double variance);
bool checkPearsonTestLogNormal(const double* quantiles, double* sequence, int num, int cellNum, int erfStepNum,
                               double mu, double quadSigma);
bool checkPearsonTestCauchy(const double* quantiles, double* sequence, int num, int cellNum, double location,
        double scale);
int* calcFrequenciesEmperic(double* sequence, int num, int cellNum, double leftBorder, double rightBorder);

bool checkKolmogorovTestNormal(double quantile, double* sequence, int num, int erfStepNum, double mean, double variance);
double calcKolmogorovDistanceNormal(double* sequence, int num, int erfStepNum, double mean, double variance);
bool checkKolmogorovTestLogNormal(double quantile, double *sequence, int num, int erfStepNum, double mu,
                                  double quadSigma);
double calcKolmogorovDistanceLogNormal(double* sequence, int num, int erfStepNum, double mu, double quadSigma);
bool checkKolmogorovTestCauchy(double quantile, double* sequence, int num, double location, double scale);
double calcKolmogorovDistanceCauchy(double* sequence, int num, double location, double scale);

double calcNormalCDF(double x, double mean, double variance, int erfStepNum);
double calcErf(double max, int stepNum);
double calcErf(double from, double to, int stepNum);

double calcLogNormalMean(double mu, double quadSigma);
double calcLogNormalVariance(double mu, double quadSigma);
double calcLogNormalCDF(double x, double mu, double quadSigma, int erfStepNum);

double calcCauchyMean(double location, double scale);
double calcCauchyVariance(double location, double scale);
double calcCauchyCDF(double x, double location, double scale);

double calcEstimationMean(const double* source, int num);
double calcEstimationVariance(const double* source, int num, double mean);

const char* printBool(bool value);

void calcSolutions()
{
    const int num = 1000;

    const long long module = 2LL << 30;
    const long long seed = 79507;
    const int multiplier = 79507;

    const double normalDistributionMean = -3;
    const double normalDistributionVariance = 16;

    const double logNormalDistributionMu = 0;
    const double logNormalDistributionQuadSigma = 4;

    const double cauchyDistributionLocation = 1;
    const double cauchyDistributionScale = 2;

    const double quantilesPearson[20] = {3.8415, 5.9915, 7.8147, 9.4877, 11.07, 12.592, 14.067, 15.507, 16.919, 18.307,
                                  19.675, 21.026, 22.362, 23.685, 24.996, 26.296, 27.587, 28.869, 30.144, 31.41};
    const double quantileKolmogorov = 1.36;
    const int cellNum = 20;
    const int erfStepNum = 200;

    PRNG* multiplicativePRNG = new MultiplicativePRNG(module, seed, multiplier);
    PRNG* normalDistributionModel =
            new NormalDistributionModel(normalDistributionMean, normalDistributionVariance, multiplicativePRNG);
    PRNG* logNormalDistributionModel =
            new LogNormalDistributionModel(logNormalDistributionMu, logNormalDistributionQuadSigma, multiplicativePRNG);
    PRNG* cauchyDistributionModel =
            new CauchyDistributionModel(cauchyDistributionLocation, cauchyDistributionScale, multiplicativePRNG);

    delete multiplicativePRNG;

    fopen_s(&outputFile, "output.txt", "w");

    auto* normalDistributionResult = new double[num];
    auto* logNormalDistributionResult = new double[num];
    auto* cauchyDistributionResult = new double[num];

    double mean;
    double variance;
    bool testResult;

    for (int i = 0; i < num; i++)
    {
        normalDistributionResult[i] = normalDistributionModel->next();
        logNormalDistributionResult[i] = logNormalDistributionModel->next();
        cauchyDistributionResult[i] = cauchyDistributionModel->next();
    }

    printf("Normal model (%lf, %lf):\n", normalDistributionMean, normalDistributionVariance);
    fprintf(outputFile, "Normal model (%lf, %lf):\n", normalDistributionMean, normalDistributionVariance);
    for (int i = 0; i < num; i++)
    {
        fprintf(outputFile, "%lf\n", normalDistributionResult[i]);
    }
    mean = calcEstimationMean(normalDistributionResult, num);
    variance = calcEstimationVariance(normalDistributionResult, num, mean);
    printf("Real mean: %f\n", mean);
    fprintf(outputFile, "Real mean: %f\n", mean);
    printf("Real variance: %f\n", variance);
    fprintf(outputFile, "Real variance: %f\n", variance);
    mean = normalDistributionMean;
    variance = normalDistributionVariance;
    printf("Theoretical mean: %f\n", mean);
    fprintf(outputFile, "Theoretical mean: %f\n", mean);
    printf("Theoretical variance: %f\n", variance);
    fprintf(outputFile, "Theoretical variance: %f\n", variance);
    testResult = checkPearsonTestNormal(quantilesPearson, normalDistributionResult, num, cellNum, erfStepNum,
            normalDistributionMean, normalDistributionVariance);
    printf("Pearson test: %s\n", printBool(testResult));
    fprintf(outputFile, "Pearson test: %s\n", printBool(testResult));
    testResult = checkKolmogorovTestNormal(quantileKolmogorov, normalDistributionResult, num, erfStepNum,
            normalDistributionMean, normalDistributionVariance);
    printf("Kolmogorov test: %s\n", printBool(testResult));
    fprintf(outputFile, "Kolmogorov test: %s\n", printBool(testResult));

    printf("\nLog-normal model (%lf, %lf):\n", logNormalDistributionMu, logNormalDistributionQuadSigma);
    fprintf(outputFile, "\nLog-normal model (%lf, %lf):\n", logNormalDistributionMu, logNormalDistributionQuadSigma);
    for (int i = 0; i < num; i++)
    {
        fprintf(outputFile, "%lf\n", logNormalDistributionResult[i]);
    }
    mean = calcEstimationMean(logNormalDistributionResult, num);
    variance = calcEstimationVariance(logNormalDistributionResult, num, mean);
    printf("Real mean: %f\n", mean);
    fprintf(outputFile, "Real mean: %f\n", mean);
    printf("Real variance: %f\n", variance);
    fprintf(outputFile, "Real variance: %f\n", variance);
    mean = calcLogNormalMean(logNormalDistributionMu, logNormalDistributionQuadSigma);
    variance = calcLogNormalVariance(logNormalDistributionMu, logNormalDistributionQuadSigma);
    printf("Theoretical mean: %f\n", mean);
    fprintf(outputFile, "Theoretical mean: %f\n", mean);
    printf("Theoretical variance: %f\n", variance);
    fprintf(outputFile, "Theoretical variance: %f\n", variance);
    testResult = checkPearsonTestLogNormal(quantilesPearson, logNormalDistributionResult, num, cellNum, erfStepNum,
            logNormalDistributionMu, logNormalDistributionQuadSigma);
    printf("Pearson test: %s\n", printBool(testResult));
    fprintf(outputFile, "Pearson test: %s\n", printBool(testResult));
    testResult = checkKolmogorovTestLogNormal(quantileKolmogorov, logNormalDistributionResult, num,  erfStepNum,
                                              logNormalDistributionMu, logNormalDistributionQuadSigma);
    printf("Kolmogorov test: %s\n", printBool(testResult));
    fprintf(outputFile, "Kolmogorov test: %s\n", printBool(testResult));

    printf("\nCauchy model (%lf, %lf):\n", cauchyDistributionLocation, cauchyDistributionScale);
    fprintf(outputFile, "\nExponential model (%lf, %lf):\n", cauchyDistributionLocation, cauchyDistributionScale);
    for (int i = 0; i < num; i++)
    {
        fprintf(outputFile, "%lf\n", cauchyDistributionResult[i]);
    }
    mean = calcEstimationMean(cauchyDistributionResult, num);
    variance = calcEstimationVariance(cauchyDistributionResult, num, mean);
    printf("Real mean: %f\n", mean);
    fprintf(outputFile, "Real mean: %f\n", mean);
    printf("Real variance: %f\n", variance);
    fprintf(outputFile, "Real variance: %f\n", variance);
    mean = calcCauchyMean(cauchyDistributionLocation, cauchyDistributionScale);
    variance = calcCauchyVariance(cauchyDistributionLocation, cauchyDistributionScale);
    printf("Theoretical mean: %f\n", mean);
    fprintf(outputFile, "Theoretical mean: %f\n", mean);
    printf("Theoretical variance: %f\n", variance);
    fprintf(outputFile, "Theoretical variance: %f\n", variance);
    testResult = checkPearsonTestCauchy(quantilesPearson, cauchyDistributionResult, num, cellNum,
            cauchyDistributionLocation, cauchyDistributionScale);
    printf("Pearson test: %s\n", printBool(testResult));
    fprintf(outputFile, "Pearson test: %s\n", printBool(testResult));
    testResult = checkKolmogorovTestCauchy(quantileKolmogorov, cauchyDistributionResult, num,
            cauchyDistributionLocation, cauchyDistributionScale);
    printf("Kolmogorov test: %s\n", printBool(testResult));
    fprintf(outputFile, "Kolmogorov test: %s\n", printBool(testResult));

    fclose(outputFile);

    delete[] normalDistributionResult;
    delete[] logNormalDistributionResult;
    delete[] cauchyDistributionResult;

    delete normalDistributionModel;
    delete logNormalDistributionModel;
    delete cauchyDistributionModel;
}

bool checkPearsonTestNormal(const double* quantiles, double* sequence, int num, int cellNum, int erfStepNum,
                            double mean, double variance)
{
    int* empericFreq;
    double chi = 0.0;
    double expectedCount;
    double step;
    double curBorder;
    double prevBorder;

    sort(sequence, &sequence[num]);
    empericFreq = calcFrequenciesEmperic(sequence, num, cellNum, sequence[0], sequence[num - 1]);

    step = (sequence[num - 1] - sequence[0]) / cellNum;
    curBorder = sequence[0];

    for(int i = 0; i < cellNum; i++)
    {
        prevBorder = curBorder;
        curBorder = (i + 1) * step;
        expectedCount = num * (calcNormalCDF(curBorder, mean, variance, erfStepNum) -
                calcNormalCDF(prevBorder, mean, variance, erfStepNum));
        chi += pow(empericFreq[i] - expectedCount, 2.0) / expectedCount;
    }

    delete[] empericFreq;

    printf("Chi: %f\n", chi);
    fprintf(outputFile, "Chi: %f\n", chi);
    printf("Quantile: %f\n", quantiles[cellNum - 2]);
    fprintf(outputFile, "Quantile: %f\n", quantiles[cellNum - 2]);

    return (chi < quantiles[cellNum - 2]);
}

bool checkPearsonTestLogNormal(const double* quantiles, double* sequence, int num, int cellNum, int erfStepNum,
                            double mu, double quadSigma)
{
    int* empericFreq;
    double chi = 0.0;
    double expectedCount;
    double step;
    double curBorder;
    double prevBorder;

    sort(sequence, &sequence[num]);
    empericFreq = calcFrequenciesEmperic(sequence, num, cellNum, sequence[0], sequence[num - 1]);

    step = (sequence[num - 1] - sequence[0]) / cellNum;
    curBorder = sequence[0];

    for(int i = 0; i < cellNum; i++)
    {
        prevBorder = curBorder;
        curBorder = (i + 1) * step;
        expectedCount = num * (calcLogNormalCDF(curBorder, mu, quadSigma, erfStepNum) -
                calcLogNormalCDF(prevBorder, mu, quadSigma, erfStepNum));
        chi += pow(empericFreq[i] - expectedCount, 2.0) / expectedCount;
    }

    delete[] empericFreq;

    printf("Chi: %f\n", chi);
    fprintf(outputFile, "Chi: %f\n", chi);
    printf("Quantile: %f\n", quantiles[cellNum - 2]);
    fprintf(outputFile, "Quantile: %f\n", quantiles[cellNum - 2]);

    return (chi < quantiles[cellNum - 2]);
}

bool checkPearsonTestCauchy(const double* quantiles, double* sequence, int num, int cellNum, double location,
        double scale)
{
    int* empericFreq;
    double chi = 0.0;
    double expectedCount;
    double step;
    double curBorder;
    double prevBorder;

    sort(sequence, &sequence[num]);
    empericFreq = calcFrequenciesEmperic(sequence, num, cellNum, sequence[0], sequence[num - 1]);

    step = (sequence[num - 1] - sequence[0]) / cellNum;
    curBorder = sequence[0];

    for(int i = 0; i < cellNum; i++)
    {
        prevBorder = curBorder;
        curBorder = (i + 1) * step;
        expectedCount = num * (calcCauchyCDF(curBorder, location, scale) - calcCauchyCDF(prevBorder, location, scale));
        chi += pow(empericFreq[i] - expectedCount, 2.0) / expectedCount;
    }

    delete[] empericFreq;

    printf("Chi: %f\n", chi);
    fprintf(outputFile, "Chi: %f\n", chi);
    printf("Quantile: %f\n", quantiles[cellNum - 2]);
    fprintf(outputFile, "Quantile: %f\n", quantiles[cellNum - 2]);

    return (chi < quantiles[cellNum - 2]);
}

int* calcFrequenciesEmperic(double* sequence, int num, int cellNum, double leftBorder, double rightBorder)
{
    int* result = new int[cellNum];
    double step = (rightBorder - leftBorder) / cellNum;
    double curBorder = step;
    int resultIndex = 0;

    for (int i = 0; i < cellNum; i++)
    {
        result[i] = 0;
    }

    for (int i = 0; i < num; i++)
    {
        if (sequence[i] < curBorder)
        {
            result[resultIndex]++;
        }
        else if (resultIndex < cellNum - 1)
        {
            resultIndex++;
            curBorder += step;
        }
    }
    printf("\nGistogramma:\n");
    for (int i = 0; i < cellNum; i++)
    {
        printf("%d\n", result[i]);
    }
    printf("\n");
    return result;
}

bool checkKolmogorovTestNormal(double quantile, double* sequence, int num, int erfStepNum, double mean, double variance)
{
    double distance;
    sort(sequence, &sequence[num]);
    distance = sqrt(num) * calcKolmogorovDistanceNormal(sequence, num, erfStepNum, mean, variance);
    printf("Kolmogorov distance: %lf\n", distance);
    fprintf(outputFile, "Kolmogorov distance: %lf\n", distance);
    printf("Quantile: %f\n", quantile);
    fprintf(outputFile, "Quantile: %f\n", quantile);
    return (distance < quantile);
}

double calcKolmogorovDistanceNormal(double* sequence, int num, int erfStepNum, double mean, double variance)
{
    double result = 0.0;
    for (int i = 0; i < num; i++)
    {
        result = max(result, abs(calcNormalCDF(sequence[i], mean, variance, erfStepNum) - (double) (i + 1) / num));
    }
    return result;
}

bool checkKolmogorovTestLogNormal(double quantile, double *sequence, int num, int erfStepNum, double mu, double quadSigma)
{
    double distance;
    sort(sequence, &sequence[num]);
    distance = sqrt(num) * calcKolmogorovDistanceLogNormal(sequence, num, erfStepNum, mu, quadSigma);
    printf("Kolmogorov distance: %lf\n", distance);
    fprintf(outputFile, "Kolmogorov distance: %lf\n", distance);
    printf("Quantile: %f\n", quantile);
    fprintf(outputFile, "Quantile: %f\n", quantile);
    return (distance < quantile);
}

double calcKolmogorovDistanceLogNormal(double* sequence, int num, int erfStepNum, double mu, double quadSigma)
{
    double result = 0.0;
    for (int i = 0; i < num; i++)
    {
        result = max(result, abs(calcLogNormalCDF(sequence[i], mu, quadSigma, erfStepNum) - (double) (i + 1) / num));
    }
    return result;
}

bool checkKolmogorovTestCauchy(double quantile, double* sequence, int num, double location, double scale)
{
    double distance;
    sort(sequence, &sequence[num]);
    distance = sqrt(num) * calcKolmogorovDistanceCauchy(sequence, num, location, scale);
    printf("Kolmogorov distance: %lf\n", distance);
    fprintf(outputFile, "Kolmogorov distance: %lf\n", distance);
    printf("Quantile: %f\n", quantile);
    fprintf(outputFile, "Quantile: %f\n", quantile);
    return (distance < quantile);
}

double calcKolmogorovDistanceCauchy(double* sequence, int num, double location, double scale)
{
    double result = 0.0;
    for (int i = 0; i < num; i++)
    {
        result = max(result, abs(calcCauchyCDF(sequence[i], location, scale) - (double) (i + 1) / num));
    }
    return result;
}

double calcNormalCDF(double x, double mean, double variance, int erfStepNum)
{
    return (1.0 + calcErf((x - mean) / sqrt(2.0 * variance), erfStepNum)) / 2.0;
}

double calcErf(double max, int stepNum)
{
    return calcErf(0.0, max, stepNum);
}

double calcErf(double from, double to, int stepNum)
{
    double step = (to - from) / stepNum;
    double sum = (1.0 + exp(-to * to)) / 2.0;

    for (int i = 1; i < stepNum; i++)
    {
        sum += exp(-pow(from + i * step, 2.0));
    }

    return 2.0 * step * sum / sqrt(M_PI);
}

double calcLogNormalMean(double mu, double quadSigma)
{
    return exp(mu + quadSigma / 2);
}

double calcLogNormalVariance(double mu, double quadSigma)
{
    return (exp(quadSigma) - 1) * exp(2 * mu + quadSigma);
}

double calcLogNormalCDF(double x, double mu, double quadSigma, int erfStepNum)
{
    return (1.0 + calcErf((log(x) - mu) / sqrt(2.0 * quadSigma), erfStepNum)) / 2.0;
}

double calcCauchyMean(double location, double scale)
{
    return NAN;
}

double calcCauchyVariance(double location, double scale)
{
    return NAN;
}

double calcCauchyCDF(double x, double location, double scale)
{
    return atan((x - location) / scale) / M_PI + 0.5;
}

double calcEstimationMean(const double* source, int num)
{
    double sum = 0.0;

    for (int i = 0; i < num; i++)
    {
        sum += source[i];
    }

    return sum / num;
}

double calcEstimationVariance(const double* source, int num, double mean)
{
    double sum = 0.0;

    for (int i = 0; i < num; i++)
    {
        sum += pow(source[i] - mean, 2.0);
    }

    return sum / (num - 1);
}

const char* printBool(bool value)
{
    return (value) ? "true" : "false";
}
