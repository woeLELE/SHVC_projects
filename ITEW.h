#pragma once
#include <vector>
#include <fstream>
#include <map>
#include <cmath>
using namespace std;

class ITEW
{
private:
    vector<double> pf;
    vector<vector<double>> wrec;
    const double a = 0.01;  // 学习率
    vector<int> rpf;
public:
    ITEW();
    ITEW(vector<double> pf_);
    double calD(int idx);
    vector<double> getW();
    double getP();
    bool meetCondition() noexcept;
};
