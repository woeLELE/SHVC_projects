#include "ITEW.h"
ofstream tt("tt.txt");
map<vector<int>, vector<double>> wp;
ITEW::ITEW()
{
    wrec.push_back({ 0, 0, 0 });
}
ITEW::ITEW(vector<double> pf_) : pf(pf_)
{
    rpf.assign(3, 0);
    for (int i = 0; i < 3; i++)
        rpf[i] = static_cast<int>(floor(pf[i] * 10));
    wrec.push_back(wp.count(rpf) ? wp[rpf] : vector<double>(3, 0));
}
double ITEW::calD(int idx)
{
    vector<double>& w = wrec.back();
    double a = 0, b = 0;
    for (int i = 0; i < 3; i++)
    {
        a += w[i] * pf[i];
        b += w[i] * (1 - pf[i]);
    }
    const double x = (2 * pf[idx] - 1) * exp(b);
    const double y = exp(a) + exp(b);
    return x / y;
}
vector<double> ITEW::getW()
{
    while (!meetCondition())
    {
        vector<double> w;
        for (int i = 0; i < 3; i++)
            w.push_back(wrec.back().at(i) + a * calD(i));
        wrec.push_back(w);
    }
    return wrec.back();
}
bool ITEW::meetCondition() noexcept
{
    if (wrec.size() < 3)
        return false;
    else if (wrec.size() >= 10000)
        return true;

    auto& w1 = *(wrec.end() - 2);
    auto& w2 = wrec.back();
    for (int i = 0; i < 3; i++)
        if (fabs(w1[i] - w2[i]) > 0.001)
            return false;
    return true;
}
double ITEW::getP()
{
    double a = 0, b = 0;
    for (int i = 0; i < 3; i++)
    {
        a += wrec.back().at(i) * pf.at(i);
        b += wrec.back().at(i) * (1 - pf.at(i));
    }
    if (!wp.count(rpf))
        wp[rpf] = wrec.back();

    tt << wrec.size() << "\t";
    for (auto n : wrec.back())
        tt << n << " ";
    tt << endl;
    return exp(a) / (exp(a) + exp(b));
}
