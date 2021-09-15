#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
using namespace std;

class GMM
{
private:
    bool stop;
    double cost;
    double r0; // ILR类的概率
    double r1; // Intra的概率
    double e0; // 期望
    double e1;
    double v0;
    double v1;
    double N0;
    double N1;
    double pre;
    double cur;
    vector<double> rec0;
    vector<double> rec1;
    vector<double> rec;
    vector<double> rec_poss;
    vector<vector<double>> member;

public:
    GMM(vector<double> rec0, vector<double> rec1, double e0, double e1, double cost) : rec0(rec0), rec1(rec1), e0(e0), e1(e1), cost(cost)
    {
        stop = false;
        for (auto n : rec0)
            rec.push_back(n);
        for (auto n : rec1)
            rec.push_back(n);
        // rec.push_back(cost);
        v0 = 0.0;
        v1 = 0.0;
        N0 = rec0.size();
        N1 = rec1.size();
        if (N0 == 0 || N1 == 0)
        {
            stop = true;
            return;
        }
        r0 = N0 / (N0 + N1);
        r1 = N1 / (N0 + N1);
        for (auto e : rec0)
            v0 += (e - e0) * (e - e0);
        for (auto e : rec1)
            v1 += (e - e1) * (e - e1);
        v0 /= N0;
        v1 /= N1;
        if (v0 == 0 || v1 == 0)
        {
            stop = true;
            return;
        }
        member.assign(rec.size(), vector<double>(2, 0.0)); // 修改
        rec_poss.push_back(e0);
    }
    void print();
    //void init();
    void train();

    double gauss(const double x, const double e, const double v);

    double get_std(const double x, const double e, const double v);
    double gauss1(const double x);
    double possibility0();
    double getResult();
};
