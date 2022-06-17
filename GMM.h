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
	vector<double> rec0;
	vector<double> rec1;
	vector<double> rec;
	vector<double> rec_poss;
	vector<vector<double>> members;
public:
	GMM(vector<double> rec0, vector<double> rec1, double e0, double e1, double cost) : rec0(rec0), rec1(rec1), e0(e0), e1(e1), cost(cost) {}
	void print();
	void init();
	void train();
	double gauss(const double x, const double e, const double v);
	double possibility0();
	double getResult();
};
