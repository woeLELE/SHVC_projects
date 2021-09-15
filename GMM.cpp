#include "GMM.h"
//ofstream out("测试输出.txt");
void GMM::print()
{
    cout << "均值：" << e0 << " : " << e1 << endl;
    cout << "方差：" << v0 << " : " << v1 << endl;
    cout << "N0和N1：" << N0 << " : " << N1 << endl;
    cout << "初始概率：" << r0 << " : " << r1 << endl;
}

void GMM::train()
{
    if (stop)
        return;
    int idx = 0;
    for (auto n : rec)
        member[idx++][0] = (r0 * gauss(n, e0, v0)) / ((r0 * gauss(n, e0, v0)) + (r1 * gauss(n, e1, v1)));
    // member[idx][0] = (r0 * gauss(cost, e0, v0)) / ((r0 * gauss(cost, e0, v0)) + (r1 * gauss(cost, e1, v1)));

    idx = 0;
    for (auto n : rec)
        member[idx++][1] = (r1 * gauss(n, e1, v1)) / ((r0 * gauss(n, e0, v0)) + (r1 * gauss(n, e1, v1)));
    // member[idx][1] = (r1 * gauss(cost, e1, v1)) / ((r0 * gauss(cost, e0, v0)) + (r1 * gauss(cost, e1, v1)));

    N0 = 0.0;
    N1 = 0.0;
    double total0 = 0.0;
    double total1 = 0.0;
    for (int i = 0; i < rec.size(); i++)
    {
        N0 += member[i][0];
        N1 += member[i][1];
        total0 += member[i][0] * rec[i];
        total1 += member[i][1] * rec[i];
    }
    if (N0 == 0 || N1 == 0)
    {
        stop = true;
        return;
    }
    // 更新期望
    e0 = total0 / N0;
    e1 = total1 / N1;

    // 更新概率
    r0 = N0 / (N0 + N1);
    r1 = N1 / (N0 + N1);

    total0 = 0.0;
    total1 = 0.0;
    for (int i = 0; i < rec.size(); i++)
    {
        total0 += member[i][0] * (rec[i] - e0) * (rec[i] * e0);
        total1 += member[i][1] * (rec[i] - e1) * (rec[i] * e1);
    }
    // 更新方差
    v0 = total0 / N0;
    v1 = total1 / N1;
    if (v0 == 0 || v1 == 0)
    {
        stop = true;
        return;
    }
    rec_poss.push_back(possibility0());
    //out << rec_poss.back() << endl;
    if (rec_poss.back() <= 0)
    {
        stop = true;
        return;
    }
    if (rec_poss.size() > 2 && rec_poss.size() < 5)
    {
        pre = *(rec_poss.end() - 2);
        cur = rec_poss.back();
        if ((fabs(pre - cur)) / pre < 0.01)
            return;
    }
    else if (rec_poss.size() >= 5)
    {
        stop = true;
        return;
    }
    train();
}

double GMM::gauss(const double x, const double e, const double v)
{
    if (stop)
        return 1;
    if (v == 0)
    {
        stop = true;
        return 1;
    }
    return 1.0 / (sqrt(2 * 3.1415926) * sqrt(v)) * exp(-0.5 * (x - e) * (x - e) / v);
}

double GMM::get_std(const double x, const double e, const double v)
{
    if (stop)
        return 1;
    if (v == 0)
    {
        stop = true;
        return 1;
    }
    return gauss1((x - e) / sqrt(v));
}

double GMM::gauss1(const double x)
{
    if (stop)
        return 1;
    return 1.0 / (sqrt(2 * 3.1415926)) * exp(-0.5 * x * x);
}

double GMM::possibility0()
{
    double y = get_std(cost, e0, v0) + get_std(cost, e1, v1);
    if (y <= 0)
    {
        stop = true;
        return 1.0;
    }
    return get_std(cost, e0, v0) / (get_std(cost, e0, v0) + get_std(cost, e1, v1));
}

double GMM::getResult()
{
    if (stop)
    {
        return rec0.size() * 1.0 / (rec0.size() + rec1.size());
    }
    double last = rec_poss.back();
    //out << last << endl;
    double r = (rec0.size() * 1.0) / rec.size();
    return last + (1 - last) * r;
}