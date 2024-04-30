#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <vector>
#include <sstream>
using namespace std;

#define a1 1.0
#define a2 1.4
#define N1 512
#define N2 512
#define dat 1033
#define cyc 10

int main()
{
    int sample, sample_constphi, count, flg = 0;
    double phiJ, t, gamma, P, sxy, G1, G2;
    double phiJ_avg, t_avg[dat], gamma_avg[dat], P_avg[dat], sxy_avg[dat], G1_avg[cyc], G2_avg[cyc];
    // int dphi[3] = {-4, -3, -2};
    // int gamma0[5] = {-5, -4, -3, -2, -1};
    double omega_sgn[10] = {1.00, 1.26, 1.58, 2.00, 2.51, 3.16, 3.98, 5.01, 6.31, 7.94};
    int omega_exp[5] = {5, 4, 3, 2, 1};
    char filename[128];

    // initialize
    sample = 0;
    sample_constphi = 0;
    phiJ_avg = 0.0;
    for (int i = 0; i < dat; i++)
    {
        gamma_avg[i] = 0.0;
        P_avg[i] = 0.0;
        sxy_avg[i] = 0.0;
    }

    // calculate and output phiJ
    ifstream ifile_p("../OS_results/AQS/phiJ_AQS_bjam_dphi1.0e-04.dat");
    if (!ifile_p)
    {
        cout << "input file (phiJ) was not found" << endl;
        return 1;
    }
    string line;
    while (true)
    {
        getline(ifile_p, line);
        ifile_p >> phiJ;
        sample++;
        sample_constphi++;
        phiJ_avg += phiJ;
        if (ifile_p.eof())
        {
            break;
        }
    }
    cout << "#sample=" << sample << endl;
    phiJ_avg /= sample;
    ifile_p.close();

    ifstream ifile("../OS_results/AQS/stress_bjam_dphi1.0e-04.dat");
    if (!ifile)
    {
        cout << "input file (data) was not found" << endl;
        return 1;
    }
    for (int i = 0; i < sample; i++)
    {
        for (int j = 0; j < dat; j++)
        {
            ifile >> gamma >> P >> sxy;
            cout << "g=" << gamma << "\t"
                 << "P=" << P << "\t"
                 << "sxy=" << sxy << endl;

            gamma_avg[j] += gamma;
            P_avg[j] += P;
            sxy_avg[j] += sxy;
        }
    }
    ifile.close();

    for (int j = 0; j < dat; j++)
    {
        gamma_avg[j] /= (double)sample;
        P_avg[j] /= (double)sample;
        sxy_avg[j] /= (double)sample;
    }

    ofstream ofile("../OS_results/avg/avg_stress_AQS_bjam_dphi1.0e-04.dat");
    if (!ofile)
    {
        cout << "output file cannot be opened" << endl;
        return 1;
    }
    for (int i = 0; i < dat; i++)
    {
        ofile << gamma_avg[i] << "\t" << P_avg[i] << "\t" << sxy_avg[i] << endl;
    }
    ofile.close();

    return 0;
}