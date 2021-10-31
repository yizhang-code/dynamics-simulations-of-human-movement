//============================================================================
// Name        : TIME CORRELATED RANDOM NUMBER GENERATOR simulation
// Author      : yi
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================
#include <iostream>
#include <fstream>
#include <random>
#include <complex>
#include <cstdlib>
#include <complex>
#include <vector>

#include "func.h"

using namespace std;

int n = 1<<22;

double dt = 0.01;
double eps = 20;
double tau = 10;
double beta = 0.33;
double w0 = 0.001;
double w1 = 0.01;
int KIDS = 2;
// const double lambda = 1000;

typedef complex<double> complexd;

void swap(complexd &datai, complexd &dataj)
{
	complexd temp_k = datai;
	datai = dataj;
	dataj = temp_k;
}

void bit_reversal(vector<complexd> &eta)
{
    unsigned int j = 1;
    for (unsigned int i = 1; i < n; i+=1)
    {
        if (j > i)
        {
            swap(eta[j-1], eta[i-1]);
        }
        unsigned int m = n >> 1;
        while (m >= 1 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

void complex_fft(vector<complexd> &eta, int isign)
{
    bit_reversal(eta);
    unsigned int mmax=2, m, i, j;
    double theta;
    while (n >= mmax)
    {
        theta = isign*(2*M_PI/mmax);
        complexd e_itheta(cos(theta), sin(theta));
        complexd twiddle_factor(1, 0);
        for (m = 0; m < mmax/2; m += 1)
        {
            for (i = m; i < n; i += mmax)
            {
                j = i + mmax/2;
                complexd temp_twillde = twiddle_factor*eta[j];
                eta[j] = eta[i]-temp_twillde;
                eta[i] += temp_twillde;
            }
            twiddle_factor *= e_itheta;
        }
        mmax <<= 1;
    }
    for (unsigned int k = 0; k < n; ++k)
    {
        eta[k] /= (n*dt);
    }
}

double power_law(unsigned int mu)
{
    double w = 2/dt*sin(mu*M_PI/n);
    return eps/pow(w*tau+w0, (1-beta));
}

size_t sysrandom(void* dst, size_t dstlen)
{
    char* buffer = reinterpret_cast<char*>(dst);
    std::ifstream stream("/dev/urandom", std::ios_base::binary | std::ios_base::in);
    stream.read(buffer, dstlen);

    return dstlen;
}


void load_eta_w(vector<complexd> &eta)
{
    unsigned int mu;
    vector2 gauss_noise = normal_distr(0, 1);
    eta[0] = sqrt(n*dt*power_law(0))*gauss_noise.val1;
    for (mu = 1; mu <= n/2-1 ; ++mu)
	{
		double stddev = sqrt(n*dt*power_law(mu)/2);
        gauss_noise = normal_distr(0, 1);
        eta[mu] = complexd(stddev*gauss_noise.val1, stddev*gauss_noise.val2);
        eta[n-mu] = conj(eta[mu]);
	}
    gauss_noise = normal_distr(0, 1);
	eta[n/2] = sqrt(n*dt*power_law(n/2))*gauss_noise.val1;
}

void load_gamma_w(vector<complexd> &eta)
{
    unsigned int mu;
    eta[0] = power_law(0);
    for (mu = 1; mu <= n/2-1 ; ++mu)
    {
        double stddev = power_law(mu);
        eta[mu] = stddev;
        eta[n-mu] = stddev;
    }
    eta[n/2] = power_law(n/2);
}

void output_corr(const vector<complexd>& eta)
{
    ofstream save("gamma-data.csv");
    int L = eta.size()/4;
    for (int i = 0; i < 5000; i ++)
    {
        double c = 0;
        int j;
        for (j = 0; i+j < n; j++)
        {
            c += eta[j].real()*eta[i+j].real();
        }
        save << i*dt << ',' << c/j << endl;
    }
}

int main()
{
    vector<complexd> eta(n);

    int dim = 2;
    for (int kid = 0; kid < dim*KIDS; ++kid)
    {
        load_eta_w(eta);
        complex_fft(eta, 1);

        std::string index; std::stringstream num;
        num << kid;
        num >> index;

        ofstream save0("eta" + index + ".csv");
        for (int i = 0; i < n; ++i)
        {
            save0 << eta[i].real() << std::endl;
        }
    }
    // output_corr(eta);
    load_gamma_w(eta);
    complex_fft(eta, 1);
    ofstream save1("corr.csv");
    for (int i = 0; i < n/4; ++i)
    {
        save1 << i*dt << ',' << eta[i].real() << std::endl;
    }
	std::cout << "completed!" << std::endl;
};
