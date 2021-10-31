#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <string.h>
#include <random>
#include<cmath>
#include <map>
#include <stdlib.h>
#include "func.h"

using namespace std;



enum room_size{ LENGTH=10, WIDTH=6, T=1, KIDS=2, M=1, I=1, K_B=1, ETA_R=1 };

static int STEPS;
static double DT;
static int UNIT_TIME;


const double D=9;
const double A=1.1;
const double R0=0.7;
// parameters for soft potential
const double EPSI=9;
const double SIGMA=0.4;
const double N=2;

const double LAMBDA = 50;
//const double DT = 0.001;
const int CHUNK = 100000;
//const int STEPS = 2000000;
const double MAG = 1.0;
//const int UNIT_TIME = 200;
const int CUTOFF=40000;

void InitStatus(vector<vector3> &locations, vector<vector<vector3> > &velocities, vector<vector3> &accels, vector<vector3> &forces)
{
	std::ifstream file("../initial_position.dat");
	std::string line;
	vector<vector3> vels(CUTOFF);
	for (int i=0; i<KIDS; i++)
	{
		std::getline(file, line);
		std::istringstream is(line);
		if (line.size())
		{
			std::vector<std::string> tokens;
			copy(std::istream_iterator<std::string>(is),
		    std::istream_iterator<std::string>(),
		    back_inserter(tokens));
		    vector3 loc = { std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2])};
		    vector3 vel = { std::stod(tokens[3]), std::stod(tokens[4]), std::stod(tokens[5])};
		    vector3 acc = { std::stod(tokens[6]), std::stod(tokens[7]), std::stod(tokens[8])};
		    vector3 force = { std::stod(tokens[9]), std::stod(tokens[10]), std::stod(tokens[11])};
		    locations[i]= loc;
		    vels[0] = vel;
		    velocities.push_back(vels);
		    accels[i] = acc;
		    forces[i] = force;

		}
	}
};

std::vector<vector3> InitZeros()
{
	std::vector<vector3> vectors;
	for (int i = 0; i < KIDS; ++i)
	{
		vector3 kid = { 0, 0, 0 };
		vectors.push_back(kid);
	}
	return vectors;
};

std::vector<vector3> InitVel()
{
	std::vector<vector3> vel;
	vector3 kid = { 0, 0, 0 };
	vel.push_back(kid);
	return vel;
};

void check_1d_bound(double& loc_val, double& vel_val, int bound)
{
	if (loc_val < 0)
	{
		loc_val = - loc_val;
		vel_val = - vel_val;
	} else if (loc_val > bound)
	{
		loc_val = 2*bound - loc_val;
		vel_val = - vel_val;
	}
};

void check_boundary(vector3& loc, vector3& vel)
{
	check_1d_bound(loc.x, vel.x, LENGTH);
	check_1d_bound(loc.y, vel.y, WIDTH);
	while (true)
	{
		if (loc.a > M_PI)
		{
			loc.a -= 2*M_PI;
		} else if (loc.a < -M_PI)
		{
			loc.a += 2*M_PI;
		} else { break; }
	}
};

double ClctKinetic(const std::vector<vector3>& velocities)
{
	double kinetic = 0;
	for (int i = 0; i < KIDS; ++i)
	{
		double vx = velocities[i].x;
		double vy = velocities[i].y;
		double va = velocities[i].a;
		kinetic += 0.5*M*(vx*vx + vy*vy) + 0.5*I*va*va;
	}
	return kinetic;
}

void update_kid(vector3& loc, vector3& force, vector3& vel, vector3& acc)
{
	loc.x += vel.x*DT + 0.5*acc.x*DT*DT;
	loc.y += vel.y*DT + 0.5*acc.y*DT*DT;
	loc.a += vel.a*DT + 0.5*acc.a*DT*DT;
	check_boundary(loc, vel);
	// std::cout << vel.y << std::endl;

	vel.x += 0.5*(acc.x+force.x/M)*DT;
	vel.y += 0.5*(acc.y+force.y/M)*DT;
	vel.a += 0.5*(acc.a+force.a/I)*DT;

	acc.x = force.x/M;
	acc.y = force.y/M;
	acc.a = force.a/I;
};

double SoftSphere(double& ds)
{
	if (ds <= SIGMA)
	{
		return -N*EPSI*pow(SIGMA/ds, N)/ds;
	} else { return 0; }
};

vector2 ClctInter(vector2& distances)
{
	double ds = distances.val1;
	double da = distances.val2;
	double grad_Vt = -2.0 * A * D * ( exp(-A*(ds-R0)) - exp(-2*A*(ds-R0)) ) * cos(da)/ds;
	double grad_Vr = -D * ( exp(-2*A*(ds-R0)) - 2*exp(-A*(ds-R0)) ) * (-sin(da));
	double sphere_grad_V = SoftSphere(ds);
	vector2 grad_vals = { grad_Vt+sphere_grad_V, grad_Vr };
	return grad_vals;
};

void ReSetVector(std::vector<vector3>& interaction)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		interaction[kid].x = 0;
		interaction[kid].y = 0;
		interaction[kid].a = 0;
	}
};

void update_interaction(std::vector<vector3>& interaction, std::vector<vector3>& locations)
{
	// std::vector<vector3> interactions = InitZeros();
	ReSetVector(interaction);
	for (int kid = 0; kid < KIDS; ++kid)
	{
		int other = kid + 1;
		while (other < KIDS)
		{
			vector2 distances = ClctDistance(locations[kid], locations[other]);
			vector2 grad_vals = ClctInter(distances);
			double grad_Vt = grad_vals.val1;
			double grad_Vr = grad_vals.val2;
			double inter_x = -grad_Vt*(locations[kid].x - locations[other].x);
			double inter_y = -grad_Vt*(locations[kid].y - locations[other].y);
			double inter_a = -grad_Vr;
			interaction[kid].x += inter_x;
			interaction[kid].y += inter_y;
			interaction[kid].a += inter_a;
			interaction[other].x -= inter_x;
			interaction[other].y -= inter_y;
			interaction[other].a -= inter_a;
			other += 1;
		}
	}
};

void update_stochastic(std::vector<vector3>& stochastic)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		vector2 xy_norm = normal_distr(0, sqrt(2 * LAMBDA * K_B * T / DT));
		vector2 a_norm = normal_distr(0, sqrt(2 * ETA_R * K_B * T / DT));
		stochastic[kid].x = xy_norm.val1;
		stochastic[kid].y = xy_norm.val2;
		stochastic[kid].a = a_norm.val1;
	}
};

void update_viscosity(std::vector<vector3>& viscosity, const std::vector<vector3>& velocities)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		viscosity[kid].x = - LAMBDA*velocities[kid].x;
		viscosity[kid].y = - LAMBDA*velocities[kid].y;
		viscosity[kid].a = - ETA_R*velocities[kid].a;
	}
};

void update_forces(std::vector<vector3>& forces, const std::vector<vector3>& interaction, const std::vector<vector3>& stochastic, const std::vector<vector3>& viscosity)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		double tot_fx = interaction[kid].x + stochastic[kid].x - viscosity[kid].x;
		double tot_fy = interaction[kid].y + stochastic[kid].y - viscosity[kid].y;
		double tot_fa = interaction[kid].a + stochastic[kid].a - viscosity[kid].a;
		forces[kid].x = tot_fx;
		forces[kid].y = tot_fy;
		forces[kid].a = tot_fa;
	}
};

void output(const std::vector<vector3>& locations, double Ek)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		std::string index;
		std::stringstream num;
		num << kid;
		num >> index;
		std::ofstream save("info_" + index + ".csv", std::ios_base::app);
		save << locations[kid].x << "," << locations[kid].y << "," << locations[kid].a << std::endl;
	}
	ofstream save1("energy.csv", ios_base::app);
	save1 << Ek << endl;
};

vector3 clct_viscosity(int kid, const vector<double> &correlation, const vector<vector3> &vels, int t, const vector3 &cur_vel, map<int, vector<double> > &sum1000vels)
{
	vector3 vis = { 0, 0, 0 };
	double corr_val; int t_eff; int lngth = sum1000vels.size()-1;
	if (t >= CUTOFF)
	{
		t_eff = CUTOFF;
		for (int i = 0; i < lngth; ++i)
		{
			double mid_val = correlation[t-i*1000-500];
			vis.x += sum1000vels[i][2*kid]*mid_val*DT;
			vis.y += sum1000vels[i][2*kid+1]*mid_val*DT;
		}
		vis.x += sum1000vels[lngth][2*kid]*correlation[t]*DT;
		vis.y += sum1000vels[lngth][2*kid+1]*correlation[t]*DT;
	}
	else { t_eff = t; }

	for (int i = 0; i <= t_eff; ++i)
	{
		corr_val = MAG*(correlation[i]+correlation[i+1])/2;
		vis.x += corr_val*vels[(t-i)%CUTOFF].x*DT;
		vis.y += corr_val*vels[(t-i)%CUTOFF].y*DT;
	}
	vis.x += LAMBDA*cur_vel.x;
	vis.y += LAMBDA*cur_vel.y;
	vis.a += ETA_R*cur_vel.a;
	return vis;
}

void output_corr(const vector<double>& eta)
{
    ofstream save("gamma-data0.csv");
    int L = eta.size()/4;
    for (int i = 0; i < 5000; i ++)
    {
        double c = 0;
        int j;
        for (j = 0; i+j < L; j++)
        {
            c += eta[j]*eta[i+j];
        }
        save << i*DT << ',' << c/j << endl;
    }
}

void load_noise(vector<vector<vector2> > &eta, vector<vector2> &temp_eta, int time)
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		LoadNoise(temp_eta, kid, time, CHUNK);
		if (time == 0)
		{
			eta.push_back(temp_eta);
		} else
		{
			eta[kid] = temp_eta;
		}
	}
}

void update_status(vector<vector3>& locations, vector<vector3>& forces, vector<vector<vector3> >& velocities, vector<vector3>& accelerations, vector<double> &correlation)
{
	vector<vector3> interaction = InitZeros(); vector<vector3> stochastic = InitZeros(); vector<vector3> viscosity = InitZeros();

	vector<vector<vector2> > eta;
	vector<vector2> sum_vels;

	vector<vector3> temp_vels(KIDS);
	vector3 vel0;
	for (int kid = 0; kid < KIDS; ++kid)
	{
		vel0 = velocities[kid][0];
		temp_vels[kid] = vel0;
	}

	vector<double> eta_test;
	map<int, vector<double> > sum1000vels;
	vector<double> v(2*KIDS);
	// ofstream check("check.dat");
	int indx;
	vector<vector2> temp_eta(CHUNK);

	for (int i=0; i<STEPS; i++)
	{
		if (i%100000 == 0) cout << i << endl;

		if (i%UNIT_TIME == 0)
		{
			double Ek = ClctKinetic(temp_vels);
			output(locations, Ek);
		}

		if (i%CHUNK == 0) {load_noise(eta, temp_eta, i);}

		if (i >= CUTOFF-1)
		{
			indx = i - CUTOFF + 1;
			if (sum1000vels.count(indx/1000) > 0)
			{
				for (int kid = 0; kid < KIDS; ++kid)
				{
					sum1000vels[indx/1000][2*kid] += velocities[kid][(i+1)%CUTOFF].x;
					sum1000vels[indx/1000][2*kid+1] += velocities[kid][(i+1)%CUTOFF].y;
				}
			}
			else
			{
				for (int kid = 0; kid < KIDS; ++kid)
				{
					v[2*kid] = velocities[kid][(i+1)%CUTOFF].x;
					v[2*kid+1] = velocities[kid][(i+1)%CUTOFF].y;
				}
				sum1000vels[indx/1000] = v;
			}
		}

		for (int kid = 0; kid < KIDS; ++kid)
		{
			viscosity[kid] = clct_viscosity(kid, correlation, velocities[kid], i, temp_vels[kid], sum1000vels);
		}

		update_interaction(interaction, locations);
		update_stochastic(stochastic);

		for (int kid = 0; kid < KIDS; ++kid)
		{
			stochastic[kid].x += eta[kid][i%CHUNK].val1; // put sqrt(k T)
			stochastic[kid].y += eta[kid][i%CHUNK].val2;
		}
		update_forces(forces, interaction, stochastic, viscosity);

		for (int kid = 0; kid < KIDS; ++kid)
		{
			update_kid(locations[kid], forces[kid], temp_vels[kid], accelerations[kid]);
			velocities[kid][(i+1)%CUTOFF] = temp_vels[kid];
		}
	}
};

std::vector<double> LoadCorrelation()
{
	std::vector<double> correlation;
	std::ifstream infile("../corr.csv");
	double t, corr; char delimiter;
	for (int i = 0; i < CUTOFF+1; ++i)
	{
		infile >> t >> delimiter >> corr;
		correlation.push_back(corr);
	}
	return correlation;
}

void InitOutput()
{
	for (int kid = 0; kid < KIDS; ++kid)
	{
		std::string index;
		std::stringstream num;
		num << kid;
		num >> index;
		std::ofstream save("info_" + index + ".csv");
	}
	ofstream save1("energy.csv");
};

 void InitParameters(char* argv[])
 {
 	STEPS=atoi(argv[1]);
 	DT = atof(argv[2]);
 	UNIT_TIME = atoi(argv[3]);
 };

int main(int argc, char* argv[])
{

	InitParameters(argv);
	InitOutput();
	std::vector<double> correlation = LoadCorrelation();

	std::vector<vector3> forces(KIDS);
	std::vector<vector3> locations(KIDS);
	std::vector<vector3> accelerations(KIDS);
	vector<vector<vector3> > velocities;
	InitStatus(locations, velocities, accelerations, forces);

	update_status(locations, forces, velocities, accelerations, correlation);
	std::cout << "Completed!" << std::endl;
};
