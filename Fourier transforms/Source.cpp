#include <vector>
#include <cmath>

#include <math.h>
#include <boost/tuple/tuple.hpp>

#include "gnuplot-iostream.h"

using namespace std;
const int SAMPLES = 1000;
double pi = 3.1415926535897;

struct datas {
public:
	double a;
	double b;
	datas(double x, double y) {
		a = x;
		b = y;
	}
};

datas getVal(double in, int n, int k, int N) {
	double bn = (2 * pi*k*n / N);
	return datas(in*(cos(-bn)),in*sin(-bn));
}

vector<double> dft(vector<double> &in) {
	vector<double> ret;
	int N = in.size();
	for (int i = 1; i < N / 2; i++) {
		double suma = 0;
		double sumb = 0;
		for (int n = 0; n < N; n++) {
			datas d = getVal(in.at(n), n, i, N);
			suma += d.a;
			sumb += d.b;
		}
		double sum = sqrt(pow(suma, 2) + pow(sumb, 2));
		ret.push_back(sum / (N/2));
	}
	return ret;
}
const int range = 20;
int main() {
	
	Gnuplot gp;

	// Gnuplot vectors (i.e. arrows) require four columns: (x,y,dx,dy)+ sin(3*x) + sin(5*x)
	gp << "set samples " << SAMPLES << "\n";
	gp << "set xrange[0:" << range << "]\n";
	gp << "plot '-' with lines title 'de'\n";
	std::vector<double> xs;
	std::vector<double> samples;
	for (int i = 0; i < SAMPLES; i++) {
		double x = 0 + ((double)range / SAMPLES)*i;
		xs.push_back(x);
		samples.push_back(sin(x) + sin(3 * x) + sin(5 * x));
	}
	gp.send1d(boost::make_tuple(xs,
		samples));
	for (int i = 0; i < SAMPLES / 2; i++) {
		xs.at(i)=i;
	}
	vector<double> dat =dft(samples);

	std::vector<double>::iterator it;
	it = dat.begin();
	dat.insert(it, 0.0);
	xs.resize(SAMPLES / 2 );

	gp << "set terminal qt 1\n";
	gp << "set xrange[0:100]\n";
	gp << "plot '-' with lines title 'ptf'\n";
	gp.send1d(boost::make_tuple(xs,
		dat));
	std::cout << "Press enter to exit" << std::endl;
	std::cin.get();
}

