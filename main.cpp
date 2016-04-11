#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cassert>


using namespace std;

const int N = 2;
const double alpha = 1;
const double beta = 0.5;
const double gama = 2;

struct point {
    vector<double> x;
    double y;
};

typedef vector<point> Simplex;

double fx(vector<double> &V) {
    return pow(100 * (V[1] - V[0] * V[0]), 2) + pow(1 - V[0], 2);
}

int cmp(const point &a, const point &b) {
    return a.y > b.y;
}


double d_random(double d_min, double d_max) {
    double d = (double) rand() / RAND_MAX;

    return d_min + d * (d_max - d_min);
}


void generate_points(vector<point> &V) {
    V.resize(N + 1);
    for (size_t i = 0; i < N + 1; i++) {
        V[i].x.resize(N);
        for (size_t j = 0; j < N; j++) {
            V[i].x[j] = d_random(-20, 20);
        }
        V[i].y = fx(V[i].x);
    }

}

void print_points(vector<point> &V) {
    for (int i = 0; i < V.size(); i++) {
        for (int j = 0; j < V[i].x.size(); j++)
            cout << V[i].x[j] << " ";
        cout << ": " << V[i].y << endl;
    }
}

point compute_center(vector<point> &V) {
    point result;
    result.x.resize(N);
    //init
    for (int i = 0; i < result.x.size(); i++) {
        result.x[i] = 0;
    }
    result.y = 0;
    for (int i = 1; i < result.x.size(); i++) {
        for (int j = 0; j < V[i].x.size(); j++)
            result.x[i] += V[i].x[j];
    }
    //divide
    for (int i = 0; i < result.x.size(); i++)
        result.x[i] /= N;
    return result;
}

point compute_mirror(point xh, point xc) {
    //multiply xc and xh
    for (int i = 0; i < xc.x.size(); i++) {
        xc.x[i] *= 1 + alpha;
        xh.x[i] *= alpha;
    }
    //substract
    for (int i = 0; i < xc.x.size(); i++) {
        xc.x[i] -= xh.x[i];
    }
    return xc;
}

point compute_xe(point xc, point xr) {
    //multiply xc and xr
    for (int i = 0; i < xc.x.size(); i++) {
        xc.x[i] *= (1 - gama);
        xr.x[i] *= gama;
    }
    //add
    for (int i = 0; i < xc.x.size(); i++) {
        xc.x[i] += xr.x[i];
    }
    cout << xc.x.size() << endl;
    return xc;
}

point compute_xs(point xh, point xc) {
    //multiply xc and xh
    for (int i = 0; i < xc.x.size(); i++) {
        xh.x[i] *= beta;
        xc.x[i] *= (1 - beta);
        xh.x[i] += xc.x[i];
    }
    return xh;
}

point shrink(point xi, point xl) {
    for (int i = 0; i < N; i++) {
        xi.x[i] = xl.x[i] + (xi.x[i]-xl.x[i])/2;
    }
    return xi;
}

int main() {

    vector<point> V;

    srand(time(NULL));
    generate_points(V);

    print_points(V);
    for(int k = 0;k<10;k++)
    {

        sort(V.begin(), V.end(), cmp);

        point xh = V[0];
        point xg = V[1];
        point xl = V.back();
        point xc = compute_center(V);
        point xr = compute_mirror(xh, xc);
        int l = V.size() - 1;

        double fr = fx(xr.x);
        double fh = fx(xh.x);
        double fg = fx(xg.x);
        double fl = fx(xl.x);

        //5
        if (fr < fl) {
            point xe = compute_xe(xc, xr);
            double fe = fx(xe.x);
            if (fe < fr) {
                xh = xe;
            }
            else if (fr < fe) {
                xh = xr;
            }
        }
        else if (fl < fr && fr < fg) {
            xh = xr;
        }
        else if (fg < fr && fr < fh) {
            swap(xr, xh);
            swap(fr, fh);
        }
//        assert(fl < fg && fg < fh && fh < fr);

        //6
        point xs = compute_xs(xh, xc);
        double fs = fx(xs.x);
        //7
        if (fs < fh) {
            xh = xs;
        }
            //8
        else if (fs > fh) {
            //shrink
            for (int j = 0; j < N; j++) {
                if (j != l) {
                    V[j] = shrink(V[j], V[l]);
                }
            }
        }

        V[0] = xh;
        V[1] = xg;
        V[V.size()-1] = xl;

        printf("iteration: %d   x=%f y=%f    x=%f y=%f    x=%f y=%f    fl=%f  fh=%f\n",k,xh.x[0],xh.x[1],xg.x[0],xg.x[1],xl.x[0],xl.x[1],fl,fh);

    }
    
    return 0;
}