#include <iostream>
#include <cmath>
#include <vector>
#include "Matrix.cpp"

using namespace std;

const int N = 2;
const double alpha = 1;
const double beta = 0.5;
const double gama = 2;
const double eps = 0.0001;

Matrix fx(Matrix &M) {
    int n = M.size()[0];
    Matrix R(n, 1);

    for (size_t i = 1; i <= n; i++) {
        double r = 100 * pow((M.get(i, 2) - M.get(i, 1) * M.get(i, 1)), 2) + pow(1 - M.get(i, 1), 2);
        R.set(i, 1, r);
    }
    return R;
}

int cmp(const vector<double> &A, const vector<double> &B) {
    int m = A.size();
    return A[m - 1] > B[m - 1];
}

bool cmpd(double a, double b,double eps)
{
    return abs(a-b) < eps;
}

int main() {
    srand(time(NULL));
    Matrix X(N + 1, N);
    Matrix y(N + 1, 1);
    X.randomize(0, 1);

    double progress = 1;
    while(progress >= 0.000000001)
    {

        y = fx(X);
        Matrix Data = X || y;


        Data.sort(cmp);
        X = Data.slice(1, N + 1, 1, N);
        y = Data.slice(1, N + 1, N + 1, -1);

        Matrix xh = X.slice(1, 1, 1, -1);
        double fh = fx(xh).to_double();
        double fhp = fh;

        Matrix xg = X.slice(2, 2, 1, -1);
        double fg = fx(xg).to_double();

        Matrix xl = X.slice(N + 1, N + 1, 1, -1);
        double fl = fx(xl).to_double();

        Matrix xc = X.slice(2, N + 1, 1, -1).col_sum() * ((double) 1 / N);
        double fc = fx(xc).to_double();

        Matrix xr = xc * (1 + alpha) - xh * alpha;
        double fr = fx(xr).to_double();

        int step = 0;
        if (fr < fl) {
            Matrix xe = xc * (1 - gama) + xr * gama;
            double fe = fx(xe).to_double();
            if (fe < fl) {
                xh = xe;
                fh = fx(xh).to_double();
            } else if (fe > fl) {
                xh = xr;
                fh = fx(xh).to_double();
            }
            step = 8;
        } else if (fl < fr && fr < fg) {
            xh = xr;
            fh = fx(xh).to_double();
            step = 8;
        } else if (fh > fr && fr > fg) {
            swap(xr, xh);
            swap(fr, fh);
            step = 5;
        } else if (fr > fh) {
            step = 5;
        }
        if (step == 5) {
            Matrix xs = xh * beta + xc * (1 - beta);
            double fs = fx(xs).to_double();

            if(fs<fh){
                xh = xs;
                fh = fx(xh).to_double();
                step = 8;
            }
            if(step != 8){
                if(fs>fh){
                    //shrink
                    for(int i=1;i<=N;++i)
                        for(int j=1; j<= N; ++j)
                            X.set(i,j,xl.get(1,j)+(X.get(i,j)-xl.get(1,j)/2));
                }
            }
        }
        X.set_row(xh,1);
        X.set_row(xg,2);
        X.set_row(xl,N+1);
        y = fx(X);
        progress = fhp-fh;
        cout<<fl<<" "<<fh<<" "<<progress<<endl;
    }
    cout<<X;
    assert(cmpd(X.get(1,1),1,eps));
    assert(cmpd(X.get(1,2),1,eps));



    return 0;
}