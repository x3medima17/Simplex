//
// Created by dumitru on 11.04.16.
//

#include <assert.h>
#include "Matrix.h"
#include <vector>
#include <string>
#include <fstream>
#include <ostream>
#include <iostream>
#include <algorithm>

Matrix::Matrix() {

}

Matrix::~Matrix() {

}

Matrix::Matrix(int n, int m) {
    this->n = n;
    this->m = m;
    data.resize(n, std::vector<double>(m, 0));
}

double Matrix::get(int i, int j) {
    --i;
    --j;
    return data[i][j];
}

void Matrix::set_all(double val) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            data[i][j] = val;
}

template <typename Comparator>
void Matrix::sort(Comparator cmp) {
    std::sort(data.begin(), data.end(),cmp);
}

std::vector<int> Matrix::size() {
    std::vector<int> tmp(2);
    tmp[0] = n;
    tmp[1] = m;
    return tmp;
}

void Matrix::set(int i, int j, double val) {
    --i;
    --j;
    //std::cout<<size()[0]<<" "<<size()[1]<<std::endl;
    this->data[i][j] = val;
}

void Matrix::set_row(Matrix &M, int row) {
    for(int i=1;i<=m;++i){
        set(row,i,M.get(1,i));
    }
}

void Matrix::load_from_file(std::string src) {
    std::ifstream fin(src.c_str());
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            fin >> data[i][j];
}

Matrix Matrix::transpose() {
    Matrix tmp(m, n);
    tmp.n = m;
    tmp.m = n;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            tmp.data[i][j] = data[j][i];
    return tmp;
}

Matrix Matrix::vectorize() {
    Matrix tmp(n * m, 1);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            tmp.data[i * m + j][0] = this->data[i][j];
    return tmp;
}

Matrix Matrix::reshape(int rows, int cols) {
    Matrix tmp(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            tmp.data[i][j] = data[i * cols + j][0];
    return tmp;
}

double Matrix::maxx() {
    int maxx = -999999;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            if (data[i][j] > maxx)
                maxx = data[i][j];
    return maxx;
}

double Matrix::sum() {
    double s = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            s += data[i][j];
    return s;
}

Matrix Matrix::col_sum(){
    double s= 0;
    Matrix R(1,m);
    R.set_all(0);
    for(int i=1;i<=m;i++){
        R.set(1,i,slice(1,-1,i,i).sum());
    }
    return R;
};

Matrix Matrix::slice(int a, int b, int c, int d) {
    if (a == -1)
        a = 1;
    if (b == -1)
        b = n;
    if (c == -1)
        c = 1;
    if (d == -1)
        d = m;

    Matrix tmp(b - a + 1, d - c + 1);
    int col = 0, row = 0;
    for (int i = a; i <= b; i++) {
        col = 0;
        for (int j = c; j <= d; j++) {
            tmp.data[row][col] = data[i - 1][j - 1];
            col++;
        }
        row++;
    }
    return tmp;
}

void Matrix::randomize() {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            data[i][j] = ((double) rand() / (RAND_MAX));
}

double Matrix::d_random(double d_min, double d_max) {
    double d = (double) rand() / RAND_MAX;
    return d_min + d * (d_max - d_min);
}

void Matrix::randomize(double d_min, double d_max) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            data[i][j] = d_random(d_min, d_max);
}

double Matrix::to_double() {
    if (size()[0] == 1 && size()[1] == 1)
        return data[0][0];
    else
        std::cout << "Alarma!\n";
}


std::ostream &operator<<(std::ostream &os, const std::vector<int> &V) {
    std::cout << V[0] << " " << V[1] << std::endl;
}

std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<double>> &V) {
    int n = V.size();
    int m = V[0].size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            std::cout << V[i][j] << ' ';
        std::cout << std::endl;
    }

}

std::ostream &operator<<(std::ostream &os, const Matrix &M) {
    int n = M.n;
    int m = M.m;
    //cout<<n<<" "<<m<<endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            std::cout << M.data[i][j] << ' ';
        std::cout << std::endl;
    }
    return os;
}


Matrix operator+(const Matrix A, const Matrix B) {
    int n = (&A)->n;
    int m = (&A)->m;
    Matrix curr(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            curr.data[i][j] = A.data[i][j] + B.data[i][j];
    return curr;
}

Matrix operator-(const Matrix A, const Matrix B) {
    int n = (&A)->n;
    int m = (&A)->m;
    Matrix curr(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            curr.data[i][j] = A.data[i][j] - B.data[i][j];
    return curr;
}

template<typename T>
Matrix operator+(const Matrix A, const T val) {
    int n = A.n;
    int m = A.m;
    Matrix tmp(n, m);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            tmp.data[i][j] = A.data[i][j] + val;
    return tmp;
}


Matrix operator*(const Matrix A, const Matrix B) {
    int n = A.n;
    int m = B.m;
    assert(A.m == B.n);
    Matrix curr(n, m);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            for (int k = 0; k < A.m; k++)
                curr.data[i][j] += A.data[i][k] * B.data[k][j];

    return curr;


}

Matrix operator*=(const Matrix A, const Matrix B) {
    int n = A.n;
    int m = A.m;

    Matrix curr(n, m);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            curr.data[i][j] = A.data[i][j] * B.data[i][j];

    return curr;
}

template<typename T>
Matrix operator*(const Matrix A, const T val) {
    int n = A.n;
    int m = A.m;
    Matrix tmp(n, m);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            tmp.data[i][j] = A.data[i][j] * val;
    return tmp;
}

Matrix operator||(const Matrix &A, const Matrix &B) {
    int an = A.n;
    int am = A.m;
    int bn = B.n;
    int bm = B.m;
    assert(an == bn);
    Matrix tmp(an, am + bm);
    for (int i = 0; i < an; i++) {
        for (int j = 0; j < am; j++) {
            tmp.data[i][j] = A.data[i][j];
        }
        int k = 0;
        for (int j = am; j < am + bm; j++) {
            tmp.data[i][j] = B.data[i][k++];
        }

    }
    return tmp;
}


