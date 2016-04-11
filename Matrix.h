//
// Created by dumitru on 11.04.16.
//

#ifndef SIMPLEX_MATRIX_H
#define SIMPLEX_MATRIX_H


#include <vector>
#include <string>

class Matrix {
public:
    int n,m;
    std::vector<std::vector<double> > data;

    Matrix();
    Matrix(int n, int m);
    void set_all(double val);
    std::vector<int> size();
    void set(int i, int j, double val);
    void load_from_file(std::string src);
    Matrix transpose();
    Matrix vectorize();
    Matrix reshape(int rows, int cols);
    double maxx();
    double sum();
    Matrix slice(int a, int b, int c, int d);
    void randomize();



};


#endif //SIMPLEX_MATRIX_H
