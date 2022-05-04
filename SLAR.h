#ifndef LAB06_SLAR_H
#define LAB06_SLAR_H

#include <iostream>
#include <cmath>
#include <fstream>

class SLAR
{
public:
    SLAR();
    SLAR(int rowAmm, int colAmm);
private:
    double **AMatrix;
    double **ATranMatrix;
    double **NMatrix;
    double **LMatrix;
    double **LTranMatrix;
    double *YMatrix;
    double *CMatrix;
    double *BMatrix;
    double *XMatrix;
    int order;
public:
    double **getAMatrix() {return AMatrix;}
    double *getBMatrix() {return BMatrix;}
    int getOrder() {return order;}
    void print();
    void readDataFromFile(std::string filepath);
    void printExtendedMatrix(double** matrA, double *matrB);
    void calcMinor(double **matrix, double **minor, int row, int col, int orderM);
    double calcDeterminant(double **matrix, int orderM);
    void swapRows(int row1, int row2, double** matrA, double *matrB);
    void addRow(int rowToAdd, int rowAdded, double mult, double** matrA, double *matrB);
    void calcTransposedMatrix(double **matrA, double **matrAT, int rowAm, int colAm);
    void multiplyMatrices(double **matr1, double **matr2, double **matr3);
    void printMatrix(double **matr, int rowAm, int colAm);
    void GaussianMethod(double **matrA, double *matrB, double *matrX);
    void SquareRootMethod();
};

#endif //LAB06_SLAR_H
