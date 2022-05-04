#include "SLAR.h"

SLAR::SLAR()
{
  order = 3;
  AMatrix = new double*[5];
  for(int i{0}; i < 5; i++)
  {
    AMatrix[i] = new double[3];
  }
  BMatrix = new double[5];
  XMatrix = new double [3];
  YMatrix = new double[3];
  CMatrix = new double[3];
  ATranMatrix = new double*[3];
  for(int i{0}; i < 3; i++)
  {
    ATranMatrix[i] = new double[5];
  }
  LMatrix = new double*[3];
  LTranMatrix = new double*[3];
  NMatrix = new double*[3];
  for(int i{0}; i < 3; i++)
  {
    LMatrix[i] = new double[3];
    LTranMatrix[i] = new double[3];
    NMatrix[i] = new double[3];
  }
}

SLAR::SLAR(int rowAmm, int colAmm)
        : order{colAmm}
{
  AMatrix = new double*[rowAmm];
  for(int i{0}; i < rowAmm; i++)
  {
    AMatrix[i] = new double[colAmm];
  }
  BMatrix = new double[rowAmm];
  XMatrix = new double[colAmm];
}

void SLAR::print()
{
  std::string divideLine(80, '-');
  std::cout << divideLine << std::endl;
  for(int i{0}; i < 5; i++)
  {
    for(int j{0}; j < order; j++)
    {
      if(AMatrix[i][j] >= 0)
      {
        std::cout << " + " << AMatrix[i][j] << "x[" << j + 1 << "] ";
      } else
      {
        std::cout << " - " << AMatrix[i][j] * (-1) << "x[" << j + 1 << "] ";
      }
    }
    std::cout << " = " << BMatrix[i] << std::endl;
  }
  std::cout << divideLine << std::endl;
}

void SLAR::readDataFromFile(std::string filepath)
{
  std::ifstream file(filepath);

  for(int i{0}; i < 5; i++)
  {
    for(int j{0}; j < order; j++)
    {
      file >> AMatrix[i][j];
    }
    file >> BMatrix[i];
  }
  file.close();
}

void SLAR::printExtendedMatrix(double** matrA, double *matrB)
{
  std::string divideLine(80, '-');
  std::cout << divideLine << std::endl;
  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      std::cout << matrA[i][j] << " ";
    }
    std::cout << "| " << matrB[i] << std::endl;
  }
  std::cout << divideLine << std::endl;
}

void SLAR::calcMinor(double **matrix, double **minor, int col, int row, int orderM)
{
  int ki, kj, di, dj;
  di = 0;
  for (ki = 0; ki < orderM - 1; ki++)
  {
    if (ki == col) di = 1;
    dj = 0;
    for (kj = 0; kj < orderM - 1; kj++)
    {
      if (kj == row) dj = 1;
      minor[ki][kj] = matrix[ki + di][kj + dj];
    }
  }
}

double SLAR::calcDeterminant(double **matrix, int orderM)
{
  int i, j, k, n;
  double det;
  double **p;
  p = new double*[orderM];
  for (i = 0; i < orderM; i++)
    p[i] = new double[orderM];
  j = 0; det = 0;
  k = 1;
  n = orderM - 1;
  if (orderM<1) std::cout << "Impossible to calc determinant" << std::endl;
  if (orderM == 1) {
    det = matrix[0][0];
    return(det);
  }
  if (orderM == 2) {
    det = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
    return(det);
  }
  if (orderM>2) {
    for (i = 0; i < orderM; i++) {
      calcMinor(matrix, p, i, 0, orderM);
      det = det + k * matrix[i][0] * calcDeterminant(p, n);
      k = -k;
    }
  }
  return(det);
}

void SLAR::swapRows(int row1, int row2, double** matrA, double *matrB)
{
  double temp[order+1];

  for(int i{0}; i < order; i++)
  {
    temp[i] = matrA[row1][i];
  }
  temp[order] = matrB[row1];
  for(int i{0}; i < order; i++)
  {
    matrA[row1][i] = matrA[row2][i];
  }
  matrB[row1] = matrB[row2];
  for(int i{0}; i < order; i++)
  {
    matrA[row2][i] = temp[i];
  }
  matrB[row2] = temp[order];
}

void SLAR::addRow(int rowToAdd, int rowAdded, double mult, double** matrA, double *matrB)
{
  for(int i{0}; i < order; i++)
  {
    matrA[rowToAdd][i] += matrA[rowAdded][i] * mult;
  }
  matrB[rowToAdd] += matrB[rowAdded] * mult;
}

void SLAR::GaussianMethod(double **matrA, double *matrB, double *matrX)
{
  std::cout << "************************* Gaussian Method **************************************"
            << std::endl;
  printExtendedMatrix(matrA, matrB);
  double biggest{0};
  int mainRow{0};
  for(int i{0}; i < order-1; i++)
  {
    biggest = NMatrix[i][i];
    mainRow = i;
    for(int j{i}; j < order; j++)
    {
      if(biggest < fabs(matrA[j][i]))
      {
        biggest = fabs(matrA[j][i]);
        mainRow = j;
      }
    }
    if(mainRow != i)
    {
      swapRows(mainRow, i, matrA, matrB);
    }
    double mult{0};
    for(int j{i+1}; j < order; j++)
    {
      if(matrA[i][i] != 0)
      {
        mult = matrA[j][i] / matrA[i][i];
        addRow(j, i, -mult, matrA, matrB);
      }
    }
    printExtendedMatrix(matrA, matrB);
  }

  matrX[order-1] = matrB[order-1]/matrA[order-1][order-1];

  for(int i{1}; i < order; i++)
  {
    double sub{0};
    for(int j{0}; j < order; j++)
    {
      sub+= matrX[order-1-j]*matrA[order-1-i][order-1-j];
    }
    matrX[order-1-i] = (matrB[order-1-i]-sub)/matrA[order-1-i][order-1-i];
  }


  std::cout << std::endl;
  for(int i{0}; i < order; i++)
  {
    std::cout << "X[" << i+1 << "] = " << matrX[i] << std::endl;
  }
}

void SLAR::calcTransposedMatrix(double **matrA, double **matrAT, int rowAm, int colAm)
{
  for(int i{0}; i < rowAm; i++)
  {
    for(int j{0}; j < colAm; j++)
    {
      matrAT[j][i] = matrA[i][j];
    }
  }
}

void SLAR::multiplyMatrices(double **matr1, double **matr2, double **matr3)
{
  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      matr3[i][j] = 0;
    }
  }

  for(int i{0}; i < order; i++)
  {
    for(int j{0}; j < order; j++)
    {
      for(int k{0}; k < 5; k++)
      {
        matr3[i][j] += matr1[i][k]*matr2[k][j];
      }
    }
  }
}

void SLAR::printMatrix(double **matr, int rowAm, int colAm)
{
  for(int i{0}; i < rowAm; i++)
  {
    for(int j{0}; j < colAm; j++)
    {
      std::cout << matr[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

void SLAR::SquareRootMethod()
{
  std::string divideLine(80, '-');
  std::cout << "Matrix A: " << std::endl;
  printMatrix(AMatrix, 5, 3);
  calcTransposedMatrix(AMatrix, ATranMatrix, 5, 3);
  std::cout << "Transposed Matrix A: " << std::endl;
  printMatrix(ATranMatrix, 3, 5);
  multiplyMatrices(ATranMatrix, AMatrix, NMatrix);
  std::cout << "Matrix N = Transposed A * A: " << std::endl;
  printMatrix(NMatrix, 3, 3);
  double NDet = calcDeterminant(NMatrix, 3);
  std::cout << "Determinant of Matrix N = " << NDet << std::endl;
  if(NDet > 0)
  {
    std::cout << "=== System is suitable ===" << std::endl;
  } else
  {
    std::cout << "~~~ System is unsuitable ~~~" << std::endl;
  }

  for(int i{0}; i < 3; i++)
  {
    for(int j{0}; j < 3; j++)
    {
      LMatrix[i][j] = 0;
    }
  }

  LMatrix[0][0] = sqrt(NMatrix[0][0]);
  LMatrix[1][0] = NMatrix[1][0] / LMatrix[0][0];
  LMatrix[2][0] = NMatrix[2][0] / LMatrix[0][0];
  LMatrix[1][1] = sqrt(NMatrix[1][1] - pow(LMatrix[1][0], 2));
  LMatrix[2][1] = (NMatrix[2][1] - LMatrix[2][0]*LMatrix[1][0])/LMatrix[1][1];
  LMatrix[2][2] = sqrt(NMatrix[2][2] - pow(LMatrix[2][0], 2) - pow(LMatrix[2][1], 2));

  calcTransposedMatrix(LMatrix, LTranMatrix, 3, 3);
  std::cout << "Matrix L:" << std::endl;
  printMatrix(LMatrix, 3, 3);
  std::cout << "Transposed Matrix L:" << std::endl;
  printMatrix(LTranMatrix, 3, 3);

  // Calculating C

  std::cout << "Matrix C: " << std::endl;

  for(int i{0}; i < 3; i++)
  {
    CMatrix[i] = 0;
  }

  for(int i{0}; i < 3; i++)
  {
      for(int k{0}; k < 5; k++)
      {
        CMatrix[i] += ATranMatrix[i][k]*BMatrix[k];
      }
      std::cout << CMatrix[i] << " ";
  }
  std::cout << std::endl;

  GaussianMethod(LMatrix, CMatrix, YMatrix);
  std::cout << divideLine << std::endl;
  for(int i{0}; i < 3; i++)
  {
    std::cout << "Y[" << i+1 << "] = " << YMatrix[i] << std::endl;
  }
  GaussianMethod(LTranMatrix, YMatrix, XMatrix);
  std::cout << divideLine << std::endl;
  std::cout << "System solution using Square Root Method: " << std::endl;
  for(int i{0}; i < 3; i++)
  {
    std::cout << "X[" << i+1 << "] = " << XMatrix[i] << std::endl;
  }

}

