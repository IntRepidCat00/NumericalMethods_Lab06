#include "SLAR.h"

int main()
{
  SLAR MySystem;

  MySystem.readDataFromFile("D:\\Work\\University\\C1S2\\Numerical Methods\\Lab06\\System");
  MySystem.print();
  MySystem.SquareRootMethod();
  return 0;
}
