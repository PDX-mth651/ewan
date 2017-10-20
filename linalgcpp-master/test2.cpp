#include <exception>
#include <random>
#include <stdio.h>
#include <assert.h>
#include <fstream>

#include "src/linalgcpp.hpp"

using namespace linalgcpp;


SparseMatrix<double> ImportEdgeList(std::string filepath)
{
  CooMatrix<double> inCoo;
  std::ifstream iFile(filepath);
  
  int x;
  int y;
  double z;

  while (iFile >> x >> y >> z)
  {
    inCoo.Add(x,y,z);
  }
  
  auto inMat = inCoo.ToSparse();
  
  return inMat;
}

/*
  DenseMatrix<double> ImportCoordinates(std::string filepath)
{
  
  

}
*/
SparseMatrix<double> ImportElementMatrix(std::string element_matrix_filepath,std::string element_node_filepath )
{
  double e11;
  double e12;
  double e13;
  double e21;
  double e22;
  double e23;
  double e31;
  double e32;
  double e33;

  std::ifstream iFile1(element_matrix_filepath);
  std::ifstream iFile2(element_node_filepath);


  while (iFile1 >> e11 >> e12 >> e13 
                >> e21 >> e22 >> e23
                >> e31 >> e32 >> e33)
  {
    
  }
  
}

void test_import()
{
  SparseMatrix<double> A = ImportEdgeList(
    "/Users/Phism/Documents/School/HSAP/C++/testMat/tridiag_3m.txt");
  
  std::vector<int> rows({100, 2003, 100002});
  std::vector<int> cols({101, 2002, 100003});
  std::vector<int> marker(A.Rows(), -1);
  
  auto submat = A.GetSubMatrix(rows, cols, marker);
  submat.PrintDense("a submatrix of A");
  
  SparseMatrix<> AA = A.Mult(A);
  auto submat2 = AA.GetSubMatrix(rows, cols, marker);
  submat2.PrintDense("a submatrix of A squared");
}


int main(int argc, char** argv)
{
    test_import();
    return EXIT_SUCCESS;
}
