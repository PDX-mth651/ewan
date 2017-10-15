///////CSRMat function implementations///////

template<class T>
void destroy(CSRMat<T> matrix)//destroys all arrays associated with the matrix.
{
  delete[]matrix.rowPoint;
  delete[]matrix.colIndex;
  if (matrix.data)
  {
    delete[]matrix.values;
  }
  if (matrix.talk)
  {
    if (matrix.data)
    {
      cout << "Freed memory for " << matrix.rows << " by " << matrix.cols << " matrix with " << matrix.nnz << " non-zero double precision entries.\n";
    }
    else
      {
        cout << "Freed memory for " << matrix.rows << " by " << matrix.cols << " binary matrix with " << matrix.nnz << " non-zero entries.\n";
      }
  }
  matrix.empty = true;
  matrix.data  = false;
  matrix.rows  = 0;
  matrix.cols  = 0;
  matrix.nnz   = 0;
	matrix.rowPoint = NULL;
	matrix.colIndex = NULL;
	matrix.values = NULL;
}
;

template<class T>
void initValues(CSRMat<T> matrix)//sets data to true. initializes values array and sets all values to zero.
{
  if (matrix.data)
  {
    delete[]matrix.values;
  }
  matrix.data = true;
  matrix.values   = new T[matrix.nnz];
  fill_n(matrix.values,matrix.nnz,0);
  if (matrix.talk)
  {
    if (matrix.data)
    {
      cout << "Allocated memory for " << matrix.rows << " by " << matrix.cols << " matrix with " << matrix.nnz << " non-zero double precision entries.\n";
    }
    else
    {
      cout << "Allocated memory for " << matrix.rows << " by " << matrix.cols << " binary matrix with " << matrix.nnz << " non-zero entries.\n";
    }
  }
  matrix.empty = false;
}
;

template<class T>
void writeEdgeList(CSRMat<T> matrix, string filepath, string option, bool includeData,bool sort)//writes an edgelist to the given (absolute) filepath. If option = "upper" or "lower", then only upper or lower triangular entries are recorded.  Otherwise, all entries are recorded. If sort, then entries of each row are sorted.
{
  if (sort)
  {
    sortRowsByCol(matrix);
  }
  ofstream outfile (filepath,ofstream::out);
  for (size_t i=0;i<matrix.rows;i++)
  {
    for (size_t j=matrix.rowPoint[i];j<matrix.rowPoint[i+1];j++)
    {
      if (
           ((option != "upper" ) && (option != "lower"))
        || ((option == "lower")&&(matrix.colIndex[j]<i))
        || ((option == "upper")&&(matrix.colIndex[j]>i))
             )
      {
        if (includeData)
        {
          outfile << i <<" "<<matrix.colIndex[j]<<" "<<matrix.values[j]<<"\r\n";
        }
        else
        {
          outfile << i <<" "<<matrix.colIndex[j]<<"\r\n";
        }
      }
    }
  }
  outfile.close();
}
;

template<class T>
void sortRowsByCol(CSRMat<T> matrix)//puts the entries of colIndex and values corresponding to each row in ascending order by collumn index
{
  for (size_t i=0; i<matrix.rows; i++)
  {
    size_t* sortedIndices = rankSort(matrix.colIndex,matrix.rowPoint[i],matrix.rowPoint[i+1]);
    permute(matrix.colIndex,sortedIndices,matrix.rowPoint[i],matrix.rowPoint[i+1]);
    if (matrix.data)
    {
      permute(matrix.values,sortedIndices,matrix.rowPoint[i],matrix.rowPoint[i+1]);
    }
  }
}
;

template<class T>
void print(CSRMat<T> matrix)//Prints matrix in compressed sparse row format in O(number nonzero entries).
{
    cout<<"Printing " << matrix.rows << " by " << matrix.cols
        << " matrix with " << matrix.nnz << " non-zero entries in CSR format.\n";

    cout << "I = ["  ;
    for (size_t j = 0; j <= matrix.rows ; j++ )
    {
      if (j == matrix.rows) {cout << matrix.rowPoint[j]<< "]\n";}
      else {cout << matrix.rowPoint[j] << ",";}
    }

    cout << "J = ["  ;
    for (size_t j = 0; j < matrix.nnz; j++ )
    {
      if (j == matrix.nnz-1) {cout << matrix.colIndex[j]<< "]\n";}
      else {cout << matrix.colIndex[j] << ",";}
    }

    cout << "D = ["  ;
    for (size_t j = 0; j < matrix.nnz; j++ )
    {
      if (j == matrix.nnz-1) {cout << matrix.getEntry(j)<< "]\n";}
      else {cout << matrix.getEntry(j) << ",";}
    }
}
;

template<class T>
void print_dense(CSRMat<T>& matrix)//Pretty prints the matrix in its full form in O(rows*cols) time. Only use "full" option for matrices under 100 by 100.
{
 cout<<"Printing " << matrix.rows << " by " << matrix.cols
		<< " matrix with " << matrix.nnz << " non-zero entries.\n";
	sortRowsByCol(matrix);
	size_t entry = 0;
	for (size_t i=0;i<matrix.rows;i++)
	{
	  entry=0;
	  for (size_t j=0;j<matrix.cols;j++)
	  {
		if (matrix.colIndex[matrix.rowPoint[i]+entry]==j && entry<(matrix.rowPoint[i+1]-matrix.rowPoint[i]))
		{
		  cout << matrix.getEntry(matrix.rowPoint[i] + entry++);
		}
		else
		{
		  cout <<"0";
		}
		if (j != matrix.cols-1)
		{
		  cout<<",";
		}
		cout << " ";
	  }
	  cout<<"\n";
	}
}
;

template<class T>
CSRMat<T> identity(CSRMat<T> matrix, string side)
{
  size_t n=0;
  if (side == "left"){n=matrix.rows;}
  else {n=matrix.cols;}
  
  CSRMat<T> id(n,n,n, false);
  for (size_t i=0; i<n; i++)
  {
    id.rowPoint[i] = i;
    id.colIndex[i] = i;
  }
  id.rowPoint[n] = n;
  return id;
}
;

template<class T, class S>
Vec<S> spMv(CSRMat<T> matrix,Vec<S> vector)// CSR Sparse Matrix Dense Vecor Product
{
  if (matrix.talk || vector.talk)
  {
    cout<< "calling spMv\n";
  }
  Vec<S> product(matrix.rows);
  fill_n(product.entry,matrix.rows,0);
  for (size_t k=0; k<matrix.rows; k++) {
    for (size_t j=matrix.rowPoint[k]; j<matrix.rowPoint[k+1] ; j++) {
      product.entry[k] += matrix.getEntry(j)*vector.entry[matrix.colIndex[j]];
    }
  }
  return product;
} 
;

template<class T>
CSRMat<T> spMt(CSRMat<T> matrix)// CSR Sparse Matrix Transpose
{
  if (matrix.talk)
  {
    cout<< "calling spMt\n";
  }
  //initialize csr matrix "transpose"
  CSRMat<T> transpose(matrix.cols, matrix.rows, matrix.nnz, matrix.data);
  
  //initialize counter   
  size_t* counter = new size_t[transpose.rows];
  fill_n(counter,transpose.rows,0);
  
  //count nnz's in transpose by rows
  for (size_t j=0; j<matrix.nnz ; j++)
  {
    transpose.rowPoint[matrix.colIndex[j]+1]++;
  }
  //Accumulate row pointers
  accumulate(transpose.rowPoint,0,transpose.rows+1);
  
  //compute collumn indices and corresponding values
  for (size_t i = 0; i<matrix.rows; i++)
  {
    for(size_t j = matrix.rowPoint[i]; j< matrix.rowPoint[i+1]; j++)
    {
      size_t l = transpose.rowPoint[matrix.colIndex[j]]+counter[matrix.colIndex[j]];
      transpose.colIndex[l] = i;
      if (matrix.data)
      {
        transpose.values[l] = matrix.values[j];
      }
      counter[matrix.colIndex[j]]++;
    }
  }
  delete[]counter;
  return transpose;
} 
;

template<class T>
CSRMat<T> spMM( CSRMat<T> matrix1,CSRMat<T> matrix2)//CSR sparse matrix product
{
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "calling spMM\n->initRowArrays\n\r";
  }
   //initialize row pointer array and flag array
  size_t* rowPoint  = new size_t[matrix1.rows+1];
  fill_n(rowPoint,matrix1.rows+1,0);
  size_t* flags     = new size_t[matrix2.cols];
  fill_n(flags,matrix2.cols,0);
  
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "initProduct\n";
  }

  size_t iCounter = 0;
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "rowPoint\n";
  }
  for (size_t i=0;i<matrix1.rows ; i++)
  {
    rowPoint[i] = iCounter;
    for (size_t jPtr = matrix1.rowPoint[i]; jPtr<matrix1.rowPoint[i+1]; jPtr++)
    {
      size_t j = matrix1.colIndex[jPtr];
      for (size_t lPtr = matrix2.rowPoint[j]; lPtr < matrix2.rowPoint[j+1]; lPtr++)
      {
        size_t l = matrix2.colIndex[lPtr];
        if (flags[l] != i+1)
        {
          iCounter++;
          flags[l] = i+1;
        }
      }
    }
  }
  rowPoint[matrix1.rows] = iCounter;
  
  //initialize CSRMat product and reinitialize flag array
  CSRMat<T> product(matrix1.rows, matrix2.cols, iCounter, true);
  product.rowPoint = rowPoint;
  fill_n(flags,product.cols,0);
  
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "init valTemp\n";
  }
  
  iCounter = 0;
  T* valTemp = new T[product.cols];
  for (size_t i=0;i<product.rows ; i++)
  {
    fill_n(valTemp,product.cols,0);
    for (size_t jPtr = matrix1.rowPoint[i]; jPtr < matrix1.rowPoint[i+1]; jPtr++)
    {
      size_t j = matrix1.colIndex[jPtr];
      T valA = matrix1.getEntry(jPtr);
      for (size_t lPtr = matrix2.rowPoint[j]; lPtr < matrix2.rowPoint[j+1]; lPtr++)
      {
        size_t l = matrix2.colIndex[lPtr];
        valTemp[l] += valA * matrix2.getEntry(lPtr);
        if (flags[l] != i+1 )
        {
          product.colIndex[iCounter] = l;
          iCounter++;
          flags[l] = i+1;
        } 
      }
    }
    for (size_t j = product.rowPoint[i]; j < product.rowPoint[i+1]; j++ )
    {
      product.values[j] = valTemp[product.colIndex[j]];
    }
  }
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "deleting flags/valTemp\n";
  }
  delete[]flags;
  delete[]valTemp;
  if (matrix1.talk || matrix2.talk)
  {
    cout<< "finished spMM\n";
  }
  return product;
} 
;

template<class T>
T dot(Vec<T> vector1,Vec<T> vector2)//dot product for Vec class
{
  size_t dim = std::min(vector1.size,vector2.size);
  double dot = 0; 
  for (size_t i = 0; i<dim; i++){
    dot += vector1.entry[i]*vector2.entry[i];
  }
  return dot;
}
;