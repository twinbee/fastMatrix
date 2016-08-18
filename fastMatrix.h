#ifndef _FASTMATRIX_H
#define _FASTMATRIX_H

/*

a fastMatrix<T> is a very simple class used for doing basic matrix operations
the matrix uses doubles, is non-templatized, and stores matrix data 
in the free store ("heap"). 

authored Matthew Bennett 02-09-2007
                 modified 04-09-07. QR factorization methods
*/       


#include "bench.h"

#include <iostream>
#include <fstream>

#include <complex> //eventually we'd like to have complex entries

#include <cmath>


std::complex <double> root_of_unity (double n, double k)
{ //tested
 std::complex <double> i = (0,1);       
 static const double _E =  2.718281828459;
 static const double _PI = 3.141592653589;
 //return pow(CONST_E, (2*CONST_PI*i));
 return std::complex <double>(cos (2*_PI*k/n), sin(2*_PI*k/n));
}


/*
// 04-09-07 matrix classification FLAGS
//several flags that will let us create matrices having special props
#define fmIDENTITY        1 //an identity matrix
#define fmDIAGONAL        2 //a diagonal matrix
#define fmPOSITIVE        4 //a matrix w positive entries
#define fmSYMMETRIC       8 //a symmetric matrix
#define fmSKEWSYMMETRIC   16 //a skew-symmetric matrix
#define fmLEXOGRAPHIC     32 //a lexographic (1 2 3) matrix
#define fmVANDERMONDE     64 //a vandermonde matrix
#define fmHERMITIAN       128 //a hermitian matrix
#define fmPERMUTATION     256 //a permutation matrix
#define fmFOURIER         512 //a permutation matrix
#define fmORTHOGONAL      1024 //an orthogonal matrix
*/

enum fmTYPE {fmIDENTITY, fmDIAGONAL, fmSYMMETRIC, fmLEXOGRAPHIC,
      fmVANDERMONDE, fmORTHOGONAL, fmRANDOM, fmZERO, fmFOURIER };
      
      //not yet implemented:
  // fmSKEWSYMMETRIC, fmHERMITIAN, fmPERMUTATION, fmCIRCULANT};

//input / output formats. enumerated since they are independant
enum  outTYPE {TEXT, PLAINTEXT, MATLAB, MAPLE, TEX, LATEX, CSV};

//ways to do QR factorization
 enum {GRAM_SCHMIDT, HOUSEHOLDER /*, GIVENS */};

template <class T=float> //default is a real matrix of 64-bit floating point
class fastMatrix
{
  public:
   //only used for creating arrays 
   fastMatrix();

   //load a random matrix      
   fastMatrix(int rows, int cols);

   //create a special matrix based upon a series of sensible flags
  fastMatrix(int rows, int cols, int type);

   //load matrix from array of Ts in row-major order
   fastMatrix(int rows, int cols, T * init);

   //load a matrix from a file
   fastMatrix(char * filename);

   //deep copy constructor
   fastMatrix(const fastMatrix<T> & rhs);

   //special copy contructor
   // takes an array of fastMatrix of the same type and creates
   // a block-diagonal form of the matrices passed in
   fastMatrix(int n, fastMatrix<T> * rhs);

   ~fastMatrix();

   //transposition operations
   void TransposeSelf();
   fastMatrix<T> & Transpose();

   //fourier transform
   fastMatrix<T> & naiveDFT();
   fastMatrix<T> & iterativeFFT();
   fastMatrix<T> & recursiveFFT();

   //inverse operations
   void InvertSelf();
   fastMatrix<T> & Inverse();   
   
   //determinant operations
   T Determinant();

   //get #rows, #cols
   inline int getR() const {return m;}
   inline int getC() const {return n;}

   //benchmarking functions
   void registe(bench * b); //register a bench class instance

   //accessors make it easy to get the r,c elem. from outside the class
   //not used   //T operator() (unsigned row, unsigned col);
   T operator() (unsigned row, unsigned col) const;

   //LU factorization
   // see en.wikipedia.org/wiki/LU_decomposition
   void LUFactor(); //destructive operation. works "in place"
   
   //matrix property test operations
   bool isIdent(); //is it the identity matrix?
   bool isOrthogonal(); //in this matrix an orthogonal one?
   bool operator==(const fastMatrix<T> &rhs); //matrix-matrix test of equivalence   
   
   //QR factorization
   // see en.wikipedia.org/wiki/QR_decomposition
   //the following QR operations are also destructive
   fastMatrix<T> * QRFactor();

   //orthogonalize using a particular mode, or just do it
   void Orthogonalize(int mode);
   void Orthogonalize();

   //orthogonalization methods
   fastMatrix<T> & Householder();   
   void GramSchmidt();        //TODO: change so it is not destructive 7-09-07
   
   //vector extraction operations
   fastMatrix<T> & getRow(unsigned k);
   fastMatrix<T> & getCol(unsigned k);

   //matrix extractions
   fastMatrix<T> & lower(); //just the lower triangualr part
   fastMatrix<T> & upper(); //just the upper triangualr part

   //row / col replacement operations
   void setRow(int k, const fastMatrix<T> & rhs);
   void setCol(int k, const fastMatrix<T> & rhs);

   //output functions
   void print();
   void write(char * filename);
   void write(char * filename, int format);
   //for possible formats: see enum OUTTYPE at beginning of this file

   //matrix-matrix operations are NON-destructive
     //instead they return a new instance of matrix
   fastMatrix<T> & operator+(const fastMatrix<T> &rhs);
   fastMatrix<T> & operator-(const fastMatrix<T> &rhs);
   fastMatrix<T> & operator*(const fastMatrix<T> &rhs);

   //cross product for vectors. nondestructive
   fastMatrix<T> & operator%(const fastMatrix<T> & rhs);
   
   //2-norm (euclidean) of a vector
   T norm2();

   //scalar arithmetic operations
   fastMatrix<T> & operator*(T a);

   //DESTRUCTIVE assign-arithmetics (in place)
   void operator*=(const fastMatrix<T> &rhs);

   void operator*=(T s);
   void operator+=(T s); 
   void operator-=(T s);
   void operator/=(T s); //division, ie scaling of vector

   //deep copy
   //void operator=(const fastMatrix<T> &rhs);
   //todo: fix implementation of copy 04-09-07

   //destructive matrix assign-arithmetics
   void operator+=(const fastMatrix<T> &rhs);
   void operator-=(const fastMatrix<T> &rhs);

  //utility functions used by other functions
   //protected:
   fastMatrix<T> & minorMatrix(int row, int col);          

   //void registe(bench *b);

  private:          

   int m, n;             //rows, cols
   T * A;                //the internal matrix data

   bool benchmark; //whether or not we are doing benchmarking for this instance
   bench * cnt;    //pointer to an instance of the bench class
};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


 //template class fastMatrix<float>;
 //template class fastMatrix<float>;

 //template class fastMatrix<std::complex<float> >; 
 //template class fastMatrix<std::complex<double> >; 

 //this ^^ is idiotic code needed by the idiotic linker ld
 
 
//mnemonic for doing loops in row-major order
// for j<m (jim) for i<n (ian) and A[j*n+i] (jani)
// remember jim, ian, and janey

//also, this is an M-row N-column matrix

template <class T>
fastMatrix<T>::~fastMatrix()
{
 delete[] A;
}

template <class T>
fastMatrix<T>::fastMatrix()
{
 benchmark = false;
 m = 0;
 n = 0;
 A = 0;
}

template <class T>
fastMatrix<T>::fastMatrix(int rows, int cols)
{ //tested
 benchmark = false;
 m = rows;
 n = cols;
 A = new T[m*n];

 for (int i = 0; i < n*m; i++) A[i] =
 (((rand()/(RAND_MAX+1.0))+rand())/(RAND_MAX+1.0))
 *1.0; //random T (0.0, 1.0]
 //*1024.0-512.0; //random T btw -512 and 512
}

//fastMatrix< std::complex<float> > ::fastMatrix(int rows, int cols, int type)
template <class T>
fastMatrix<T>::fastMatrix(int rows, int cols, int type)
{
 //typedef std::complex <float> T;
 
 benchmark = false;
 m = rows;
 n = cols;
 A = new T[m*n];
 
 switch (type)
 { 
 case fmIDENTITY: //identity matrix is one of the most popular flavors

  for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
          A[j*n + i] = ((i == j) ? T(1.0) : T(0.0));
 break;

 case fmZERO: //identity matrix is one of the most popular flavors

  for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
          A[j*n + i] = T(0.0);
 break;

 case fmDIAGONAL: //midway flavor

  for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
       A[j*n + i] = ((i == j) ? T((((rand()/(RAND_MAX+1.0))+rand())/(RAND_MAX+1.0))) : T(0.0));
  break;

 case fmLEXOGRAPHIC: //counting flavor

  for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
          A[j*n + i] = i+j*n+1;
 break;

 case fmVANDERMONDE: //if you gotta ask, you aint never gonna know
 //this is the vandermonde of the lex vector by default
  for (int j = 0; j < m; j++)
  {
   T a = T(j+1);
   //T a = T((((rand()/(RAND_MAX+1.0))+rand())/(RAND_MAX+1.0)));
      for (int i = 0; i < n; i++)
          A[j*n + i] = pow(a,i);
  }
 break;

 case fmSYMMETRIC: //random symmetric matrix
  for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
       A[i*m + j] =
       A[j*n + i] = T((((rand()/(RAND_MAX+1.0))+rand())/(RAND_MAX+1.0)));
 break;

 case fmFOURIER: //matrix used for discrete fourier transform
  if (n != m) 
  {
   std::cerr << "Fourier matrix must be of square dimension\n";      
   m = n = 0;
   return;
  }
  for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
       A[i*m + j] = T( root_of_unity(n, i*j) / sqrt(n));
       
 break;


 default: //random matrix
  for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
       A[j*n + i] = T((((rand()/(RAND_MAX+1.0))+rand())/(RAND_MAX+1.0)));         
 break;
}//end switch

}


template <class T>
fastMatrix<T>::fastMatrix(int rows, int cols, T * s)
{ //tested
 //the input array s is in row major order                          
                           
 benchmark = false;
 m = rows;
 n = cols;
 A = new T[m*n];
  
 for (int i = 0; i < m*n; i++) A[i] = s[i];
}

template <class T>
fastMatrix<T>::fastMatrix(int numMat, fastMatrix<T> * rhs)
{ //tested

 bool debug = false;

    //special copy contructor
   // takes an array of fastMatrix of the same type and creates
   // a block-diagonal form of the matrices passed in

m = 0;
n = 0;

int rowoffset = 0;
int coloffset = 0;

for (int k = 0; k < numMat; k++ )
{
 m += rhs[k].getR();
 n += rhs[k].getC(); 
}

A = new T[m*n];

if (debug) std::cerr << "m: " << m << " n: " << n << std::endl;

for (int i = 0; i < m*n; i++) A[i] = T(0.0);


for (int k = 0; k < numMat; k++)
{
 int mm = rhs[k].getR();
 int nn = rhs[k].getC();
 
 for (int j = 0; j < nn; j++)
 {
  for (int i = 0; i < mm; i++)
  {
   A[(j+coloffset)*n + (i+rowoffset)] = rhs[k](j,i);
  }
 }
 
 rowoffset += mm;
 coloffset += nn;

}

}

template <class T> bool fastMatrix<T>::isIdent()
{ //tested
 //this function is very picky. needs to be rewritten to allow for some error

 if (n != m) return false; //identity must be square

 for (int j = 0; j < m; j++)
 { 
  for (int i = 0; i < n; i++) 
  {
   if (i != j && 
   A[j*n + i] != T(0.0) ) return false;
   if (i == j && 
   A[j*n + i] != T(1.0) )return false;
  }
 }

 return true;
}



template <class T>
fastMatrix<T>::fastMatrix(char* filename)
// the array must be in row major order
{ //tested
 benchmark = false;
 std::ifstream f;

 f.open(filename);

 f >> m; //num rows
 f >> n; //num cols
 
 A = new T[m*n]; //the matrix to decompose

 for (int i = 0; i < m*n; i++) f >> A[i]; //loads in row major order

 f.close();
}

template <class T>
void fastMatrix<T>::print()
{ //tested
 bool nline = false;
 bool print2 = true;
 
 if (nline) for (int i = 0; i < m*n; i++) std::cout << A[i] << " ";

 else if (print2) //test the accessor and print
 {
  for (int i = 0; i < m; i++)
  {
   for (int j = 0; j < n; j++)
   {
    std::cout << this->operator()(i,j) << "\t"; //row major order
   }
   std::cout << "\n";
  } 
 }
 
 else //original print routine
 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   std::cout << A[j*n+i] << "\t"; //row major order
  }
  std::cout << "\n";
 }
 
 
 
 return;
}

template <class T>
void fastMatrix<T>::write(char * filename)
{ //tested
 std::ofstream f;

 f.open(filename);

 f << n << " " << m << "\n";

 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   f << A[j*n+i] << "\t"; //row major order
  }
  f << "\n";
 }

 f.close();
 return;
}

template <class T>
void fastMatrix<T>::write(char * filename, int format)
{ 
 static char var = 'a';
 
 std::ofstream f;

 f.open(filename, std::ios::app);

 switch (format)
 {
  case MATLAB: //tested
   f << var++ << " = [";      
 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   f << A[j*n+i] << " ";
  }
  f << ";";
 }    
 f << "];\n";
  break;

  case MAPLE: //tested
   f << var++ << " := matrix([";      
 for (int j = 0; j < m; j++)
 {
  f << "[";
  for (int i = 0; i < n; i++)
  {
   f << A[j*n+i];
   if (i < n-1) f << ",";
  }
  if (j < m-1) f << "],";
  else f << "]";
 }    
 f << "]);\n";
  break;

  case LATEX: //tested
  case TEX:       
 f << "$$" << var++ << " = \\left( \\begin{array}{";
 for (int i = 0; i < n; i++) f << "c";
 f << "}";
 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   f << A[j*n+i];
   if (i < n-1) f << " & ";
  }
  f << "\\\\\n";
 }    
 f << "\\end{array} \\right) $$\n";  
  break;

  case CSV: //tested
 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   f << A[j*n+i] << ",";
  }
  f << "\n";
 }
  break;  
  
  default:
  case TEXT: //tested
  case PLAINTEXT:
 f << n << " " << m << "\n";

 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   f << A[j*n+i] << "\t"; //row major order
  }
  f << "\n";
 }
  break;
 } // end switch on format

 f.close();
 return;
}


template <class T>
void fastMatrix<T>::LUFactor()
{
 if (m != n)
 {
  std::cerr << "Cannot LU-Factor a non-square matrix\n";
  return; //don't factor non-square matrices
 }

 for (int k = 0; k < n-1; k++)
 {
  for (int j = k+1; j < n; j++)
  {
   A[n*j + k] /= A[n*k + k];
   if (benchmark) cnt->mul();
   for (int i = k+1; i < n; i++)   
   {
    A[n*j + i] -= A[n*k + i]*A[n*j + k];
    if (benchmark)
    {
     cnt->mul();
     cnt->add();    
    } //count the benchmarked operations
   } //end for i
  } //end for j
 } //end for k
 
 return;
}

template <class T>
T fastMatrix<T>::operator() (unsigned r, unsigned c) const
{ //tested
 if (c >= n)
 {
  std::cerr << "Matrix col index out of bounds on accessor\n";
  return 0.0f;        
 }
 
 if (r >= m)
 {
  std::cerr << "Matrix row index out of bounds on accessor\n";
  return 0.0f;        
 }

 return A[r*n + c];
}

template <class T>
fastMatrix<T> & fastMatrix<T>::upper()
{ //tested
 T * data = new T[n*m];

 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   data[j*n + i] = ((j <= i) ? A[j*n + i] : 0 );
  }
 }

 fastMatrix * result = new fastMatrix(m, n, data);
 return *result;
}

template <class T>
fastMatrix<T> & fastMatrix<T>::lower()
{ //tested
 T * data = new T[n*m];
 for (int i = 0; i < n*m; i++) data[i] = 0;

 for (int j = 0; j < m; j++)
 {
  for (int i = 0; i < n; i++)
  {
   data[j*n + i] = ((j > i) ? A[j*n + i] : ((i == j) ? 1 : 0 ));
  }
 }

 fastMatrix * result = new fastMatrix(m, n, data);
 return *result;
}

template <class T>
fastMatrix<T> & fastMatrix<T>::operator+(const fastMatrix<T> &rhs)
{
 T * data = new T[n*m];

 for (int j = 0; j < n; j++)
 {
  for (int i = 0; i < m; i++)
  {
   data[n*j+i] = A[n*j+i] + rhs(i,j);
    if (benchmark)
    {
     cnt->add();    
    }
  }
 }

 fastMatrix * result = new fastMatrix(m, n, data);
 return *result;
}

template <class T>
fastMatrix<T> & fastMatrix<T>::operator-(const fastMatrix<T> &rhs)
{
 T * data = new T[n*m];

 for (int j = 0; j < n; j++)
 {
  for (int i = 0; i < m; i++)
  {
   data[n*j+i] = A[n*j+i] - rhs(i,j);
    if (benchmark)
    {
     cnt->add();    
    }
  }
 }

 fastMatrix * result = new fastMatrix(m, n, data);
 return *result;
}

template <class T>
fastMatrix<T> & fastMatrix<T>::operator%(const fastMatrix &rhs)
{
 // vector cross product returns a row vector
 // not yet implemented
 /*
 T data [n];
 
 if (this->getR() == 1 && rhs.getC() == 1)
 {      
  
 }

 else if (this->getC() == 1 && rhs.getR() == 1)
 {
      
      
 }
 
 else 
 {
  std::cerr << "Vector cross product is only defined on a 1row x 1col combo\n";
  return *(new fastMatrix<T> (1,n));
 }
 */
}

template <class T>
fastMatrix<T> & fastMatrix<T>::operator*(const fastMatrix &rhs)
{ //tested
 //04-09-07 version 0.2 now does arbitrary-sized multiplications
 // so long as the inner dimensions match

 bool debug = false;

 int ln = this->getC();
 int rn = rhs.getC();

 int lm = this->getR();
 int rm = rhs.getR();
           
 if (ln != rm) //if left rows and right cols dont match
 {
  std::cerr << "Matrix dimensions mismatch on multiplication\n";      
  std::cerr << m << "x" << n <<" times "<<rhs.getR() << "x" << rhs.getC() <<"\n";
  //return;
 }
 
 T * data = new T[lm*rn];

 if (debug) std::cerr << lm << "x" << rn << "created\n";

 //clear the matrix
 for (int i = 0; i < lm*rn; i++) data[i] = 0.0f;

 for (int i = 0; i < lm; i++) //left rows
 {
  for (int j = 0; j < rn; j++) //right cols
  {

    for (int k = 0; k < ln; k++)
    {    
     data[i*lm+j] += (*this)(i,k)*rhs(k,j);
    }
 }
}
 
 fastMatrix * result = new fastMatrix(lm, rn, data);
 return *result;
}

template <class T>
void fastMatrix<T>::operator-=(T s)
{ //tested
 for (int i = 0; i < m*n; i++) A[i] -= s;
}

template <class T>
void fastMatrix<T>::operator+=(T s)
{ //tested
 for (int i = 0; i < m*n; i++) A[i] += s;
}

template <class T>
void fastMatrix<T>::operator+=(const fastMatrix<T> &rhs)
{//tested
 if (m != rhs.getR() || n != rhs.getC())
 {
  std::cerr << "Matrix dimension mismatch on assign-addition. Aborting";
  return;
 }
 
 for (int j = 0; j < n; j++)
 {
  for (int i = 0; i < m; i++)
  {
   A[i*n + j] += rhs(i,j);
  }
 }
 
}

template <class T>
void fastMatrix<T>::operator-=(const fastMatrix<T> &rhs)
{//tested
 if (m != rhs.getR() || n != rhs.getC())
 {
  std::cerr << "Matrix dimension mismatch on assign-subtraction. Aborting";
  return;
 }
 
 for (int j = 0; j < n; j++)
 {
  for (int i = 0; i < m; i++)
  {
   A[i*n + j] -= rhs(i,j);
  }
 }
 
}

template <class T>
fastMatrix<T> & fastMatrix<T>::minorMatrix(int row, int col)
{ //tested
 bool debug = false;
 //returns the submatrix below i,j.
 //used primarily to assist the house holder transform but might have other uses            
 if (row >= m || col >=n)
 {
  std::cerr << "Index out of bounds on matrix Minor call.\n";
 }

 T * data = new T[(m-row)*(n-col)];

 for (int j = row; j < m; j++)
 {
  for (int i = col; i < n; i++)     
  {
   T as = A[j*n + i];
   data[(j-row)*(n-col) + (i-col)] = as;
   if (debug) std::cerr << (j-row)*(n-col) + (i-col) << std::endl;
  }
   if (debug)    std::cerr << "\n";  
 }

 fastMatrix * result = new fastMatrix((m-row), (n-col), data);

 return (*result);
}

/*
template <class T>
void fastMatrix<T>::operator=(const fastMatrix<T> &rhs)
{
 //deep copy operator
 delete[] A;
 
 m = rhs.getR();
 n = rhs.getC();
 
 A = new T[m*n];
 
 for (int j = 0; j < n; j++)
 {
  for (int i = 0; i < m; i++)
  {
   A[i*n + j] = rhs(i,j);
  }
 }
 
}
*/


template <class T>
fastMatrix<T>::fastMatrix(const fastMatrix<T> &rhs)
{
 bool debug = false;

 if (debug) std::cerr << "copy constructor bein called\n";

 //deep copy constructor
 //delete[] A;
 
 m = rhs.getR();
 n = rhs.getC();
 
 A = new T[m*n];
 
 for (int j = 0; j < m; j++) //row
 {
  for (int i = 0; i < n; i++) //col
  {
   A[j*n + i] = rhs(j,i);
  }
 }
 
}


template <class T>
void fastMatrix<T>::operator*=(const fastMatrix &rhs)
{ 
 //this multiply and assign may change the dimensions of the matrix
  // destructive. use with extreme caution    
  bool debug = false;
  
 int ln = this->getC();
 int rn = rhs.getC();

 int lm = this->getR();
 int rm = rhs.getR(); 
     
 if (ln != rm) //if left rows and right cols dont match
 {
  std::cerr << "Matrix dimensions mismatch on assign-mult. Aborting.\n";      
  std::cerr << m << "x" << n <<" times "<<rhs.getR() << "x" << rhs.getC() <<"\n";
  return;
 }
 
 T * data = new T[lm*rn];

 if (debug) std::cerr << lm << "x" << rn << "created\n";

 for (int i = 0; i < lm*rn; i++) data[i] = 0.0f;

 //I DO NOT HAVE ANY RATIONAL FOR NAMING THE VARIABLE NAMES AS THEY ARE
 // IT JUST WORKS ok? 04-09-07

 for (int j = 0; j < lm; j++) //left rows
 {
  for (int i = 0; i < rn; i++) //right cols
  {
   for (int r = 0; r < ln; r++) //iterations
   {
    if (debug) std::cerr <<  " + " << A[j*n+r] << "x" << rhs(r,i);
    data[j*lm+i] += A[j*n+r]*rhs(r,i);
   }
   if (debug) std::cerr << "\n";
  }
 }

 m = lm;
 n = rn;
 
 delete[] A;
 
 A = new T[lm*rn];
 
 for (int i = 0; i < m*n; i++) A[i] = data[i]; 
 
}

//scalar assign-multiplication 
template <class T>
void fastMatrix<T>::operator*=(T a)
{ //tested
 //this multiply and assign is
  // destructive. use with extreme caution    

 for (int i = 0; i < m*n; i++) A[i] *= a;
} 

template <class T>
void fastMatrix<T>::operator/=(T a)
{ //tested
 //this multiply and assign is
  // destructive. use with extreme caution    

 if (a == T(0.0))
 {
  std::cerr << "Scale by zero. Aborted.\n";
  return;
 }
 
 for (int i = 0; i < m*n; i++) A[i] /= a;
} 

//scalar multiplication
template <class T>
fastMatrix<T> & fastMatrix<T>::operator*(T a)
{ //tested
 //this multiply and assign is
  // destructive. use with extreme caution    
 T * data = new T[n*m];

 for (int i = 0; i < m*n; i++) data[i] = A[i] * a;
 
 fastMatrix * result = new fastMatrix(m,n,data);
 return *result;
} 

template <class T> 
void fastMatrix<T>::registe(bench *b)
{
 benchmark = true;
 cnt = b;
}

template <class T>
bool fastMatrix<T>::operator==(const fastMatrix<T> &rhs)
{ //tested
 bool debug = false;

 if (rhs.getR() != this->getR() || rhs.getC() != this->getC() ) 
 {
  std::cerr << "Dimension mismatch on equality\n";
  return false;
 }

 for (int j = 0; j < m; j++)
 {
 for (int i = 0; i < n; i++)
  {
   if (A[i*n + j] != rhs(i,j))
   {
    if (debug) std::cerr << A[i*n + j] << " (" << i << "," << j << ") " << rhs(i,j) << "\n";
    return false;
   }
  }
 }

 return true;
}

template <class T>
fastMatrix<T> & fastMatrix<T>::Transpose()
{ //tested
 //this transpose may change the dimensions of the matrix
 bool debug = false;
 
 T * data = new T[n*m];

 if (debug) std::cerr << n << "x" << m << "created\n";

 for (int j = 0; j < m; j++) //left rows
 {
  for (int i = 0; i < n; i++) //right cols
  {
   data[i*m + j] = A[j*n + i];
  }
 }
 
 fastMatrix * result = new fastMatrix (n,m,data);
 
 return *result;
}

template <class T>
void fastMatrix<T>::TransposeSelf()
{ //tested
 //this transpose may change the dimensions of the matrix
 bool debug = false;
 
 T * data = new T[n*m];

 if (debug) std::cerr << n << "x" << m << "created\n";

 for (int j = 0; j < m; j++) //left rows
 {
  for (int i = 0; i < n; i++) //right cols
  {
   data[i*m + j] = A[j*n + i];
  }
 }
 
 int temp = m;
 m = n;
 n = temp;
 
 for (int i = 0; i < m*n; i++) A[i] = data[i]; 
 
 //return *this; 
 
 //not appropriate! vulgar! obscene!
}

template <class T>
fastMatrix<T> * fastMatrix<T>::QRFactor()
{
 //returns an array of fastMatrix
 // first is Q (orthogonal)
 // second is R (upper triangular)


}

template <class T>
fastMatrix<T> & fastMatrix<T>::Householder()
{
 //idea taken from wikipedia
  // http://en.wikipedia.org/wiki/QR_decomposition#Connection_to_a_determinant_or_a_product_of_eigenvalues
 
 bool debug = false;
 bool naive = true;
 
 
if (naive) //reference implementation. extremely inefficient
{
 int nn = (m < n) ? m : n;

 fastMatrix<T>  I(m, n, fmIDENTITY);

 fastMatrix<T>  v(1,n); //vector used in each step

 fastMatrix<T>  Q[nn];

 fastMatrix<T>  current(*this);
 
 for (int i = 0; i < nn-1; i++)
 {
  v = (current.getCol(0)).Transpose();
  
  fastMatrix<T> ident(nn-i,nn-i, fmIDENTITY); 
 
  v -= ident.getRow(0)*v.norm2();
  v /= v.norm2();
  
  if (debug)
  {
   std::cout << "v[" << i << "]\n";
   v.print();   
  }

 fastMatrix<T> mini_ident(i+1,i+1, fmZERO); 

 fastMatrix<T> H = (v.Transpose()*v*T(2.0));

 fastMatrix<T> join[2] = {mini_ident, H};
 fastMatrix<T> hottt(2,join); //create a block diagonal so we can get Q[i]

  Q[i] = I - hottt;
  

  if (debug)
  {
   std::cout << "Q[" << i << "]:\n";
   Q[i].print();
  }
  
  if (debug)
  {
   std::cout << "Q[" << i << "]A:\n";
   (Q[i]*(*this)).print();
  }
  
  current = (Q[i]*(*this)).minorMatrix(i+1,i+1);
  //if(debug) current.print();
 }
 
  fastMatrix<T> *result = new fastMatrix(I);

  for (int i = 0; i < nn-1; i++)
  {
   *result *= Q[i]; //create the R matrix
  }
  
  *result *= (*this);
  
  return *result;
} //end naive reference implementation


}


template <class T>
fastMatrix<T> & fastMatrix<T>::getCol(unsigned k)
{ //tested
 if (k >= n)
 {
  std::cerr << "Column index out of range on vector extraction\n";      
 }

 T * data = new T[m*1];
  
 for (int i = 0; i < m; i++)
 {
  data[i] = A[k + i*n];
 }

 fastMatrix<T> * C = new fastMatrix<T>(m,1,data);
 return *C;
}

template <class T>
fastMatrix<T> & fastMatrix<T>::getRow(unsigned k)
{ //tested
 if (k >= m)
 {
  std::cerr << "Row index out of range on vector extraction\n";      
 }

 T * data = new T[n*1];
  
 for (int i = 0; i < n; i++)
 {
  data[i] = A[k*m + i];
 }

 fastMatrix<T> * C = new fastMatrix<T>(1,n,data);
 return *C;
}

template <class T>
void fastMatrix<T>::setRow(int k, const fastMatrix<T> & rhs)
{
 if (rhs.getR() != 1)
 {
  std::cerr << "Cannot set row for a non-vector rhs\n";               
  return;
 }
 
 for (int i = 0; i < n; i++)
 {
  A[k*m + i] = rhs(0,i);
 } 
 
}

template <class T>
void fastMatrix<T>::setCol(int k, const fastMatrix<T> & rhs)
{
 if (rhs.getC() != 1)
 {
  std::cerr << "Cannot set col for a non-vector rhs\n";               
  return;
 }
 
 for (int i = 0; i < m; i++)
 {
  A[i*n + k] = rhs(i,0);
 } 
 
}

template <class T>
T fastMatrix<T>::norm2()
{
 T summand = T(0.0);
 
 if (m == 1)
 {
  for (int i = 0; i < n; i++) summand += A[i]*A[i];
  return sqrt(summand);
 }
 else if (n == 1)
 {
  for (int i = 0; i < m; i++) summand += A[i]*A[i];
  return sqrt(summand);
 }
 else
 {
  std::cerr << "norm2() is meaningless for a non-vector\n";
  return T(0.0);
 }
}

template <class T>
void fastMatrix<T>::GramSchmidt()
{
 bool naive = true;
 
 if (naive)
 {
 //first a naive and nonoptimal way
 //destructively performs gram-schmidt in place
  fastMatrix<T> vi (1,n);  
  fastMatrix<T> vj (1,n);
  
 for (int j = 0; j < n; j++) //column
 {
  vj = getCol(j);
  
  for (int i = 0; i <= j; i++) //also column
  {
   //subtract off part that is not orthogonal
   vi = getCol(i);
   vj -= vi*((vj.Transpose()*vi)(0,0));
   setCol(j, vj);
  }
  
  //normalize
  vj /= vj.norm2();
  setCol(j, vj);
 }     
 }
 
}

template <class T>
void fastMatrix<T>::Orthogonalize(int mode = GRAM_SCHMIDT)
{
 switch (mode)
 {
 case HOUSEHOLDER:
      Householder();
 break;
 
 default:
 case GRAM_SCHMIDT:
      GramSchmidt();
      break;
 }
}


template <class T>
bool fastMatrix<T>::isOrthogonal()
{
 //this part takes advantage of the fact that for an orthogonal matrix A, 
  // A*A.transpose() = I
  bool property1 = true;
 
 if (property1)
 {
  fastMatrix AA = this->Transpose();
  if (((*this)*AA).isIdent() ) return true;
  else return false;
 } 
     
}

template <class T>
fastMatrix<T> & fastMatrix<T>::recursiveFFT()
{
 //rfft is in roots_of_unity.h and is the meat and bones of this operation
 //return rfft(*this);              
}


template <class T>
fastMatrix<T> & fastMatrix<T>::iterativeFFT()
{
 //not implemented. see recursiveFFT   
              
}


template <class T>
fastMatrix<T> & fastMatrix<T>::naiveDFT()
{
 int n = getC();
 fastMatrix<T> fourier (n,n,fmFOURIER);
 
 
 return (*this)*fourier;
}

#endif
