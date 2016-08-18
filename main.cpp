 // quick -n- dirty QR decomposition
 // this is just a demo program for the matrix class

#include <iostream>

#include <cmath>      //for float log(float)

#include <complex> //incase you want complex entries

#include <vector>

#include "fastMatrix.h"

using namespace std;

 typedef std::complex <float> cd;
         
typedef fastMatrix < std::complex<float> > m;

//typedef fastMatrix <double> m;


fastMatrix <cd> & rfft(fastMatrix <cd> & a)
{

 bool debug = false;
 //takes ROW vector, returns ROW vector
 int n = a.getC();


 if (n == 1)
 {
  fastMatrix <cd> * b = new fastMatrix <cd> (a);
  return *b;
 }
 
 static const double _PI = 3.141592653589;
 cd omega_n = cd(cos (2*_PI/n), sin(2*_PI/n));
 
 cd omega = 1.0; 
 
 cd * data = new cd [n];

 for (int i = 0; i < n/2; i++) data[i] = a(0,2*i);
 fastMatrix<cd> a0(1,n/2,data);

 for (int i = 0; i < n/2; i++) data[i] = a(0,2*i+1);
 fastMatrix<cd> a1(1,n/2,data);

 //a0.print();
 //a1.print();

 fastMatrix<cd> y0 = rfft(a0);
 fastMatrix<cd> y1 = rfft(a1);
 
 for (int k = 0; k < n/2; k++)
 {
  data[k] = y0(0,k) + omega*y1(0,k);
  data[k+(n/2)] = y0(0,k) - omega*y1(0,k);
  omega *= omega_n;
 }

  fastMatrix<cd> * y = new fastMatrix <cd> (1,n,data);
  
  //cout << y->getC() << " " << y0.getC() << " " << y1.getC() << endl;
  
  return *y;             
}

int main()
{
 cout.setf(ios::fixed);    

 cout << "DFT test times:\n";

 int mode = 1;

 for (int i = 0; i < 14; i++) //up to a gigabyte
 {
  int j = int(pow(2.0,i));

  if (mode == 1) cout << i << " Vector Size: " << j << " DFT time: ";
  if (mode == 2) cout << " " << log(j) << ", ";  

  m * X = new m(1,j, fmLEXOGRAPHIC);

  bench c("DFT");
  (X->naiveDFT());

  c.clockit();
  if (mode == 1) cout << log(c.get_time()) << endl;
  if (mode == 2) cout << log(c.get_time()) << "; " << endl;
  delete X;
 }

 cout << "FFT test times:\n";

 for (int i = 0; i < 20; i++) //up to a gigabyte
 {
  int j = int(pow(2.0,i));

  if (mode == 1) cout << i << " Vector Size: " << j << " FFT time: ";
  if (mode == 2) cout << " " << log(j) << ", ";  

  m * X = new m(1,j, fmLEXOGRAPHIC);

  bench c("FFT");
  (rfft(*X));

  c.clockit();
  if (mode == 1) cout << log(c.get_time()) << endl;
  if (mode == 2) cout << log(c.get_time()) << "; " << endl;
  delete X;
 }


 return 0;
}
