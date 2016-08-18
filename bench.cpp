#include "bench.h"
#include <iostream>
#include <fstream>

#include <cstring>

#include <ctime>

using namespace std;

bench::bench(char * d)
{
 strcpy(descr, d);
 muls = adds = 0;
 wall = 0.0;
 last_time = clock();
}

bench::bench(){ muls = adds = 0; descr[0] = '\0';}

void bench::mul(){ muls++; }

void bench::add(){ adds++; }

void bench::add(unsigned a) { adds += a;}
void bench::mul(unsigned a) { muls += a;}

void bench::clockit()
{ 
 wall += (static_cast<double>(clock() - last_time))/CLOCKS_PER_SEC;
 last_time = clock();
}

unsigned bench::gMul() { return muls;}

unsigned bench::gAdd() { return adds;}

void bench::write()
{
 std::cout << descr << "\n";
 std::cout << "DOUBLE PRECISION ADDS: " << adds << "\n";
 std::cout << "DOUBLE PRECISION MULT: " << muls << "\n";
 std::cout << "WALL TIME:             " << wall << " seconds\n";
}

double bench::get_time()
{
 return wall;
}

void bench::write(char* filename)
{
 std::ofstream f;
 f.open(filename);
 f << descr << "\n";
 f <<   "DOUBLE PRECISION ADDITIONS:       \t" << adds;
 f << "\nDOUBLE PRECISION MULTIPLICATIONS: \t" << muls;
 f.close();
}
