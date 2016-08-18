/*

bench just keeps track of the number of multiplications and additions
 it is dumb and driven by its user class

02-10-2007

*/

#include <ctime>

class bench
{
 public:
  bench();
  bench(char * description);

  void mul();
  void add();

  void mul(unsigned);
  void add(unsigned);
  void clockit();

  unsigned gMul();
  unsigned gAdd();

  void write(char* filename);
  void write();
  double get_time();

 private:
  unsigned muls; //num multiplies
  unsigned adds; //num adds

  double wall; //total time
  clock_t last_time; //stop watch

  char descr[100]; //description
};
