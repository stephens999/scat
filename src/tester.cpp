#include <map>
#include "tester.hpp"

using namespace std;

int main() {
  map<string,string> mymap;
  mymap["foo"] = string("bar");
  mymap["baz"] = string("yoiks");
  exit(0);
}
