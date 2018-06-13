#include <iostream>
#include <sstream>

int main() {
  std::stringstream ss;

  ss << "Garbage" << std::end;

  std::cout << ss.str() << std::endl;
}
