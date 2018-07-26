#include <iostream>
#include <sstream>

int main() {
  std::stringstream ss;

  ss << "Garbage" << std::endl;

  std::cout << ss.str() << std::endl;
}
