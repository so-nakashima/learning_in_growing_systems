#include <colormap-shaders/include/colormap/colormap.h>
#include <iomanip>
#include <iostream>
#include <boost/format.hpp>

int main()
{
    std::cout << boost::format("%02x") % 0
    << std::endl;
}