#include <iostream>
#include "tensors.h"

int main()
{
    tensor4D<double> a(2,5,7,9);

    a(0,0,0,0)=1;
    a(1,4,6,8)=3;
    
    
    std::cout<< a(1,4,6,8) << " "  << a.n1 <<std::endl;
    tensor4D<double> b;
    b=a;
    std::cout<< b(1,4,6,8) << " "  << b.n1 <<std::endl;    
    std::cout<< " kk" <<  sizeof(double) << " " << sizeof(long double) << std::endl;

    a(4,6,8,1)=5;
} 
