//
//  main.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 13/03/2021.
//

#include <iostream>
#include "stdafx.hpp"
#include "testing.hpp"
#include "PolynomialBasis.hpp"
#include "Matrix.hpp"
#include "matlib.hpp"

int main(int argc, const char * argv[])
{
    // insert code here...
    std::cout << "Hello there!\n";
    testMatrix();
    testPolynomialBasis();
    testMatlib();
    return 0;
}
