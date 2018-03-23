#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>

/*int  main(int argc, char *argv[])
{
    if(argc != 3) {
        std::cerr << "Wrong arguments passed" << std::endl;
        return -1;
    }

    long long a = std::atoll(argv[1]);
    long long b = std::atoll(argv[2]);

    //bool cmp1 = *argv[1] == argv[1];

    std::cout << argv[1] << std::endl;
    //std::cout << cmp1 << std::endl;
    std::cout << argv[2] << std::endl;
    //std::cout << int(*argv[2]) << std::endl;

    long long c;
    c = a + b;

    std::cout << "c = " << c << " a = " << a << " b = " << b << std::endl;

    std::cout << "EOF" << std::endl;
}*/

int main(int argc, char *argv[])
{
    std::ifstream inFile(argv[1]);
    int sum = 0;
    std::string a;
    char sep = ',';
    //char c;
    //int cmp = int(std::strncmp(*c, ",", 1));
    //inFile.open(argv[1]);
    //while ((inFile >> a >> sep) && (sep == ","))
    while (inFile.good())
    {
        std::getline(inFile, a, sep);
        sum = sum + std::stoi(a);
        std::cout << sum << std::endl;
        //inFile.getline(inFile, b, ",");
        //std::getline(inFile, a, ",");
        //sum = sum + std::stoi(a);  // + std::stoi(b);
        //std::cout << std::stoi(a) << std::endl;
    }

    //while (inFile >> a >> "," >> b)
    //{
    //    sum = sum + int(a) + int(b);
    //}
    inFile.close();

    //std::cout << std::stoi(a) << std::endl;
    //std::cout << b << std::endl;
    std::cout << "Final value: " << sum << std::endl;

    /*float a = std::atof(argv[1]);
    double c = std::atof(argv[1]);
    int b = std::atoi(argv[2]);

    std::cout << "&argv[1]: " << &argv[1] << " &argv[2]: " << &argv[2] << std::endl;
    std::cout << "*argv[1]: " << *argv[1] << " *argv[2]: " << *argv[2] << std::endl;
    std::cout << "argv[1]: " << argv[1] << " argv[2]: " << argv[2] << std::endl;

    float sum = 0;
    double sumd = 0;

    for(int i = 0; i < b; i++)
    {
        sum = sum + a;
        sumd = sumd + c;
    }

    std::cout << std::setprecision (15) << sum << " " << std::setprecision (15) << sumd << std::endl;
    */

    std::cout << "EOF" << std::endl;
}