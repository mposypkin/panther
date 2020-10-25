/*
 * File: testSimulatedAnnealing.cpp
 * Author: Maksim Galynchik
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "sannestand.hpp"
#include "sannecomponents.hpp"


using namespace std;

//input function(for example)
double func(const double* const x) {
    return (x[0] - 10)*(x[0] - 10); // global min 0 at x = 10
}

//input function(for example)
double func2(const double* const x) {
    return 5 * sin(2 * x[0]) + x[0] * x[0] ; // global min -4.4393 at x = -0.71378
}

//input function(for example)
double func3(const double* const x) {
    return (x[0] * x[0] - 1) / (x[0] * x[0] + 1); // global min -1 at x = 0
}

using namespace panther;

int main() {
    const int n = 1;
    double lowerBound[n], upperBound[n], startPoint[n];

    double delta = 0.025;
    RandomCandidate<double> D(delta);

    Metropolis<double> A;

    QuickCooling<double> Temp;

    unsigned int maxIter = 6000;
    double accuracy = 0.025;
    unsigned int stopingIter = 10;
    bool printable = true;
    StandartStoping<double> Stop(maxIter, accuracy, stopingIter, printable);

    StandartSimulatedAnnealing<double> SA(D, A, Temp, Stop);
    std::cout << SA.about();

    std::fill(lowerBound, lowerBound + n, 5);
    std::fill(upperBound, upperBound + n, 15);
    std::fill(startPoint, startPoint + n, 7);
    unsigned int start_time = clock();
    std::cout << "find " << SA.search(n, startPoint, lowerBound, upperBound, func) << std::endl;
    unsigned int end_time = clock();
    std::cout << "runtime = " << (end_time - start_time) / 1000.0 << std::endl;
    std::cout << "at point " << startPoint[0] << std::endl;

    std::fill(lowerBound, lowerBound + n, -3);
    std::fill(upperBound, upperBound + n, 7);
    std::fill(startPoint, startPoint + n, 2);
    start_time = clock();
    std::cout << "find " << SA.search(n, startPoint, lowerBound, upperBound, func2) << std::endl;
    end_time = clock();
    std::cout << "runtime = " << (end_time - start_time) / 1000.0 << std::endl;
    std::cout << "at point " << startPoint[0] << std::endl;

    std::fill(lowerBound, lowerBound + n, -4);
    std::fill(upperBound, upperBound + n, 4);
    std::fill(startPoint, startPoint + n, 2);
    start_time = clock();
    std::cout << "find " << SA.search(n, startPoint, lowerBound, upperBound, func3) << std::endl;
    end_time = clock();
    std::cout << "runtime = " << (end_time - start_time) / 1000.0 << std::endl;
    std::cout << "at point " << startPoint[0] << std::endl;
    return 0;

}