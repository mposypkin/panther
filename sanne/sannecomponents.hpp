/*
 * File: sannecomponents.hpp
 * Author: Maksim Galynchik, Romanova Karina
 */

#ifndef SANNECOMPONENTS_HPP
#define SANNECOMPONENTS_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include "sannestand.hpp"
#include "sanneutilities.hpp"
#include <string>
#include <deque>

using namespace std;

template<class T> class RandomCandidate : public panther::NextCandidateDistribution<T> {
public:
    struct Options {
        /**
         * Step of the NextCandidateDistribution proportional to the function domain
         */
        T delta;
    };
    RandomCandidate(T step = 0.025) {
        mOpt.delta = step;
    }
    void nextCandidate(int n, T* point, const T* lowerBound, const T* upperBound) const override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(-1, 1);
        for (unsigned int i = 0; i < n; i++) {
            T num;
            do {
                num = point[i] + dis(gen) * (upperBound[i] - lowerBound[i]) * mOpt.delta;
            } while (num >= upperBound[i] || num <= lowerBound[i]);
            point[i] = num;
        }
    }

    std::string about() override {
        std::ostringstream options;
        options << "NextCandidateDistribution : RandomCandidate\n";
        options << "options:\n";
        options << "Step of the NextCandidateDistribution proportional to the function domain " << mOpt.delta << "\n";
        return options.str();
    }

    Options& getOptions() {
        return mOpt;
    }

private:
    Options mOpt;
};

template<class T> class Metropolis : public panther::AcceptanceFunction<T> {
public:
    bool acceptance(T state1, T state2, T t) const override {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        return (state1 > state2) ? true :
            (exp((state1 - state2) * pow(t, -1)) >= dis(gen));
    }

    std::string about() override {
        std::ostringstream options;
        options << "AcceptanceFunction : Metropolis\n";
        return options.str();
    }
};

template<class T> class QuickCooling : public panther::CoolingSchedule<T> {
public:

    struct Options {
        /*
         * Start temperature
         */
        T temp0;
    };

    QuickCooling(T temperature = 70) {
        mOpt.temp0 = temperature;
    }
    T coolingSchedule(unsigned int iteration) const override {
        return mOpt.temp0 / (exp(iteration));
    }

    std::string about() override {
        std::ostringstream options;
        options << "CoolingSchedule :  QuickCooling\n";
        return options.str();
    }

    Options& getOptions() {
        return mOpt;
    }

private:
    Options mOpt;
};

template<class T> class  BoltzmanCooling : public panther::CoolingSchedule<T> {
public:
    struct Options {
        /*
         * Start temperature
         */
        T temp0;
    };

    BoltzmanCooling(T temperature = 70) {
        mOpt.temp0 = temperature;
    }

    T coolingSchedule(unsigned int iteration) const override {
        return mOpt.temp0 / (T)log(iteration + exp(1));
    }

    std::string about() override {
        std::ostringstream options;
        options << "CoolingSchedule : BoltzmannCooling\n";
        return options.str();
    }

    Options& getOptions() {
        return mOpt;
    }

private:
    Options mOpt;
};

template<class T> class StandartStoping : public panther::StopingCriterion<T> {
public:
    struct Options {
        /**
         * Maximum number of iterations of Simulated Annealing
         */
        unsigned int maxIter;
        /**
         * Minimum difference between the value of the previous point and the new
         */
        T accuracy;
        /**
         * Maximum number of iterations at which the difference between points is less than accuracy
         */
        unsigned int stopingIter;
        /**
         * Is stoping information printable
         */
        bool printable;
    };

    StandartStoping(unsigned int maxI = 3000, T acc = 0.01, unsigned int stopIt = 1, bool printable = false) {
        mOpt.maxIter = maxI;
        mOpt.accuracy = acc;
        mOpt.stopingIter = stopIt;
        mOpt.printable = printable;
    }

    bool stoping(unsigned int iter, int numberStop, std::deque<T>& infLastPoints) const override {
        int nowIter = 0;
        for (int j = 0; j < infLastPoints.size() - 1; j++) {
            if (infLastPoints[j] - infLastPoints[j + 1] <= mOpt.accuracy && infLastPoints[j] - infLastPoints[j + 1] >= 0) {
                nowIter++;
            } else {
                nowIter = 0;
            }
        }
        if (nowIter >= mOpt.stopingIter) {
            if (mOpt.printable) std::cout << "Stoping with low value difference on iteration " << iter << std::endl;
            return 0;
        }
        if (iter == mOpt.maxIter) {
            if (mOpt.printable) std::cout << "Stoping with iteration limit." << std::endl;
            return 0;
        }
        return 1;
    }

    std::string about() override {
        std::ostringstream options;
        options << "StopingCriterion : StandartStoping\n";
        options << "options:\n";
        options << "number of steps = " << mOpt.maxIter << "\n";
        options << "Minimum difference between the value of the previous point and the new " << mOpt.accuracy << "\n";
        options << "Maximum number of iterations at which the difference between points is less than accuracy " << mOpt.stopingIter << "\n";
        return options.str();
    }
  
    Options& getOptions() {
        return mOpt;
    }

private:
    Options mOpt;
};

#endif
