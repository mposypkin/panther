/*
 * File:   sannestand.hpp
 * Author: Maksim Galynchik, Romanova Karina
 */

#ifndef SANNESTAND_HPP
#define SANNESTAND_HPP

#include <algorithm>
#include <cmath>
#include <common/bbsolver.hpp>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <omp.h>
#include <random>
#include "sanneutilities.hpp"
#include <sstream>


namespace panther {
    /**
     * Standart Simulated Annealing method
     * Discription P.M. Pardalos and H.E. Rom.eijn (eds.), Handbook of Global Optimization, Volume 2, 179-229. © 2002 Kluwer Academic Publishers.
     */
    template<class T> class StandartSimulatedAnnealing : public BlackBoxSolver<T> {
    public:
	      /**
	       * Construct Standart Simulated Annealing class
	       * @param next class that defines the next candidate point
	       * @param accept class that determines whether to go to a next candidate point
	       * @param t class that determines the temperature at the current iteration
	       * @param stop class that determines whether to stop the algorithm at the current iteration
           * @param numberI number of points to remember
	       */
        StandartSimulatedAnnealing(NextCandidateDistribution<T>& next, AcceptanceFunction<T>& accept,
            CoolingSchedule<T>& t, StopingCriterion<T>& stop, int numberI = 20) : mD(next), mA(accept), mTemp(t), mStop(stop), numberStopP(numberI) {}

        /**
	       * Performs search
	       * @param n space dimension
	       * @param currentPoint start point and result
	       * @param lowerBound lower bound of domain
         * @param upperBound upper bound of domain
	       * @param f target function
	       * @return value of the result point
         */
        T search(int n, T* currentPoint, const T* lowerBound, const T* upperBound, const std::function<T(const T*)>& f) override {
            T* newPoint = new T[n];
            std::copy(currentPoint, currentPoint + n, newPoint);
            T fOldPoint = f(currentPoint);
            std::deque<T> informationLastPoints;
            informationLastPoints.push_back(fOldPoint);
            unsigned int i = 0;
            do {
                mD.nextCandidate(n, newPoint, lowerBound, upperBound);
                if (mA.acceptance(f(currentPoint), f(newPoint), mTemp.coolingSchedule(i))) {
                    fOldPoint = f(currentPoint);
                    std::copy(newPoint, newPoint + n, currentPoint);
                    if (informationLastPoints.size() == numberStopP)
                        informationLastPoints.pop_front();
                    informationLastPoints.push_back(f(currentPoint));
                }
                i++;
            } while (mStop.stoping(i - 1, informationLastPoints.size(), informationLastPoints));
            delete[] newPoint;
            return f(currentPoint);
        }


        std::string about() {
            std::ostringstream options;
            options << "Simulated Annealing\n";
            options << mD.about();
            options << mA.about();
            options << mTemp.about();
            options << mStop.about();
            return options.str();
        }

    private:
        NextCandidateDistribution<T>& mD;
        AcceptanceFunction<T>& mA;
        CoolingSchedule<T>& mTemp;
        StopingCriterion<T>& mStop;
        int numberStopP;
    };
}

#endif
