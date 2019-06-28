/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   advancedcoordescent.hpp
 * Author: posypkin
 *
 * Created on June 28, 2019, 4:59 PM
 */

#ifndef ADVANCEDCOORDESCENT_HPP
#define ADVANCEDCOORDESCENT_HPP

#include <sstream>
#include <vector>
#include <functional>
#include <memory>
#include <common/bbsolver.hpp>
#include <common/vec.hpp>

namespace panther {

    /**
     * An adaptive advanced coordinate descent solver
     */
    template <class T> class AdvancedCoorDescent : public BlackBoxSolver <T> {
    public:

        T search(int n, T* x, const T * const a, const T * const b, const std::function<T(const T * const)> &f) override {
            return 0;
        }

    private:
        

    };
}



#endif /* ADVANCEDCOORDESCENT_HPP */

