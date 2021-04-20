#include <sstream>
#include <vector>
#include <functional>
#include <memory>
#include <common/bbsolver.hpp>
#include <common/vec.hpp>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <iomanip>
#pragma once

namespace panther {
    template <typename T> class CoordinateDescentMethod : public BlackBoxSolver<T> {
    public:
        struct Options {
            T delta = 1e-3;
            T mMinStep = 1e-4;
            T mDec = 0.5;
        } mOptions;

        bool check_border(int N, const T* left_border, const T* right_border) {
            //left >= right
            for (int i = 0; i < N; ++i) {
                if (left_border[i] > right_border[i]) {
                    return false;
                }
            }
            return true;
        }

        void add_scalar(int N, T delta, const T* first, T* second) {
            for (int i = 0; i < N; ++i) {
                second[i] = first[i] + delta;
            }
        }
        //покоординатный спуск
        T search(int N, T* start, const T* left_border, const T* right_border, const std::function<T(const T * const)> &f) override {
            //проверка коррректности введенных данных:
            if (!check_border(N, left_border, right_border)) {
                throw std::invalid_argument("Границы интервала введены некорректно\n");
            }
            if (is_in_box(N, start, left_border, right_border) != 0) {
                throw std::invalid_argument("Точка не принадлежит интервалу\n");
            }
            //данные корректны

            //заполняем единичные векторы //надо пофиксить будет(...)
            T unit_vectors[N][N];
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    unit_vectors[i][j] = i == j ? 1 : 0;
                }
            }
            //ищем минимальную дельту
            std::vector<T> deltas_array(N, mOptions.delta);
            auto maxStep = [&deltas_array, N]() {
                T max_delta = 0;
                for (int i = 0; i < N; i++) {
                    max_delta = std::max(max_delta, deltas_array[i]);
                }
                return max_delta;
            };

            while (maxStep() > mOptions.mMinStep) { // уменьшаем погрешность, пока это возможно
                for (int i = 0; i < N; ++i) {    // идём по всем направлениям
                    T along_the_axis[N], against_the_axis[N];//вдоль оси, против оси
                    snowgoose::VecUtils vec;

                    vec.vecSaxpy(N, start, unit_vectors[i], deltas_array[i], along_the_axis);//точка вдоль оси
                    vec.vecSaxpy(N, start, unit_vectors[i], -deltas_array[i], against_the_axis);//точка против оси

                    //если мы сделали шаг вправо и вышли за правую границу
                    if (is_in_box(N, along_the_axis, left_border, right_border) == 1) {
                        add_scalar(N, -deltas_array[i], right_border, along_the_axis);
                    }
                    //если мы сделали шаг влево и вышли за левую границу
                    if (is_in_box(N, against_the_axis, left_border, right_border) == 2) {
                        add_scalar(N, deltas_array[i], left_border, against_the_axis);
                    }

                    //по дефолту скажем, что правильное направление вдоль оси
                    T right_direction[N];
                    vec.vecCopy(N, along_the_axis, right_direction);


                    //сравниваем значения функции вдоль оси и против оси
                    if (f(along_the_axis) > f(against_the_axis)) {
                        vec.vecCopy(N, against_the_axis, right_direction);
                        deltas_array[i] *= -1;
                    }
                    //сравниваем наименьшее значение со значением в нашей точке
                    if (f(right_direction) >= f(start)) {
                        continue;
                    }
                    //идем в выбранном направлении, пока функция убывает
                    while (f(right_direction) < f(start)) {
                        vec.vecCopy(N, right_direction, start);
                        vec.vecSaxpy(N, right_direction, unit_vectors[i], deltas_array[i], right_direction);
                        if (is_in_box(N, right_direction, left_border, right_border) != 0) {
                            break;
                        }
                    }
                    //уменьшаем дельту, если не нашли куда идти
                    deltas_array[i] = i == N - 1 ? deltas_array[i] * mOptions.mDec : deltas_array[i];
                }
            }
            return f(start);
        }


        int is_in_box(int n, const T* vec_for_check, const T* left_border, const T* right_border) {
            for (int i = 0; i < n; i++) {
                if (vec_for_check[i] >= right_border[i]) {
                    return 1;
                }
                if (vec_for_check[i] <= left_border[i]) {
                    return 2;
                }
            }
            return 0;
        }
    };
}
