#include "gtest/gtest.h"
#include "CoordinateDescentMethod/CoordinateDescentMethod.h"
#include "brute/bruteforce.hpp"
#include <vector>
#include "advcoordesc/advancedcoordescent.hpp"
#include <cmath>

const int N = 2;
double x[N];    //начальная точка поиска минимума
double k[N];    //вектор коэфф
double left[N]; //левая граница поиска
double right[N];//правая граница
double b = 0;   //смещение

void print(double* x, double f) {
    std::cout << "\nFOUND\t" << std::fixed << std::setprecision(7) << f << "\n";
    std::cout << "At\t";
    std::copy(x, x + N , std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n\n";
}

void fill(double* x, double* k, double* left, double* right) {
    std::fill(left, left + N, -2);
    std::fill(right, right + N, -1);
    std::fill(x, x + N, -1.5);
    std::fill(k, k + N, 1);
}

class Function { // родительский класс функций
public:
    int N;
    std::vector<double> a;
    double b;
    Function(int n = 0, double* x = {}, double b1 = 0) : N(n) {
        for (int i = 0; i < n; ++i) {
            a.push_back(x[i]);
        }
        b = b1;
    }
};

class Line: public Function { // класс линейных функций
public:
    using Function::Function;
    double operator() (const double* x) {
        double v;
        for (int i = 0; i < N; i++) {
            v += x[i] * a[i];
        }
        return v + b;
    }
};

class Parabola: public Function { // класс парабол
public:
    using Function::Function;
    double operator() (const double* x) {
        double v;
        for (int i = 0; i < N; i++) {
            v += x[i] * x[i] * a[i] * pow(-1, i); // седло
        }
        return v + b;
    }
};

class Trig: public Function { // класс тригонометрических функций
public:
    using Function::Function;
    double operator() (const double* x) {
        double v;
        for (int i = 0; i < N; i++) {
            v += sin(x[i]) * a[i];
        }
        return v + b;
    }
};

class TestLine : public ::testing::Test {
protected:
    void SetUp()
    {
        fill(x, k, left, right);
        line = new Line(N, k, b);
        res = adv.search(N, x, left, right, std::ref(*line)); // наш минимум
        res_adv_coord_d = adv_coord_d.search(N, x, left, right, std::ref(*line)); // минимум рабочего метода
    }
    void TearDown()
    {
        delete line;
    }
    panther::CoordinateDescentMethod<double> adv;
    panther::AdvancedCoorDescent<double> adv_coord_d;
    double res;
    double res_adv_coord_d;
    Line *line;
};

class TestParabola : public ::testing::Test {
protected:
    void SetUp()
    {
        fill(x, k, left, right);
        parabola = new Parabola(N, k, b);
        res1 = adv1.search(N, x, left, right, std::ref(*parabola));
        res_adv_coord_d1 = adv_coord_d1.search(N, x, left, right, std::ref(*parabola));
    }
    void TearDown()
    {
        delete parabola;
    }
    panther::CoordinateDescentMethod<double> adv1;
    panther::AdvancedCoorDescent<double> adv_coord_d1;
    Parabola *parabola;
    double res1;
    double res_adv_coord_d1;
};

class TestTrig : public ::testing::Test {
protected:
    void SetUp()
    {
        fill(x, k, left, right);
        trig = new Trig(N, k, b);
        res2 = adv2.search(N, x, left, right, std::ref(*trig));
        res_adv_coord_d2 = adv_coord_d2.search(N, x, left, right, std::ref(*trig));
    }
    void TearDown()
    {
        delete trig;
    }
    panther::CoordinateDescentMethod<double> adv2;
    panther::AdvancedCoorDescent<double> adv_coord_d2;
    Trig *trig;
    double res2;
    double res_adv_coord_d2;
};

TEST_F(TestLine, test1) {
double result = abs(res_adv_coord_d - res);
print(x, res);
ASSERT_LE(result, 1e-1); // проверка с заданной точностью
}

TEST_F(TestParabola, test2) {
double result = abs(res_adv_coord_d1 - res1);
print(x, res1);
ASSERT_LE(result, 1e-1);
}

TEST_F(TestTrig, test3) {
double result = abs(res_adv_coord_d2 - res2);
print(x, res2);
ASSERT_LE(result, 1e-1);
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}