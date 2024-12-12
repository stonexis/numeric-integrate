#pragma once
#include <string>
#include <utility>

namespace Task_const {
    /// Редактируемые параметры
    inline constexpr long double A = -5.5312; ///Концы отрезка
    inline constexpr long double B = 3.32; ///Концы отрезка
    inline constexpr std::size_t K = 39; ///Количество узлов сетки

    /// Нередактируемые параметры 
    inline const long double H = std::abs(B-A)/(K-1); ///Шаг равномерной сетки
    inline const long double STEP_H_2 = Task_const::H / 2; /// Шаг сетки h/2 
}

template <typename T>
const T* gen_grid_func_and_analyt_integrate(
                                    T& analytical_integral,
                                    std::size_t& count_nodes_out,
                                    const T* func_rare=nullptr, 
                                    const std::size_t count_nodes_init=Task_const::K,
                                    const std::size_t ratio=1, 
                                    const T a=Task_const::A, const T b=Task_const::B 
                                    );
template <typename T> 
const T* gen_uniform_grid(
                const T step,
                const std::size_t count_nodes=Task_const::K, 
                const T a=Task_const::A, 
                const T b=Task_const::B
                );
template <typename T>
const T* calculate_numerical_integrals(const T* func, const std::size_t count_nodes, const T step);

template <typename T>
const T* calculate_errors(const T analytical, const T* methods);

#include "numerical_integration.tpp"