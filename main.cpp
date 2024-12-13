#include <iostream>
#include "numerical_integration.hpp"
using namespace std;

int main(){
    long double analytic_integral; //Переменная для хранения аналитического интеграла
    std::size_t count_nodes_h; //Количетсво узлов сетки с шагом h

    const auto grid_f_in_h = gen_grid_func_and_analyt_integrate(analytic_integral, count_nodes_h); //Сетка значений функции f с шагом h
    const auto integrals_in_h = calculate_numerical_integrals(grid_f_in_h, count_nodes_h, Task_const::H); // Массив значений интегралов вычисленных соответсвующим методом на сетке h
    const auto errors_in_h = calculate_errors(analytic_integral, integrals_in_h); // Массив значений относительных ошибок на сетке h

    std::size_t count_nodes_h_2; //Количество узлов сетки с шагом h/2
    const auto grid_f_in_h_2 = gen_grid_func_and_analyt_integrate(analytic_integral, count_nodes_h_2, grid_f_in_h, count_nodes_h, 2); //Измельчаем сетку
    const auto integrals_in_h_2 = calculate_numerical_integrals(grid_f_in_h_2, count_nodes_h_2, Task_const::STEP_H_2);
    const auto errors_in_h_2 = calculate_errors(analytic_integral, integrals_in_h_2);

    print_error_table(errors_in_h, errors_in_h_2);
    
    delete[] grid_f_in_h;
    delete[] grid_f_in_h_2;
    delete[] integrals_in_h;
    delete[] integrals_in_h_2;
    delete[] errors_in_h;
    delete[] errors_in_h_2;

    return 0;
}