#include <memory>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <functional>

namespace Method {
    enum Type : std::size_t { // Не вызывает никаких накладных расходов, поскольку на этапе компиляции преобразуется в числа
        Rectangles,  // 0 Метод прямоугольников
        Trapeze,     // 1 Метод трапеций
        Simpson,     // 2 Метод Симпсона
        NewtonCotes, // 3 Метод Ньютона-Котеса
        Gauss,       // 4 Метод Гаусса
        Count        // 5 Количество методов
    };
}
/**
 * @brief Функция для генерации или измельчения массива значений заданной функции на отрезке и вычисления аналитического интеграла
 * @tparam T Тип данных (float, double, long double).
 * @param[out] analytical_integral Значение аналитического интеграла (Выходной параметр)
 * @param[out] count_nodes_out Количество узлов получившейся сетки (Выходной параметр)
 * @param func_rare Старый массив сетки, на основе которого строится измельчение (По умолчанию nullptr)
 * @param count_nodes_init Количество узлов для инициализации сетки, при существовании более редкой сетки ДОЛЖЕН быть равен количеству узлов в ней (По умолчанию Task_const::K)
 * @param ratio Во сколько раз необходимо измельчить сетку (По умолчанию = 1)
 * @param a Левый конец отрезка (По умолчанию Task_const::A)
 * @param b Правый конец отрезка (По умолчанию Task_const::B)
 * @return Пара массивов pair(func, antiderivative)
 * @note Каждый массив имеет размер count_nodes_out, накладные расходы по памяти - размер пары 2*размер указателя = 16 байт
 */
template <typename T>
const T* gen_grid_func_and_analyt_integrate(
                                    T& analytical_integral,
                                    std::size_t& count_nodes_out,
                                    const T* func_rare, 
                                    const std::size_t count_nodes_init,
                                    const std::size_t ratio, 
                                    const T a, const T b
                                    ){
    if (count_nodes_init < 2) throw std::invalid_argument("Invalid count_nodes_start values");                                    
    if (std::abs(b - a) < std::numeric_limits<T>::epsilon()) throw std::invalid_argument("Invalid a, b values");
    if (a > b) throw std::invalid_argument("Invalid a, b values");
    if (ratio < 1) throw std::invalid_argument("Invalid ratio");
    if (func_rare == nullptr && ratio != 1) throw std::invalid_argument("Incorrect initialization");

    T* arr_func = nullptr; // Массив значений функции на сетке

    static auto func = [](T x) -> T {return std::sin(x); }; // Заданная функция
    static auto ant = [](T x) -> T {return -1.0 * std::cos(x); }; // Аналитическая первообразная

    static auto calc_analyt = [](T a, T b, std::function<T(T)> ant) -> T {
        static T prev_a = a; static T prev_b = b; // Тут значения вычисляются только при первом запуске функции,
        static T integral = ant(b) - ant(a); // в дальнейшем, чтобы их изменить, нужно их переприсвоить
        //Если a и b не изменились, то аналитический интеграл пересчитывать не нужно
        if (abs(a - prev_a) < std::numeric_limits<T>::epsilon() || 
            abs(b - prev_b) < std::numeric_limits<T>::epsilon()){
            integral = ant(b) - ant(a);
            prev_a = a; prev_b = b;
        }
        return integral;
    };
    analytical_integral = calc_analyt(a, b, ant); // Поскольку ant не захватывает переменные, то при передаче не создается ее копия, в функцию передается ее ссылка

    if (func_rare == nullptr){ // Массив другой сетки не существует, это первый запуск функции
        count_nodes_out = count_nodes_init;
        T step = std::abs(b-a) / (count_nodes_out - 1);
        const T* grid_x = gen_uniform_grid(step, count_nodes_out, a, b); // Создаем сетку на оси x
        arr_func = new T[count_nodes_out]{};
    
        for (std::size_t i = 0; i < count_nodes_out; i++)
            arr_func[i] = func(grid_x[i]); //заданная функция

        delete[] grid_x;
    }
    else { // Массив старой сетки существует, нужно измельчить сетку функции
        count_nodes_out = count_nodes_init * ratio - 1; //-1 проверяется на листке бумаги
        T step = std::abs(b-a) / (count_nodes_out - 1);
        const T* grid_x = gen_uniform_grid(step, count_nodes_out, a, b); //Пересоздаем сетку на оси x(По производительности уступает на 15%, однако поскольку в дальнейшем не используется существенная экономия памяти)
        arr_func = new T[count_nodes_out]{};

        for (std::size_t i = 0; i < count_nodes_init - 1; i++){ //Не доходим до последнего элемента, чтобы не выходить за границы массива
            arr_func[ratio * i] = func_rare[i]; //Каждый nй (2-й) элемент новой сетки это элемент старой сетки
            arr_func[ratio * i + 1] = func(grid_x[ratio * i + 1]); // Пересчитываем значение
        }
        arr_func[count_nodes_out - 1] = func_rare[count_nodes_init - 1]; // Последние элементы совпадают
        delete[] grid_x;
    }
    
    return arr_func;
}

/**
 * @brief Функция для генерации равномерной сетки на отрезке [a,b]
 * @tparam T Тип данных (float, double, long double).
 * @param step - Шаг равномерной сетки.
 * @param count_nodes - Количество узлов сетки. (По умолчанию Task_const::K)
 * @param a Начало отрезка. (По умолчанию Task_const::A)
 * @param b Конец отрезка. (По умолчанию Task_const::B)
 * @return T* Указатель на массив равномерной сетки.
 * @note Массив имеет размер count_nodes. (По умолчанию Task_const::K)
 */
template <typename T>
const T* gen_uniform_grid(const T step, const std::size_t count_nodes, const T a, const T b) {
    if (count_nodes < 2) throw std::invalid_argument("Invalid count_nodes values");
    if (abs(b - a) < std::numeric_limits<T>::epsilon()) throw std::invalid_argument("Invalid a, b values");
    if (a > b) throw std::invalid_argument("Invalid a, b values");
    T* array = new T[count_nodes]{}; 
    for (std::size_t i = 0; i < count_nodes; i++) 
        array[i] = a + step * i; // Заполняем значения, включая последний узел, равный b
    if (array[count_nodes - 1] != b)
        array[count_nodes - 1] = b;
    
    return array;
}
/**
 * @brief Функция для численного вычисления интеграла от функции заданной на равномерной сетке на отрезке.
 * @tparam T Тип данных (float, double, long double).
 * @param func Массив сетки функции.
 * @param count_nodes Количество узлов сетки. 
 * @param step Шаг равномерной сетки
 * @return T* Указатель на массив значений используемых методов (Имеет размер Methods::Count=5) 
 */
template <typename T>
const T* calculate_numerical_integrals(const T* func, const std::size_t count_nodes, const T step){
    if (func == nullptr) throw std::invalid_argument("func is null");
    if (count_nodes < 5) throw std::invalid_argument("Invalid count_nodes values");

    T rectangles = T(); T trapeze = T(); T simpson = T(); T newton_cotes = T(); T gauss = T();
    T* methods = new T[Method::Count]{}; //Создаем массив размера count для хранения вычисленных значений

    // Формула средних прямоугольников
    static auto rectangles_calc = [](const T* func, const T step, const std::size_t i) -> T { 
        T f = func[i]; 
        return 2.0 * step * f; // (x_2 - x_0)*f(x_1) = 2h*f(x_1)
    };

    // Формула трапеций
    static auto trapeze_calc = [](const T* func, const T step, const std::size_t i) -> T {
        T f = func[i] + func[i + 1];
        return (step / 2.0) * f; // ((x_1 - x_0) / 2) * (f(x_0) + f(x_1))
    };

    // Формула Симпсона
    static auto simpson_calc = [](const T* func, const T step, const std::size_t i) -> T {
        constexpr T weigth[3] = {1.0, 4.0, 1.0};
        T f = weigth[0] * func[i] + weigth[1] * func[i + 1] + weigth[2] * func[i + 2]; 
        return (step / 3.0)* f; // ((x_2 - x_0) / 6) * (f(x_0) + 4f(x_1) + f(x_2)) 
    };

    // Формула Ньютона-Котеса для 5 узлов 
    static auto newton_cotes_calc = [](const T* func, const T step, const std::size_t i) -> T {
        constexpr T weigth[5] = {7, 32, 12, 32, 7};
        T f = weigth[0] * func[i] + weigth[1] * func[i + 1] + weigth[2] * func[i + 2] + weigth[3] * func[i + 3] + weigth[4] * func[i + 4];
        return (2.0 / 45.0) * step * f; // 2/45 * h * (7f_0 + 32f_1 + 12f_2 + 32f_3 + 7f_4)
    };

    // Формула Гаусса для 3х узлов
    static auto gauss_calc = [](const T* func, const T step, const std::size_t i) -> T {
        //Поскольку по трем узлам(средння точка центральная), то подотрезок будет длинны 2*step, веса метода гаусса опредлены на [-1,1], отображаем на [-h, h]
        const T weigth[3] = {(5.0 / 9.0) * step, (8.0 / 9.0) * step, (5.0 / 9.0) * step};
        return weigth[0] * func[i] + weigth[1] * func[i + 1] + weigth[2] * func[i + 2];
    };
    // Занимаемый обьем памяти при создании каждой лямбда функции - размер указателя на функцию, те 8 байт

    for (std::size_t i = 0; i < count_nodes; i++) { // Не по каждому индексу можно итерироваться в формулах, поскольку могут возникнуть наложения
        // Прямоугольники (центральные узлы, исключаем первый и последний, идем по 3 точки, тоесть по всем нечетным)
        if (i % 2 != 0 && i != count_nodes - 1)
            rectangles += rectangles_calc(func, step, i);
        // Трапеции (каждый узел, исключая последний)
        if (i < count_nodes - 1) 
            trapeze += trapeze_calc(func, step, i);
        // Симпсон (каждые 3 узла)
        if (i % 2 == 0 && i < count_nodes - 2) // Идем по каждому четному, проверяется на бумаге
            simpson += simpson_calc(func, step, i);
        // Ньютон-Котес (каждые 5 узлов)
        if (i % 4 == 0 && i < count_nodes - 4) // Идем по каждому четвертому, проверяется на бумаге
            newton_cotes += newton_cotes_calc(func, step, i);
        // Гаусс (каждые 3 узла)
        if (i % 2 == 0 && i < count_nodes - 2)
            gauss += gauss_calc(func, step, i);
    }
    methods[Method::Rectangles] = rectangles; //Преобразование типа enum class к size_t, для заполнения массива по индексу
    methods[Method::Trapeze] = trapeze;
    methods[Method::Simpson] = simpson;
    methods[Method::NewtonCotes] = newton_cotes;
    methods[Method::Gauss] = gauss;

    return methods;
}
/**  
* @brief Функция для вычисления относительных погрешностей вычисления интегралов 
* @tparam T Тип данных (float, double, long double).
* @param analytical Аналитическое значение интеграла.
* @param methods Массив значений интегралов вычесленных соответсвующими методами 
* @return T* Указатель на массив с соответсвующими относительными ошибками вычисления
* @note Массив имеет размер Methods::Count(=5)  
*/ 
template <typename T>
const T* calculate_errors(const T analytical, const T* methods){
    T* errors = new T[Method::Count]{}; // Массив относительных ошибок соответсвующих методов
    static auto error_calc = [](const T analytical, const T numerical) -> T {
        return (std::abs(analytical - numerical) / analytical);
    };
    for (std::size_t i = Method::Rectangles; i < Method::Count; i++)
        errors[i] = error_calc(analytical, methods[i]);

    return errors;
}
/**
 * @brief Функция вывода значений абсолютной и относительной погрешностей в формате таблицы
 * @tparam T Тип данных (float, double, long double).
 * @param errors_h Массив ошибок на сетке h
 * @param errors_h_2 Массив ошибок на сетке h_2
 */
template <typename T>
void print_error_table(const T* errors_h, const T* errors_h_2){
    if (errors_h == nullptr || errors_h_2 == nullptr) throw std::invalid_argument("Input arrays cannot be null");
    // Функция для вывода строки таблицы
    auto print_row = [](const std::string& label, const T data_1, const T data_2){
            std::cout << std::left // Выравнивание текста влево
                      << std::setw(18) << label // Устанавливает фиксированную ширину вывода
                      << std::scientific << std::setprecision(6) //  Устанавливает точность
                      << std::setw(15) << data_1
                      << std::setw(15) << data_2 << "\n";
        };
    // Заголовок таблицы
    std::cout << std::left
              << std::setw(18) << " " // Пропуск для выравнивания
              << std::setw(15) << "h"
              << std::setw(15) << "h/2" << "\n";

    std::cout << std::string(61, '-') << "\n";
    print_row("Rectangles", errors_h[Method::Rectangles], errors_h_2[Method::Rectangles]);
    std::cout << std::string(61, '-') << "\n";
    print_row("Trapeze", errors_h[Method::Trapeze], errors_h_2[Method::Trapeze]);
    std::cout << std::string(61, '-') << "\n";
    print_row("Simpson", errors_h[Method::Simpson], errors_h_2[Method::Simpson]);
    std::cout << std::string(61, '-') << "\n";
    print_row("NewtonCotes", errors_h[Method::NewtonCotes], errors_h_2[Method::NewtonCotes]);
    std::cout << std::string(61, '-') << "\n";
    print_row("Gauss", errors_h[Method::Gauss], errors_h_2[Method::Gauss]);
    std::cout << std::string(61, '-') << "\n";
}
