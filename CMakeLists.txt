cmake_minimum_required(VERSION 3.10)
project(NumericalIntegration)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "C:/MinGW/mingw64/bin/g++.exe")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0") # -g: включить отладочную информацию, -O0: отключить оптимизации
set(CMAKE_CXX_FLAGS_RELEASE "-O3")  # -O3: максимальная оптимизация для Release режима

# Установить режим сборки по умолчанию (Debug)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
# Исходный файл
set(SOURCES
    main.cpp
)
# Заголовочные файлы
set(HEADERS
    numerical_integration.hpp
)
# Подключение `.tpp` файла в зависимостях
set(TEMPLATES
    numerical_integration.tpp
)

# Создание исполняемого файла
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} ${TEMPLATES})
# Добавляем директории заголовков
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})