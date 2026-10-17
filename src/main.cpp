#include "fluido.h"

#include <exception>
#include <iostream>

int main() {
    try {
        return ejecutarFluido();
    } catch (const std::exception& error) {
        std::cerr << "Error durante la simulacion: " << error.what() << '\n';
        return 1;
    }
}
