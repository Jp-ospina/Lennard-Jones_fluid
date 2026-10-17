#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <string>
#include <numeric>
#include "../include/fluido.h"


class fluido { // U(r) = 4ε[(σ/r)^^12 − (σ/r)^6]        Walter Greiner - Statistical Mechanics Pág. 401
    private:

        struct Particula {
            double x;
            double y;
            double vx;
            double vy;
        };

        struct Vector2 {
            double x;
            double y;
        };

        double epsilon; // Profundidad del pozo de potencial
        double sigma;   // Distancia a la que el potencial es cero
        double rad = 0.005; // Radio de la particula
        double m = 1.0; // Masa de la particula
        std::vector<Particula> posiciones;
        std::vector<std::vector<Particula>> historial;
        std::vector<std::vector<Vector2>> historialPosicionesDesenvueltas;

        double dt, l, tf, t = 0.0, ax = 0.0, ay = 0.0, g = 0.0;
        const double PI = 3.14159265358979323846;
        int N;
        unsigned int seed = 654321;
        std::mt19937 generador{seed};
        int metodoSeleccionado = 0;
        double anchoBinDistribucion = 0.0;
        double inicioAjusteDifusion = 0.0;
        const std::string directorioResultados = "data/";

        unsigned long long comprobacionesUbicacion = 0;
        mutable unsigned long long evaluacionesParesFuerza = 0;
        mutable unsigned long long evaluacionesComponentesFuerza = 0;
        mutable unsigned long long evaluacionesParesHamiltoniano = 0;
        unsigned long long evaluacionesParesDistribucion = 0;

        double diferenciaMinima(double diferencia) const {
            return diferencia - l * std::round(diferencia / l);
        }

        double envolver(double coordenada) const {
            const double resultado = std::fmod(coordenada, l);
            return resultado < 0.0 ? resultado + l : resultado;
        }
        
        
    public:

        void pedirDatos(){
            std::cout << "Ingrese el valor de epsilon: "; std::cin >> epsilon;
            std::cout << "Ingrese el valor de sigma: "; std::cin >> sigma;
            std::cout << "Ingrese el paso de tiempo dt: "; std::cin >> dt;
            std::cout << "Ingrese el tiempo final tf: "; std::cin >> tf;
            std::cout << "Ingrese el numero de particulas: "; std::cin >> N;
            std::cout << "Ingrese el tamanho del sistema l: "; std::cin >> l;
            std::cout << std::endl;
            std::cout << std::endl;
        };

        void verificar(){
            while (epsilon <= 0 || sigma <= 0 || dt <= 0 || tf <= 0 || N <= 0 || l <= 2.0 * rad) {
                std::cerr << "Error: epsilon, sigma, dt, tf y N deben ser positivos; l debe ser mayor que el diámetro de una partícula." << std::endl;
                pedirDatos();
            }
        };

        fluido(){
            pedirDatos();
            verificar();
            if (!ubicar()) {
                throw std::runtime_error("No se pudieron ubicar todas las partículas.");
            }
        };

        void mostrarDatos(){
            std::cout << "Parametros de la simulacion:\n"
                      << "  epsilon: " << epsilon << '\n'
                      << "  sigma: " << sigma << '\n'
                      << "  dt: " << dt << '\n'
                      << "  tf: " << tf << '\n'
                      << "  pasos: "
                      << static_cast<unsigned long long>(std::ceil(tf / dt)) << '\n'
                      << "  N: " << N << '\n'
                      << "  caja cuadrada: " << l << " x " << l << '\n'
                      << "  masa: " << m << '\n'
                      << "  radio de particula: " << rad << '\n'
                      << "  velocidades iniciales: uniforme en [-1, 1]\n"
                      << "  separacion minima inicial: "
                      << std::max(2.0 * rad, 0.8 * sigma) << '\n'
                      << "  semilla: " << seed << std::endl;
        };

        Particula generarParticula() {
            std::uniform_real_distribution<double> distribucionPosicion(rad, l - rad);
            std::uniform_real_distribution<double> distribucionVelocidad(-1.0, 1.0);

            return {
                distribucionPosicion(generador),
                distribucionPosicion(generador),
                distribucionVelocidad(generador),
                distribucionVelocidad(generador)
            };
        };

        bool fueraDeCaja(const Particula& part) const {
            return part.x < rad || part.x > l - rad || part.y < rad || part.y > l - rad;
        };

        bool haySolapamiento(const Particula& part1, const Particula& part2) {
            ++comprobacionesUbicacion;
            const double dx = diferenciaMinima(part1.x - part2.x);
            const double dy = diferenciaMinima(part1.y - part2.y);
            const double distanciaMinima = std::hypot(dx, dy);
            // Evita iniciar con pares demasiado próximos para el potencial LJ.
            const double separacionMinima = std::max(2.0 * rad, 0.8 * sigma);
            return distanciaMinima < separacionMinima;
        };

        bool ubicar() {
            posiciones.clear();
            posiciones.reserve(N);
            const unsigned long long maxIntentos =
                static_cast<unsigned long long>(N) * 10000ULL;
            unsigned long long intentos = 0;

            for (int i = 0; i < N; ) {
                if (intentos >= maxIntentos) {
                    std::cerr << "Error: no se pudieron ubicar todas las partículas sin solapamiento. "
                                 "Revise N y el tamaño del sistema l."
                              << std::endl;
                    return false;
                }
                ++intentos;

                Particula part = generarParticula();

                bool invalida = false;

                // Verifica que no salga de la caja
                if (fueraDeCaja(part)) {
                    invalida = true;
                }

                // Verifica solapamiento con otras partículas ya ubicadas
                for (const Particula& part2 : posiciones) {
                    if (haySolapamiento(part, part2)) {
                        invalida = true;
                        break;
                    }
                }

                // Si es inválida, vuelve a generar esa misma partícula
                if (invalida) {
                    continue;
                }

                // Si pasa todas las pruebas, la guardamos
                posiciones.push_back(part);
                ++i;
            }

            const std::string rutaArchivo = directorioResultados + "particulas.csv";
            std::ofstream archivo(rutaArchivo.c_str());
            if (!archivo) {
                std::cerr << "Error: no se pudo crear " << rutaArchivo << "." << std::endl;
                return false;
            }

            archivo << "x,y,vx,vy\n";
            archivo << std::setprecision(17);

            for (const Particula& part : posiciones) {
                archivo << part.x << ',' << part.y << ','
                        << part.vx << ',' << part.vy << '\n';
            }
            archivo.close();
            if (!archivo) {
                std::cerr << "Error: ocurrió un problema al escribir "
                          << rutaArchivo << "." << std::endl;
                return false;
            }
            return true;
        };

        double distancia(const Particula& part1, const Particula& part2) const {
            const double dx = diferenciaMinima(part1.x - part2.x);
            const double dy = diferenciaMinima(part1.y - part2.y);
            return std::hypot(dx, dy);
        };

        double calcularHamiltoniano() const {
            double T = 0.0;
            double U = 0.0;

            for (std::size_t i = 0; i < posiciones.size(); ++i) {
                const Particula& part = posiciones[i];
                T += 0.5 * m * (part.vx * part.vx + part.vy * part.vy);

                for (std::size_t j = i + 1; j < posiciones.size(); ++j) {
                    ++evaluacionesParesHamiltoniano;
                    const double r = distancia(part, posiciones[j]);
                    if (r == 0.0 || !std::isfinite(r)) {
                        throw std::runtime_error(
                            "Distancia inválida al calcular el potencial Lennard-Jones.");
                    }

                    const double sr6 = std::pow(sigma / r, 6);
                    U += 4.0 * epsilon * (sr6 * sr6 - sr6);
                }
            }

            return T + U;
        };

        double calcularFuerza(const Particula& part1, const Particula& part2, std::string xy) const {
            ++evaluacionesComponentesFuerza;
            const double dx = diferenciaMinima(part1.x - part2.x);
            const double dy = diferenciaMinima(part1.y - part2.y);
            const double r2 = dx * dx + dy * dy;
            double fx = 0.0;
            double fy = 0.0;
            if (!(r2 > 0.0) || !std::isfinite(r2)) {
                throw std::runtime_error(
                    "Distancia inválida al calcular la fuerza Lennard-Jones.");
            }

            const double sigma2SobreR2 = sigma * sigma / r2;
            const double sr6 = sigma2SobreR2 * sigma2SobreR2 * sigma2SobreR2;
            const double factor = (24.0 * epsilon / r2)
                                * (2.0 * sr6 * sr6 - sr6);
            fx = factor * dx;
            fy = factor * dy;

            if (xy == "x") {
                return fx;
            } else if (xy == "y") {
                return fy;
            } else {
                throw std::invalid_argument("Parámetro 'xy' debe ser 'x' o 'y'.");
            }
        };

        std::vector<Vector2> calcularAceleraciones() const {
            std::vector<Vector2> aceleraciones(posiciones.size(), {0.0, 0.0});
            for (std::size_t i = 0; i < posiciones.size(); ++i) {
                for (std::size_t j = i + 1; j < posiciones.size(); ++j) {
                    ++evaluacionesParesFuerza;
                    const double fx = calcularFuerza(posiciones[i], posiciones[j], "x") / m;
                    const double fy = calcularFuerza(posiciones[i], posiciones[j], "y") / m;
                    aceleraciones[i].x += fx;
                    aceleraciones[i].y += fy;
                    aceleraciones[j].x -= fx;
                    aceleraciones[j].y -= fy;
                }
            }
            return aceleraciones;
        }

        void funcionDistribucion() {
            double dr;
            std::cout << "Ingrese el ancho de bin para g(r): ";
            if (!(std::cin >> dr) || !(dr > 0.0) || !std::isfinite(dr)) {
                throw std::invalid_argument(
                    "El ancho de bin debe ser un número positivo y finito.");
            }
            anchoBinDistribucion = dr;
            std::cout << "Ancho de bin usado para g(r): "
                      << anchoBinDistribucion << std::endl;

            const double radioMaximo = l / 2.0;
            const double binsReales = std::ceil(radioMaximo / dr);

            const std::size_t cantidadBins = std::max(
                static_cast<std::size_t>(1), static_cast<std::size_t>(binsReales));
            std::vector<double> conteos(cantidadBins, 0.0);

            for (const std::vector<Particula>& configuracion : historial) {
                for (std::size_t i = 0; i < configuracion.size(); ++i) {
                    for (std::size_t j = i + 1; j < configuracion.size(); ++j) {
                        ++evaluacionesParesDistribucion;
                        const double r = distancia(configuracion[i], configuracion[j]);
                        if (r >= radioMaximo) {
                            continue;
                        }
                        const std::size_t bin =
                            static_cast<std::size_t>(r / dr);
                        if (bin < cantidadBins) {
                            conteos[bin] += 2.0;
                        }
                    }
                }
            }

            const double area = l * l;
            const double n = static_cast<double>(N);
            const double muestras = static_cast<double>(historial.size());
            const std::string rutaArchivo =
                directorioResultados + "funcion_distribucion_radial.csv";
            std::ofstream archivo(rutaArchivo.c_str());
            if (!archivo) {
                throw std::runtime_error("No se pudo crear " + rutaArchivo + ".");
            }

            
            archivo << "r,g(r),pares_promedio\n" << std::setprecision(17);

            for (std::size_t bin = 0; bin < cantidadBins; ++bin) {
                const double rMin = static_cast<double>(bin) * dr;
                const double rMax = std::min(rMin + dr, radioMaximo);
                const double deltaR = rMax - rMin;
                const double rMedio = 0.5 * (rMin + rMax);
                const double paresPromedio = conteos[bin] / muestras;
                const double g = (rMedio > 0.0 && deltaR > 0.0)
                    ? area * paresPromedio
                      / (n * n * 2.0 * PI * rMedio * deltaR)
                    : 0.0;
                archivo << rMedio << ',' << g << ',' << paresPromedio << '\n';
            }

            archivo.close();
            if (!archivo) {
                throw std::runtime_error("Error al escribir " + rutaArchivo + ".");
            }

            std::cout << "g(r) calculada sobre " << historial.size()
                      << " configuraciones en " << rutaArchivo << "."
                      << std::endl;
        }

        void calcularCoste() const {
            const unsigned long long totalPares =
                comprobacionesUbicacion
                + evaluacionesParesFuerza
                + evaluacionesParesHamiltoniano
                + evaluacionesParesDistribucion;
            const unsigned long long trabajoContabilizado =
                totalPares + evaluacionesComponentesFuerza;
            const unsigned long long paresPorConfiguracion =
                static_cast<unsigned long long>(N) * (N - 1ULL) / 2ULL;

            std::cout << "\nCosto computacional por evaluaciones de pares\n"
                      << "Pares posibles por configuracion: "
                      << paresPorConfiguracion << '\n'
                      << "Comparaciones al ubicar particulas: "
                      << comprobacionesUbicacion << '\n'
                      << "Visitas a pares para calcular fuerzas: "
                      << evaluacionesParesFuerza << '\n'
                      << "Calculos de componentes de fuerza (x/y): "
                      << evaluacionesComponentesFuerza << '\n'
                      << "Visitas a pares para Hamiltoniano: "
                      << evaluacionesParesHamiltoniano << '\n'
                      << "Visitas a pares para g(r): "
                      << evaluacionesParesDistribucion << '\n'
                      << "Total de visitas/comparaciones de pares: "
                      << totalPares << '\n'
                      << "Total incluyendo componentes de fuerza: "
                      << trabajoContabilizado << '\n'
                      << "Orden dominante: O(S N^2), con S = "
                      << historial.size() << " configuraciones."
                      << std::endl;
        }

        void guardarResumenCSV(double hamiltonianoInicial,
                               double hamiltonianoFinal,
                               double coeficienteDifusion) const {
            const unsigned long long totalPares =
                comprobacionesUbicacion
                + evaluacionesParesFuerza
                + evaluacionesParesHamiltoniano
                + evaluacionesParesDistribucion;
            const unsigned long long trabajoContabilizado =
                totalPares + evaluacionesComponentesFuerza;
            const std::string rutaArchivo =
                directorioResultados + "resumen_simulacion.csv";
            std::ofstream archivo(rutaArchivo.c_str());

            archivo << "metrica,valor\n" << std::setprecision(17)
                    << "hamiltoniano_inicial," << hamiltonianoInicial << '\n'
                    << "hamiltoniano_final," << hamiltonianoFinal << '\n'
                    << "coeficiente_difusion_2d," << coeficienteDifusion << '\n'
                    << "metodo_calculo_difusion,"
                       "pendiente_lineal_MSD_sobre_4_ultimos_50_por_ciento\n"
                    << "inicio_ajuste_difusion," << inicioAjusteDifusion << '\n'
                    << "fin_ajuste_difusion," << tf << '\n'
                    << "epsilon," << epsilon << '\n'
                    << "sigma," << sigma << '\n'
                    << "dt," << dt << '\n'
                    << "tf," << tf << '\n'
                    << "N," << N << '\n'
                    << "l," << l << '\n'
                    << "masa," << m << '\n'
                    << "velocidad_inicial_min,-1\n"
                    << "velocidad_inicial_max,1\n"
                    << "separacion_minima_inicial,"
                    << std::max(2.0 * rad, 0.8 * sigma) << '\n'
                    << "pasos," << (historial.size() - 1) << '\n'
                    << "semilla," << seed << '\n'
                    << "metodo,"
                    << (metodoSeleccionado == 1 ? "Euler"
                        : metodoSeleccionado == 2 ? "Verlet" : "Velocity-Verlet")
                    << '\n'
                    << "ancho_bin_gr," << anchoBinDistribucion << '\n'
                    << "particulas," << N << '\n'
                    << "configuraciones," << historial.size() << '\n'
                    << "comparaciones_ubicacion," << comprobacionesUbicacion << '\n'
                    << "evaluaciones_pares_fuerza," << evaluacionesParesFuerza << '\n'
                    << "evaluaciones_componentes_fuerza,"
                    << evaluacionesComponentesFuerza << '\n'
                    << "evaluaciones_pares_hamiltoniano,"
                    << evaluacionesParesHamiltoniano << '\n'
                    << "evaluaciones_pares_distribucion,"
                    << evaluacionesParesDistribucion << '\n'
                    << "total_visitas_pares," << totalPares << '\n'
                    << "trabajo_incluyendo_componentes,"
                    << trabajoContabilizado << '\n';
            archivo.close();
            if (!archivo) {
                throw std::runtime_error("Error al escribir " + rutaArchivo + ".");
            }
        }

        void mover() {
            int metodo;
            std::cout << "Seleccione el metodo: 1) Euler, 2) Verlet, 3) Velocity-Verlet: ";
            if (!(std::cin >> metodo) || metodo < 1 || metodo > 3) {
                throw std::invalid_argument("Seleccione un método entre 1 y 3.");
            }
            metodoSeleccionado = metodo;
            std::cout << "Metodo de integracion seleccionado: "
                      << (metodo == 1 ? "Euler"
                          : metodo == 2 ? "Verlet" : "Velocity-Verlet")
                      << std::endl;

            const double numeroPasosReal = std::ceil(tf / dt);

            
            const std::size_t numeroPasos = std::max(
                static_cast<std::size_t>(1),
                static_cast<std::size_t>(numeroPasosReal));

            t = 0.0;
            historial.clear();
            historial.push_back(posiciones);
            std::vector<double> xActual(posiciones.size());
            std::vector<double> yActual(posiciones.size());
            std::vector<double> xPrevio(posiciones.size());
            std::vector<double> yPrevio(posiciones.size());
            std::vector<double> xNuevo(posiciones.size());
            std::vector<double> yNuevo(posiciones.size());
            for (std::size_t i = 0; i < posiciones.size(); ++i) {
                xActual[i] = posiciones[i].x;
                yActual[i] = posiciones[i].y;
            }

            std::vector<Vector2> aceleraciones = calcularAceleraciones();
            historialPosicionesDesenvueltas.clear();
            historialPosicionesDesenvueltas.emplace_back(posiciones.size());
            for (std::size_t i = 0; i < posiciones.size(); ++i) {
                historialPosicionesDesenvueltas.back()[i] =
                    {xActual[i], yActual[i]};
            }
            if (metodo == 2) {
                const double hInicial = std::min(dt, tf);
                for (std::size_t i = 0; i < posiciones.size(); ++i) {
                    xPrevio[i] = xActual[i] - posiciones[i].vx * hInicial
                               + 0.5 * aceleraciones[i].x * hInicial * hInicial;
                    yPrevio[i] = yActual[i] - posiciones[i].vy * hInicial
                               + 0.5 * aceleraciones[i].y * hInicial * hInicial;
                }
            }

            const std::string rutaArchivo =
                directorioResultados + "simulacionParticulas.csv";
            std::ofstream archivo2(rutaArchivo.c_str());
            if (!archivo2) {
                throw std::runtime_error("No se pudo crear " + rutaArchivo + ".");
            }

            
            archivo2 << std::setprecision(17) << "t,x,y,vx,vy\n";
            for (const Particula& part : posiciones) {
                archivo2 << 0.0 << ',' << part.x << ',' << part.y << ','
                         << part.vx << ',' << part.vy << '\n';
            }

            double hAnterior = std::min(dt, tf);
            for (std::size_t paso = 0; paso < numeroPasos; ++paso) {
                const double h = std::min(dt, tf - t);

                std::vector<Vector2> aceleracionesNuevas;
                switch (metodo) {
                    case 1:
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            xNuevo[i] = xActual[i] + posiciones[i].vx * h;
                            yNuevo[i] = yActual[i] + posiciones[i].vy * h;
                            posiciones[i].vx += aceleraciones[i].x * h;
                            posiciones[i].vy += aceleraciones[i].y * h;
                        }
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            xActual[i] = xNuevo[i];
                            yActual[i] = yNuevo[i];
                            posiciones[i].x = envolver(xActual[i]);
                            posiciones[i].y = envolver(yActual[i]);
                        }
                        aceleraciones = calcularAceleraciones();
                        break;

                    case 2:
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            xNuevo[i] = xActual[i]
                                + (h / hAnterior) * (xActual[i] - xPrevio[i])
                                + 0.5 * aceleraciones[i].x * h * (h + hAnterior);
                            yNuevo[i] = yActual[i]
                                + (h / hAnterior) * (yActual[i] - yPrevio[i])
                                + 0.5 * aceleraciones[i].y * h * (h + hAnterior);
                        }
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            xPrevio[i] = xActual[i];
                            yPrevio[i] = yActual[i];
                            xActual[i] = xNuevo[i];
                            yActual[i] = yNuevo[i];
                            posiciones[i].x = envolver(xActual[i]);
                            posiciones[i].y = envolver(yActual[i]);
                        }
                        aceleracionesNuevas = calcularAceleraciones();
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            posiciones[i].vx = (xActual[i] - xPrevio[i]) / h
                                + 0.5 * aceleracionesNuevas[i].x * h;
                            posiciones[i].vy = (yActual[i] - yPrevio[i]) / h
                                + 0.5 * aceleracionesNuevas[i].y * h;
                        }
                        aceleraciones.swap(aceleracionesNuevas);
                        hAnterior = h;
                        break;

                    case 3:
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            xNuevo[i] = xActual[i] + posiciones[i].vx * h
                                      + 0.5 * aceleraciones[i].x * h * h;
                            yNuevo[i] = yActual[i] + posiciones[i].vy * h
                                      + 0.5 * aceleraciones[i].y * h * h;
                        }
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            xActual[i] = xNuevo[i];
                            yActual[i] = yNuevo[i];
                            posiciones[i].x = envolver(xActual[i]);
                            posiciones[i].y = envolver(yActual[i]);
                        }
                        aceleracionesNuevas = calcularAceleraciones();
                        for (std::size_t i = 0; i < posiciones.size(); ++i) {
                            posiciones[i].vx += 0.5
                                * (aceleraciones[i].x + aceleracionesNuevas[i].x) * h;
                            posiciones[i].vy += 0.5
                                * (aceleraciones[i].y + aceleracionesNuevas[i].y) * h;
                        }
                        aceleraciones.swap(aceleracionesNuevas);
                        break;
                }

                for (const Particula& part : posiciones) {
                    if (!std::isfinite(part.x) || !std::isfinite(part.y)
                        || !std::isfinite(part.vx) || !std::isfinite(part.vy)) {
                        throw std::runtime_error(
                            "La simulación se volvió numéricamente inestable. Reduzca dt.");
                    }
                }
                t = paso + 1 == numeroPasos ? tf : std::min(tf, t + h);
                historial.push_back(posiciones);
                historialPosicionesDesenvueltas.emplace_back(posiciones.size());
                for (std::size_t i = 0; i < posiciones.size(); ++i) {
                    historialPosicionesDesenvueltas.back()[i] =
                        {xActual[i], yActual[i]};
                }
                for (const Particula& part : posiciones) {
                    archivo2 << t << ',' << part.x << ',' << part.y << ','
                             << part.vx << ',' << part.vy << '\n';
                }
            }
            archivo2.close();
            if (!archivo2) {
                throw std::runtime_error("Error al escribir " + rutaArchivo + ".");
            }
        };

        double calcularDifusion() {
            if (historialPosicionesDesenvueltas.size() < 2
                || historialPosicionesDesenvueltas.size() != historial.size()) {
                throw std::runtime_error(
                    "Se requiere una trayectoria completa para calcular la difusión.");
            }

            const std::string rutaArchivo =
                directorioResultados + "desplazamiento_cuadratico_medio.csv";
            std::ofstream archivo(rutaArchivo.c_str());
            if (!archivo) {
                throw std::runtime_error("No se pudo crear " + rutaArchivo + ".");
            }

            std::vector<double> tiempos;
            std::vector<double> msd;
            tiempos.reserve(historialPosicionesDesenvueltas.size());
            msd.reserve(historialPosicionesDesenvueltas.size());
            archivo << "t,msd,d_inst\n" << std::setprecision(17);

            const auto& inicial = historialPosicionesDesenvueltas.front();
            for (std::size_t frame = 0;
                 frame < historialPosicionesDesenvueltas.size(); ++frame) {
                double sumaDesplazamientos = 0.0;
                for (std::size_t i = 0; i < inicial.size(); ++i) {
                    const double dx =
                        historialPosicionesDesenvueltas[frame][i].x - inicial[i].x;
                    const double dy =
                        historialPosicionesDesenvueltas[frame][i].y - inicial[i].y;
                    sumaDesplazamientos += dx * dx + dy * dy;
                }
                const double tiempo = frame == 0
                    ? 0.0
                    : (frame + 1 == historialPosicionesDesenvueltas.size()
                        ? tf
                        : std::min(tf, static_cast<double>(frame) * dt));
                const double msdFrame =
                    sumaDesplazamientos / static_cast<double>(inicial.size());
                const double difusionInstantanea = tiempo > 0.0
                    ? msdFrame / (4.0 * tiempo)
                    : 0.0;
                tiempos.push_back(tiempo);
                msd.push_back(msdFrame);
                archivo << tiempo << ',' << msdFrame << ','
                        << difusionInstantanea << '\n';
            }

            archivo.close();
            if (!archivo) {
                throw std::runtime_error("Error al escribir " + rutaArchivo + ".");
            }

            const std::size_t inicioAjuste = tiempos.size() / 2;
            inicioAjusteDifusion = tiempos[inicioAjuste];
            const double cantidadAjuste =
                static_cast<double>(tiempos.size() - inicioAjuste);
            const double promedioTiempo =
                std::accumulate(tiempos.begin() + inicioAjuste,
                                tiempos.end(), 0.0) / cantidadAjuste;
            const double promedioMsd =
                std::accumulate(msd.begin() + inicioAjuste,
                                msd.end(), 0.0) / cantidadAjuste;
            double numerador = 0.0;
            double denominador = 0.0;
            for (std::size_t i = inicioAjuste; i < tiempos.size(); ++i) {
                const double deltaTiempo = tiempos[i] - promedioTiempo;
                numerador += deltaTiempo * (msd[i] - promedioMsd);
                denominador += deltaTiempo * deltaTiempo;
            }
            if (!(denominador > 0.0)) {
                throw std::runtime_error(
                    "No se puede ajustar la pendiente del MSD en el tiempo.");
            }
            const double pendienteMsd = numerador / denominador;
            const double coeficienteDifusion = pendienteMsd / 4.0;
            std::cout << "MSD guardado en " << rutaArchivo << '\n'
                      << "Ventana del ajuste lineal: t = "
                      << inicioAjusteDifusion << " a " << tf << '\n'
                      << "Coeficiente de difusión 2D (pendiente MSD/4): "
                      << coeficienteDifusion << std::endl;
            return coeficienteDifusion;
        }


};

int ejecutarFluido() {
    std::cout << "Bienvenido al programa de simulacion de un fluido con potencial de Lennard-Jones. \n" << std::endl;
    
    fluido obj;
    obj.mostrarDatos();
    const double hamiltonianoInicial = obj.calcularHamiltoniano();
    std::cout << std::setprecision(17)
                << "Hamiltoniano inicial: " << hamiltonianoInicial << std::endl;
    obj.mover();
    obj.funcionDistribucion();
    const double hamiltonianoFinal = obj.calcularHamiltoniano();
    std::cout << "Hamiltoniano final: " << hamiltonianoFinal << std::endl;
    const double coeficienteDifusion = obj.calcularDifusion();
    obj.calcularCoste();
    obj.guardarResumenCSV(
        hamiltonianoInicial, hamiltonianoFinal, coeficienteDifusion);
    return 0;
}
