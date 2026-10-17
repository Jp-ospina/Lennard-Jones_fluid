# Fluido de Lennard-Jones en dos dimensiones

Simulacion de dinamica molecular de un fluido bidimensional con interacciones
de Lennard-Jones y condiciones periodicas. El programa permite elegir entre
Euler, Verlet y Velocity-Verlet, y calcula el Hamiltoniano, la funcion de
distribucion radial `g(r)`, el desplazamiento cuadratico medio (MSD) y el
coeficiente de difusion.

## Marco Teórico

Para `N` particulas de masa `m`, posiciones `r_i` y velocidades `v_i`, el
Hamiltoniano es

$$
H = sum_i (1/2) m |v_i|^2 + sum_{i<j} U(r_ij)
$$

La interaccion entre cada par de particulas se describe mediante el potencial
de Lennard-Jones:

$$
U(r) = 4 epsilon [(sigma/r)^12 - (sigma/r)^6]
$$

`epsilon` fija la profundidad del pozo de potencial y `sigma` la distancia a
la que el potencial se anula. La fuerza correspondiente es el gradiente
negativo del potencial:

$$
F_ij = (24 epsilon/r_ij^2)
       [2(sigma/r_ij)^12 - (sigma/r_ij)^6] r_ij
$$

En una caja cuadrada de lado `l` con condiciones periodicas, cada componente
del desplazamiento entre particulas se calcula usando la convencion de imagen
minima:

$$
Delta x = x_i - x_j - l round((x_i - x_j)/l)
Delta y = y_i - y_j - l round((y_i - y_j)/l)
$$

La segunda ley de Newton determina las aceleraciones:

$$
dr_i/dt = v_i
dv_i/dt = a_i = (1/m) sum_{j != i} F_ij
$$

El integrador Velocity-Verlet actualiza primero las posiciones y luego las
velocidades:

$$
r_(n+1) = r_n + v_n dt + (1/2) a_n dt^2
v_(n+1) = v_n + (1/2) (a_n + a_(n+1)) dt
$$

Ademas de las trayectorias y energias, se calculan dos observables:

$$
g(r) = l^2 <n(r)> / (N^2 2 pi r Delta r)
MSD(t) = (1/N) sum_i |r_i(t) - r_i(0)|^2
MSD(t) ~= 4 D t       (difusion en dos dimensiones)
$$

`g(r)` mide la distribucion radial de particulas respecto a un gas ideal; `D`
se estima como la pendiente de la region lineal del MSD dividida entre cuatro.