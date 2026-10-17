# Fluido de Lennard-Jones en dos dimensiones

Simulacion de dinamica molecular de un fluido bidimensional con interacciones
de Lennard-Jones y condiciones periodicas. El programa permite elegir entre
Euler, Verlet y Velocity-Verlet, y calcula el Hamiltoniano, la funcion de
distribucion radial `g(r)`, el desplazamiento cuadratico medio (MSD) y el
coeficiente de difusion.

## Marco Teórico

Para `N` particulas de masa `m`, posiciones `r_i` y velocidades `v_i`, el
Hamiltoniano es

$$H = \sum_i \frac{1}{2} m |v_i|^2 + \sum_{i < j} U(r_{ij})$$

La interaccion entre cada par de particulas se describe mediante el potencial
de Lennard-Jones:

$$
U(r) = 4 \epsilon \Bigg[\Bigg(\frac{\sigma}{r}\Bigg)^{12} - \Bigg(\frac{\sigma}{r}\Bigg)^6 \Bigg]
$$

`epsilon` fija la profundidad del pozo de potencial y `sigma` la distancia a
la que el potencial se anula. La fuerza correspondiente es el gradiente
negativo del potencial:

$$F_{ij} = \frac{24 \epsilon}{r_{ij}^2} \Bigg[ 2\Bigg(\frac{\sigma}{r_{ij}}\Bigg)^{12} - \Bigg(\frac{\sigma}{r_{ij}}\Bigg)^{6} \Bigg] r_{ij}$$

En una caja cuadrada de lado `l` con condiciones periodicas, cada componente
del desplazamiento entre particulas se calcula usando la convencion de imagen
minima:

$$\Delta x = x_i - x_j - l \hspace{5pt} round \Bigg(\frac{x_i - x_j}{l} \Bigg)$$

$$\Delta y = y_i - y_j - l \hspace{5pt} round\Bigg(\frac{y_i - y_j}{l} \Bigg)$$

La segunda ley de Newton determina las aceleraciones:

$$\dot{r}_i = v_i$$

$$\dot{v}_i = a_i = \frac{1}{m} \sum_{j \neq i} F_{ij}$$

El integrador Velocity-Verlet actualiza primero las posiciones y luego las
velocidades:

$$r_{n+1} = r_n + v_n dt + \frac{1}{2} a_n dt^2$$
$$v_{n+1} = v_n + \frac{1}{2} \Big(a_n + a_{n+1}\Big) dt$$

Ademas de las trayectorias y energias, se calculan dos observables:

$$g(r) = \frac{l^2}{N^2} \frac{\langle n(r) \rangle}{2 \pi r \Delta r}$$

$$MSD(t) = \frac{1}{N} \sum_i |r_i(t) - r_i(0)|^2$$

$$MSD(t) ~= 4 D t \hspace{10pt} \text{(difusion en dos dimensiones)}$$

`g(r)` mide la distribucion radial de particulas respecto a un gas ideal.
