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
H = \sum_i \frac{1}{2} m |v_i|^2 + \sum_{i<j} U(r_{ij})
$$

La interaccion entre cada par de particulas se describe mediante el potencial
de Lennard-Jones:

$$
U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}
- \left(\frac{\sigma}{r}\right)^6\right]
$$

$$
\mathbf{F}_{ij} =
\frac{24\epsilon}{r_{ij}^2}
\left[2\left(\frac{\sigma}{r_{ij}}\right)^{12}
- \left(\frac{\sigma}{r_{ij}}\right)^6\right]\mathbf{r}_{ij}
$$

En una caja cuadrada de lado `l` con condiciones periodicas, cada componente
del desplazamiento entre particulas se calcula usando la convencion de imagen
minima:

$$
\Delta x = x_i - x_j - l\,\operatorname{round}\left(\frac{x_i-x_j}{l}\right),
\qquad
\Delta y = y_i - y_j - l\,\operatorname{round}\left(\frac{y_i-y_j}{l}\right)
$$

$$
\frac{d\mathbf{r}_i}{dt} = \mathbf{v}_i,
\qquad
\frac{d\mathbf{v}_i}{dt} = \mathbf{a}_i =
\frac{1}{m}\sum_{j\ne i}\mathbf{F}_{ij}
$$

El integrador Velocity-Verlet actualiza primero las posiciones y luego las
velocidades:

$$
\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_n\,dt
+ \frac{1}{2}\mathbf{a}_n\,dt^2
$$

$$
\mathbf{v}_{n+1} = \mathbf{v}_n
+ \frac{1}{2}(\mathbf{a}_n+\mathbf{a}_{n+1})\,dt
$$

Ademas de las trayectorias y energias, se calculan dos observables:

$$
g(r) = \frac{l^2\langle n(r)\rangle}{N^2\,2\pi r\,\Delta r},
\qquad
\operatorname{MSD}(t) = \frac{1}{N}\sum_i
|\mathbf{r}_i(t)-\mathbf{r}_i(0)|^2
$$

$$
\operatorname{MSD}(t) \approx 4Dt
\qquad \text{(difusión en dos dimensiones)}
$$

`g(r)` mide la distribucion radial de particulas respecto a un gas ideal; `D`
se estima como la pendiente de la region lineal del MSD dividida entre cuatro.