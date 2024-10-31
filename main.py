# Начальные условия
x = 0.75 * 10**-4
y = 10**-15
w = 0

# Параметры
r = 1  # Замените на нужное значение r
s = 0.1  # Замените на нужное значение s
alpha = 0.01  # Замените на нужное значение alpha
H = 0.5  # Замените на нужное значение H
dt = 0.1  # Шаг интегрирования
m1 = 1e-14  # Задайте значение m1
t = 0  # Начальное время
integral = 0  # Начальное значение интеграла

# Определение функции f1
def f1(x, y, w):
    return x * y + w

# Определение производных
def dxdt(r, s, x, y, w):
    return -1 * r * x * y 

def dydt(r, alpha, H, x, y):
    return -1 * r * x * y 

def dwdt(r, x, y):
    return r * x * y

# Реализация метода Рунге-Кутты 4-го порядка
def rk4_step(x, y, w, r, s, alpha, H, dt):
    k1x = dxdt(r, s, x, y, w)
    k1y = dydt(r, alpha, H, x, y)
    k1w = dwdt(r, x, y)

    k2x = dxdt(r, s, x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, w + 0.5 * dt * k1w)
    k2y = dydt(r, alpha, H, x + 0.5 * dt * k1x, y + 0.5 * dt * k1y)
    k2w = dwdt(r, x + 0.5 * dt * k1x, y + 0.5 * dt * k1y)

    k3x = dxdt(r, s, x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, w + 0.5 * dt * k2w)
    k3y = dydt(r, alpha, H, x + 0.5 * dt * k2x, y + 0.5 * dt * k2y)
    k3w = dwdt(r, x + 0.5 * dt * k2x, y + 0.5 * dt * k2y)

    k4x = dxdt(r, s, x + dt * k3x, y + dt * k3y, w + dt * k3w)
    k4y = dydt(r, alpha, H, x + dt * k3x, y + dt * k3y)
    k4w = dwdt(r, x + dt * k3x, y + dt * k3y)

    x += (dt / 6) * (k1x + 2 * k2x + 2 * k3x + k4x)
    y += (dt / 6) * (k1y + 2 * k2y + 2 * k3y + k4y)
    w += (dt / 6) * (k1w + 2 * k2w + 2 * k3w + k4w)

    return x, y, w

# Цикл для нахождения t_0
while integral < m1:
    # Обновление значений x, y, w с помощью RK4
    x, y, w = rk4_step(x, y, w, r, s, alpha, H, dt)

    # Вычисление текущего значения функции f1
    current_f1 = f1(x, y, w)

    # Обновление интеграла (используем простое приближение для интеграла)
    integral += current_f1 * dt
    
    # Обновление времени
    t += dt

    # Вывод текущих значений для мониторинга
    print(f"t: {t:.2f}, x: {x:.6e}, y: {y:.6e}, w: {w:.6e}, f1: {current_f1:.6e}, integral: {integral:.6e}")

# t является искомым t_0
print(f"Final t0: {t:.6e}, Final integral value: {integral:.6e}")
