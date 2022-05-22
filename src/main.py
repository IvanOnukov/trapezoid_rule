import matplotlib.pyplot as plt
import numpy as np
import sympy
from numpy import zeros
from sympy import parse_expr


def f(x_in, form="sin(x)"):
    """конвертирование из строгого формата в численный формат"""
    expr = parse_expr(form)
    x = sympy.symbols('x')
    fun = sympy.lambdify(x, expr, 'numpy')
    return float(fun(x_in))

def range_f(x_in, form):
    """вычисление значений функции в заданных точках """
    expr = parse_expr(form)
    x = sympy.symbols('x')
    f = sympy.lambdify(x, expr, 'numpy')
    result = []
    for i in x_in:
        result.append(float(f(i)))
    return result


def trapezoid_method(func, a, b, nseg, chart=False):
    """
        Формула трапеций

        In = h( (f(a)+f(b)) / 2 + sum from i=1 to n-1 f(a + i*h) )

        h =( b - a ) / n

        Входные параметры
        func - строковое предстовление формулы
        a - нижняя граница
        b - верхняя граница
        nseg - число отрезков, на которые разбивается [a;b]
        chart - флаг отображения графика
    """
    dx = (float(b) - float(a)) / float(nseg)
    sum_area = 0.5 * (f(a, func) + f(b, func))
    for i in range(1, nseg):
        sum_area += f(a + i * dx, func)

    answer = sum_area * dx

    print("integral " + form + "  " + str(a) + " < x < " + str(b) + "\t=\t" +
          str(answer)
          )

    if chart:
        graph(a, b, form, 0.1, nseg)

    return answer


def _rectangle_rule(func, a, b, nseg, frac):
    """Обобщённое правило прямоугольников."""
    dx = 1.0 * (b - a) / nseg
    sum = 0.0
    xstart = a + frac * dx  # 0 <= frac <= 1 задаёт долю смещения точки,
    # в которой вычисляется функция,
    # от левого края отрезка dx

    for i in range(nseg):
        sum += f(xstart + i * dx, func)

    return sum * dx


def midpoint_rectangle_rule(func, a, b, nseg):
    """Правило прямоугольников со средней точкой"""
    return _rectangle_rule(func, a, b, nseg, 0.5)


def trapezoid_rule(func, a, b, rtol=0.01, nseg0=100):
    """Правило трапеций
       rtol - желаемая относительная точность вычислений
       nseg0 - начальное число отрезков разбиения"""
    nseg = nseg0
    old_ans = 0.0
    dx = 1.0 * (b - a) / nseg
    ans = 0.5 * (f(a, func) + f(b, func))
    for i in range(1, nseg):
        ans += f(a + i * dx, func)

    ans *= dx
    err_est = max(1, abs(ans))
    while err_est > abs(rtol * ans):
        old_ans = ans
        ans = 0.5 * (ans + midpoint_rectangle_rule(func, a, b, nseg))  # новые точки для уточнения интеграла
        # добавляются ровно в середины предыдущих отрезков
        nseg *= 2
        err_est = abs(ans - old_ans)

    print("integral " + func + "  " + str(a) + " < x < " + str(b) + "\t=\t" +
          str(ans)
          )

    return ans


def graph(a, b, form="sin(x)", delt=0.5, nseg=50):
    x = np.linspace(a - delt, b + delt, nseg)
    y = range_f(x, form)
    plt.plot(x, y, c='r')

    x_0 = a
    x_1 = b

    y_0 = f(x_0, form)
    y_1 = f(x_1, form)

    plt.title(form, fontsize=10)
    plt.fill_between([x_0, x_1], [y_0, y_1])

    plt.xlim([a - delt, b + delt])
    plt.ylim([f(a, form), f(b, form) + delt])

    plt.show()


def PrintTriangular(A, i):
    #print(' ', end=' ')
    for l in range(len(A)):
        print('p={0:<4d}'.format(p + l * q), end=' ')
    print()
    for m in range(len(A)):
        print('s={0:<2d}'.format(m), end=' ')
        for l in range(m + 1 - i):
            print('{0:7.4f}'.format(A[m, l]), end=' ')
        print()
    print()


if __name__ == "__main__":

    path = "testData/form3"

    file = open(path, 'r')
    try:
        form = file.readline().strip()
        a = float(eval(file.readline().strip()))
        b = float(eval(file.readline().strip()))
        segments = int(file.readline().strip())
        exactness = float(file.readline().strip())
    finally:
        file.close()

    # trapezoid_method(form, a, b, segments)
    # trapezoid_rule(form, a, b, exactness, segments)

    r = 3   # коифициент сгушениея
    S = 4   # количество сеток на которых будет проводится вычисления
    p = 2   #порядок точности исходной формулы
    q = 2   #для метода трапеций равен 2

    U = zeros((S, S))
    R = zeros((S, S))

    for s in range(S):
        U[s, 0] = trapezoid_method(form, a, b, (r ** s) * segments)


    for s in range(1, S):
        for l in range(s):
            R[s, l] = (U[s, l] - U[s - 1, l]) / (r ** (p + l * q) - 1)
            U[s, l + 1] = U[s, l] + R[s, l]
    print('Таблица приближённых значений интеграла:')
    PrintTriangular(U, 0)
    print('Таблица оценок ошибок:')
    PrintTriangular(R, 1)
