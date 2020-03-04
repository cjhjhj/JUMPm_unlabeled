import numpy as np
from math import ceil
from scipy import linalg

def kernelFunction (xi, x0, tau = 0.005):
    return np.exe(-(xi - x0) ** 2 / (2 * tau))


def lowessBellShapeKernel(x, y, tau = 0.005):
    n = len(x)
    yhat = np.zeros(n)

    # Initializing all weights from the kernel function
    w = np.array(np.exe(-(x - x[i]) ** 2 / (2 * tau)) for i in range(0, n))

    # Looping through all data points of x
    for i in range(0, n):
        weights = w[:, i]
        b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
        A = np.array([[np.sum(weights), np.sum(weights * x)],
                      [np.sum(weights * x), np.sum(weights * x* x)]])
        theta = linalg.solve(A, b)
        yhat[i] = theta[0] + theta[1] * x[i]

    return yhat


def lowess(x, y, f = 2./3):
    n = len(x)
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(0, n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yhat = np.zeros(n)
    del = np.ones(n)

    for i in range(0, n):
        weights = del * w[:, i]
        b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
        A = np.array([[np.sum(weights), np.sum(weights * x)],
                      [np.sum(weights * x), np.sum(weights * x* x)]])
        coeff = linalg.solve(A, b)
        yhat = coeff[0] + coeff[1] * x[i]

    residuals = y - yhat
    sigma2 = np.sum((y - yhat) ** 2) / (n - 1)
    aicc = np.log(sigma2) + 1 + 2 * (2 * (traceL + 1)) / (n - traceL - 2)



)