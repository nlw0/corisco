from attest import TestBase, test, Assert

import numpy as np

from filtersqp import filterSQP

## "Flat" plane
def val_ff(x):
    return x[0] + 2.0 * x[1]
def grad_ff(x):
    return np.array([1.0, 2.0])
def hess_ff(x):
    return np.array([[0.0, 0.0], [0.0, 0.0]])

## Parabolic sheet
def val_fp(x):
    return x[0] ** 2 + 2 * x[0] * x[1] + x[1] ** 2
def grad_fp(x):
    return np.array([2 * x[0] + 2 * x[1], 2 * x[0] + 2 * x[1]])
def hess_fp(x):
    return np.array([[2.0, 2.0], [2.0, 2.0]])

## Conic constraint
def val_c(x):
    return np.linalg.norm(x)
def grad_c(x):
    return x/np.linalg.norm(x)
def hess_c(x):
    nx = np.linalg.norm(x)
    return (nx ** 2 * np.identity(x.shape[0]) - np.outer(x, x)) / nx ** 3

class Optimization(TestBase):

    @test
    def circle_parabolic_sheet(self):
        x0 = np.array([1.0, 2.0])
        lam0 = 1.0
        rho0 = 0.1
        funcs = (val_c, grad_c, hess_c, val_fp, grad_fp, hess_fp)
        args_f = ()
        x_opt, val_opt, iterations, lam_opt, rho_opt = filterSQP(x0, lam0, rho0, funcs, args_f)
        
        print 'x_opt', x_opt, 'val_opt', val_opt

        Assert(abs(x_opt[0] / - 0.5 ** 0.5 - 1)) < 1.2e-15
        Assert(abs(x_opt[1] / + 0.5 ** 0.5 - 1)) < 1.2e-15

    @test
    def circle_plane(self):
        x0 = np.array([1.0, 2.0])
        lam0 = 1.0
        rho0 = 0.1
        funcs = (val_c, grad_c, hess_c, val_ff, grad_ff, hess_ff)
        args_f = ()
        x_opt, val_opt, iterations, lam_opt, rho_opt = filterSQP(x0, lam0, rho0, funcs, args_f)

        print 'x_opt', x_opt, 'val_opt', val_opt

        Assert(abs(x_opt[0] / - 0.2 ** 0.5 - 1)) < 1.2e-15
        Assert(abs(x_opt[1] / - 0.8 ** 0.5 - 1)) < 1.2e-15

if __name__ == '__main__':
    from attest import Tests
    suite = Tests([Optimization()])
    suite.run()
