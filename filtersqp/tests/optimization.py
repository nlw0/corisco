from attest import TestBase, test, Assert

import numpy as np

from filtersqp import filterSQP

def val_f(x):
    return x[0]**2 + 2*x[0]*x[1] + x[1]**2
def grad_f(x):
    return np.array([2*x[0] + 2*x[1], 2*x[0] + 2*x[1]])
def hess_f(x):
    return np.array([[2.0,2],[2,2]])

## Conic constraint
def val_c(x):
    return np.linalg.norm(x)
def grad_c(x):
    return x/np.linalg.norm(x)
def hess_c(x):
    nx = np.linalg.norm(x)
    return (nx**2 * np.identity(x.shape[0]) - np.outer(x,x)) / nx**3

class Optimization(TestBase):

    @test
    def circle_plane(self):

        x0 = np.array([1.0, 2.0])
        lam0 = 1.0
        rho0 = 0.1
        funcs = (val_c, grad_c, hess_c, val_f, grad_f, hess_f)
        args_f = ()
        x_opt, val_opt, iterations, lam_opt, rho_opt = filterSQP(x0, lam0, rho0, funcs, args_f)
        
        print x_opt, val_opt

        Assert(x_opt[0]) == -np.sqrt(0.5)
        Assert(x_opt[1]) == np.sqrt(0.5)

if __name__ == '__main__':
    from attest import Tests
    suite = Tests([Optimization()])
    suite.run()
