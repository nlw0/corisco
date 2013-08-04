from attest import TestBase, test, Assert

from itertools import product
from numpy import array, ones

from filtersqp.trust_rootfind import find_step_size, calculate_distance

class Rootfind(TestBase):

    @test
    def find_root(self):
        L_alpha = [1.0, 10.0]
        L_lam = [[1.0, 2.0, 3.0],
                 [1.0, 1.0],
                 [0.0],
                 [-1.0, 1.0],
                 [-1.0, 0.0, 1.0],
                 [1.0, -2.0, -3.0],
                 [-3.0, -4.0, -5.0, -6.0]]
        L_rho = [1e-2, 1e-1, 1e0, 1e1]
        L_tol = [1e-3, 1e-5, 1e-10, 1e-14]

        for alpha_r, lam, rho, tol in product(L_alpha, L_lam, L_rho, L_tol):
            alpha = ones(len(lam)) * alpha_r
            lam = array(lam)
            tol = rho * tol
            if min(lam) > 0:
                if calculate_distance(alpha, lam, 0.0) < rho:
                    continue

            print 'alpha', alpha, 'lam', lam, 'rho', rho, 'tol', tol
            nu = find_step_size(alpha, lam, rho, tol)
            f_nu = calculate_distance(alpha, lam, nu)
            print 'nu', nu, 'f(nu)', f_nu

            abs(rho - calculate_distance(alpha, lam, nu)) < Assert(tol)
            #abs(rho - calculate_distance(alpha, lam, nu)) == 0

if __name__ == '__main__':
    from attest import Tests
    suite = Tests([Rootfind()])
    suite.run()
