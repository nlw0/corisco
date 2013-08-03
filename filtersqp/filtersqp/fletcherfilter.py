#!/usr/bin/python
#coding:utf-8

# Copyright 2012 Nicolau Leal Werneck, Anna Helena Reali Costa and
# Universidade de São Paulo
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

## This module contains a class that implements the filter that is
## used in the FilterSQP constrained optimization algorithm.

import numpy as np

class FletcherFilter:
    def __init__(self):
        max_pt = 2000
        self.values = np.zeros((max_pt,2))
        self.q = np.zeros(max_pt)
        self.mu = np.zeros(max_pt)
        self.valid = np.zeros(max_pt) > 0  ## All false at the beginning
        self.beta = 99e-2
        self.alpha1 = 25e-2
        self.alpha2 = 1e-4

    def dominated(self, pt):
        ## Test if point is acceptable and return the opposite
        testh = pt[1] <= (self.beta * self.values[self.valid, 1])
        testf1 = pt[0] <= (self.values[self.valid, 0]
                           - self.alpha1 * self.q[self.valid])
        testf2 = pt[0] <= (self.values[self.valid, 0]
                           - self.alpha2 * self.values[self.valid, 1] * self.mu[self.valid])

        ## The north-west rule.
        nw_rule = True
        if self.valid.any():
            leftmost = np.nonzero(self.valid)[0][np.argmin(self.values[self.valid,0])]
            lmu = 1e3 * self.mu[leftmost]
            nw_rule = ((pt[0] + lmu * pt[1]) <=
                       (self.values[leftmost,0] + lmu * self.values[leftmost, 1]))

        return not nw_rule * (testh + testf1 * testf2).all()

    def add(self, pt, delta_q, lam):
        if self.dominated(pt):
            return

        new_index = self.get_new_index()
        self.values[new_index] = pt
        self.valid[new_index] = True
        self.q[new_index] = delta_q
        self.mu[new_index] = np.clip(1e-6, np.max(np.abs(lam)), 1e6)

        ## Invalidate points that are dominated by the new point. Not
        ## strictly necessary, might perhaps improve performance.
        dom = (pt < self.values)
        dom = (dom[:,0] * dom[:,1]) * self.valid
        self.valid[dom] = False

    def get_new_index(self):
        '''Pick the first invalid, "empty" position on the point memory to save a new point.'''
        emp = np.nonzero(1-self.valid)[0]
        assert emp != []
        return  emp[0]
