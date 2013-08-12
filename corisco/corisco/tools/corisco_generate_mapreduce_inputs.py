#!/usr/bin/python
#coding:utf-8

# Copyright 2011 Nicolau Leal Werneck, Anna Helena Reali Costa and
# Universidade de São Paulo
#
# Licensed under theglim Apache License, Version 2.0 (the "License");
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

import corisco.process

import argparse
import code
import json
from corisco.quaternion import Quat

def main():
    ## Command-line argument parsing
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', type=str)
    parser.add_argument('frames', type=int)

    args = parser.parse_args()

    with open(args.input_file) as finput:
        job_params = json.load(finput)

    intrinsic_parameters = job_params['intrinsic_parameters']

    output = [
        {'image_filename': job_params['filename_format'] % fn,
         'proc_args':
             {'grid': {'size': gs,
                       'lim': 20.0,
                       'method': 'shigeru'
                       },
              'ransac_itr': rit,
              'smooth': 0.0,
              'op_tol': 1e-12,
              'err_func' : [3.0, 1.0, 0.15],
              'int_parm': intrinsic_parameters,
              }
         }
        for rit in [200, 1000, 10000]
        for gs in [128, 64, 32]
        for fn in range(args.frames)
        ]

    for kk in output:
        print json.dumps(kk)
