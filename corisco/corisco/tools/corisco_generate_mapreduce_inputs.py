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
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--output_file', type=str)
    parser.add_argument('--frames', type=int)
    parser.add_argument('--set_name', type=str)
    parser.add_argument('--subset_offset', type=int, default=0)
    parser.add_argument('--ransac_iterations', type=int, nargs='+')
    parser.add_argument('--grid_size', type=int, nargs='+')


    args = parser.parse_args()

    with open(args.input_file) as finput:
        job_params = json.load(finput)

    intrinsic_parameters = job_params['intrinsic_parameters']

    output = [
        {'img_filename': job_params['filename_format'] % fn,
         'ransac_iterations': args.ransac_iterations,
         'grid_size': args.grid_size,
         'img_id': {'set_name': args.set_name, 'index': fn + args.subset_offset},
         'proc_args':
             {'grid': {'size': -1,
                       'lim': 20.0,
                       'method': 'shigeru'
                       },
              'ransac_itr': -1,
              'smooth': 0.0,
              'op_tol': 1e-12,
              'err_func' : [3.0, 1.0, 0.15],
              'int_parm': intrinsic_parameters,
              }
         }
        for fn in range(args.frames)
        ]

    with open(args.output_file, 'a') as fp:
        for kk in output:
            fp.write(json.dumps(kk))
            fp.write('\n')
