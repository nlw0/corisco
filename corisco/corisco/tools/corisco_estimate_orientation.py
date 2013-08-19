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

###############################################################################
## Analyze a single image. Read file name format and intrinsic
## parameters from an input file, and frame number from command
## line. Also reads all possible parameters from command line, then
## assemble a dict with everything the "process" method needs, and
## call it. You can optionally plot the results.

import corisco.process

import argparse
import code
import json
from corisco.quaternion import Quat

def main():
    ## Command-line argument parsing
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', type=str)
    parser.add_argument('frame_number', type=int)
    parser.add_argument('--do_plot', default=False, action='store_true')

    parser.add_argument('--edge_detection_method', type=str, default='shigeru')
    parser.add_argument('--grid_size', type=int, default=4)
    parser.add_argument('--gradient_threshold', type=float, default=20.0)
    parser.add_argument('--ransac_iterations', type=int, default=1000)
    parser.add_argument('--optimization_tolerance', type=float, default=1e-12)
    parser.add_argument('--smooth', type=float)
    parser.add_argument('--error_function_parameters', type=float, nargs=3, default=[3.0, 1.0, 0.15])

    parser.add_argument('--directions_plot_spacing', type=int, default=47)
    parser.add_argument('--save_graphics', default=False, action='store_true')

    args = parser.parse_args()

    with open(args.input_file) as finput:
        job_params = json.load(finput)

    image_filename = job_params['filename_format'] % args.frame_number
    image_file = get_image(image_filename)
    intrinsic_parameters = job_params['intrinsic_parameters']

    process_args = {'grid': {'size': args.grid_size,
                             'lim': args.gradient_threshold,
                             'method': args.edge_detection_method
                             },
                    'ransac_itr': args.ransac_iterations,
                    'smooth': args.smooth,
                    'op_tol': args.optimization_tolerance,
                    'err_func' : args.error_function_parameters,
                    'int_parm': intrinsic_parameters,
                    }

    output_data, pic = corisco.process.estimate_orientation(
        process_args,
        image_file)

    print json.dumps(output_data)

    if args.do_plot or args.save_graphics:
        import matplotlib
        if args.save_graphics:
            matplotlib.use('Agg') 

        import pylab as plt

        if not args.save_graphics:
            plt.ion()
        plt.figure(1,figsize=(6,4.5))
        plt.title('Extracted edgels')
        plt.imshow(.25+pic.frame*.75/255.)
        pic.plot_edgels(plt.gca(), output_data['input']['proc_args']['grid']['size'] * 0.35)
        plt.axis([0, pic.Iwidth, pic.Iheight, 0])
        if args.save_graphics:
            plt.savefig('edg.pdf',dpi=300)

        plt.figure(11,figsize=(6,4.5))
        plt.title(u'Estimated orientation')
        plt.imshow(.2 + pic.frame * .8 / 255.)
        pic.plot_vdirs(plt.gca(), args.directions_plot_spacing, output_data['ori_est'])
        plt.axis([0, pic.Iwidth, pic.Iheight, 0])
        
        if args.save_graphics:
            plt.savefig('ori.pdf',dpi=300)

        if not args.save_graphics:
            code.interact()

def get_image(filename):
    if filename[:5] == 's3://':
        import boto

        conn = boto.connect_s3()

        split_path = filename.split('/')
        bucket_name = split_path[2]
        remaining_path = '/'.join(split_path[3:])

        k = boto.s3.key.Key(conn.get_bucket(bucket_name))
        k.key = remaining_path
        return k.get_contents_as_string()

    else:
        return open(filename).read()
