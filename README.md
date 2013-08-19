corisco
=======

Camera orientation from monocular images based on edgels

This is a rewrite from the project hosted at https://code.google.com/p/corisco/

How to use
==========

We recommend the use of virtualenv to install this experimental
code. Corisco can then be installed as usual for any python package,
with:

  $ python setup.py develop

We notice that the package contains some dependencies such as Cython and PIL.

The corisco command line tool needs to read a json file that defines
the intrinsic paramenters from the camere used to capture an image
set, and also defines the filenames. An example files and a few xample
images are packaged along with corisco. Once corisco is installed you
can test it with these images by running:

  $ cd examples/
  $ corisco set-apasample.json 0 --do_plot

You should see two graphics be plotted, and also an output such as

  {"Nedgels": 16249, "ransac_ori_est": [0.9364874945887336, -0.03621275256917897, 0.34785885021235524, -0.02596207540275075], "sqp_itr": 54, "ori_est": [0.9258680293148183, -0.032798202485475124, 0.3760959820129496, -0.015635936816569097], "time": {"ransac": 1.5499999999999998, "total": 2.89, "edgel": 0.16, "sqp": 1.1500000000000001}, "input": {"img_md5": "cf6b240cb4630e6a2aa36138077dec3b", "proc_args": {"ransac_itr": 1000, "int_parm": {"projection_model": "harris", "principal_point": [500, 375.0], "focal_distance": 847.667796003, "distortion_coefficient": 111.897839472}, "err_func": [3.0, 1.0, 0.15], "op_tol": 1e-12, "smooth": null, "grid": {"lim": 20.0, "method": "shigeru", "size": 4}}}}

