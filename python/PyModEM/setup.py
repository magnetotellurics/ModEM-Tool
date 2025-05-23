#!/usr/bin/env python

from distutils.core import setup

requirements = [
    "numpy"
]

setup(
    author="Miles Curry",
    version='0.1',
    install_requires=requirements,
    long_description_content_type="text/markdown",
    include_package_data=True,
    name="PyModEM",
    packages=["PyModEM"],
    scripts=["scripts/make_mesh",
            "scripts/make_modem_data",
            "scripts/modem_mesh_dump",
            "scripts/plot_mesh",
            "scripts/plot_mesh_2d",
            "scripts/plot_mesh_3d",
            "scripts/plot_mesh_3d_plotly",
            "scripts/read_mesh",
            "scripts/read_modem_data",
            "scripts/sort_report",
            "scripts/synth_error",
            "scripts/csem_covariance"
            ]
)
