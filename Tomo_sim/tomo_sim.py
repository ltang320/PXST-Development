"""Three-dimensional Cone-beam Tomography simulation.
The struture study from Nikolay.st_sim.py

Generate intensity frames based on Fresnel diffraction
theory. :class:`TomoSim` does the heavy lifting of calculating
the wavefront propagation to the detector plane.
:class:`STConverter` exports simulated data to a `CXI`_ format
file accordingly to the provided :class:`Protocol` object and
saves the protocol and experimental parameters to the same folder.

Logs of the simulation process are saved in `logs` folder.

.. _CXI: https://www.cxidb.org/cxi.html
"""
from __future__ import division
from __future__ import absolute_import

import os
import logging
import datetime
import argparse

import numpy as np
import scipy as sp
import h5py

from ..protocol import cxi_protocol, ROOT_PATH
from ..data_processing import TomoData
from .tomo_sim_param import TomoParams, parameters
from ..bin import aperture_wp, barcode_steps, barcode_profile, lens_wp
from ..bin import make_frames, make_whitefield
from ..bin import fraunhofer_1d, fraunhofer_1d_scan

class tomoSim:
    """Three-dimensional Cone-beam Tomography simulation class.
    Generates circle's transmission profile and lens' abberated
    wavefront, whereupon propagates the wavefront to the detector's
    plane. Logs all the steps to `logs` folder.

    Parameters
    ----------
    tomo_params : TomoParams
        Experimental parameters.
    bsteps : numpy.ndarray, optional (should change)
        Array of barcode's bar coordinates. Generates the array
        automatically if it's not provided.

    See Also
    --------
    st_sim_param : Full list of experimental parameters.!!
    """
    # generate the log file during calculation
    log_dir = os.path.join(ROOT_PATH, '../logs')

    def __init__(self,tomo_sim_params, bsteps=None):
        self.parameters = tomo_sim_params
        self._init_logging()
        self._init_coord()
        self._init_lens()
        self._init_barcode(bsteps)
        self._init_detector()

    def __getattr__(self,attr):
        """
        Get the value of initial parameters form the file
        """
        if attr in self.parameters:
            return self.parameters.__getattr__(attr)

    def _init_logging(self):
        """
        Initial the log file 
        """
        os.makedirs(self.log_dir, exist_ok=True)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.level = logging.INFO
        filename = os.path.join(
            self.log_dir,
            datetime.datetime.now().strftime('%d-%m-%Y_%H-%M-%S.log'))
        self.logger.addHandler(logging.FileHandler(filename))
        if self.verbose:
            self.logger.addHandler(logging.StreamHandler(stdout))
        for handler in self.logger.handlers:
            handler.setFormatter(
                logging.Formatter(
                    '%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        self.logger.info('Initializing')
        self.logger.info('Current parameters')

    def _init_coord(self):
        # Initializing coordinate parameters
        xx_span = self.fs_size * self.pix_size
        yy_span = self.ss_size * self.pix_size

        x_span = 1.6 * self.ap_x / self.focus * self.defocus
        y_span = 1.2 * self.ap_y
        n_x = int(x_span * xx_span / self.wl / self.det_dist)
        n_y = int(y_span * yy_span / self.wl / self.det_dist)

        # Initializing coordinate arrays
        self.logger.info(
            "Initializing coordinate arrays at the sample's plane")
        self.logger.info('Number of points in x axis: {:d}'.format(n_x))
        self.logger.info('Number of points in y axis: {:d}'.format(n_y))

        self.x_arr = np.linspace(-x_span / 2, x_span / 2, n_x)
        self.y_arr = np.linspace(-y_span / 2, y_span / 2, n_y)
        # Space disperse
        self.xx_arr = np.linspace(-xx_span / 2,
                                  xx_span / 2,
                                  self.fs_size,
                                  endpoint=False)
        self.yy_arr = np.linspace(-yy_span / 2,
                                  yy_span / 2,
                                  self.ss_size,
                                  endpoint=False)


def _init_lens(self):
    #Initializing wavefields at the sample's plane
    self.logger.info("Generating wavefields at the sample's plane")
    self.wf0_x = lens_wp(x_arr=self.x_arr,
                         wl=self.wl,
                         ap=self.ap_x,
                         focus=self.focus,
                         defoc=self.defocus,
                         alpha=self.alpha,
                         xc=(self.x0 - 0.5) * self.ap_x)
    self.wf0_y = aperture_wp(x_arr=self.y_arr,
                             z=self.focus + self.defocus,
                             wl=self.wl,
                             ap=self.ap_y)
    self.i0 = self.p0 / self.ap_x / self.ap_y
    self.smp_c = 1 / self.wl / (self.focus + self.defocus)
    self.logger.info("The wavefields have been generated")
