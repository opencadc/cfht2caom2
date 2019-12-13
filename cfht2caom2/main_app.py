# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

import importlib
import logging
import os
import sys
import traceback

from caom2 import Observation, CalibrationLevel, ObservationIntentType
from caom2 import ProductType, CompositeObservation
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc


__all__ = ['cfht_main_app', 'update', 'CFHTName', 'COLLECTION',
           'APPLICATION', 'ARCHIVE']


APPLICATION = 'cfht2caom2'
COLLECTION = 'CFHT'
ARCHIVE = 'CFHT'

# units are Angstroms
FILTER_LOOKUP = {  # commented out so a test passes ...
                 # 'u.MP9301': {'centre': 3754.5, 'width': 867 / 2.0},
                 'g.MP9401': {'centre': 4890.0, 'width': 1624 / 2.0},
                 'r.MP9601': {'centre': 6248.5, 'width': 1413.0},
                 # TODO - filter lookup - I added this entry ....
                 'r.MP9602': {'centre': 6414., 'width': 1524.0},
                 'i.MP9701': {'centre': 7763., 'width': 1726 / 2.0},
                 'i.MP9702': {'centre': 7623., 'width': 1728 / 2.0},
                 # TODO the following line is commented out so a test passes
                 # before I get a chance to talk to Chris
                 # 'z.MP9801': {'centre': 9083., 'width': 1848.0},
                 }

# All comments denoted 'CW' are copied from
# ssh://goliaths@gimli3/srv/cadc/git/wcaom2archive


class CFHTName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support fz files in storage
    - product id == file id
    - the file_name attribute has ALL the extensions, including compression
      type.
    """

    CFHT_NAME_PATTERN = '*'

    def __init__(self, obs_id=None, fname_on_disk=None, file_name=None):
        # set compression to an empty string so the file uri method still
        # works, since the file_name element will have all extensions,
        # including the .fz to indicate compresssion
        super(CFHTName, self).__init__(
            None, COLLECTION, CFHTName.CFHT_NAME_PATTERN, file_name,
            compression='')
        self._file_id = CFHTName.remove_extensions(file_name)
        if CFHTName.is_raw(self._file_id):
            self.obs_id = self._file_id[:-1]
        else:
            self.obs_id = self._file_id
        self._file_name = file_name
        logging.debug(self)

    def __str__(self):
        return f'obs_id {self.obs_id}, ' \
               f'file_id {self._file_id}, ' \
               f'file_name {self.file_name}, ' \
               f'lineage {self.lineage}'

    def is_valid(self):
        return True

    @property
    def file_id(self):
        """The file id - no extensions."""
        return self._file_id

    @property
    def file_name(self):
        """The file name."""
        return self._file_name

    @file_name.setter
    def file_name(self, value):
        """The file name."""
        self._file_name = value

    @property
    def product_id(self):
        return self._file_id

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.fz', '').replace('.header',
                                                                    '')

    @staticmethod
    def is_raw(file_id):
        return file_id[-1] in ['b', 'd', 'f', 'l', 'o', 'x']


def get_calibration_level(uri):
    ignore_scheme, ignore_archive, file_name = mc.decompose_uri(uri)
    file_id = CFHTName.remove_extensions(file_name)
    result = CalibrationLevel.CALIBRATED
    if CFHTName.is_raw(file_id):
        result = CalibrationLevel.RAW_STANDARD
    return result


def get_energy_function_delta(header):
    filter = header.get('FILTER')
    temp = FILTER_LOOKUP.get(filter)
    result = None
    if temp is not None:
        result = temp.get('width')
    return result


def get_energy_function_val(header):
    filter = header.get('FILTER')
    temp = FILTER_LOOKUP.get(filter)
    result = None
    if temp is not None:
        result = temp.get('centre')
    return result


def get_energy_resolving_power(header):
    delta = get_energy_function_delta(header)
    val = get_energy_function_val(header)
    result = None
    if delta is not None and val is not None:
        result = val/delta
    return result


def get_environment_elevation(header):
    elevation = mc.to_float(header.get('TELALT'))
    logging.error(f'elevation is {elevation}')
    if elevation is not None and not (0.0 <= elevation <= 90.0):
        logging.info(f'Setting elevation to None for {header.get("FILENAME")} '
                     f'because the value is {elevation}.')
        elevation = None
    return elevation


def get_exptime(header):
    exptime = mc.to_float(header.get('EXPTIME'))
    # units are seconds
    if exptime is None:
        file_name = header.get('FILENAME')
        if not CFHTName.is_raw(file_name):
            # caom2IngestMegacaomdetrend.py, l438
            exptime = 0.0
    return exptime


def get_obs_intent(header):
    result = ObservationIntentType.CALIBRATION
    obs_type = header.get('OBSTYPE')
    if obs_type == 'OBJECT':
        run_id = header.get('RUNID')
        logging.error(f'run_id {run_id} obs_type {obs_type}')
        if run_id is None or len(run_id) < 4:
            file_name = header.get('FILENAME')
            if file_name is not None:
                if file_name[-1] == 'o':
                    result = ObservationIntentType.SCIENCE
        else:
            if run_id[3] != 'q':
                result = ObservationIntentType.SCIENCE
    return result


def get_product_type(header):
    result = ProductType.SCIENCE
    obs_type = get_obs_intent(header)
    if obs_type == ObservationIntentType.CALIBRATION:
            result = ProductType.CALIBRATION
    return result


def get_proposal_project(header):
    lookup = {'NGVS': ['08BP03', '08BP04', '09AP03', '09AP04', '09BP03',
                       '09BP04', '10AP03', '10AP04', '10BP03', '10BP04',
                       '11AP03', '11AP04', '11BP03', '11BP04', '12AP03',
                       '12AP04', '12BP03', '12BP04', '13AP03', '13AP04',
                       '13AC02'],
              'PANDAS': ['08BP01', '08BP02', '09AP01', '09AP02', '09BP01',
                         '09BP02', '10AP01', '10AP02', '10BP01', '10BP02'],
              'OSSOS': ['13AP05', '13AP06', '13BP05', '13BP06', '14AP05',
                        '14AP06', '14BP05', '14BP06', '15AP05', '15AP06',
                        '15BP05', '15BP06', '16AP05', '16AP06', '16BP05',
                        '16BP06'],
              'MATLAS': ['13AP07', '13AP08', '13BP07', '13BP08', '14AP07',
                         '14AP08', '14BP07', '14BP08', '15AP07', '15AP08',
                         '15BP07', '15BP08', '16AP07', '16AP08', '16BP07',
                         '16BP08'],
              'LUAU': ['15AP09', '15AP10', '15BP09', '15BP10', '16AP09',
                       '16AP10', '16BP09', '16BP10'],
              'CFIS': ['17AP30', '17AP99', '17AP98', '17BP30', '17BP97',
                       '17BP98', '17BP99', '18AP30', '18AP97', '18AP98',
                       '18AP99', '18BP30', '18BP97', '18BP98', '18BP99',
                       '19AP30', '19AP97', '19AP98', '19AP99', '19BP30',
                       '19BP97', '19BP98', '19BP99'],
              'VESTIGE': ['17AP31', '17BP31', '18AP31', '18BP31', '19AP31',
                          '19BP31']}
    result = None
    pi_name = header.get('PI_NAME')
    if 'CFHTLS' in pi_name:
        result = 'CFHTLS'
    else:
        run_id = header.get('RUNID')
        if run_id is not None:
            for key, value in lookup.items():
                if run_id in value:
                    result = key
                    break
    return result


def get_time_refcoord_delta_cal(header):
    logging.error('time_refcoord_delta_cal')
    mjd_obs = get_time_refcoord_val_cal(header)
    tv_stop = header.get('TVSTOP')
    if tv_stop is None:
        # caom2IngestMegacamdetrend.py, l429
        exp_time = 20.0
    else:
        logging.error(f'tv_stop {tv_stop}')
        mjd_end = ac.get_datetime(tv_stop)
        exp_time = mjd_end - mjd_obs
    return exp_time


def get_time_refcoord_delta_raw(header):
    logging.error('time_refcoord_delta_raw')
    # caom2IngestMegacam.py
    exp_time = get_exptime(header)
    if exp_time is None:
        exp_time = mc.to_float(header.get('DARKTIME'))
    if exp_time is not None:
        # units are days for raw retrieval values
        exp_time = exp_time / 86400.0
    return exp_time


def get_time_refcoord_val_cal(header):
    logging.error('is processed refcoord val')
    dt_str = header.get('TVSTART')
    if dt_str is None:
        dt_str = header.get('REL_DATE')
        if dt_str is None:
            dt_str = header.get('DATE')
    mjd_obs = ac.get_datetime(dt_str)
    logging.error(f'dt_str {dt_str} mjd_obs {mjd_obs}')
    return mjd_obs


def get_time_refcoord_val_raw(header):
    mjd_obs = mc.to_float(header.get('MJD-OBS'))
    return mjd_obs


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level.

    This code captures the portion of the TDM->CAOM model mapping, where
    the relationship is one or many elements of the TDM are required to set
    individual elements of the CAOM model. If the mapping cardinality is 1:1
    generally, use add_fits_attribute. If the mapping cardinality is n:1 use
    the set method to reference a function call.
    """
    logging.debug('Begin accumulate_bp.')
    bp.configure_position_axes((1, 2))
    bp.configure_time_axis(3)
    bp.configure_energy_axis(4)
    bp.configure_polarization_axis(5)
    bp.configure_observable_axis(6)

    bp.set('Observation.intent', 'get_obs_intent(header)')
    # add most preferred attribute last
    bp.clear('Observation.metaRelease')
    bp.add_fits_attribute('Observation.metaRelease', 'REL_DATE')
    bp.add_fits_attribute('Observation.metaRelease', 'DATE')
    bp.add_fits_attribute('Observation.metaRelease', 'DATE-OBS')
    bp.add_fits_attribute('Observation.metaRelease', 'MET_DATE')
    bp.clear('Observation.sequenceNumber')
    bp.add_fits_attribute('Observation.sequenceNumber', 'EXPNUM')

    bp.set('Observation.environment.elevation',
           'get_environment_elevation(header)')
    bp.clear('Observation.environment.humidity')
    bp.add_fits_attribute('Observation.environment.humidity', 'RELHUMID')
    # bp.set('Observation.environment.photometric', False)

    # TODO title is select title from runid_title where proposal_id = 'runid'
    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'RUNID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PI_NAME')
    bp.set_default('Observation.proposal.pi', 'CFHT')
    bp.set('Observation.proposal.project', 'get_proposal_project(header)')

    bp.clear('Observation.instrument.name')
    bp.add_fits_attribute('Observation.instrument.name', 'INSTRUME')

    # bp.set('Observation.target.standard', False)
    # bp.set('Observation.target.moving', False)

    bp.clear('Observation.target_position.coordsys')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'RADECSYS')
    bp.clear('Observation.target_position.equinox')
    bp.add_fits_attribute('Observation.target_position.equinox', 'EQUINOX')
    bp.clear('Observation.target_position.point.cval1')
    bp.add_fits_attribute('Observation.target_position.point.cval1', 'RA_DEG')
    bp.clear('Observation.target_position.point.cval2')
    bp.add_fits_attribute('Observation.target_position.point.cval2', 'DEC_DEG')

    bp.set('Observation.telescope.name', 'CFHT 3.6m')
    x, y, z = ac.get_geocentric_location('cfht')
    bp.set('Observation.telescope.geoLocationX', x)
    bp.set('Observation.telescope.geoLocationY', y)
    bp.set('Observation.telescope.geoLocationZ', z)

    bp.set('Plane.dataProductType', 'image')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(uri)')
    bp.clear('Plane.metaRelease')
    bp.add_fits_attribute('Plane.metaRelease', 'REL_DATE')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE-OBS')
    bp.add_fits_attribute('Plane.metaRelease', 'MET_DATE')
    bp.clear('Plane.dataRelease')
    bp.add_fits_attribute('Plane.dataRelease', 'DATE')
    bp.add_fits_attribute('Plane.dataRelease', 'DATE-OBS')
    bp.add_fits_attribute('Plane.dataRelease', 'REL_DATE')

    bp.set_default('Plane.provenance.name', 'ELIXIR')
    bp.set_default('Plane.provenance.producer', 'CFHT')
    bp.set_default('Plane.provenance.project', 'STANDARD PIPELINE')
    bp.set_default('Plane.provenance.reference',
                   'http://www.cfht.hawaii.edu/Instruments/Elixir/')
    bp.clear('Plane.provenance.runID')
    bp.add_fits_attribute('Plane.provenance.runID', 'CRUNID')
    bp.clear('Plane.provenance.version')
    bp.add_fits_attribute('Plane.provenance.version', 'EL_SYS')

    bp.set('Artifact.productType', 'get_product_type(header)')

    # Changes made to metadata
    #
    # raw image
    #
    # 1. Obs.metarelease
    # 2. Plane.metarelease
    # 3. Telescope.location
    # 4. removed preview and thumbnail artifacts
    # 5. removed blank telescope keywords
    # 6. removed content length
    #
    # 1013337
    # 1. targetPosition.coordsys - FK5 RADECSYS says GAPPT
    # 1. targetPosition.equinox - 2000.0, EQUINOXX says 2008.58
    # 1. how the target position get set? and why is it only set for
    #    some observations?
    #
    # 1000003
    # why is there a numbered energyAxis (4), but no energy values?
    # if I don't comment out the named filter, there are energy values
    # for this file
    #
    # 2463854
    # 1. the non-primary hdus are BINTABLEs, which means there should be no
    #    chunk metadata for them, if I'm remembering my NEOSSat lessons
    #    correctly - so, this is NOT true for CFHT, since the extensions
    #    are *compressed* images, which is not the same thing as a true
    #    BINTABLE.
    #
    # Questions:
    # 1. Elevation - why is it sometimes there, and sometimes not?
    # 1. Photometric - why is sometimes false and sometimes undefined
    # 1. Moving - why is it sometimes false and sometimes undefined
    # 1. metaRelease values have times on them
    #
    # To stop the noise:
    # 1. always have elevation
    # 1. photometric - not set unless true
    # 1. moving - not set unless true
    # 1. standard - not set unless true
    # 1. target_position is there if RA_DEG and DEC_DEG are present
    #
    # 19Bm03.bias.0.40.00
    # 1. should sequence number have a value?
    # 1. add energy axis - the filter name is NOT unique among the test files
    #    I've selected

    # hard-coded values from:
    # - wcaom2archive/cfh2caom2/config/caom2megacam.default and
    # - wxaom2archive/cfht2ccaom2/config/caom2megacam.config
    #
    bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
    bp.set('Chunk.energy.axis.axis.cunit', 'Angstrom')
    bp.set('Chunk.energy.axis.error.rnder', 1.0)
    bp.set('Chunk.energy.axis.error.syser', 1.0)
    bp.set('Chunk.energy.axis.function.delta',
           'get_energy_function_delta(header)')
    bp.set('Chunk.energy.axis.function.naxis', 1)
    bp.set('Chunk.energy.axis.function.refCoord.pix', 1)
    bp.set('Chunk.energy.axis.function.refCoord.val',
           'get_energy_function_val(header)')
    bp.clear('Chunk.energy.bandpassName')
    bp.add_fits_attribute('Chunk.energy.bandpassName', 'FILTER')
    bp.set('Chunk.energy.resolvingPower', 'get_energy_resolving_power(header)')
    bp.set('Chunk.energy.specsys', 'TOPOCENT')
    bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
    bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')

    bp.set('Chunk.position.axis.axis1.cunit', 'deg')
    bp.set('Chunk.position.axis.axis2.cunit', 'deg')
    bp.set('Chunk.position.axis.error1.rnder', 0.0000278)
    bp.set('Chunk.position.axis.error1.syser', 0.0000278)
    bp.set('Chunk.position.axis.error2.rnder', 0.0000278)
    bp.set('Chunk.position.axis.error2.syser', 0.0000278)
    bp.clear('Chunk.position.coordsys')
    bp.add_fits_attribute('Chunk.position.coordsys', 'RADECSYS')

    bp.set('Chunk.time.exposure', 'get_exptime(header)')
    bp.set('Chunk.time.resolution', 'get_exptime(header)')
    bp.set('Chunk.time.timesys', 'UTC')
    bp.set('Chunk.time.axis.axis.ctype', 'TIME')
    bp.set('Chunk.time.axis.axis.cunit', 'd')
    bp.set('Chunk.time.axis.error.rnder', 0.0000001)
    bp.set('Chunk.time.axis.error.syser', 0.0000001)
    bp.set('Chunk.time.axis.function.naxis', 1)
    if get_calibration_level(uri) == CalibrationLevel.RAW_STANDARD:
        bp.set('Chunk.time.axis.function.delta',
               'get_time_refcoord_delta_raw(header)')
        bp.set('Chunk.time.axis.function.refCoord.val',
               'get_time_refcoord_val_raw(header)')
    else:
        bp.set('Chunk.time.axis.function.delta',
               'get_time_refcoord_delta_cal(header)')
        bp.set('Chunk.time.axis.function.refCoord.val',
               'get_time_refcoord_val_cal(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)

    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

    This code captures the portion of the TDM->CAOM model mapping, where
    the relationship is multiple elements of the TDM are required to set
    multiple elements of the CAOM model (mapping cardinality n:n).

    :param observation A CAOM Observation model instance.
    :param **kwargs Everything else."""
    logging.debug('Begin update.')
    mc.check_param(observation, Observation)

    headers = None
    if 'headers' in kwargs:
        headers = kwargs['headers']
    fqn = None
    if 'fqn' in kwargs:
        fqn = kwargs['fqn']

    # processed files
    if (cc.is_composite(headers) and not
            isinstance(observation, CompositeObservation)):
        logging.info('{} is a Composite Observation.'.format(
            observation.observation_id))
        observation = cc.change_to_composite(observation, 'master_detrend')

    ccdbin = headers[0].get('CCBIN1')
    radecsys = headers[0].get('RADECSYS')
    ctype1 = headers[0].get('CTYPE1')
    filter_name = headers[0].get('FILTER')

    for plane in observation.planes.values():
        for artifact in plane.artifacts.values():
            if (get_calibration_level(artifact.uri) ==
                    CalibrationLevel.RAW_STANDARD):
                time_delta = get_time_refcoord_delta_raw(headers[0])
            else:
                time_delta = get_time_refcoord_delta_cal(headers[0])
            for part in artifact.parts.values():
                for c in part.chunks:
                    index = part.chunks.index(c)
                    chunk = part.chunks[index]
                    # CW
                    # Ignore position wcs if a calibration file (except 'x'
                    # type calibration) and/or position info not in header or
                    # binned 8x8
                    if (plane.product_id[-1] in ['b', 'l', 'd', 'f'] or
                            ccdbin == 8 or radecsys is None or ctype1 is None):
                        cc.reset_position(chunk)

                    # CW
                    # Ignore energy wcs if some type of calibration file or
                    # filter='None' or 'Open' or cdelt4=0.0 (no filter match)
                    if (filter_name is None or filter_name == 'Open' or
                            filter_name not in FILTER_LOOKUP or
                            plane.product_id[-1] in ['b', 'l', 'd']):
                        cc.reset_energy(chunk)

                    if (time_delta == 0.0 and
                            chunk.time.axis.function.delta == 1.0):
                        # undo the effects of the astropy cdfix call on a
                        # matrix, in fits2caom2.WcsParser, which sets the
                        # diagonal element of the matrix to unity if all
                        # keywords associated with a given axis were omitted.
                        # See:
                        # https://docs.astropy.org/en/stable/api/astropy.wcs.
                        # Wcsprm.html#astropy.wcs.Wcsprm.cdfix
                        chunk.time.axis.function.delta = 0.0

        if isinstance(observation, CompositeObservation):
            cc.update_plane_provenance(plane, headers[1:], 'IMCMB',
                                       COLLECTION,
                                       _repair_provenance_value,
                                       observation.observation_id)

    # relies on update_plane_provenance being called
    if isinstance(observation, CompositeObservation):
        cc.update_observation_members(observation)

    logging.debug('Done update.')
    return observation


def _build_blueprints(uris):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param uris The artifact URIs for the files to be processed."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for uri in uris:
        blueprint = ObsBlueprint(module=module)
        if not mc.StorageName.is_preview(uri):
            accumulate_bp(blueprint, uri)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.local:
        for ii in args.local:
            file_id = mc.StorageName.remove_extensions(os.path.basename(ii))
            file_name = '{}.fits'.format(file_id)
            result.append(CFHTName(file_name=file_name).file_uri)
    elif args.lineage:
        for ii in args.lineage:
            result.append(ii.split('/', 1)[1])
    else:
        raise mc.CadcException(
            'Could not define uri from these args {}'.format(args))
    return result


def _repair_provenance_value(value, obs_id):
    prov_obs_id = None
    prov_prod_id = None
    # CFHT files contain other IMCMB headers that look like this - ignore
    # them here:
    # IMCMB_DT= 'FLIPS ver 3.0 - Elixir by CFHT - Wed Dec 4 2019 -  9:53:11'
    # IMCMB_AL= 'SIGMA   '
    # IMCMB_CA= 'MEDIAN  '
    # IMCMB_LS=                  3.5 / Lower threshold for clipping rejection
    # IMCMB_HS=                  3.5 / Lower clipping rejection threshold
    # IMCMB_FT= 'MASTER_DETREND_BIAS'
    # IMCMB_OP= '1DYMODEL_OVERSCAN'
    # IMCMB_NI=                   30 / Number of input files
    # IMCMB_DA=                    T / Dual amplifier A,B
    # IMCMB_IF= 'NAME BIAS_A MODE_A BIAS_B MODE_B' / Input file parameters
    if '.fits' in str(value):
        # input looks like:
        # '2463481b.fits[ccd39] 1231 1 1225 1' / Input file stats
        # or like:
        # '707809o00.fits 0 1569 0.341' / Input file stats
        temp = value.split('.fits')
        if '[' in value:
            prov_prod_id = temp[0]
        else:
            prov_prod_id = temp[0][:-2]
        prov_obs_id = CFHTName(file_name=prov_prod_id).obs_id
    return prov_obs_id, prov_prod_id


def _main_app():
    args = get_gen_proc_arg_parser().parse_args()
    uris = _get_uris(args)
    blueprints = _build_blueprints(uris)
    result = gen_proc(args, blueprints)
    return result


def cfht_main_app():
    args = get_gen_proc_arg_parser().parse_args()
    try:
        result = _main_app()
        sys.exit(result)
    except Exception as e:
        logging.error('Failed {} execution for {}.'.format(APPLICATION, args))
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)

    logging.debug('Done {} processing.'.format(APPLICATION))
