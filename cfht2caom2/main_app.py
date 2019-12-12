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
from caom2 import ProductType
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc


__all__ = ['cfht_main_app', 'update', 'CFHTName', 'COLLECTION',
           'APPLICATION', 'ARCHIVE']


APPLICATION = 'cfht2caom2'
COLLECTION = 'CFHT'
ARCHIVE = 'CFHT'

FILTER_LOOKUP = {'u.MP9301': {'centre': 3754.5, 'width': 867 / 2.0},
                 'g.MP9401': {'centre': 4890.0, 'width': 1624 / 2.0},
                 'r.MP9601': {'centre': 6248.5, 'width': 1413.0},
                 'i.MP9701': {'centre': 7763., 'width': 1726 / 2.0},
                 'i.MP9702': {'centre': 7623., 'width': 1728 / 2.0},
                 # 'z.MP9801': {'centre': 9083., 'width': 1848.0},
                 }

# All comments denoted 'CW' are copied from
# ssh://goliaths@gimli3/srv/cadc/git/wcaom2archive


class CFHTName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support fz files in storage
    - product id == file id
    """

    CFHT_NAME_PATTERN = '*'

    def __init__(self, obs_id=None, fname_on_disk=None, file_name=None):
        super(CFHTName, self).__init__(
            None, COLLECTION, CFHTName.CFHT_NAME_PATTERN, file_name,
            compression='.fz')
        self._file_id = CFHTName.remove_extensions(file_name)
        self.obs_id = self._file_id[:-1]
        logging.debug(self)

    def __str__(self):
        return f'obs_id {self.obs_id}, ' \
               f'file_id {self._file_id}, ' \
               f'file_name {self.file_name}, ' \
               f'lineage {self.lineage}, ' \
               f'external urls {self.external_urls}'

    def is_valid(self):
        return True

    @property
    def product_id(self):
        return self._file_id

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        return name.replace('.fits', '').replace('.fz', '').replace('.header',
                                                                    '')


def get_calibration_level(uri):
    ignore_scheme, ignore_archive, file_name = mc.decompose_uri(uri)
    file_id = CFHTName.remove_extensions(file_name)
    logging.error(file_id)
    result = None
    if file_id[-1] in ['b', 'd', 'f', 'l', 'o', 'x']:
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


def get_exptime(header):
    exptime = mc.to_float(header.get('EXPTIME'))
    # units are seconds
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


def get_time_delta(header):
    exptime = get_exptime(header)
    if exptime is None:
        exptime = mc.to_float(header.get('DARKTIME'))
    if exptime is not None:
        # units are days
        exptime = exptime / 86400.0
    return exptime


def get_time_refcoord_val(header):
    mjd_obs = mc.to_float(header.get('MJD-OBS'))
    if mjd_obs is None:
        pass
    return mjd_obs


def accumulate_bp(bp, uri):
    """Configure the telescope-specific ObsBlueprint at the CAOM model 
    Observation level."""
    logging.debug('Begin accumulate_bp.')
    bp.configure_position_axes((1, 2))
    bp.configure_time_axis(3)
    bp.configure_energy_axis(4)
    bp.configure_polarization_axis(5)
    bp.configure_observable_axis(6)

    bp.set('Observation.intent', 'get_obs_intent(header)')
    bp.clear('Observation.sequenceNumber')
    bp.add_fits_attribute('Observation.sequenceNumber', 'EXPNUM')

    bp.clear('Observation.environment.elevation')
    bp.add_fits_attribute('Observation.environment.elevation', 'TELALT')
    bp.clear('Observation.environment.humidity')
    bp.add_fits_attribute('Observation.environment.humidity', 'RELHUMID')
    # bp.set('Observation.environment.photometric', False)

    # TODO title is select title from runid_title where proposal_id = 'runid'
    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'RUNID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PI_NAME')

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
    bp.add_fits_attribute('Plane.metaRelease', 'MET_DATE')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE-OBS')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE')
    bp.clear('Plane.dataRelease')
    bp.add_fits_attribute('Plane.dataRelease', 'REL_DATE')
    bp.add_fits_attribute('Plane.dataRelease', 'DATE-OBS')
    bp.add_fits_attribute('Plane.dataRelease', 'DATE')

    bp.set_default('Plane.provenance.name', 'ELIXIR')
    bp.set_default('Plane.provenance.producer', 'CFHT')
    bp.set_default('Plane.provenance.project', 'STANDARD PIPELINE')
    bp.set_default('Plane.provenance.reference',
                   'http://www.cfht.hawaii.edu/Instruments/Elixir/')
    bp.clear('Plane.provenance.runID')
    bp.add_fits_attribute('Plane.provenance.runID', 'CRUNID')

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
    bp.set('Chunk.time.axis.function.delta', 'get_time_delta(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)
    bp.set('Chunk.time.axis.function.refCoord.val',
           'get_time_refcoord_val(header)')

    logging.debug('Done accumulate_bp.')


def update(observation, **kwargs):
    """Called to fill multiple CAOM model elements and/or attributes, must
    have this signature for import_module loading and execution.

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

    ccdbin = headers[0].get('CCBIN1')
    radecsys = headers[0].get('RADECSYS')
    ctype1 = headers[0].get('CTYPE1')
    filter_name = headers[0].get('FILTER')
    time_delta = get_time_delta(headers[0])

    for plane in observation.planes.values():
        for artifact in plane.artifacts.values():
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
