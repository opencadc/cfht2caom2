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

"""
CFHT Cardinality:

CW - 02-01-20
The CADC page describing file types:
http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/cfht/extensions.html

For Espadons, Sitelle and Spirou there are composite observations made up from
several exposures. These composite observations are given the ‘p’ ending to
distinguish them from the simple observations that are based on a single
exposure.

For WIRcam, based on the CAOM model the o and g are included as separate
artifacts of the same plane because they correspond to different groups of
raw pixels on the detector. This has caused much trouble for the need to
connect them in etransfer. For SPIRou we did some combining of file types
within a plane e.g. ‘v’ within ‘e’, and CFHT say they would prefer to keep
them separate so we should do that.

Conversation with CW/SF 02-01-19:
- SITELLE - has hdf5 files, but there is currently no way of extracting WCS
            information from them
          - there will be (eventually) a fits file for every hdf5 file
          - copy the fits wcs to the hdf5 wcs, depending on the order of
            arrival of the files
- SPIROU - no hdf5 files, but there are two file types that have binary table
           extensions that won't support cut-outs
- WIRCam - no hdf5 files
         - 'g' are cubes of guide window(s)
- 'p' files:
    - MegaCam/WIRCam - 'p' is a processed file that is an additional plane to
            a SimpleObservation
    - SITELLE - 'p' is processed + Composite - a different observation
    - SPIROU/Espadons - 'p' is polarized + Composite (TBC on Composite), a
      different observation
    - there are other processed files for single exposures
    - users want 'p' files

- conclusion - one plane / file, because users want to see one row / file in
  the results tab
"""
import copy
import importlib
import logging
import math
import os
import sys
import traceback

from caom2 import Observation, CalibrationLevel, ObservationIntentType
from caom2 import ProductType, CompositeObservation, TypedList, Chunk
from caom2 import CoordRange2D, CoordAxis2D, Axis, Coord2D, RefCoord
from caom2 import SpatialWCS, DataProductType, ObservationURI, PlaneURI
from caom2 import CoordAxis1D, CoordRange1D, SpectralWCS, Slice
from caom2 import ObservableAxis
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2utils import FitsParser, WcsParser, get_cadc_headers
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc

from cfht2caom2 import metadata as md


__all__ = ['cfht_main_app', 'update', 'CFHTName', 'COLLECTION',
           'APPLICATION', 'ARCHIVE']


APPLICATION = 'cfht2caom2'
COLLECTION = 'CFHT'
ARCHIVE = 'CFHT'

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

    def __init__(self, obs_id=None, fname_on_disk=None, file_name=None,
                 instrument=None, ad_uri=None):
        # set compression to an empty string so the file uri method still
        # works, since the file_name element will have all extensions,
        # including the .fz to indicate compression  # TODO no longer true
        super(CFHTName, self).__init__(
            None, COLLECTION, CFHTName.CFHT_NAME_PATTERN, file_name,
            compression='')
        self._instrument = md.Inst(instrument)
        if ad_uri is not None and file_name is None:
            file_name = mc.CaomName(ad_uri).file_name
        self._file_id = CFHTName.remove_extensions(file_name)
        self._suffix = self._file_id[-1]
        if self.is_simple and not self.is_master_cal:
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
        result = self._file_id
        # TODO this has changed
        if (self._file_id[-1] in ['g', 'o'] and
                self._instrument is md.Inst.WIRCAM):
            result = f'{self._file_id[:-1]}og'
        return result

    @property
    def is_master_cal(self):
        return ('weight' in self._file_id or 'master' in self._file_id or
                    'hotpix' in self._file_id or 'badpix' in self._file_id or
                    'deadpix' in self._file_id or 'dark' in self._file_id)

    @property
    def has_polarization(self):
        return self._suffix in ['p'] and self._instrument is md.Inst.ESPADONS

    @property
    def is_simple(self):
        result = False
        if (self._suffix in ['a', 'b', 'c', 'd', 'f', 'g', 'l', 'm', 'o', 'x']
                or self.simple_by_suffix or self.is_master_cal):
            result = True
        return result

    @property
    def simple_by_suffix(self):
        return ((self._suffix == 'p' and
                 self._instrument in [md.Inst.MEGACAM, md.Inst.WIRCAM]) or
                (self._suffix == 'i' and self._instrument is md.Inst.ESPADONS))

    @property
    def suffix(self):
        return self._suffix

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        # ESPaDOnS files have a .gz extension ;)
        return name.replace('.fits', '').replace('.fz', '').replace(
            '.header', '').replace('.gz', '')


def get_bandpass_name(header):
    instrument = _get_instrument(header)
    result = header.get('FILTER')
    if instrument is md.Inst.WIRCAM:
        wheel_a = header.get('WHEELADE')
        wheel_b = header.get('WHEELBDE')
        if wheel_a == 'Open' and wheel_b != 'Open':
            result = wheel_b
        elif wheel_b == 'Open' and wheel_a != 'Open':
            result = wheel_a
        elif wheel_a == 'Open' and wheel_b == 'Open':
            result = 'Open'
    return result


def get_calibration_level(params):
    header, suffix, uri = _decompose_params(params)
    instrument = _get_instrument(header)
    cfht_name = CFHTName(ad_uri=uri, instrument=instrument)
    result = CalibrationLevel.CALIBRATED
    if (cfht_name.is_simple and not cfht_name.simple_by_suffix and not
            cfht_name.is_master_cal):
        result = CalibrationLevel.RAW_STANDARD
    return result


def get_energy_ctype(header):
    result = None
    if _has_energy(header):
        result = header.get('CTYPE3')
        if result is None:
            result = 'WAVE'
    return result


def get_energy_cunit(header):
    result = None
    if _has_energy(header):
        result = header.get('CUNIT3')
        if result is None:
            result = 'Angstrom'
    return result


def get_energy_function_delta(params):
    result = None
    header, suffix, uri = _decompose_params(params)
    if _has_energy(header):
        if _is_espadons_energy(header, uri):
            # caom2IngestEspadons.py l639
            result = 0.0031764
        elif _is_sitelle_energy(header, uri):
            result = header.get('CDELT3')
            if result is None:
                # units in file are nm, units in blueprint are Angstroms
                result = 10.0 * mc.to_float(header.get('FILTERBW'))
        else:
            filter_name = header.get('FILTER')
            instrument = _get_instrument(header)
            temp = md.filter_cache.get_svo_filter(instrument, filter_name)
            result = ac.FilterMetadataCache.get_fwhm(temp)
    return result


def get_energy_function_naxis(params):
    result = 1.0
    header, suffix, uri = _decompose_params(params)
    if _is_espadons_energy(header, uri):
        # caom2IngestEspadons.py l636
        result = 213542
    elif _is_sitelle_energy(header, uri):
        result = header.get('NAXIS3')
        if result is None:
            result = 1.0
    return result


def get_energy_function_pix(params):
    result = None
    header, suffix, uri = _decompose_params(params)
    if _has_energy(header):
        result = 1.0
        if _is_espadons_energy(header, uri):
            # caom2IngestEspadons.py l637
            result = 0.5
        elif _is_sitelle_energy(header, uri):
            result = header.get('CRPIX3')
            if result is None:
                result = 0.5
    return result


def get_energy_function_val(params):
    result = None
    header, suffix, uri = _decompose_params(params)
    if _has_energy(header):
        if _is_espadons_energy(header, uri):
            # caom2IngestEspadons.py l638
            result = 370.0
        elif _is_sitelle_energy(header, uri):
            result = header.get('CRVAL3')
            if result is None:
                # units in file are nm, units in blueprint are Angstroms
                result = 10.0 * mc.to_float(header.get('FILTERLB'))
        else:
            filter_name = header.get('FILTER')
            instrument = _get_instrument(header)
            temp = md.filter_cache.get_svo_filter(instrument, filter_name)
            result = ac.FilterMetadataCache.get_central_wavelength(temp)
    return result


def get_energy_resolving_power(params):
    result = None
    header, suffix, uri = _decompose_params(params)
    if _has_energy(header):
        delta = get_energy_function_delta(params)
        val = get_energy_function_val(params)
        result = None
        if delta is not None and val is not None:
            result = val/delta
        instrument = _get_instrument(header)
        if instrument is md.Inst.ESPADONS:
            instmode = header.get('INSTMODE')
            if instmode is None or 'R=' not in instmode:
                # CW - Default if resolving power value not in header
                # caom2IngestEspadons.py, l377
                result = 65000.0
            else:
                # CW - This string is already in instrument keywords but also
                # need to extract resolving power from it:
                # 'Spectroscopy, star only, R=80,000'
                temp = instmode.split('R=')
                values = temp[1].split(',')
                if len(values) == 1:
                    result = values[0]
                else:
                    result = f'{values[0]}{values[1]}'
                result = mc.to_float(result)
        elif instrument is md.Inst.SITELLE:
            sitresol = header.get('SITRESOL')
            if sitresol is not None and sitresol > 0.0:
                result = sitresol
            if result is None:
                result = 1.0
    return result


def get_environment_elevation(header):
    elevation = mc.to_float(header.get('TELALT'))
    if elevation is not None and not (0.0 <= elevation <= 90.0):
        logging.info(f'Setting elevation to None for {_get_filename(header)} '
                     f'because the value is {elevation}.')
        elevation = None
    return elevation


def get_exptime(params):
    header, suffix, uri = _decompose_params(params)
    exptime = mc.to_float(header.get('EXPTIME'))
    instrument = _get_instrument(header)
    if instrument is md.Inst.ESPADONS and suffix == 'p':
        # caom2IngestEspadons.py, l406
        exptime = 0.0
        polar_seq = mc.to_int(header.get('POLARSEQ'))
        for ii in range(1, polar_seq + 1):
            exptime += mc.to_float(header.get(f'EXPTIME{ii}'))
    elif instrument is md.Inst.SITELLE:
        if suffix == 'p':
            num_steps = header.get('STEPNB')
            exptime = exptime * num_steps
    # units are seconds
    if exptime is None:
        cfht_name = CFHTName(ad_uri=uri, instrument=instrument)
        if cfht_name.is_simple:
            # caom2IngestMegacaomdetrend.py, l438
            exptime = 0.0
    return exptime


def get_instrument_keywords(header):
    inst_mode = header.get('INSTMODE')
    temp = header.get('SITSTEP')
    sit_step = None
    if temp is not None:
        sit_step = f'SITSTEP={temp}'
    temp = header.get('SITSTEPS')
    sit_steps = None
    if temp is not None:
        sit_steps = f'SITSTEPS={temp}'
    result = ','.join(filter(None, (inst_mode, sit_step, sit_steps)))
    if 'Unknown' in result:
        result = 'Unknown'
    return result


def get_obs_intent(header):
    # CW
    # Determine Observation.intent = obs.intent = "science" or "calibration"
    # phot & astr std & acquisitions/align are calibration.
    result = ObservationIntentType.CALIBRATION
    obs_type = _get_obstype(header)
    if obs_type == 'OBJECT':
        run_id = header.get('RUNID')
        if run_id is None or len(run_id) < 4:
            file_name = _get_filename(header)
            if file_name is not None:
                if file_name[-1] == 'o':
                    result = ObservationIntentType.SCIENCE
        else:
            if run_id[3] != 'q':
                result = ObservationIntentType.SCIENCE
    return result


def get_obs_sequence_number(params):
    header, suffix, uri = _decompose_params(params)
    instrument = _get_instrument(header)
    cfht_name = CFHTName(ad_uri=uri, instrument=instrument)
    result = None
    exp_num = header.get('EXPNUM')
    if cfht_name.is_simple and not cfht_name.is_master_cal:
        result = exp_num
    elif (instrument in [md.Inst.ESPADONS, md.Inst.SITELLE] and
            suffix == 'p' and exp_num is None):
        result = cfht_name.file_id[:-1]
    return result


def get_obs_type(params):
    header, suffix, uri = _decompose_params(params)
    result = header.get('OBSTYPE')
    instrument = _get_instrument(header)
    if instrument is md.Inst.WIRCAM:
        # caom2IngestWircamdetrend.py, l369
        if 'weight' in uri:
            result = 'WEIGHT'
        elif 'badpix' in uri or 'hotpix' in uri or 'deadpix' in uri:
            result = 'BPM'
        elif suffix == 'g' and result is None:
            result = 'GUIDE'
    return result


def get_plane_data_product_type(header):
    instrument = _get_instrument(header)
    result = DataProductType.IMAGE
    if instrument is md.Inst.ESPADONS:
        result = DataProductType.SPECTRUM
    return result


def get_polarization_function_val(header):
    lookup = {'I': 1,
              'Q': 2,
              'U': 3,
              'V': 4,
              'W': 5}
    result = 6
    temp = header.get('CMMTSEQ')
    if temp is not None:
        result = lookup.get(temp[0], result)
    return result


def get_product_type(params):
    header, suffix, ignore = _decompose_params(params)
    result = ProductType.SCIENCE
    obs_type = get_obs_intent(header)
    if obs_type == ObservationIntentType.CALIBRATION:
        result = ProductType.CALIBRATION
    if suffix in ['m', 'w', 'y']:
        result = ProductType.AUXILIARY
    return result


def get_proposal_project(header):
    # TODO - put in a configuration file, cause it changes by semester
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
                          '19BP31'],
              'SIGNALS': ['18BP41', '19AP41', '19BP41', '20AP41', '20BP41',
                          '21AP41', '21BP41', '22AP41'],
              'WIRDS': ['06AC01', '06AC99', '06BC31', '06BF97', '07AC20',
                        '07AF34', '07BC99', '07BF98', '07BF97', '08AC98',
                        '08AC99', '08AF97', '08BC99', '08BF99'],
              'IPMOS': ['06AF01', '06BF23', '07AF22', '07AF99', '07BF99',
                        '08AF98', '08BF98'],
              'TETrEs': ['10BP21', '11BP21', '12AP21', '12BP21'],
              'CFHQSIR': ['10BP22', '11AP22', '11BP22', '11BP23', '12AD95',
                          '12AP22', '12AP23', '12BP22', '12BP23', '13AF20'
                          '13AC14'],
              'CIPP': ['17AP32', '17BP32', '18AP32', '18BP32', '19AP32',
                       '19BP32'],
              'MAPP': ['08BP11', '08BP12', '09AP11', '09AP12', '09BP11',
                       '09BP12', '10AP11', '10AP12', '10BP11', '10BP12',
                       '11AP11', '11AP12', '11BP11', '11BP12', '12AP11',
                       '12AP12', '12BP11', '12BP12'],
              'MIMES': ['08BP13', '08BP14', '09AP13', '09AP14', '09BP13',
                        '09BP14', '10AP13', '10AP14', '10BP13', '10BP14',
                        '11AP13', '11AP14', '11BP13', '11BP14', '12AP13',
                        '12AP14', '12BP13', '12BP14'],
              'BINAMICS': ['13AP15', '13AP16', '13BP15', '13BP16', '14AP15',
                           '14AP16', '14BP15', '14BP16', '15AP15', '15AP16',
                           '15BP15', '15BP16', '16AP15', '16AP16', '16BP15',
                           '16BP16'],
              'MATYSSE': ['13AP17', '13AP18', '13BP17', '13BP18', '14AP17',
                          '14AP18', '14BP17', '14BP18', '15AP17', '15AP18',
                          '15BP17', '15BP18', '16AP17', '16AP18', '16BP17',
                          '16BP18'],
              'HMS': ['15AP19', '15AP20', '15BP19', '15BP20', '16AP19',
                      '16AP20', '16BP19', '16BP20']}
    result = None
    pi_name = header.get('PI_NAME')
    if pi_name is not None and 'CFHTLS' in pi_name:
        result = 'CFHTLS'
    else:
        run_id = header.get('RUNID')
        if run_id is not None:
            for key, value in lookup.items():
                if run_id in value:
                    result = key
                    break
    return result


def get_provenance_keywords(params):
    header, suffix, ignore = _decompose_params(params)
    instrument = _get_instrument(header)
    result = None
    if instrument is md.Inst.WIRCAM and suffix in ['p', 's']:
        # caom2IngestWircam.py, l1063
        if suffix == 'p':
            result = 'skysubtraction=yes'
        else:
            result = 'skysubtraction=no'
    elif instrument is md.Inst.ESPADONS and suffix in ['i', 'p']:
        temp = header.get('REDUCTIO')
        if temp is not None:
            result = f'reduction={temp}'
    return result


def get_provenance_last_executed(header):
    result = None
    instrument = _get_instrument(header)
    if instrument is md.Inst.ESPADONS:
        comments = header.get('COMMENT')
        if comments is not None:
            for comment in comments:
                if 'Upena processing date:' in comment:
                    result = comment.split('Upena processing date: ')[1]
                    # format like Fri Mar 13 22:51:55 HST 2009, which default
                    # code doesn't understand
                    result = mc.make_time(result)
                    break
                elif 'opera-' in comment:
                    result = comment.split('opera-')[1].split(' build date')[0]
                    break
    else:
        result = header.get('PROCDATE')
        if result is not None:
            # format like 2018-06-05HST17:21:20, which default code doesn't
            # understand
            result = mc.make_time(result)
    return result


def get_provenance_name(header):
    result = 'TCS'  # ESPaDOnS
    instrument = _get_instrument(header)
    if instrument is md.Inst.ESPADONS:
        comments = header.get('COMMENT')
        if comments is not None:
            for comment in comments:
                if 'Upena' in comment:
                    result = 'UPENA'
                    break
                elif 'opera-' in comment:
                    result = 'OPERA'
                    break
    elif instrument is md.Inst.MEGAPRIME:
        result = 'ELIXIR'
    elif instrument is md.Inst.SITELLE:
        result = 'ORBS'
    elif instrument is md.Inst.WIRCAM:
        result = 'IIWI'
    return result


def get_provenance_project(header):
    result = 'STANDARD PIPELINE'
    if get_provenance_name(header) == 'TCS':
        result = None
    return result


def get_provenance_version(header):
    result = None
    instrument = _get_instrument(header)
    if instrument is md.Inst.ESPADONS:
        comments = header.get('COMMENT')
        if comments is not None:
            for comment in comments:
                if 'Upena version' in comment:
                    result = comment.split('Upena version')[1]
                    break
                elif 'opera-' in comment and 'build date' in comment:
                    result = comment.split(' build date')[0]
                    break
    else:
        result = header.get('IIWIVER')
        if result is None:
            result = header.get('ORBSVER')
            if result is None:
                result = header.get('EL_SYS')
    return result


def get_provenance_reference(header):
    lookup = {
        md.Inst.ESPADONS:
            'http://www.cfht.hawaii.edu/Instruments/Spectroscopy/Espadons/',
        md.Inst.MEGAPRIME: 'http://www.cfht.hawaii.edu/Instruments/Elixir/',
        md.Inst.SITELLE: 'http://ascl.net/1409.007',
        md.Inst.WIRCAM:
            'http://www.cfht.hawaii.edu/Instruments/Imaging/WIRCam'}
    instrument = _get_instrument(header)
    return lookup.get(instrument)


def get_target_position_cval1(header):
    ra, ignore_dec = _get_ra_dec(header)
    if ra is None:
        instrument = _get_instrument(header)
        if instrument is md.Inst.ESPADONS:
            ra = header.get('RA_DEG')
    return ra


def get_target_position_cval2(header):
    ignore_ra, dec = _get_ra_dec(header)
    if dec is None:
        instrument = _get_instrument(header)
        if instrument is md.Inst.ESPADONS:
            dec = header.get('DEC_DEG')
    return dec


def get_target_standard(header):
    obs_type = _get_obstype(header)
    run_id = _get_run_id(header)
    result = None
    if run_id is not None:
        run_id_type = run_id[3].lower()
        if run_id_type == 'q' and obs_type == 'OBJECT':
            obj_name = header.get('OBJECT').lower()
            instrument = _get_instrument(header)
            if instrument is md.Inst.SITELLE:
                if 'std' in obj_name:
                    result = True
                else:
                    result = False
            else:
                if ('flat' in obj_name or 'focus' in obj_name or
                        'zenith' in obj_name):
                    result = False
                else:
                    result = True
        else:
            result = False
    return result


def get_time_refcoord_delta_composite(params):
    header, suffix, uri = _decompose_params(params)
    instrument = _get_instrument(header)
    if instrument is md.Inst.ESPADONS and suffix == 'p':
        exp_time = get_exptime(params) / 86400.0  # units are d
    elif instrument is md.Inst.SITELLE:
        mjd_start = _get_mjd_start(header)
        mjd_end = mc.to_float(header.get('MJDEND'))
        if mjd_end is None:
            exp_time = header.get('EXPTIME')
            if exp_time is None:
                exp_time = mjd_start
        else:
            # caom2IngestSitelle.py, l704
            exp_time = mjd_end - mjd_start
            logging.error(f'exp_time {exp_time}')
    else:
        mjd_obs = get_time_refcoord_val_composite(header)
        tv_stop = header.get('TVSTOP')
        if tv_stop is None:
            # caom2IngestMegacamdetrend.py, l429
            # caom2IngestWircamdetrend.py, l422
            exp_time = 20.0
        else:
            mjd_end = ac.get_datetime(tv_stop)
            exp_time = mjd_end - mjd_obs
    return exp_time


def get_time_refcoord_delta_simple(params):
    # caom2IngestMegacam.py
    exp_time = get_exptime(params)
    if exp_time is None:
        header, suffix, uri = _decompose_params(params)
        exp_time = mc.to_float(header.get('DARKTIME'))
    if exp_time is not None:
        # units are days for raw retrieval values
        exp_time = exp_time / 86400.0
    return exp_time


def get_time_refcoord_val_composite(header):
    instrument = _get_instrument(header)
    mjd_start1 = header.get('MJDSTART1')
    mjd_date1 = header.get('MJDATE1')
    if (instrument is md.Inst.ESPADONS and
            (mjd_start1 is not None or mjd_date1 is not None)):
        # caom2IngestEspadons.py, l406
        if mjd_start1 is not None:
            mjd_obs = mjd_start1
        else:
            mjd_obs = mjd_date1
    elif instrument is md.Inst.SITELLE:
        mjd_obs = _get_mjd_start(header)
    else:
        # CW
        # caom2IngestWircamdetrend.py, l388
        # Time - set exptime as time of one image, start and stop dates
        # as one pixel so this means crval3 is not equal to exptime
        # if TVSTART  not defined, use release_date as mjdstart
        dt_str = header.get('TVSTART')
        if dt_str is None:
            dt_str = header.get('REL_DATE')
            if dt_str is None:
                dt_str = header.get('DATE')
        mjd_obs = ac.get_datetime(dt_str)
    return mjd_obs


def get_time_refcoord_val_simple(header):
    mjd_obs = _get_mjd_obs(header)
    if mjd_obs is None:
        instrument = _get_instrument(header)
        if instrument is md.Inst.ESPADONS:
            date_obs = header.get('DATE-OBS')
            time_obs = header.get('TIME-OBS')
            if (date_obs is None or time_obs is None or
                    date_obs == '1970-01-01' or date_obs == '1970-00-01'):
                hst_time = header.get('HSTTIME)')
                # fmt 'Mon Nov 27 15:58:17 HST 2006'
                mjd_obs = ac.get_datetime(hst_time)
            else:
                mjd_obs_str = f'{date_obs}T{time_obs}'
                mjd_obs = ac.get_datetime(mjd_obs_str)
    return mjd_obs


def _decompose_params(params):
    header = params.get('header')
    uri = params.get('uri')
    instrument = _get_instrument(header)
    suffix = CFHTName(ad_uri=uri, instrument=instrument)._suffix
    return header, suffix, uri


def _get_filename(header):
    return header.get('FILENAME')


def _get_instrument(header):
    return md.Inst(header.get('INSTRUME'))


def _get_mjd_obs(header):
    return mc.to_float(header.get('MJD-OBS'))


def _get_mjd_start(header):
    mjd_obs = _get_mjd_obs(header)
    if mjd_obs is None:
        date_str = header.get('DATE-OBS')
        if date_str is None:
            dt_str = header.get('DATE')
            mjd_obs = ac.get_datetime(dt_str)
        else:
            time_str = header.get('TIME-OBS')
            date_obs = ac.get_datetime(date_str)
            time_obs = ac.get_datetime(time_str)
            mjd_obs = date_obs + time_obs
    return mjd_obs


def _get_obstype(header):
    return header.get('OBSTYPE')


def _get_ra_dec(header):
    obj_ra = header.get('OBJRA')
    obj_dec = header.get('OBJDEC')
    obj_ra_dec = header.get('OBJRADEC')
    if obj_ra_dec is not None:
        obj_ra_dec = obj_ra_dec.lower()
    ra = None
    dec = None
    if obj_ra is not None and obj_dec is not None and obj_ra_dec is not None:
        if obj_ra_dec == 'gappt':
            # SF 18-12-19
            # seb 4:01 PM
            # this is a flat. i have the impression in this case you can
            # ignore the ra/dec stuff
            logging.warning(f'OBSRADEC is GAPPT for {_get_filename(header)}')
        else:
            ra, dec = ac.build_ra_dec_as_deg(obj_ra, obj_dec, obj_ra_dec)
    return ra, dec


def _get_run_id(header):
    run_id = header.get('RUNID')
    if (run_id is not None and (len(run_id) < 3 or len(run_id) > 9 or
                                run_id == 'CFHT')):
        # a well-known default value that indicates the past
        run_id = '17BE'
    return run_id


def _get_types(params):
    header, suffix, ignore = _decompose_params(params)
    dp_result = DataProductType.IMAGE
    pt_result = ProductType.SCIENCE
    obs_type = get_obs_intent(header)
    if obs_type == ObservationIntentType.CALIBRATION:
        pt_result = ProductType.CALIBRATION
    if suffix in ['m', 'w', 'y']:
        dp_result = DataProductType.AUXILIARY
        pt_result = ProductType.AUXILIARY
    return dp_result, pt_result


def _has_energy(header):
    obs_type = _get_obstype(header)
    # from conversation with CW, SF
    # also from caom2IngestEspadons.py, l393, despite an existing example
    # with energy information
    return obs_type not in ['BIAS', 'DARK']


def _is_espadons_energy(header, uri):
    instrument = _get_instrument(header)
    result = False
    if instrument is md.Inst.ESPADONS:
        cfht_name = CFHTName(ad_uri=uri, instrument=instrument)
        if cfht_name.suffix in ['a', 'b', 'c', 'd', 'f', 'i', 'o', 'p', 'x']:
            result = True
    return result


def _is_sitelle_energy(header, uri):
    instrument = _get_instrument(header)
    result = False
    if instrument is md.Inst.SITELLE:
        cfht_name = CFHTName(ad_uri=uri, instrument=instrument)
        if cfht_name.suffix in ['a', 'c', 'f', 'o', 'p', 'x']:
            result = True
    return result


def accumulate_bp(bp, uri, instrument):
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
    bp.configure_observable_axis(6)

    bp.set('Observation.intent', 'get_obs_intent(header)')
    # add most preferred attribute last
    bp.clear('Observation.metaRelease')
    bp.add_fits_attribute('Observation.metaRelease', 'REL_DATE')
    bp.add_fits_attribute('Observation.metaRelease', 'DATE')
    # caom2IngestEspadons.py, l625
    bp.add_fits_attribute('Observation.metaRelease', 'DATE-OB1')
    bp.add_fits_attribute('Observation.metaRelease', 'DATE-OBS')
    bp.add_fits_attribute('Observation.metaRelease', 'MET_DATE')

    bp.set('Observation.sequenceNumber', 'get_obs_sequence_number(params)')
    bp.set('Observation.type', 'get_obs_type(params)')

    bp.set('Observation.environment.elevation',
           'get_environment_elevation(header)')
    bp.clear('Observation.environment.humidity')
    bp.add_fits_attribute('Observation.environment.humidity', 'RELHUMID')

    # TODO title is select title from runid_title where proposal_id = 'runid'
    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'RUNID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PI_NAME')
    bp.set_default('Observation.proposal.pi', 'CFHT')
    bp.set('Observation.proposal.project', 'get_proposal_project(header)')

    bp.clear('Observation.instrument.name')
    bp.add_fits_attribute('Observation.instrument.name', 'INSTRUME')
    bp.set('Observation.instrument.keywords',
           'get_instrument_keywords(header)')

    bp.set('Observation.target.standard', 'get_target_standard(header)')

    bp.clear('Observation.target_position.coordsys')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'OBJRADEC')
    bp.clear('Observation.target_position.equinox')
    bp.add_fits_attribute('Observation.target_position.equinox', 'OBJEQN')
    bp.add_fits_attribute('Observation.target_position.equinox', 'OBJEQUIN')
    bp.set('Observation.target_position.point.cval1',
           'get_target_position_cval1(header)')
    bp.set('Observation.target_position.point.cval2',
           'get_target_position_cval2(header)')

    bp.set('Observation.telescope.name', 'CFHT 3.6m')
    x, y, z = ac.get_geocentric_location('cfht')
    bp.set('Observation.telescope.geoLocationX', x)
    bp.set('Observation.telescope.geoLocationY', y)
    bp.set('Observation.telescope.geoLocationZ', z)

    bp.set('Plane.dataProductType', 'get_plane_data_product_type(header)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(params)')
    bp.clear('Plane.metaRelease')
    bp.add_fits_attribute('Plane.metaRelease', 'REL_DATE')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE')
    bp.add_fits_attribute('Plane.metaRelease', 'DATE-OBS')
    bp.add_fits_attribute('Plane.metaRelease', 'MET_DATE')
    bp.clear('Plane.dataRelease')
    bp.add_fits_attribute('Plane.dataRelease', 'DATE')
    bp.add_fits_attribute('Plane.dataRelease', 'DATE-OBS')
    bp.add_fits_attribute('Plane.dataRelease', 'REL_DATE')

    bp.set('Plane.provenance.keywords', 'get_provenance_keywords(params)')
    bp.set('Plane.provenance.lastExecuted',
           'get_provenance_last_executed(header)')
    bp.set('Plane.provenance.name', 'get_provenance_name(header)')
    bp.set_default('Plane.provenance.producer', 'CFHT')
    bp.set('Plane.provenance.project', 'get_provenance_project(header)')
    bp.set('Plane.provenance.reference', 'get_provenance_reference(header)')
    bp.clear('Plane.provenance.runID')
    bp.add_fits_attribute('Plane.provenance.runID', 'CRUNID')
    bp.set('Plane.provenance.version', 'get_provenance_version(header)')

    bp.set('Artifact.productType', 'get_product_type(params)')

    # hard-coded values from:
    # - wcaom2archive/cfh2caom2/config/caom2megacam.default and
    # - wxaom2archive/cfht2ccaom2/config/caom2megacam.config
    #
    bp.set('Chunk.energy.axis.axis.ctype', 'get_energy_ctype(header)')
    bp.set('Chunk.energy.axis.axis.cunit', 'get_energy_cunit(header)')
    bp.set('Chunk.energy.axis.error.rnder', 1.0)
    bp.set('Chunk.energy.axis.error.syser', 1.0)
    bp.set('Chunk.energy.axis.function.delta',
           'get_energy_function_delta(params)')
    bp.set('Chunk.energy.axis.function.naxis',
           'get_energy_function_naxis(params)')
    bp.set('Chunk.energy.axis.function.refCoord.pix',
           'get_energy_function_pix(params)')
    bp.set('Chunk.energy.axis.function.refCoord.val',
           'get_energy_function_val(params)')
    bp.set('Chunk.energy.bandpassName', 'get_bandpass_name(header)')
    bp.set('Chunk.energy.resolvingPower', 'get_energy_resolving_power(params)')
    bp.set('Chunk.energy.specsys', 'TOPOCENT')
    bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
    bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')

    bp.set('Chunk.position.axis.axis1.cunit', 'deg')
    bp.set('Chunk.position.axis.axis2.cunit', 'deg')

    if instrument is md.Inst.ESPADONS:
        # constants from caom2espadons.config
        bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN')
        bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN')
        bp.set('Chunk.position.axis.function.dimension.naxis1', 1)
        bp.set('Chunk.position.axis.function.dimension.naxis2', 1)
        bp.set('Chunk.position.axis.function.refCoord.coord1.pix', 1.0)
        bp.clear('Chunk.position.axis.function.refCoord.coord1.val')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.refCoord.coord1.val', 'RA_DEG')
        bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 1.0)
        bp.clear('Chunk.position.axis.function.refCoord.coord2.val')
        bp.add_fits_attribute(
            'Chunk.position.axis.function.refCoord.coord2.val', 'DEC_DEG')
        # CW
        # Fibre size is 1.6", i.e. 0.000444 deg
        bp.set('Chunk.position.axis.function.cd11', -0.000444)
        bp.set('Chunk.position.axis.function.cd12', 0.0)
        bp.set('Chunk.position.axis.function.cd21', 0.0)
        bp.set('Chunk.position.axis.function.cd22', 0.000444)
    bp.set('Chunk.position.axis.error1.rnder', 0.0000278)
    bp.set('Chunk.position.axis.error1.syser', 0.0000278)
    bp.set('Chunk.position.axis.error2.rnder', 0.0000278)
    bp.set('Chunk.position.axis.error2.syser', 0.0000278)

    bp.clear('Chunk.position.coordsys')
    bp.add_fits_attribute('Chunk.position.coordsys', 'RADECSYS')

    bp.set('Chunk.time.exposure', 'get_exptime(params)')
    bp.set('Chunk.time.resolution', 'get_exptime(params)')
    bp.set('Chunk.time.timesys', 'UTC')
    bp.set('Chunk.time.axis.axis.ctype', 'TIME')
    bp.set('Chunk.time.axis.axis.cunit', 'd')
    bp.set('Chunk.time.axis.error.rnder', 0.0000001)
    bp.set('Chunk.time.axis.error.syser', 0.0000001)
    bp.set('Chunk.time.axis.function.naxis', 1)
    cfht_name = CFHTName(ad_uri=uri, instrument=instrument)
    # TODO - this is really really wrong that is_simple is not sufficient
    # to make the distinction between the appropriate implementations.
    if cfht_name.is_simple and not cfht_name.is_master_cal:
        bp.set('Chunk.time.axis.function.delta',
               'get_time_refcoord_delta_simple(params)')
        bp.set('Chunk.time.axis.function.refCoord.val',
               'get_time_refcoord_val_simple(header)')
    else:
        bp.set('Chunk.time.axis.function.delta',
               'get_time_refcoord_delta_composite(params)')
        bp.set('Chunk.time.axis.function.refCoord.val',
               'get_time_refcoord_val_composite(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)

    if cfht_name.has_polarization:
        bp.configure_polarization_axis(6)
        # caom2IngestEspadons.py, l209, lTODO
        bp.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
        bp.set('Chunk.polarization.axis.function.delta', 1)
        bp.set('Chunk.polarization.axis.function.naxis', 1)
        bp.set('Chunk.polarization.axis.function.refCoord.pix', 1)
        bp.set('Chunk.polarization.axis.function.refCoord.val',
               'get_polarization_function_val(header)')

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

    idx = _update_observation_metadata(observation, headers, fqn)
    ccdbin = headers[idx].get('CCBIN1')
    radecsys = headers[idx].get('RADECSYS')
    ctype1 = headers[idx].get('CTYPE1')
    filter_name = headers[idx].get('FILTER')
    instrument = _get_instrument(headers[idx])

    cfht_name = CFHTName(file_name=os.path.basename(fqn), instrument=instrument)
    is_composite, composite_type = _is_composite(
        headers, cfht_name, observation.observation_id)
    if is_composite and not isinstance(observation, CompositeObservation):
        logging.info('{} will be changed to a Composite Observation.'.format(
            observation.observation_id))
        algorithm_name = 'master_detrend'
        if observation.observation_id[-1] == 'p':
            if cfht_name.has_polarization:
                algorithm_name = 'polarization'
            else:
                algorithm_name = 'scan'
        observation = cc.change_to_composite(observation, algorithm_name)

    for plane in observation.planes.values():
        if observation.algorithm.name == 'scan':
            plane.data_product_type = DataProductType.CUBE
            if plane.provenance is not None:
                plane.provenance.last_executed = mc.make_time(
                    headers[idx].get('DATE'))
        for artifact in plane.artifacts.values():
            params = {'uri': artifact.uri,
                      'header': headers[idx]}
            if (get_calibration_level(params) ==
                    CalibrationLevel.RAW_STANDARD):
                time_delta = get_time_refcoord_delta_simple(params)
            else:
                time_delta = get_time_refcoord_delta_composite(params)
            for part in artifact.parts.values():
                for c in part.chunks:
                    chunk_idx = part.chunks.index(c)
                    chunk = part.chunks[chunk_idx]

                    if (time_delta == 0.0 and
                            chunk.time.axis.function.delta == 1.0):
                        # undo the effects of the astropy cdfix call on a
                        # matrix, in fits2caom2.WcsParser, which sets the
                        # diagonal element of the matrix to unity if all
                        # keywords associated with a given axis were
                        # omitted.
                        # See:
                        # https://docs.astropy.org/en/stable/api/astropy.
                        # wc.Wcsprm.html#astropy.wcs.Wcsprm.cdfix
                        chunk.time.axis.function.delta = 0.0

                    if (chunk.energy is not None and
                            chunk.energy.bandpass_name in ['NONE', 'Open']):
                        # CW
                        # If no or "open" filter then set filter name to
                        # null
                        chunk.energy.bandpass_name = None

                    if instrument is md.Inst.ESPADONS:
                        if chunk.time is not None:
                            # consistent with caom2IngestEspadons.py
                            chunk.time_axis = 5
                        _update_energy_espadons(
                            chunk, plane.product_id[-1], headers, idx,
                            artifact.uri, fqn, observation.observation_id)

                        if chunk.position is not None:
                            # CW - Ignore position wcs if a calibration file
                            # suffix list from caom2IngestEspadons.py, l389
                            # 'b', 'd', 'c', 'f', 'x'
                            # with missing spatial indicator keywords
                            if (plane.product_id[-1] in ['a', 'i', 'o', 'p']
                                    and radecsys.lower() != 'null' and
                                    headers[idx].get('RA_DEG') is not None and
                                    headers[idx].get('DEC_DEG') is not None):
                                chunk.position_axis_1 = 3
                                chunk.position_axis_2 = 4
                            else:
                                cc.reset_position(chunk)

                        if plane.product_id[-1] in ['i', 'p']:
                            _update_observable(part, chunk, cfht_name.suffix,
                                               observation.observation_id)

                    elif instrument is md.Inst.MEGAPRIME:
                        # CW
                        # Ignore position wcs if a calibration file (except 'x'
                        # type calibration) and/or position info not in header
                        # or binned 8x8
                        if (plane.product_id[-1] in ['b', 'l', 'd', 'f'] or
                                ccdbin == 8 or radecsys is None or
                                ctype1 is None):
                            cc.reset_position(chunk)

                        # CW
                        # Ignore energy wcs if some type of calibration file
                        # or filter='None' or 'Open' or cdelt4=0.0 (no filter
                        # match)
                        cdelt_value = -1.0
                        if (chunk.energy is not None and
                                chunk.energy.axis is not None and
                                chunk.energy.axis.function is not None):
                            cdelt_value = chunk.energy.axis.function.delta
                        if (filter_name is None or filter_name == 'Open' or
                                cdelt_value == 0.0 or
                                plane.product_id[-1] in ['b', 'l', 'd']):
                            cc.reset_energy(chunk)

                    elif instrument is md.Inst.SITELLE:
                        if chunk.energy_axis is not None:
                            # CW, SF 18-12-19 - SITELLE biases and darks have
                            # no energy, all other observation types do
                            if observation.type in ['BIAS', 'DARK']:
                                cc.reset_energy(chunk)
                            else:
                                chunk.energy_axis = 3
                        if chunk.time_axis is not None:
                            chunk.time_axis = 4

                        if (plane.product_id[-1] in ['a', 'o', 'x'] and
                                chunk.position is None):
                            _update_position_sitelle(chunk, headers[idx],
                                                     observation.observation_id)

                        if (plane.product_id[-1] == 'p' and
                                chunk.position.axis.function is None):
                            _update_position_function_sitelle(
                                chunk, headers[idx],
                                observation.observation_id, idx)

                    elif instrument is md.Inst.WIRCAM:
                        if (chunk.time is not None and
                                observation.type not in
                                ['BPM', 'DARK', 'FLAT', 'WEIGHT']):
                            _update_wircam_time(
                                chunk, headers[idx],
                                observation.observation_id)

                        if (plane.product_id[-1] in ['f'] or
                                observation.type in
                                ['BPM', 'DARK', 'FLAT', 'WEIGHT']):
                            cc.reset_position(chunk)

        if isinstance(observation, CompositeObservation):
            if composite_type == 'IMCMB':
                cc.update_plane_provenance(plane, headers[1:],
                                           composite_type, COLLECTION,
                                           _repair_imcmb_provenance_value,
                                           observation.observation_id)
            elif composite_type == 'COMMON':
                cc.update_plane_provenance_single(
                    plane, headers, composite_type, COLLECTION,
                    _repair_comment_provenance_value,
                    observation.observation_id)
            else:
                cc.update_plane_provenance(plane, headers, composite_type,
                                           COLLECTION,
                                           _repair_filename_provenance_value,
                                           observation.observation_id)
        if instrument is md.Inst.WIRCAM and plane.product_id[-1] == 'p':
            # caom2IngestWircam.py, l193
            # TODO change this from the existing behaviour
            _update_plane_provenance_p(plane, observation.observation_id, 'og')
        elif instrument is md.Inst.ESPADONS and plane.product_id[-1] == 'i':
            # caom2IngestEspadons.py, l714
            _update_plane_provenance_p(plane, observation.observation_id, 'o')

    # relies on update_plane_provenance being called
    if isinstance(observation, CompositeObservation):
        cc.update_observation_members(observation)

    logging.debug('Done update.')
    return observation


def _is_composite(headers, cfht_name, obs_id):
    result = False
    composite_type = ''
    if cc.is_composite(headers):
        result = True
        composite_type = 'IMCMB'
    else:
        file_type = headers[0].get('FILETYPE')
        if file_type is not None and 'alibrat' in file_type:
            logging.info(
                f'Treating {obs_id} with filetype {file_type} as composite. ')
            result = True
            composite_type = 'COMMENT'
    if not result and not cfht_name.is_simple:
        logging.error(f'changing is composite for {cfht_name.suffix} {cfht_name.file_name}')
        result = True
        composite_type = 'FILENAM'
    return result, composite_type


def _update_energy_espadons(chunk, suffix, headers, idx, uri, fqn, obs_id):
    logging.debug(f'Begin _update_energy_espadons for {obs_id}')
    cfht_name = CFHTName(file_name=os.path.basename(fqn),
                         instrument=md.Inst.ESPADONS)
    if cfht_name.suffix == suffix:
        axis = Axis('WAVE', 'nm')
        params = {'header': headers[idx],
                  'uri': uri}
        resolving_power = get_energy_resolving_power(params)
        if suffix in ['a', 'c', 'f', 'o', 'x']:
            naxis1 = get_energy_function_naxis(params)
            cdelt1 = get_energy_function_delta(params)
            crval1 = get_energy_function_val(params)
            ref_coord_1 = RefCoord(0.5, crval1)
            ref_coord_2 = RefCoord(1.5, crval1 + float(naxis1)*cdelt1)
            coord_range = CoordRange1D(ref_coord_1, ref_coord_2)
            coord_axis = CoordAxis1D(axis=axis, range=coord_range)
        elif suffix in ['i', 'p']:
            # PD slack 08-01-20
            # espadons is a special case because using bounds allows one to
            # define "tiles" and then the SODA cutout service can extract the
            # subset of tiles that overlap the desired region. That's the best
            # that can be done because it is not possible to create a
            # CoordFunction1D to say what the wavelength of each pixel is
            #
            # If the coverage had significant gaps (eg SCUBA or SCUBA2 from
            # JCMT)  then the extra detail in bounds would enable better
            # discovery (as the gaps would be captured in the plane metadata).
            # In the case of espadons I don't think the gaps  are significant
            # (iirc, espadons is an eschelle spectrograph but I don't recall
            # whether the discontinuity between eschelle was a small gap or an
            # overlap)
            #
            # So: bounds provides more detail and it can in principle improve
            # data discovery (if gaps) and enable extraction of subsections of
            # the spectrum via the SODA service. Espadons was one of the use
            # cases that justified having bounds there

            # SF slack 08-01-20
            # We need the information that is contained in bounds. Gaps need
            # to be captured. So keep bounds. If you decide to remove range,
            # then advanced users would have to dig in the info to understand
            # range is first and last bounds.

            # read in the complete fits file, including the data
            logging.info(f'Reading data from {fqn}.')
            hdus = ac.read_fits_data(fqn)
            wave = hdus[idx].data[0, :]
            coord_bounds = ac.build_chunk_energy_bounds(wave, axis)
            coord_axis = CoordAxis1D(axis=axis, bounds=coord_bounds)
            hdus.close()
        chunk.energy = SpectralWCS(coord_axis,
                                   specsys='TOPOCENT',
                                   ssysobs='TOPOCENT',
                                   ssyssrc='TOPOCENT',
                                   resolving_power=resolving_power)
        chunk.energy_axis = 1
    logging.debug(f'End _update_energy_espadons for {obs_id}')


def _update_observable(part, chunk, suffix, obs_id):
    logging.debug(f'Begin _update_observable for {obs_id}')
    # caom2IngestEspadons.py, l828
    # CW Set up observable axes, row 1 is wavelength, row 2 is normalized
    # flux, row 3 ('p' only) is Stokes spectrum

    # this check is here, because it's quite difficult to find the 'right'
    # chunk, and blind updating generally causes both chunks to have the
    # same metadata values.
    if chunk.observable is None:
        chunk.observable_axis = 2
        independent_axis = Axis('WAVE', 'nm')
        independent = Slice(independent_axis, 1)
        dependent_axis = Axis('flux', 'counts')
        dependent = Slice(dependent_axis, 2)
        chunk.observable = ObservableAxis(dependent, independent)

        if suffix == 'p' and len(part.chunks) == 1:
            logging.error('yes or no?')
            # caom2IngestEspadons.py, l863
            dependent_axis = Axis('polarized flux', 'percent')
            dependent = Slice(dependent_axis, 3)
            new_chunk = copy.deepcopy(chunk)
            new_chunk.observable = ObservableAxis(dependent, independent)
            new_chunk._id = Chunk._gen_id()
            part.chunks.append(new_chunk)
    logging.debug(f'End _update_observable for {obs_id}')


def _update_observation_metadata(obs, headers, fqn):
    """
    Why this method exists:

    There are CFHT files that have almost no metadata in the primary HDU, but
    all the needed metadata in subsequent HDUs.

    It's not possible to apply extension
    numbers for non-chunk blueprint entries, so that means that to use the
    information captured in the blueprint, the header that's provided
    must be manipulated instead. There is only access to the header
    information in this extension of the fitscaom2 module (i.e. this file)
    during the execution of the 'update' step of fits2caom2 execution.
    Hence the self-referential implementation. Maybe it will result in an
    infinite loop and I'll be sad.
    """

    # notes to myself
    # - can only set blueprint extension stuff for chunk bits
    # - need to do the extension/header replacement here for the
    #   Observation, but the chunk level blueprint has been hosed by the
    #   previous settings
    # - the correct values are set in extension 1 for the chunk
    #   and need to be accessible as extension 0 for the observation, plane,
    #   artifact
    # - so then - do I need the blueprint before it's been modified?
    #   maybe instead of storing the blueprint I just call accumulatee_bp?

    # and I'm back trying to figure out how to undo the doing ...
    # the first time through, the CD4_4 value is set from CDELT4,
    # the second time through, it is not, because

    # check for files with primary headers that have NO information
    # - e.g. 2445848a
    idx = 0
    run_id = headers[0].get('RUNID')
    if run_id is None:
        idx = 1

        logging.warning(f'Resetting the header/blueprint '
                        f'relationship for {obs.observation_id}')

        # use the fqn to define the URI
        # TODO - leaking name structure here
        product_id = CFHTName.remove_extensions(os.path.basename(fqn))
        extension = '.fz'
        instrument = _get_instrument(headers[idx])
        if instrument is md.Inst.ESPADONS:
            extension = '.gz'
        uri = mc.build_uri(ARCHIVE,
                           os.path.basename(fqn).replace('.header',
                                                         extension))
        module = importlib.import_module(__name__)
        bp = ObsBlueprint(module=module)
        accumulate_bp(bp, uri, instrument)

        # re-read the headers from disk, because the first pass through
        # caom2gen will have modified the header content based on the
        # original blueprint  # TODO this is not a long-term implementation
        # the execution goes through and sets the CDELT4 value, then from
        # that sets the CD4_4 value. Then the second time through, the
        # CD4_4 value update is expressly NOT done, because the CD4_4 value
        # is already set from the first pass-through - need to figure out a
        # way to fix this .... sigh
        unmodified_headers = get_cadc_headers(f'file://{fqn}')

        # make a list of length 1, because then
        # fits2caom2.FitsParser.augment_observation won't try to read the
        # headers from a file on disk
        # print('{!r}'.format(headers[1]))
        # parser = FitsParser([headers[1]], bp, uri)
        parser = FitsParser([unmodified_headers[1]], bp, uri)
        parser.augment_observation(obs, uri, product_id)

        # re-home the chunk information to be consistent with accepted CAOM
        # patterns of part/chunk relationship - i.e. part 0 never has chunks
        for plane in obs.planes.values():
            for artifact in plane.artifacts.values():
                if artifact.uri == uri:
                    part0 = artifact.parts['0']
                    part1 = artifact.parts['1']
                    part1.chunks[0] = part0.chunks[0]
                    part0.chunks = TypedList(Chunk, )

    return idx


def _update_plane_provenance_p(plane, obs_id, suffix):
    logging.debug(f'Begin _update_plane_provenance_p for {obs_id}')
    obs_member_str = mc.CaomName.make_obs_uri_from_obs_id(COLLECTION,
                                                          obs_id)
    obs_member = ObservationURI(obs_member_str)
    plane_uri = PlaneURI.get_plane_uri(obs_member, f'{obs_id}{suffix}')
    plane.provenance.inputs.add(plane_uri)
    logging.debug(f'End _update_plane_provenance_p for {obs_id}')


def _update_position_sitelle(chunk, header, obs_id):
    logging.debug(f'Begin _update_position_sitelle for {obs_id}')
    # from caom2IngestSitelle.py l894
    obs_ra = header.get('RA_DEG')
    obs_dec = header.get('DEC_DEG')
    if obs_ra is None or obs_dec is None:
        logging.error(
            'RA_DEG {obs_ra} DEC_DEG {obs_dec} for {obs_id} are not set.')
        return

    pix_scale2 = mc.to_float(header.get('PIXSCAL2'))
    delta_dec = (1024.0 * pix_scale2) / 3600.0
    obs_dec_bl = obs_dec - delta_dec
    obs_dec_tr = obs_dec + delta_dec

    pix_scale1 = mc.to_float(header.get('PIXSCAL1'))
    # obs_dec_tr == obs_dec_tl
    delta_ra_top = (1024.0 * pix_scale1) / (
                3600.0 * (math.cos(obs_dec_tr * 3.14159 / 180.0)))
    delta_ra_bot = (1024.0 * pix_scale1) / (
                3600.0 * (math.cos(obs_dec_bl * 3.14159 / 180.0)))
    obs_ra_bl = obs_ra + delta_ra_bot
    obs_ra_tr = obs_ra - delta_ra_top

    axis = CoordAxis2D(Axis('RA', 'deg'), Axis('DEC', 'deg'))
    axis.range = CoordRange2D(Coord2D(RefCoord(0.5, obs_ra_bl),
                                      RefCoord(0.5, obs_dec_bl)),
                              Coord2D(RefCoord(2048.5, obs_ra_tr),
                                      RefCoord(2048.5, obs_dec_tr)))
    position = SpatialWCS(axis)
    position.coordsys = 'FK5'
    position.equinox = 2000.0
    chunk.position = position
    logging.debug(f'End _update_position_sitelle for {obs_id}')


def _update_position_function_sitelle(chunk, header, obs_id, extension):
    logging.debug(f'Begin _update_position_function_sitelle for {obs_id}')
    cd1_1 = header.get('CD1_1')
    # caom2IngestSitelle.py, l745
    # CW
    # Be able to handle any of the 3 wcs systems used
    if cd1_1 is None:
        pc1_1 = header.get('PC1_1')
        if pc1_1 is not None:
            cdelt1 = mc.to_float(header.get('CDELT1'))
            if cdelt1 is None:
                cd1_1 = mc.to_float(header.get('PC1_1'))
                cd1_2 = mc.to_float(header.get('PC1_2'))
                cd2_1 = mc.to_float(header.get('PC2_1'))
                cd2_2 = mc.to_float(header.get('PC2_2'))
            else:
                cdelt2 = mc.to_float(header.get('CDELT2'))
                cd1_1 = cdelt1 * mc.to_float(header.get('PC1_1'))
                cd1_2 = cdelt1 * mc.to_float(header.get('PC1_2'))
                cd2_1 = cdelt2 * mc.to_float(header.get('PC2_1'))
                cd2_2 = cdelt2 * mc.to_float(header.get('PC2_2'))
            header['CD1_1'] = cd1_1
            header['CD1_2'] = cd1_2
            header['CD2_1'] = cd2_1
            header['CD2_2'] = cd2_2

    wcs_parser = WcsParser(header, obs_id, extension)
    if chunk is None:
        chunk = Chunk()
    wcs_parser.augment_position(chunk)
    logging.debug(f'End _update_position_function_sitelle for {obs_id}')


def _update_wircam_time(chunk, header, obs_id):
    logging.debug(f'Begin _update_wircam_time for {obs_id}')
    n_exp = header.get('NEXP')
    if (n_exp is not None and chunk.time.axis is not None and
            chunk.time.axis.function is not None):
        # caom2IngestWircam.py, l843
        chunk.time.axis.function.naxis = mc.to_int(n_exp)
    logging.debug(f'End _update_wircam_time for {obs_id}')


def _identify_instrument(uri):
    # TODO - make this a header read when I get to it ...
    result = md.Inst.SITELLE  # 1944968p, 2445397p
    if '2463796o' in uri:
        result = md.Inst.MEGACAM
    elif '2281792p' in uri or '2157095o' in uri or 'weight' in uri:
        result = md.Inst.WIRCAM
    elif ('1001063b' in uri or '1001836x' in uri or '1003681' in uri or
            '1219059' in uri or '1883829c' in uri or '2460602a' in uri or
            '760296f' in uri or '881162d' in uri or '979339' in uri or
            '2460503p' in uri):
        result = md.Inst.ESPADONS
    return result


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
        instrument = _identify_instrument(uri)
        blueprint = ObsBlueprint(module=module)
        if not mc.StorageName.is_preview(uri):
            accumulate_bp(blueprint, uri, instrument)
        blueprints[uri] = blueprint
    return blueprints


def _get_uris(args):
    result = []
    if args.lineage:
        for ii in args.lineage:
            result.append(ii.split('/', 1)[1])
    elif args.local:
        for ii in args.local:
            # TODO hack that leaks naming format - sigh ... :(
            file_uri = mc.build_uri(ARCHIVE, os.path.basename(ii))
            result.append(file_uri)
    else:
        raise mc.CadcException(
            'Could not define uri from these args {}'.format(args))
    return result


def _repair_comment_provenance_value(value, obs_id):
    results = []
    # COMMENT headers with provenance:
    # COMMENT Scan member=2445653o st=174 iq=1.2200 bk=5.5214 ex=0.024000 ...
    # COMMENT Flat member=2445211f
    # COMMENT Standard member=2445849o
    if 'member=' in str(value):
        for entry in value:
            if 'Scan member' in entry:
                temp = str(entry).split('member=')
                prov_prod_id = temp[1].split()[0]
                prov_obs_id = CFHTName(file_name=prov_prod_id).obs_id
                # 0 - observation
                # 1 - plane
                results.append([prov_obs_id, prov_prod_id])
    return results


def _repair_filename_provenance_value(value, obs_id):
    # values require no repairing, because they look like:
    # FILENAME= '2460503p'
    # FILENAM1= '2460503o'           / Base filename at acquisition
    # FILENAM2= '2460504o'           / Base filename at acquisition
    # FILENAM3= '2460505o'           / Base filename at acquisition
    # FILENAM4= '2460506o'           / Base filename at acquisition
    prov_prod_id = None
    prov_obs_id = None
    if value != obs_id:
        prov_prod_id = value
        prov_obs_id = value[:-1]
    return prov_obs_id, prov_prod_id


def _repair_imcmb_provenance_value(value, obs_id):
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
        logging.debug('Done {} processing.'.format(APPLICATION))
        sys.exit(result)
    except Exception as e:
        logging.error('Failed {} execution for {}.'.format(APPLICATION, args))
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
