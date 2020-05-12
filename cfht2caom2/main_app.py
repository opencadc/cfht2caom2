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

For Espadons, Sitelle and Spirou there are derived observations made up from
several exposures. These derived observations are given the ‘p’ ending to
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
    - SITELLE - 'p' is processed + Derived - a different observation
    - SPIROU/Espadons - 'p' is polarized + Derived (TBC on Derived), a
      different observation
    - there are other processed files for single exposures
    - users want 'p' files

- conclusion - one plane / file, because users want to see one row / file in
  the query results

  JJK - slack - 01-04-20
  CFHT files are independent, in that, for example, a user does not require a
  'p' file to understand the content of a 'b' file. Given this independence, it
  is acceptable to map one plane / file.

  SF - slack - 02-04-20
  - MegaCam - the logic should be probably be 2 planes: p and o for science.
            - all cfht exposures are sorted by EXPNUM if i understand their
            data acquisition. b,f,d,x should be 1 plane observations.
            - my assumption is that the b,f,d,x have no reason to have a
            processed equivalent.

  Typical science with the actual data - more or less what Elixir does:
  - do something like: p = (o - <b>) / ( <f - <b> > )
  - where <x> is the average of a bunch of x frames


CFHT Energy:

slack - 08-01-20
- PD - Well, espadons is a special case because using bounds allows one to
define "tiles" and then the SODA cutout service can extract the subset of
tiles that overlap the desired region. That's the best that can be done
because it is not possible to create a CoordFunction1D to say what the
wavelength of each pixel is
- PD - Bounds provides more detail than range and enables a crude tile-based
cutout operation later. If the covergage had significant gaps (eg SCUBA or
SCUBA2 from JCMT) then the extra detail in bounds would enable better
discovery (as the gaps would be captured in the plane metadata). In the case
of espadons I don't think the gaps are significant (iirc, espadons is an
eschelle spetrograph but I don't recall whether the discontinuity between
eschelle was a small gap or an overlap)
- PD - So: bounds provides more detail and it can in principle improve data
discovery (if gaps) and enable extraction of subsections of the spectrum via
the SODA service. Espadons was one of the use cases that justified having
bounds there
- SF - ok now that is clearer, i think we need the information that is
contained in bounds, gaps need to be captured. so keep bounds. if you decide
to remove range, then advanced users would have to dig in the info to
understand range is first and last bounds if i understand correctly.


CFHT WCS:
- CW - 28-04-20
- These raw SITELLE observations have a weird data format where 2048x2048
  images from the two interferometer arms (which are images of the same single
  field) are stitched together into a single 2048x4100ish image. So there is
  no contiguous wcs describing the raw data and in any case nobody wants to
  cut out of it.

"""

import copy
import importlib
import logging
import math
import os
import sys
import traceback

from enum import Enum

from caom2 import Observation, CalibrationLevel, ObservationIntentType
from caom2 import ProductType, DerivedObservation, TypedList, Chunk
from caom2 import CoordRange2D, CoordAxis2D, Axis, Coord2D, RefCoord
from caom2 import SpatialWCS, DataProductType, ObservationURI, PlaneURI
from caom2 import CoordAxis1D, CoordRange1D, SpectralWCS, Slice
from caom2 import ObservableAxis, CoordError
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2utils import FitsParser, WcsParser, get_cadc_headers
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from caom2pipe import translate_composable as tc

from cfht2caom2 import cfht_builder as cb
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import metadata as md


__all__ = ['cfht_main_app', 'to_caom2', 'update', 'APPLICATION']


APPLICATION = 'cfht2caom2'

# All comments denoted 'CW' are copied from
# ssh://goliaths@gimli3/srv/cadc/git/wcaom2archive.git
# from /cfht2caom2/scripts


class ProvenanceType(Enum):
    """The different types of header values that identify provenance
    information. Used to specify different functions for
    execution. """
    COMMENT = 'COMMENT'
    FILENAME = 'FILENAM'
    IMCMB = 'IMCMB'
    UNDEFINED = 'UNDEFINED'  # hdf5 file support


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
    cfht_name = cn.CFHTName(ad_uri=uri, instrument=instrument)
    bp.configure_position_axes((1, 2))
    bp.configure_time_axis(3)
    if cfht_name.has_energy:
        bp.configure_energy_axis(4)
    bp.configure_observable_axis(6)

    bp.set('Observation.intent', 'get_obs_intent(header)')

    meta_producer = mc.get_version(APPLICATION)
    bp.set('Observation.metaProducer', meta_producer)
    bp.set('Observation.metaRelease', 'get_meta_release(header)')

    bp.set('Observation.sequenceNumber', 'get_obs_sequence_number(params)')
    bp.set('Observation.type', 'get_obs_type(header)')

    bp.set('Observation.environment.elevation',
           'get_environment_elevation(header)')
    bp.set('Observation.environment.humidity',
           'get_obs_environment_humidity(header)')

    # TODO title is select title from runid_title where proposal_id = 'runid'
    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'RUNID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PI_NAME')
    bp.set_default('Observation.proposal.pi', 'CFHT')
    bp.set('Observation.proposal.project', 'get_proposal_project(header)')
    bp.set('Observation.proposal.title', 'get_proposal_title(header)')

    bp.set('Observation.instrument.name', instrument.value)
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
    bp.set('Plane.dataRelease', 'get_plane_data_release(header)')
    bp.set('Plane.metaRelease', 'get_meta_release(header)')

    bp.set('Plane.provenance.lastExecuted',
           'get_provenance_last_executed(header)')
    bp.set('Plane.metaProducer', meta_producer)
    bp.set_default('Plane.provenance.producer', 'CFHT')
    bp.set('Plane.provenance.project', 'STANDARD PIPELINE')
    bp.clear('Plane.provenance.runID')
    bp.add_fits_attribute('Plane.provenance.runID', 'CRUNID')
    bp.set('Plane.provenance.version', 'get_provenance_version(header)')

    bp.set('Artifact.metaProducer', meta_producer)
    bp.set('Artifact.productType', 'get_product_type(params)')
    bp.set('Artifact.releaseType', 'data')

    bp.set('Chunk.metaProducer', meta_producer)
    # hard-coded values from:
    # - wcaom2archive/cfh2caom2/config/caom2megacam.default and
    # - wxaom2archive/cfht2ccaom2/config/caom2megacam.config
    #
    # Gemini is all range, make Mega* range too, and WIRCam too
    if instrument not in [md.Inst.MEGACAM, md.Inst.MEGAPRIME,
                          md.Inst.WIRCAM]:
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
        bp.clear('Chunk.energy.bandpassName')
        bp.add_fits_attribute('Chunk.energy.bandpassName', 'FILTER')
        bp.set('Chunk.energy.resolvingPower',
               'get_energy_resolving_power(params)')
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

    if cfht_name.suffix != 'g':
        bp.set('Chunk.time.exposure', 'get_exptime(params)')
        bp.set('Chunk.time.resolution', 'get_exptime(params)')
        bp.set('Chunk.time.timesys', 'UTC')
        bp.set('Chunk.time.axis.axis.ctype', 'TIME')
        bp.set('Chunk.time.axis.axis.cunit', 'd')
        bp.set('Chunk.time.axis.error.rnder', 0.0000001)
        bp.set('Chunk.time.axis.error.syser', 0.0000001)
        bp.set('Chunk.time.axis.function.naxis', 1)

    # TODO - remove this declaration
    cfht_name = cn.CFHTName(ad_uri=uri, instrument=instrument)
    # TODO - this is really really wrong that is_simple is not sufficient
    # to make the distinction between the appropriate implementations.
    if cfht_name.is_simple and not cfht_name.is_master_cal:
        bp.set('Chunk.time.axis.function.delta',
               'get_time_refcoord_delta_simple(params)')
        bp.set('Chunk.time.axis.function.refCoord.val',
               'get_time_refcoord_val_simple(header)')
    else:
        bp.set('Chunk.time.axis.function.delta',
               'get_time_refcoord_delta_derived(header)')
        bp.set('Chunk.time.axis.function.refCoord.val',
               'get_time_refcoord_val_derived(header)')
    bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)

    if instrument is md.Inst.ESPADONS:
        _accumulate_espadons_bp(bp, cfht_name)
    elif instrument is md.Inst.MEGACAM or instrument is md.Inst.MEGAPRIME:
        _accumulate_mega_bp(bp, uri, cfht_name)
    elif instrument is md.Inst.SITELLE:
        _accumulate_sitelle_bp(bp, uri, cfht_name)
    elif instrument is md.Inst.SPIROU:
        _accumulate_spirou_bp(bp, uri, cfht_name)
    elif instrument is md.Inst.WIRCAM:
        _accumulate_wircam_bp(bp, uri, cfht_name)

    logging.debug('Done accumulate_bp.')


def _accumulate_espadons_bp(bp, cfht_name):
    """Configure the ESPaDOnS-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_espadons_bp.')

    bp.add_fits_attribute('Observation.target_position.coordsys', 'RADECSYS')

    bp.set('Plane.provenance.keywords',
           'get_espadons_provenance_keywords(params)')
    bp.set('Plane.provenance.lastExecuted',
           'get_espadons_provenance_last_executed(header)')
    bp.set('Plane.provenance.name', 'get_espadons_provenance_name(header)')
    bp.set('Plane.provenance.project',
           'get_espadons_provenance_project(header)')
    bp.set_default('Plane.provenance.reference',
                   'http://www.cfht.hawaii.edu/Instruments/Spectroscopy/'
                   'Espadons/')
    bp.set('Plane.provenance.version',
           'get_espadons_provenance_version(header)')

    bp.set('Chunk.energy.resolvingPower',
           'get_espadons_energy_resolving_power(params)')

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

    bp.set('Chunk.time.axis.function.delta',
           'get_espadons_time_refcoord_delta(params)')
    bp.set('Chunk.time.axis.function.refCoord.val',
           'get_espadons_time_refcoord_val(params)')

    bp.set('Chunk.time.exposure', 'get_espadons_exptime(params)')
    bp.set('Chunk.time.resolution', 'get_espadons_exptime(params)')

    if cfht_name.suffix == 'p':
        bp.configure_polarization_axis(6)
        # caom2IngestEspadons.py, l209, lTODO
        bp.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
        bp.set('Chunk.polarization.axis.function.delta', 1)
        bp.set('Chunk.polarization.axis.function.naxis', 1)
        bp.set('Chunk.polarization.axis.function.refCoord.pix', 1)
        bp.set('Chunk.polarization.axis.function.refCoord.val',
               'get_polarization_function_val(header)')

    logging.debug('Done _accumulate_espadons_bp.')


def _accumulate_mega_bp(bp, uri, cfht_name):
    """Configure the MegaCam/MegaPrime-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_mega_bp.')

    bp.set('Plane.provenance.lastExecuted',
           'get_mega_provenance_last_executed(header)')
    bp.set_default('Plane.provenance.name', 'ELIXIR')
    bp.set_default('Plane.provenance.reference',
                   'http://www.cfht.hawaii.edu/Instruments/Elixir/')

    logging.debug('Done _accumulate_mega_bp.')


def _accumulate_sitelle_bp(bp, uri, cfht_name):
    """Configure the Sitelle-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_sitelle_bp.')
    bp.set('Plane.dataProductType', 'get_sitelle_plane_data_product_type(uri)')
    bp.set_default('Plane.provenance.name', 'ORBS')
    bp.set_default('Plane.provenance.reference', 'http://ascl.net/1409.007')

    bp.set('Chunk.energy.resolvingPower',
           'get_sitelle_energy_resolving_power(params)')

    bp.set('Chunk.time.axis.function.delta',
           'get_sitelle_time_refcoord_delta(params)')
    bp.set('Chunk.time.axis.function.refCoord.val',
           '_get_mjd_start(header)')

    logging.debug('End _accumulate_sitelle_bp.')


def _accumulate_spirou_bp(bp, uri, cfht_name):
    """Configure the SPIRou-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    if cfht_name.suffix == 'r':
        pass
    elif cfht_name.suffix in ['a', 'c', 'd', 'f', 'o', 'x']:
        bp.set('Plane.provenance.name', 'get_spirou_provenance_name(header)')
        bp.set('Plane.provenance.reference',
               'http://www.cfht.hawaii.edu/Instruments/SPIRou/')
        bp.set('Plane.provenance.version',
               'get_spirou_provenance_version(header)')
    else:
        bp.set('Plane.provenance.name', 'DRS')
        bp.set('Plane.provenance.reference',
               'https://www.cfht.hawaii.edu/Instruments/SPIRou/'
               'SPIRou_pipeline.php')
        bp.clear('Plane.provenance.version')
        bp.add_fits_attribute('Plane.provenance.version', 'VERSION')

    # from caom2IngestSpirou.py, l654
    # CW - Do energy stuff for raw and calibrated data
    # Initial numbers to be updated once have some reduced spectra
    # This was based on 50 orders of 4000 pix each covering 0.98 to 2.35
    # microns. New numbers based on reduced spectra 0.955 to 2.515 microns
    # TODO - should update old ones some time
    # naxis6 = 1
    # crpix6 = 0.5
    # crval6 = 955.0
    # cdelt6 = 1560.0
    #
    # other values from caom2spirou.default
    bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
    bp.set('Chunk.energy.axis.axis.cunit', 'nm')
    bp.set('Chunk.energy.axis.function.delta', 1560.0)
    bp.set('Chunk.energy.axis.function.naxis', 1)
    bp.set('Chunk.energy.axis.function.refCoord.pix', 0.5)
    bp.set('Chunk.energy.axis.function.refCoord.val', 955.0)
    bp.set('Chunk.energy.axis.error.rnder', 0.001)
    bp.set('Chunk.energy.axis.error.syser', 0.001)
    bp.set('Chunk.energy.resolvingPower', 73000.0)

    # values from caom2spirou.default
    bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN')
    bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN')
    bp.set('Chunk.position.axis.function.dimension.naxis1', 1)
    bp.set('Chunk.position.axis.function.dimension.naxis2', 1)
    bp.set('Chunk.position.axis.function.refCoord.coord1.pix', 1.0)
    bp.set('Chunk.position.axis.function.refCoord.coord1.val',
           'get_ra_deg_from_0th_header(header)')
    bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 1.0)
    bp.set('Chunk.position.axis.function.refCoord.coord2.val',
           'get_dec_deg_from_0th_header(header)')
    bp.set('Chunk.position.axis.function.cd11', -0.00035833)
    bp.set('Chunk.position.axis.function.cd12', 0.0)
    bp.set('Chunk.position.axis.function.cd21', 0.0)
    bp.set('Chunk.position.axis.function.cd22', 0.00035833)

    bp.set('Chunk.position.coordsys',
           'get_position_coordsys_from_0th_header(header)')
    bp.set('Chunk.position.equinox',
           'get_position_equinox_from_0th_header(header)')

    bp.set('Chunk.time.axis.function.delta',
           'get_spirou_time_refcoord_delta(params)')
    bp.set('Chunk.time.axis.function.naxis',
           'get_spirou_time_refcoord_naxis(params)')
    bp.set('Chunk.time.exposure', 'get_spirou_exptime(params)')
    bp.set('Chunk.time.resolution', 'get_spirou_resolution(params)')


def _accumulate_wircam_bp(bp, uri, cfht_name):
    """Configure the WIRCam-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_wircam_bp.')
    bp.set('Observation.type', 'get_wircam_obs_type(params)')

    bp.set('Plane.provenance.keywords', 'get_wircam_provenance_keywords(uri)')
    bp.set_default('Plane.provenance.name', 'IIWI')
    bp.set_default('Plane.provenance.reference',
                   'http://www.cfht.hawaii.edu/Instruments/Imaging/WIRCam')

    bp.set('Chunk.energy.bandpassName', 'get_wircam_bandpass_name(header)')

    logging.debug('Done _accumulate_wircam_bp.')


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
    uri = None
    if 'uri' in kwargs:
        uri = kwargs['uri']

    ingesting_hdf5 = False

    # TODO - use the instrument name in the observation? might not work
    # for a redo
    if uri is None:
        instrument = md.Inst(observation.instrument.name)
    else:
        suffix = cn.CFHTName.remove_extensions(uri)[-1]
        if suffix == 'z':
            instrument = md.Inst.SITELLE
            ingesting_hdf5 = True
            logging.info(
                f'Ingesting the hdf5 plane for {observation.observation_id}')
        else:
            instrument = cb.CFHTBuilder.get_instrument(headers, uri)

    if instrument is md.Inst.MEGACAM:
        # need the 'megacam' for the filter lookup at SVO, but there is only
        # a 'MegaPrime' instrument in the CAOM collection at CADC
        # see e.g. 2003A.frpts.z.36.00
        observation.instrument = cc.copy_instrument(observation.instrument,
                                                    md.Inst.MEGAPRIME.value)

    # fqn is not always defined, ffs
    if uri is not None:
        ignore_scheme, ignore_archive, f_name = mc.decompose_uri(uri)
        cfht_name = cn.CFHTName(
            file_name=f_name, instrument=instrument)
    elif fqn is not None:
        cfht_name = cn.CFHTName(
            file_name=os.path.basename(fqn), instrument=instrument)
    else:
        raise mc.CadcException(f'Cannot define a CFHTName instance for '
                               f'{observation.observation_id}')

    if ingesting_hdf5:
        # avoid all the code that references undefined headers variable
        if not isinstance(observation, DerivedObservation):
            observation = cc.change_to_composite(
                observation, 'scan', cn.COLLECTION)
        _update_sitelle_plane(observation, uri)
        logging.debug('Done hdf5 update.')
        return observation

    is_derived, derived_type = _is_derived(
        headers, cfht_name, observation.observation_id)
    if is_derived and not isinstance(observation, DerivedObservation):
        logging.info(f'{observation.observation_id} will be changed to a '
                     f'Derived Observation.')
        algorithm_name = 'master_detrend'
        if observation.observation_id[-1] == 'p':
            if cfht_name.has_polarization:
                algorithm_name = 'polarization'
            else:
                algorithm_name = 'scan'
        observation = cc.change_to_composite(
            observation, algorithm_name, cn.COLLECTION)

    idx = _update_observation_metadata(observation, headers, cfht_name,
                                       fqn, uri)
    ccdbin = headers[idx].get('CCBIN1')
    radecsys = headers[idx].get('RADECSYS')
    ctype1 = headers[idx].get('CTYPE1')
    filter_name = headers[idx].get('FILTER')
    if filter_name is None and len(headers) > idx + 1:
        filter_name = headers[idx+1].get('FILTER')

    for plane in observation.planes.values():
        if plane.product_id != cfht_name.product_id:
            # do only the work for the applicable plane
            continue

        if instrument is md.Inst.SPIROU and not cfht_name.suffix == 'r':
            cc.rename_parts(observation, headers)

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
                time_delta = get_time_refcoord_delta_derived(headers[idx])
            for part in artifact.parts.values():

                for c in part.chunks:
                    chunk_idx = part.chunks.index(c)
                    chunk = part.chunks[chunk_idx]

                    cc.undo_astropy_cdfix_call(chunk, time_delta)

                    if chunk.energy is not None:
                        if chunk.energy.bandpass_name in ['NONE', 'Open']:
                            # CW
                            # If no or "open" filter then set filter name to
                            # null
                            chunk.energy.bandpass_name = None

                        if (chunk.energy.axis is not None and
                                chunk.energy.axis.axis is not None and
                                chunk.energy.axis.axis.ctype is not None):
                            # PD 08-04-20
                            # the correct way to express "inverse meter" is
                            # either  m**-1 or m^-1
                            #
                            # we support both exponentiations but convert ^
                            # to ** so I guess at that time we thought ** was
                            # the more common style.
                            if chunk.energy.axis.axis.cunit == '1 / m':
                                chunk.energy.axis.axis.cunit = 'm**-1'

                    if instrument is md.Inst.ESPADONS:
                        # if chunk.time is not None:
                        #     # consistent with caom2IngestEspadons.py
                        #     chunk.time_axis = 5
                        # TODO - use cfht_name as a parameter here
                        _update_energy_espadons(
                            chunk, plane.product_id[-1], headers, idx,
                            artifact.uri, fqn, observation.observation_id)

                        if chunk.position is not None:
                            # conform to stricter WCS validation
                            chunk.position_axis_1 = None
                            chunk.position_axis_2 = None
                            chunk.naxis = None
                            # CW - Ignore position wcs if a calibration file
                            # suffix list from caom2IngestEspadons.py, l389
                            # 'b', 'd', 'c', 'f', 'x'
                            # with missing spatial indicator keywords
                            if not (cfht_name.suffix in ['a', 'i', 'o', 'p']
                                    and radecsys.lower() != 'null' and
                                    headers[idx].get('RA_DEG') is not None and
                                    headers[idx].get('DEC_DEG') is not None):
                                cc.reset_position(chunk)

                        if cfht_name.suffix in ['i', 'p']:
                            _update_observable(part, chunk, cfht_name.suffix,
                                               observation.observation_id)

                        if chunk.energy is not None:
                            chunk.energy_axis = None
                        if chunk.time is not None:
                            chunk.time_axis = None
                    elif instrument in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
                        # CW
                        # Ignore position wcs if a calibration file (except 'x'
                        # type calibration) and/or position info not in header
                        # or binned 8x8
                        if (cfht_name.suffix in ['b', 'l', 'd', 'f'] or
                                ccdbin == 8 or radecsys is None or
                                ctype1 is None):
                            cc.reset_position(chunk)
                            chunk.naxis = None

                        # CW
                        # Ignore energy wcs if some type of calibration file
                        # or filter='None' or 'Open' or there is no filter
                        # match
                        #
                        # SGo - use range for energy with filter information
                        filter_md = _get_filter_md(
                            instrument, filter_name, uri)
                        if (filter_name is None or
                                filter_name in ['Open', 'NONE'] or
                                ac.FilterMetadataCache.get_fwhm(
                                    filter_md) is None or
                                cfht_name.suffix in ['b', 'l', 'd'] or
                                observation.type in ['DARK']):
                            cc.reset_energy(chunk)
                        else:
                            _update_energy_range(chunk, filter_name, filter_md)

                        # PD - in general, do not assign, unless the wcs
                        # metadata is in the FITS header
                        if chunk.time_axis is not None:
                            chunk.time_axis = None

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

                        if (cfht_name.suffix in ['a', 'o', 'x'] and
                                chunk.position is None):
                            _update_position_sitelle(chunk, headers[idx],
                                                     observation.observation_id)

                        if cfht_name.suffix == 'p':
                            if chunk.position.axis.function is None:
                                _update_position_function_sitelle(
                                    chunk, headers[idx],
                                    observation.observation_id, idx)
                            _update_sitelle_plane(observation, uri)

                    elif instrument is md.Inst.SPIROU:
                        _update_position_spirou(
                            chunk, headers[idx], observation.observation_id)

                        if cfht_name.suffix == 's':
                            part.chunks = TypedList(Chunk,)
                        # stricter WCS validation
                        chunk.naxis = None
                        if chunk.energy is not None:
                            chunk.energy_axis = None
                            if observation.type == 'DARK':
                                # caom2IngestSpirou.py, l514
                                chunk.energy = None
                        if chunk.time is not None:
                            chunk.time_axis = None
                        if chunk.position is not None:
                            chunk.position_axis_1 = None
                            chunk.position_axis_2 = None

                    elif instrument is md.Inst.WIRCAM:
                        _update_wircam_time(
                            part, chunk, headers, idx, cfht_name,
                            observation.type,
                            observation.observation_id)

                        if (cfht_name.suffix in ['f'] or
                                observation.type in
                                ['BPM', 'DARK', 'FLAT', 'WEIGHT']):
                            cc.reset_position(chunk)
                            chunk.naxis = None

                        if cfht_name.suffix == 'g':
                            _update_wircam_position_g(
                                part, chunk, headers, idx,
                                observation.observation_id)
                            temp_bandpass_name = headers[idx].get('FILTER')
                            if temp_bandpass_name == 'FakeBlank':
                                cc.reset_energy(chunk)

                        if cfht_name.suffix == 'o':
                            _update_wircam_position_o(
                                part, chunk, headers, idx,
                                observation.observation_id)

                        # position axis check is to determine if naxis should
                        # be set
                        if (cfht_name.suffix in ['d', 'f', 'g'] and
                                chunk.position_axis_1 is None):
                            # PD - 17-01-20
                            #  This is a FLAT field exposure so the position
                            #  is not relevant and only the energy and time is
                            #  relevant for discovery. The easiest correct
                            #  thing to do is to leave naxis, energyAxis, and
                            #  timeAxis all null: the same plane metadata
                            #  should be generated and that should be valid.
                            #  The most correct thing to do would be to set
                            #  naxis=2, positionAxis1=1, positionAxis2=2 (to
                            #  indicate image) and then use a suitable
                            #  coordinate system description that meant "this
                            #  patch of the inside of the dome" or maybe some
                            #  description of the pixel coordinate system
                            #  (because wcs kind of treats the sky and the
                            #  pixels as two different systems)... I don't
                            #  know how to do that and it adds very minimal
                            #  value (it allows Plane.position.dimension to be
                            #  assigned a value).
                            chunk.naxis = None
                            chunk.energy_axis = None
                            chunk.time_axis = None

                        filter_md = _get_filter_md(
                            instrument, filter_name, uri)
                        _update_energy_range(chunk, filter_name, filter_md)

                        if chunk.naxis == 2:
                            chunk.time_axis = None
        if isinstance(observation, DerivedObservation):
            if derived_type is ProvenanceType.IMCMB:
                cc.update_plane_provenance(plane, headers[1:],
                                           derived_type.value, cn.COLLECTION,
                                           _repair_imcmb_provenance_value,
                                           observation.observation_id)
            elif derived_type is ProvenanceType.COMMENT:
                cc.update_plane_provenance_single(
                    plane, headers, derived_type.value, cn.COLLECTION,
                    _repair_comment_provenance_value,
                    observation.observation_id)
            else:
                cc.update_plane_provenance(plane, headers, derived_type.value,
                                           cn.COLLECTION,
                                           _repair_filename_provenance_value,
                                           observation.observation_id)
        if instrument is md.Inst.WIRCAM and cfht_name.suffix in ['p', 's']:
            # caom2IngestWircam.py, l193
            # CW 09-01-20
            # Only the 'o' is input
            _update_plane_provenance_p(plane, observation.observation_id, 'o')
        elif instrument is md.Inst.ESPADONS and cfht_name.suffix == 'i':
            # caom2IngestEspadons.py, l714
            _update_plane_provenance_p(plane, observation.observation_id, 'o')
        elif (instrument in [md.Inst.MEGAPRIME, md.Inst.MEGACAM] and
              cfht_name.suffix == 'p'):
            # caom2IngestMegacam.py, l142
            _update_plane_provenance_p(plane, observation.observation_id, 'o')
        elif (instrument is md.Inst.SPIROU and
              cfht_name.suffix in ['e', 's', 't']):
            # caom2IngestSpirou.py, l584
            _update_plane_provenance_p(plane, observation.observation_id, 'o')

    # relies on update_plane_provenance being called
    if isinstance(observation, DerivedObservation):
        cc.update_observation_members(observation)

    logging.debug('Done update.')
    return observation


def get_calibration_level(params):
    header, suffix, uri = _decompose_params(params)
    instrument = _get_instrument(header)
    cfht_name = cn.CFHTName(ad_uri=uri, instrument=instrument)
    result = CalibrationLevel.CALIBRATED
    if (cfht_name.is_simple and not cfht_name.simple_by_suffix and not
            cfht_name.is_master_cal):
        result = CalibrationLevel.RAW_STANDARD
    return result


def get_dec_deg_from_0th_header(header):
    return header.get('DEC_DEG')


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
            temp = _get_filter_md(instrument, filter_name, uri)
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
            temp = _get_filter_md(instrument, filter_name, uri)
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
    if instrument is md.Inst.SITELLE:
        if suffix == 'p':
            num_steps = header.get('STEPNB')
            exptime = exptime * num_steps
    # units are seconds
    if exptime is None:
        cfht_name = cn.CFHTName(ad_uri=uri, instrument=instrument)
        if cfht_name.is_simple:
            # caom2IngestMegacaomdetrend.py, l438
            exptime = 0.0
    return exptime


def get_espadons_energy_resolving_power(params):
    result = None
    header, suffix, uri = _decompose_params(params)
    if _has_energy(header):
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
    return result


def get_espadons_exptime(params):
    header, suffix, uri = _decompose_params(params)
    exptime = mc.to_float(header.get('EXPTIME'))
    if suffix == 'p':
        # caom2IngestEspadons.py, l406
        exptime = 0.0
        polar_seq = mc.to_int(header.get('POLARSEQ'))
        for ii in range(1, polar_seq + 1):
            exptime += mc.to_float(header.get(f'EXPTIME{ii}'))
    # units are seconds
    if exptime is None:
        cfht_name = cn.CFHTName(ad_uri=uri, instrument=md.Inst.ESPADONS)
        if cfht_name.is_simple:
            # caom2IngestMegacaomdetrend.py, l438
            exptime = 0.0
    return exptime


def get_espadons_provenance_keywords(params):
    header, suffix, ignore = _decompose_params(params)
    result = None
    if suffix in ['i', 'p']:
        temp = header.get('REDUCTIO')
        if temp is not None:
            result = f'reduction={temp}'
    return result


def get_espadons_provenance_last_executed(header):
    result = None
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
    return result


def get_espadons_provenance_name(header):
    result = 'TCS'  # ESPaDOnS
    comments = header.get('COMMENT')
    if comments is not None:
        for comment in comments:
            if 'Upena' in comment:
                result = 'UPENA'
                break
            elif 'opera-' in comment:
                result = 'OPERA'
                break
    return result


def get_espadons_provenance_project(header):
    result = 'STANDARD PIPELINE'
    if get_espadons_provenance_name(header) == 'TCS':
        result = None
    return result


def get_espadons_provenance_version(header):
    result = None
    comments = header.get('COMMENT')
    if comments is not None:
        for comment in comments:
            if 'Upena version' in comment:
                result = comment.split('Upena version')[1]
                break
            elif 'opera-' in comment and 'build date' in comment:
                result = comment.split(' build date')[0]
                break
    return result


def get_espadons_time_refcoord_delta(params):
    exptime = get_espadons_exptime(params)
    return exptime / 86400.0  # units are d


def get_espadons_time_refcoord_val(params):
    header, suffix, uri = _decompose_params(params)
    if suffix == 'p':
        mjd_start1 = header.get('MJDSTART1')
        mjd_date1 = header.get('MJDATE1')
        mjd_obs = None
        if mjd_start1 is not None or mjd_date1 is not None:
            # caom2IngestEspadons.py, l406
            if mjd_start1 is not None:
                mjd_obs = mjd_start1
            else:
                mjd_obs = mjd_date1
    else:
        mjd_obs = _get_mjd_obs(header)
        if mjd_obs is None:
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


def get_instrument_keywords(header):
    inst_mode = header.get('INSTMODE')
    if inst_mode is not None:
        inst_mode = f'INSTMODE={inst_mode}'
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


def get_mega_provenance_last_executed(header):
    result = get_provenance_last_executed(header)
    if result is None:
        result = header.get('DATEPROC')
        if result is not None:
            result = mc.make_time(result)
    return result


def get_meta_release(header):
    # order set from:
    # caom2IngestWircam.py, l777
    result = header.get('MET_DATE')
    if result is None:
        result = header.get('DATE-OBS')
        if result is None:
            # caom2IngestEspadons.py, l625
            result = header.get('DATE-OB1')
            if result is None:
                result = header.get('DATE')
                if result is None:
                    result = header.get('REL_DATE')
                if result is None:
                    # caom2IngestMegadetrend.py, l445
                    result = header.get('TVSTART')
    return result


def get_obs_environment_humidity(header):
    result = header.get('RELHUMID')
    if result is not None and result < 0.0:
        logging.warning(f'RELHUMID invalid value {result}.')
        result = None
    return result


def get_obs_intent(header):
    # CW
    # Determine Observation.intent = obs.intent = "science" or "calibration"
    # phot & astr std & acquisitions/align are calibration.
    # from caom2IngestWircam.py, l731
    result = ObservationIntentType.CALIBRATION
    obs_type = _get_obstype(header)
    if obs_type == 'OBJECT':
        run_id = _get_run_id(header)
        if run_id[3].lower() != 'q':
            result = ObservationIntentType.SCIENCE
    return result


def get_obs_sequence_number(params):
    header, suffix, uri = _decompose_params(params)
    instrument = _get_instrument(header)
    cfht_name = cn.CFHTName(ad_uri=uri, instrument=instrument)
    result = None
    # SF 09-01-20
    # *y files are produced from other files, I am guessing the sky
    # subtraction software at CFHT copies the header from one of the exposure
    # and does not update the EXPNUM.
    #
    # SGo - because of this, use the file name to find the sequence number,
    # not the 'EXPNUM' keyword as in the originating caom2Ingest*.py scripts.
    if ((cfht_name.is_simple and not cfht_name.is_master_cal) or (
            instrument in [md.Inst.ESPADONS,
                           md.Inst.SITELLE] and suffix == 'p')):
        result = cfht_name.file_id[:-1]
    return result


def get_obs_type(header):
    result = header.get('OBSTYPE')
    if result is not None:
        if result == 'FRPTS':
            result = 'FRINGE'
        elif result == 'scatter':
            result = 'FLAT'
    return result


def get_plane_data_product_type(header):
    instrument = _get_instrument(header)
    result = DataProductType.IMAGE
    # caom2spirou.default
    if instrument in [md.Inst.ESPADONS, md.Inst.SPIROU]:
        result = DataProductType.SPECTRUM
    return result


def get_plane_data_release(header):
    # order set from:
    # caom2IngestWircam.py, l756
    #
    # from http://www.cfht.hawaii.edu/en/science/QSO/
    #
    # "The proprietary period of QSO data extends by default to 1 year + 1
    # month starting at the end of the QSO semester. For instance, data taken
    # for the 2009B semester (August 1 - January 31) will have a default
    # release date set to 02/28/2011. The extra month is allowed because of
    # possible delays in the data reduction distribution of observations
    # carried out near the end of a semester. If an extension is requested
    # during the Phase 1 period and is approved by TAC, a new date will be
    # set for this program through the QSO system. This release date for the
    # QSO data is indicated in the fits headers by the keyword REL_DATE."

    result = header.get('REL_DATE')
    if result is None:
        date_obs = header.get('DATE-OBS')
        run_id = _get_run_id(header)
        if run_id is not None:
            if run_id == 'SMEARING':
                result = header.get('DATE')
            elif ((run_id[3].lower() == 'e' or run_id[3].lower() == 'q')
                  and date_obs is not None):
                result = f'{date_obs}T00:00:00'
            else:
                obs_intent = get_obs_intent(header)
                if obs_intent == ObservationIntentType.CALIBRATION:
                    # from caom2IngestMegacamdetrend.py, l445
                    result = header.get('DATE')
                    if result is None:
                        result = header.get('TVSTART')
                if result is None:
                    logging.warning(
                        f'REL_DATE not in header. Derive from RUNID {run_id}.')
                    semester = mc.to_int(run_id[0:2])
                    rel_year = 2000 + semester + 1
                    if run_id[2] == 'A':
                        result = f'{rel_year}-08-31T00:00:00'
                    else:
                        rel_year += 1
                        result = f'{rel_year}-02-28T00:00:00'
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


def get_position_coordsys_from_0th_header(header):
    return header.get('RADECSYS')


def get_position_equinox_from_0th_header(header):
    return header.get('EQUINOX')


def get_product_type(params):
    header, suffix, ignore = _decompose_params(params)
    result = ProductType.SCIENCE
    obs_type = get_obs_intent(header)
    if obs_type == ObservationIntentType.CALIBRATION:
        result = ProductType.CALIBRATION
    if suffix in ['g', 'm', 'w', 'y']:
        result = ProductType.CALIBRATION

    # The goal is to make all file types easily findable by archive users,
    # which means having each file type show as a row in the search results.
    # With the search results limitation, a single file must be owned by a
    # plane, as that is how it gets displayed in the search results. The
    # planes must also contain plane-level metadata for the search to find.
    # Plane-level metadata is only calculated for science or calibration
    # artifacts, so any file types that might conceivably be auxiliary
    # product types are labeled as calibration, so that plane-level metadata
    # is calculated.
    #
    # Confirm the goal is find-ability in conversation with CW, SF on
    # 27-01-20.

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


def get_proposal_title(header):
    result = None
    run_id = header.get('RUNID')
    if run_id is not None:
        result = md.cache.get_title(run_id)
    return result


def get_provenance_last_executed(header):
    result = header.get('PROCDATE')
    if result is not None:
        # format like 2018-06-05HST17:21:20, which default code doesn't
        # understand
        result = mc.make_time(result)
    return result


def get_provenance_version(header):
    result = header.get('IIWIVER')
    if result is None:
        result = header.get('ORBSVER')
        if result is None:
            result = header.get('EL_SYS')
    return result


def get_ra_deg_from_0th_header(header):
    return header.get('RA_DEG')


def get_sitelle_energy_resolving_power(params):
    result = None
    header, suffix, uri = _decompose_params(params)
    if _has_energy(header):
        # from caom2IngestSitelle.py, l555+
        sitresol = header.get('SITRESOL')
        if sitresol is not None and sitresol > 0.0:
            result = sitresol
        if result is None:
            result = 1.0
            if suffix in ['a', 'c', 'f', 'o', 'x']:
                # from caom2IngestSitelle.py, l596
                crval3 = mc.to_float(header.get('FILTERLB'))
                cdelt3 = mc.to_float(header.get('FILTERBW'))
                if crval3 is not None and cdelt3 is not None:
                    result = crval3 / cdelt3
    return result


def get_sitelle_plane_data_product_type(uri):
    cfht_name = cn.CFHTName(ad_uri=uri, instrument=md.Inst.SITELLE)
    result = DataProductType.IMAGE
    if cfht_name.is_derived_sitelle:
        result = DataProductType.CUBE
    return result


def get_sitelle_time_refcoord_delta(params):
    header, suffix, uri = _decompose_params(params)
    if suffix == 'p':
        delta = None
        mjd_start = _get_mjd_start(header)
        mjd_end = mc.to_float(header.get('MJDEND'))
        # caom2IngestSitelle.py, l704
        if mjd_start is not None and mjd_end is not None:
            delta = mjd_end - mjd_start
        else:
            exp_time = header.get('EXPTIME')
            if exp_time is None:
                delta = mjd_start
    else:
        exp_time = mc.to_float(header.get('EXPTIME'))
        if exp_time is None:
            exp_time = mc.to_float(header.get('DARKTIME'))
        delta = exp_time / 86400.0
    return delta


def get_spirou_exptime(params):
    # caom2IngestSpirou.py, l530+
    header, suffix, uri = _decompose_params(params)
    result = None
    if suffix in ['a', 'c', 'd', 'f', 'o', 'r', 'x']:
        result = header.get('EXPTIME')
    elif suffix == 'p':
        result = header.get('TOTETIME')
    else:
        result = header.get('DARKTIME')
    if result is None:
        logging.warning(f'No Time WCS refcoord.delta value for {uri}.')
    return result


def get_spirou_provenance_name(header):
    result = None
    temp = header.get('RAMPSWV')
    if temp is not None:
        result = temp.split(' v')[0]
    return result


def get_spirou_provenance_version(header):
    result = None
    temp = header.get('RAMPSWV')
    if temp is not None:
        result = temp.split(' v')[1]
    return result


def get_spirou_resolution(params):
    # caom2IngestSpirou.py, l530+
    header, suffix, uri = _decompose_params(params)
    result = get_spirou_time_refcoord_delta(params)
    if suffix == 'r':
        result = result * (24.0 * 3600.0)
    else:
        result = get_spirou_exptime(params)
    if result is None:
        logging.warning(f'No Time WCS resolution value for {uri}.')
    return result


def get_spirou_time_refcoord_delta(params):
    # caom2IngestSpirou.py, l530+
    header, suffix, uri = _decompose_params(params)
    result = None
    if suffix == 'r':
        temp = header.get('FRMTIME')
    elif suffix == 'p':
        temp = header.get('TOTETIME')
    else:
        temp = header.get('DARKTIME')
    if temp is None:
        logging.warning(f'No Time WCS refcoord.delta value for {uri}.')
    else:
        result = temp / (24.0 * 3600.0)
    return result


def get_spirou_time_refcoord_naxis(params):
    # caom2IngestSpirou.py, l557
    header, suffix, uri = _decompose_params(params)
    result = 1.0
    if suffix == 'r':
        result = header.get('NREADS')
    if result is None:
        logging.warning(f'No Time WCS refcoord.naxis value for {uri}.')
    return result


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


def get_time_refcoord_delta_derived(header):
    mjd_obs = get_time_refcoord_val_derived(header)
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


def get_time_refcoord_val_derived(header):
    # CW
    # caom2IngestWircamdetrend.py, l388
    # Time - set exptime as time of one image, start and stop dates
    # as one pixel so this means crval3 is not equal to exptime
    # if TVSTART not defined, use release_date as mjdstart
    dt_str = header.get('TVSTART')
    if dt_str is None:
        dt_str = header.get('REL_DATE')
        if dt_str is None:
            dt_str = header.get('DATE')
    mjd_obs = ac.get_datetime(dt_str)
    if mjd_obs is None:
        logging.warning(f'Chunk.time.axis.function.refCoord.val is None for '
                        f'{_get_filename(header)}')
    return mjd_obs


def get_time_refcoord_val_simple(header):
    result = _get_mjd_obs(header)
    if result is None:
        temp = header.get('DATE-OBS')
        if temp is None:
            # from caom2IngestMegacam.py, l549
            temp = header.get('DATE')
        result = ac.get_datetime(temp)
    return result


def get_wircam_bandpass_name(header):
    wheel_a = header.get('WHEELADE')
    wheel_b = header.get('WHEELBDE')
    result = None
    if wheel_a == 'Open' and wheel_b != 'Open':
        result = wheel_b
    elif wheel_b == 'Open' and wheel_a != 'Open':
        result = wheel_a
    elif wheel_a == 'Open' and wheel_b == 'Open':
        result = 'Open'
    return result


def get_wircam_obs_type(params):
    header, suffix, uri = _decompose_params(params)
    result = header.get('OBSTYPE')
    # caom2IngestWircamdetrend.py, l369
    if 'weight' in uri:
        result = 'WEIGHT'
    elif 'badpix' in uri or 'hotpix' in uri or 'deadpix' in uri:
        result = 'BPM'
    elif suffix == 'g' and result is None:
        result = 'GUIDE'
    return result


def get_wircam_provenance_keywords(uri):
    suffix = cn.CFHTName(ad_uri=uri, instrument=md.Inst.WIRCAM)._suffix
    result = None
    if suffix in ['p', 's']:
        # caom2IngestWircam.py, l1063
        if suffix == 'p':
            result = 'skysubtraction=yes'
        else:
            result = 'skysubtraction=no'
    return result


def _decompose_params(params):
    header = params.get('header')
    uri = params.get('uri')
    instrument = _get_instrument(header)
    suffix = cn.CFHTName(ad_uri=uri, instrument=instrument)._suffix
    return header, suffix, uri


def _get_filename(header):
    return header.get('FILENAME')


def _get_instrument(header):
    try:
        result = md.Inst(header.get('INSTRUME'))
    except ValueError:
        # set the entry parameter to nothing, as it's only used to check
        # for hdf5
        result = cb.CFHTBuilder.get_instrument([header, {}], '')
    return result


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
    if run_id is None:
        run_id = header.get('CRUNID')
    if run_id is not None:
        if len(run_id) < 3 or len(run_id) > 9 or run_id == 'CFHT':
            # a well-known default value that indicates the past, as
            # specified in
            # caom2IngestMegacam.py, l392
            # caom2IngestWircamdetrend.py, l314
            # caom2IngestEspadons.py, l522
            logging.warning(f'Setting RUNID to default 17BE for '
                            f'{header.get("FILENAME")}.')
            run_id = '17BE'
        else:
            run_id = run_id.strip()
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
        cfht_name = cn.CFHTName(ad_uri=uri, instrument=instrument)
        if cfht_name.suffix in ['a', 'b', 'c', 'd', 'f', 'i', 'o', 'p', 'x']:
            result = True
    return result


def _is_sitelle_energy(header, uri):
    instrument = _get_instrument(header)
    result = False
    if instrument is md.Inst.SITELLE:
        cfht_name = cn.CFHTName(ad_uri=uri, instrument=instrument)
        if cfht_name.suffix in ['a', 'c', 'f', 'o', 'p', 'x']:
            result = True
    return result


def _get_filter_md(instrument, filter_name, entry):
    temp = md.filter_cache.get_svo_filter(instrument, filter_name)
    if not md.filter_cache.is_cached(instrument, filter_name):
        # want to stop ingestion if the filter name is not expected
        raise mc.CadcException(f'Could not find filter metadata for '
                               f'{filter_name} in {entry}.')
    return temp


def _is_derived(headers, cfht_name, obs_id):
    result = False
    derived_type = ''
    if cc.is_composite(headers):
        result = True
        derived_type = ProvenanceType.IMCMB
    else:
        file_type = headers[0].get('FILETYPE')
        if file_type is not None and 'alibrat' in file_type:
            logging.info(
                f'Treating {obs_id} with filetype {file_type} as derived. ')
            result = True
            derived_type = ProvenanceType.COMMENT
    if not result and not cfht_name.is_simple:
        result = True
        derived_type = ProvenanceType.FILENAME
    if cfht_name.is_derived_sitelle and cfht_name.suffix == 'z':
        result = True
        derived_type = ProvenanceType.UNDEFINED
    if cfht_name.instrument is md.Inst.MEGAPRIME and cfht_name.suffix == 'p':
        # 'p' files are processed and do have IMCMB inputs, but they are
        # additional planes on SimpleObservations, not Derived. See header
        # discussion.
        result = False
    return result, derived_type


def _update_energy_espadons(chunk, suffix, headers, idx, uri, fqn, obs_id):
    logging.debug(f'Begin _update_energy_espadons for {obs_id} {uri} {fqn}')
    if fqn is None:
        cfht_name = cn.CFHTName(ad_uri=uri,
                                instrument=md.Inst.ESPADONS)
    else:
        cfht_name = cn.CFHTName(file_name=os.path.basename(fqn),
                                instrument=md.Inst.ESPADONS)
    if cfht_name.suffix == suffix:
        axis = Axis('WAVE', 'nm')
        params = {'header': headers[idx],
                  'uri': uri}
        resolving_power = get_espadons_energy_resolving_power(params)
        coord_axis = None
        if suffix in ['a', 'c', 'f', 'o', 'x']:
            naxis1 = get_energy_function_naxis(params)
            cdelt1 = get_energy_function_delta(params)
            crval1 = get_energy_function_val(params)
            ref_coord_1 = RefCoord(0.5, crval1)
            ref_coord_2 = RefCoord(1.5, crval1 + float(naxis1)*cdelt1)
            coord_range = CoordRange1D(ref_coord_1, ref_coord_2)
            coord_axis = CoordAxis1D(axis=axis, range=coord_range)
            ssysobs = 'TOPOCENT'
            if suffix == 'o':
                ssysobs = None
        elif suffix in ['b', 'd', 'i', 'p']:
            # i, p, are done in the espadons energy data visitor, and b, d are
            # not done at all
            # CW caom2IngestEspadons.py, l393
            # Ignore energy wcs if some type of calibration file
            return
        chunk.energy = SpectralWCS(coord_axis,
                                   specsys='TOPOCENT',
                                   ssysobs=ssysobs,
                                   ssyssrc='TOPOCENT',
                                   resolving_power=resolving_power)
        chunk.energy_axis = 1
    logging.debug(f'End _update_energy_espadons for {obs_id}')


def _update_energy_range(chunk, filter_name, filter_md):
    """
    Make the MegaCam/MegaPrime energy metadata look more like the
    Gemini metadata.
    """
    cc.build_chunk_energy_range(chunk, filter_name, filter_md)
    if chunk.energy is not None:
        chunk.energy.ssysobs = 'TOPOCENT'
        chunk.energy.ssyssrc = 'TOPOCENT'
        # values from caom2megacam.default, caom2megacamdetrend.default
        chunk.energy.axis.error = CoordError(1.0, 1.0)


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
            # caom2IngestEspadons.py, l863
            dependent_axis = Axis('polarized flux', 'percent')
            dependent = Slice(dependent_axis, 3)
            new_chunk = copy.deepcopy(chunk)
            new_chunk.observable = ObservableAxis(dependent, independent)
            new_chunk._id = Chunk._gen_id()
            part.chunks.append(new_chunk)
    logging.debug(f'End _update_observable for {obs_id}')


def _update_observation_metadata(obs, headers, cfht_name, fqn, uri):
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
    #   maybe instead of storing the blueprint I just call accumulate_bp?

    # and I'm back trying to figure out how to undo the doing ...
    # the first time through, the CD4_4 value is set from CDELT4,
    # the second time through, it is not, because

    # more notes to myself - why does this not work within the pipeline?
    # it should???

    # check for files with primary headers that have NO information
    # - e.g. 2445848a
    logging.debug(f'Begin _update_observation_metadata for '
                  f'{cfht_name.file_name}')
    idx = 0
    run_id = headers[0].get('RUNID')
    if run_id is None:
        run_id = headers[0].get('CRUNID')
        if run_id is None:
            idx = 1

            logging.warning(f'Resetting the header/blueprint relationship for '
                            f'{cfht_name.file_name} in {obs.observation_id}')
            instrument = _get_instrument(headers[idx])
            if fqn is not None:
                # use the fqn to define the URI
                # TODO - leaking name structure here
                extension = '.fz'
                if instrument is md.Inst.ESPADONS:
                    extension = '.gz'
                uri = mc.build_uri(cn.ARCHIVE,
                                   os.path.basename(fqn).replace('.header',
                                                                 extension))
                # this is the fits2caom2 implementation, which returns
                # a list structure
                unmodified_headers = get_cadc_headers(f'file://{fqn}')
            elif uri is not None:
                # this is the fits2caom2 implementation, which returns
                # a list structure
                unmodified_headers = get_cadc_headers(uri)
            else:
                raise mc.CadcException(f'Cannot retrieve un-modified headers '
                                       f'for {uri}')
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

            tc.add_headers_to_obs_by_blueprint(
                obs, unmodified_headers[1:], bp, uri, cfht_name.product_id)

    logging.debug(f'End _update_observation_metadata.')
    return idx


def _update_plane_provenance_p(plane, obs_id, suffix):
    logging.debug(f'Begin _update_plane_provenance_p for {obs_id}')
    obs_member_str = mc.CaomName.make_obs_uri_from_obs_id(cn.COLLECTION,
                                                          obs_id)
    obs_member = ObservationURI(obs_member_str)
    plane_uri = PlaneURI.get_plane_uri(obs_member, f'{obs_id}{suffix}')
    plane.provenance.inputs.add(plane_uri)
    logging.debug(f'End _update_plane_provenance_p for {obs_id}')


def _update_position_spirou(chunk, header, obs_id):
    logging.debug(f'Begin _update_position_spirou for {obs_id}')
    # from caom2IngestSpirou.py, l499+
    # CW - ignore position wcs if a calibration file
    obs_type = _get_obstype(header)
    ra_deg = header.get('RA_DEG')
    dec_deg = header.get('DEC_DEG')
    ra_dec_sys = header.get('RADECSYS')
    if (obs_type not in ['OBJECT', 'ALIGN'] or
            (ra_deg is None and dec_deg is None and
             (ra_dec_sys is None or ra_dec_sys.lower() == 'null'))):
        cc.reset_position(chunk)
    logging.debug(f'Begin _update_position_spirou for {obs_id}')


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


def _update_sitelle_plane(observation, uri):
    logging.debug(
        f'Begin _update_sitelle_plane for {observation.observation_id}')
    # if the 'p' plane exists, the observation id is the same as the plane id,
    # so copy the metadata to the 'z' plane
    z_plane_key = observation.observation_id.replace('p', 'z')
    temp_z_uri = uri.replace('p', 'z')
    z_artifact_key = f'{cn.CFHTName.remove_extensions(temp_z_uri)}.hdf5'

    # fix the plane-level information for the z plane
    if z_plane_key in observation.planes.keys():
        z_plane = observation.planes[z_plane_key]
        z_plane.data_product_type = DataProductType.CUBE
        z_plane.calibration_level = CalibrationLevel.CALIBRATED
        z_plane.meta_producer = mc.get_version(APPLICATION)
        observation.meta_producer = z_plane.meta_producer
        z_plane.artifacts[z_artifact_key].meta_producer = z_plane.meta_producer

        if observation.observation_id in observation.planes.keys():
            # replicate the plane-level information from the p plane to the
            # z plane
            p_plane = observation.planes[observation.observation_id]
            p_artifact_key = uri.replace('z', 'p', 1).replace('.hdf5', '.fits.fz')
            if p_artifact_key not in p_plane.artifacts.keys():
                p_artifact_key = uri.replace(
                    'z', 'p').replace('.hdf5', '.fits')
                if p_artifact_key not in p_plane.artifacts.keys():
                    p_artifact_key = uri.replace('z', 'p', 1).replace(
                        '.hdf5', '.fits.gz')
                    if p_artifact_key not in p_plane.artifacts.keys():
                        p_artifact_key = uri.replace('z', 'p', 1).replace(
                            '.hdf5', '.fits.header')
                        if p_artifact_key not in p_plane.artifacts.keys():
                            raise mc.CadcException(
                                f'Unexpected extension name pattern for '
                                f'artifact URI {p_artifact_key} in '
                                f'{observation.observation_id}.')
            features = mc.Features()
            features.supports_latest_caom = True
            for part in p_plane.artifacts[p_artifact_key].parts.values():
                z_plane.artifacts[z_artifact_key].parts.add(cc.copy_part(part))
                for chunk in part.chunks:
                    z_plane.artifacts[z_artifact_key].parts[part.name].chunks.\
                        append(cc.copy_chunk(chunk, features))
            z_plane.artifacts[z_artifact_key].meta_producer = \
                p_plane.artifacts[p_artifact_key].meta_producer
            z_plane.provenance = p_plane.provenance
            z_plane.calibration_level = p_plane.calibration_level
            z_plane.data_product_type = p_plane.data_product_type
            z_plane.data_release = p_plane.data_release
            z_plane.meta_producer = p_plane.meta_producer
            z_plane.meta_release = p_plane.meta_release

    logging.debug('End _update_sitelle_plane')


def _update_wircam_position_o(part, chunk, headers, idx, obs_id):
    logging.debug(f'Begin _update_wircam_position_o for {obs_id}')
    part_index = mc.to_int(part.name)
    header = headers[part_index]
    ra_deg = header.get('RA_DEG')
    dec_deg = header.get('DEC_DEG')
    if chunk.position is None and ra_deg is not None and dec_deg is not None:
        logging.error(f'Adding position information for {obs_id}')
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CRVAL1'] = ra_deg
        header['CRVAL2'] = dec_deg
        wcs_parser = WcsParser(header, obs_id, idx)
        if chunk is None:
            chunk = Chunk()
        wcs_parser.augment_position(chunk)
        chunk.position_axis_1 = 1
        chunk.position_axis_2 = 2
    logging.debug(f'End _update_wircam_position_o')


def _update_wircam_position_g(part, chunk, headers, idx, obs_id):
    """'g' file position handling, which is quite unique."""
    logging.debug(f'Begin _update_wircam_position_g for {obs_id}')
    header = headers[idx]

    obj_name = header.get('OBJNAME')
    part_num = mc.to_int(part.name)
    if obj_name == 'zenith' or part_num >= 5:
        # SGo - the second clause is here, because there are only four sets
        # of position information in headers (RA/DEC of guide start on
        # arrays 1 2 3 4), and that's the only thing that is calculated
        # in the original code. Values for HDU 5+ are not written as part of
        # the override file, and thus default to 0, which fails ingestion
        # to the service.
        logging.warning(f'obj_name is {obj_name}. part_num is {part_num} No '
                        f'position for {obs_id}')
        cc.reset_position(chunk)
        return

    part_index = mc.to_int(part.name)
    header = headers[part_index]
    cd1_1 = None
    cd2_2 = None
    if header.get('CRVAL2') is not None:
        cd1_1 = mc.to_float(header.get('PIXSCAL1')) / 3600.0
        cd2_2 = mc.to_float(header.get('PIXSCAL2')) / 3600.0

    if cd1_1 is None or cd2_2 is None:
        cc.reset_position(chunk)
        return

    naxis_1 = header.get('NAXIS1')
    if naxis_1 is None:
        naxis_1 = header.get('ZNAXIS1')
    naxis_2 = header.get('NAXIS2')
    if naxis_2 is None:
        naxis_2 = header.get('ZNAXIS2')

    wcgd_ra = header.get(f'WCGDRA{part.name}')
    wcgd_dec = header.get(f'WCGDDEC{part.name}')
    cr_val1, cr_val2 = ac.build_ra_dec_as_deg(wcgd_ra, wcgd_dec, frame='fk5')
    if mc.to_float(obs_id) < 980000:
        cr_val2 *= 15.0

    # caom2IngestWircam.py, l367
    header['NAXIS1'] = naxis_1
    header['NAXIS2'] = naxis_2
    header['CRPIX1'] = naxis_1 / 2.0
    header['CRPIX2'] = naxis_2 / 2.0
    header['CRVAL1'] = cr_val1
    header['CRVAL2'] = cr_val2
    header['CD1_1'] = -1.0 * cd1_1
    header['CD1_2'] = 0.0
    header['CD2_1'] = 0.0
    header['CD2_2'] = cd2_2

    wcs_parser = WcsParser(header, obs_id, idx)
    if chunk is None:
        chunk = Chunk()
    wcs_parser.augment_position(chunk)
    chunk.position_axis_1 = 1
    chunk.position_axis_2 = 2
    logging.debug(f'End _update_wircam_position_g for {obs_id}')


def _update_wircam_time(part, chunk, headers, idx, cfht_name, obs_type,
                        obs_id):
    logging.debug(f'Begin _update_wircam_time for {obs_id}')
    if cfht_name.suffix == 'g':
        # construct TemporalWCS for 'g' files from the CAOM2 pieces
        # because the structure of 'g' files is so varied, it's not possible
        # to hand over even part of the construction to the blueprint.

        # SF - 07-05-20
        # so NAXIS here (where 'here' is a 'g' file) is ZNAXIS=3: time sequence
        # of images of the guiding camera => this means try ZNAXIS* keyword
        # values before trying NAXIS*, hence the header lookup code

        ref_coord_val = headers[0].get('MJD-OBS')
        if ref_coord_val is None:
            ref_coord_val = headers[1].get('MJD-OBS')
        part_index = mc.to_int(part.name)
        part_header = headers[part_index]

        if chunk.time is None:
            from caom2 import TemporalWCS
            chunk.time = TemporalWCS(CoordAxis1D(Axis('TIME', 'd')),
                                     timesys='UTC')

        if chunk.time.axis is None:
            chunk.time.axis = CoordAxis1D(axis=Axis('TIME', 'd'),
                                          error=None,
                                          range=None,
                                          bounds=None,
                                          function=None)

        if chunk.time.axis.error is None:
            from caom2 import CoordError
            chunk.time.axis.error = CoordError(rnder=0.0000001,
                                               syser=0.0000001)

        if chunk.time.axis.function is None:
            from caom2 import CoordFunction1D
            ref_coord = RefCoord(pix=0.5,
                                 val=mc.to_float(ref_coord_val))

            time_index = part_header.get('ZNAXIS')
            if time_index is None:
                time_index = part_header.get('NAXIS')
            naxis_key = f'ZNAXIS{time_index}'
            time_naxis = part_header.get(naxis_key)
            if time_naxis is None:
                naxis_key = f'NAXIS{time_index}'
                time_naxis = part_header.get(naxis_key)

            if (time_naxis is not None and time_index is not None and
                    time_index == 3):
                # caom2.4 wcs validation conformance
                chunk.time_axis = 3

            # CW
            # Define time samples for guidecube data
            # Guiding time doesn't seem to match up very well, so just say that
            # use MJD-OBS gnaxis3 and WCPERIOD
            #
            # code from caom2IngestWircam.py, l876+
            wcgdra1_0 = headers[0].get('WCGDRA1')
            wc_period_0 = headers[0].get('WCPERIOD')
            wc_period = None
            if wcgdra1_0 is not None and wc_period_0 is not None:
                wc_period = wc_period_0
            else:
                wcgdra1_1 = headers[1].get('WCGDRA1')
                wc_period_1 = headers[1].get('WCPERIOD')
                if wcgdra1_1 is not None and wc_period_1 is not None:
                    wc_period = wc_period_1

            time_delta = None
            if wc_period is not None:
                if wc_period < 0.0:
                    wc_period = 100.0
                # caom2IngestWircam.py, l375
                chunk.time.exposure = wc_period / 1000.0
                chunk.time.resolution = chunk.time.exposure
                time_delta = chunk.time.exposure / 86400.0
            chunk.time.axis.function = CoordFunction1D(naxis=time_naxis,
                                                       delta=time_delta,
                                                       ref_coord=ref_coord)

    # fits2caom2 prefers ZNAXIS to NAXIS, but the originating scripts
    # seem to prefer NAXIS, so odd as this placement seems, do not rely
    # on function execution, because it affects NAXIS, not ZNAXIS - sigh
    if (obs_type not in ['BPM', 'DARK', 'FLAT', 'WEIGHT'] and
            cfht_name.suffix != 'g'):
        n_exp = headers[idx].get('NEXP')
        if n_exp is not None:
            # caom2IngestWircam.py, l843
            chunk.time.axis.function.naxis = mc.to_int(n_exp)

    logging.debug(f'End _update_wircam_time for {obs_id}')


def _identify_instrument(uri):
    logging.debug(f'Begin _identify_instrument for uri {uri}.')
    config = mc.Config()
    config.get_executors()
    cfht_builder = cb.CFHTBuilder(config)
    storage_name = cfht_builder.build(mc.CaomName(uri).file_name)
    logging.debug(
        f'End _identify_instrument for uri {uri} instrument is '
        f'{storage_name.instrument}.')
    return storage_name.instrument


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
        if ((not mc.StorageName.is_preview(uri)) and
                (not mc.StorageName.is_hdf5(uri))):
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
            file_uri = mc.build_uri(cn.ARCHIVE, os.path.basename(ii))
            result.append(file_uri)
    else:
        raise mc.CadcException(
            f'Could not define uri from these args {args}')
    return result


def _repair_comment_provenance_value(value, obs_id):
    logging.debug(f'Begin _repair_comment_provenance_value for {obs_id}')
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
                prov_obs_id = cn.CFHTName(file_name=prov_prod_id).obs_id
                # 0 - observation
                # 1 - plane
                results.append([prov_obs_id, prov_prod_id])
    logging.debug(f'End _repair_comment_provenance_value')
    return results


def _repair_filename_provenance_value(value, obs_id):
    logging.debug(f'Begin _repair_filename_provenance_value for {obs_id}')
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
    logging.debug(f'End _repair_filename_provenance_value')
    return prov_obs_id, prov_prod_id


def _repair_imcmb_provenance_value(value, obs_id):
    logging.debug(f'Begin _repair_imcmb_provenance_value for {obs_id}')
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
        prov_obs_id = cn.CFHTName(file_name=prov_prod_id).obs_id
    logging.debug(f'End _repair_imcmb_provenance_value')
    return prov_obs_id, prov_prod_id


def _cfht_args_parser():
    args = get_gen_proc_arg_parser().parse_args()
    if args.not_connected:
        # stop attempts to connect to SVO
        md.filter_cache.connected = False
    return args


def to_caom2():
    """This function is called by pipeline execution. It must have this name.
    """
    args = _cfht_args_parser()
    uris = _get_uris(args)
    blueprints = _build_blueprints(uris)
    result = gen_proc(args, blueprints)
    return result


def cfht_main_app():
    """External application."""
    args = _cfht_args_parser()
    try:
        result = to_caom2()
        logging.debug(f'Done {APPLICATION} processing.')
        sys.exit(result)
    except Exception as e:
        logging.error(f'Failed {APPLICATION} execution for {args}.')
        tb = traceback.format_exc()
        logging.error(tb)
        sys.exit(-1)
