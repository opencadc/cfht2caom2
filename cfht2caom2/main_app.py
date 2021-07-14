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
    - SPIROU/Espadons - 'p' is polarized + Derived, a different observation
    - there are other processed files for single exposures
    - users want 'p' files

- conclusion - one plane / file, because users want to see one row / file in
  the query results

  JJK - slack - 01-04-20
  CFHT files are independent, in that, for example, a user does not require a
  'p' file to understand the content of a 'b' file. Given this independence,
  it is acceptable to map one plane / file.

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
cutout operation later. If the coverage had significant gaps (eg SCUBA or
SCUBA2 from JCMT) then the extra detail in bounds would enable better
discovery (as the gaps would be captured in the plane metadata). In the case
of espadons I don't think the gaps are significant (iirc, espadons is an
eschelle spectrograph but I don't recall whether the discontinuity between
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

SITELLE 'v' files:
- SF - 02-07-20
- the v files are not raw, they were derived from others, and the data
  reduction pipeline - based on orbs - probably did not do proper header
  copy/paste
- so ok to assume sequence number (RUN_ID, which is not in 'v' headers)
  == observationID minus suffix

SF 22-12-20
- as a general rule, fix typos from header metadata in CAOM2 records

SF 12-04-21
- SPIRou 'r' files are 'ramp' files, and should be 'RAW' for caom2:
  https://www.cfht.hawaii.edu/Instruments/SPIRou/FileStructure.php

"""

import importlib
import logging
import os
import sys
import traceback

from datetime import timedelta
from enum import Enum

from caom2 import Observation, CalibrationLevel, ObservationIntentType
from caom2 import ProductType, DerivedObservation, TypedList, Chunk
from caom2 import DataProductType
from caom2utils import ObsBlueprint, get_gen_proc_arg_parser, gen_proc
from caom2utils import get_cadc_headers
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import client_composable
from caom2pipe import manage_composable as mc
from caom2pipe import translate_composable as tc

from cfht2caom2 import cfht_builder as cb
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import instruments
from cfht2caom2 import metadata as md


__all__ = ['cfht_main_app', 'to_caom2', 'update', 'APPLICATION']


APPLICATION = 'cfht2caom2'

# All comments denoted 'CW' are copied from
# ssh://goliaths@gimli3/srv/cadc/git/wcaom2archive.git
# from /cfht2caom2/scripts


class ProvenanceType(Enum):
    """The different types of header values that identify provenance
    information. Used to specify different functions for
    execution."""

    COMMENT = 'COMMENT'
    FILENAME = 'FILENAM'
    IMCMB = 'IMCMB'
    UNDEFINED = 'UNDEFINED'  # hdf5 file support


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
    cfht_name = cn.cfht_names.get(uri)
    bp.configure_position_axes((1, 2))
    if not (
        cfht_name.suffix == 'p' and cfht_name.instrument == md.Inst.SPIROU
    ):
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

    bp.set_default('Observation.algorithm.name', None)

    bp.set(
        'Observation.environment.elevation',
        'get_environment_elevation(header)',
    )
    bp.set(
        'Observation.environment.humidity',
        'get_obs_environment_humidity(header)',
    )

    # title is select title from runid_title where proposal_id = 'runid'
    # this obtained from cache.yml now
    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'RUNID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PI_NAME')
    bp.set('Observation.proposal.project', 'get_proposal_project(header)')
    bp.set('Observation.proposal.title', 'get_proposal_title(header)')

    bp.set('Observation.instrument.name', cfht_name.instrument.value)
    bp.set(
        'Observation.instrument.keywords',
        'get_instrument_keywords(header)',
    )

    bp.set('Observation.target.standard', 'get_target_standard(params)')

    bp.clear('Observation.target_position.coordsys')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'OBJRADEC')
    bp.clear('Observation.target_position.equinox')
    bp.add_fits_attribute('Observation.target_position.equinox', 'OBJEQN')
    bp.add_fits_attribute('Observation.target_position.equinox', 'OBJEQUIN')
    bp.set(
        'Observation.target_position.point.cval1',
        'get_target_position_cval1(params)',
    )
    bp.set(
        'Observation.target_position.point.cval2',
        'get_target_position_cval2(params)',
    )

    bp.set('Observation.telescope.name', 'CFHT 3.6m')
    x, y, z = ac.get_geocentric_location('cfht')
    bp.set('Observation.telescope.geoLocationX', x)
    bp.set('Observation.telescope.geoLocationY', y)
    bp.set('Observation.telescope.geoLocationZ', z)

    bp.set('Plane.dataProductType', 'get_plane_data_product_type(params)')
    bp.set('Plane.calibrationLevel', 'get_calibration_level(params)')
    bp.set('Plane.dataRelease', 'get_plane_data_release(header)')
    bp.set('Plane.metaRelease', 'get_meta_release(header)')

    bp.set(
        'Plane.provenance.lastExecuted',
        'get_provenance_last_executed(header)',
    )
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
    if cfht_name.instrument not in [
        md.Inst.MEGACAM, md.Inst.MEGAPRIME, md.Inst.WIRCAM
    ]:
        bp.set('Chunk.energy.axis.axis.ctype', 'get_energy_ctype(header)')
        bp.set('Chunk.energy.axis.axis.cunit', 'get_energy_cunit(header)')
        bp.set('Chunk.energy.axis.error.rnder', 1.0)
        bp.set('Chunk.energy.axis.error.syser', 1.0)
        bp.set(
            'Chunk.energy.axis.function.delta',
            'get_energy_function_delta(params)',
        )
        bp.set(
            'Chunk.energy.axis.function.naxis',
            'get_energy_function_naxis(params)',
        )
        bp.set(
            'Chunk.energy.axis.function.refCoord.pix',
            'get_energy_function_pix(params)',
        )
        bp.set(
            'Chunk.energy.axis.function.refCoord.val',
            'get_energy_function_val(params)',
        )
        bp.clear('Chunk.energy.bandpassName')
        bp.add_fits_attribute('Chunk.energy.bandpassName', 'FILTER')
        bp.set(
            'Chunk.energy.resolvingPower',
            'get_energy_resolving_power(params)',
        )
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

    # TODO - this is really really wrong that is_simple is not sufficient
    # to make the distinction between the appropriate implementations.
    if cfht_name.is_simple and not cfht_name.is_master_cal:
        bp.set(
            'Chunk.time.axis.function.delta',
            'get_time_refcoord_delta_simple(params)',
        )
        bp.set(
            'Chunk.time.axis.function.refCoord.val',
            'get_time_refcoord_val_simple(header)',
        )
    else:
        bp.set(
            'Chunk.time.axis.function.delta',
            'get_time_refcoord_delta_derived(header)',
        )
        bp.set(
            'Chunk.time.axis.function.refCoord.val',
            'get_time_refcoord_val_derived(header)',
        )
    bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)

    if cfht_name.instrument is md.Inst.ESPADONS:
        _accumulate_espadons_bp(bp, cfht_name)
    elif cfht_name.instrument in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
        _accumulate_mega_bp(bp, uri, cfht_name)
    elif cfht_name.instrument is md.Inst.SITELLE:
        _accumulate_sitelle_bp(bp, uri, cfht_name)
    elif cfht_name.instrument is md.Inst.SPIROU:
        _accumulate_spirou_bp(bp, uri, cfht_name)
    elif cfht_name.instrument is md.Inst.WIRCAM:
        _accumulate_wircam_bp(bp, uri, cfht_name)

    logging.debug('Done accumulate_bp.')


def _accumulate_espadons_bp(bp, cfht_name):
    """Configure the ESPaDOnS-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_espadons_bp.')

    # bp.set('Observation.target.targetID', '_get_gaia_target_id(header)')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'RADECSYS')

    bp.set(
        'Plane.provenance.keywords',
        'get_espadons_provenance_keywords(params)',
    )
    bp.set(
        'Plane.provenance.lastExecuted',
        'get_espadons_provenance_last_executed(header)',
    )
    bp.set('Plane.provenance.name', 'get_espadons_provenance_name(header)')
    bp.set(
        'Plane.provenance.project',
        'get_espadons_provenance_project(header)',
    )
    bp.set(
        'Plane.provenance.reference',
        'get_espadons_provenance_reference(header)',
    )
    bp.set(
        'Plane.provenance.version',
        'get_espadons_provenance_version(header)',
    )

    bp.set(
        'Chunk.energy.resolvingPower',
        'get_espadons_energy_resolving_power(params)',
    )

    # constants from caom2espadons.config
    bp.set('Chunk.position.axis.axis1.ctype', 'RA---TAN')
    bp.set('Chunk.position.axis.axis2.ctype', 'DEC--TAN')
    bp.set('Chunk.position.axis.function.dimension.naxis1', 1)
    bp.set('Chunk.position.axis.function.dimension.naxis2', 1)
    bp.set('Chunk.position.axis.function.refCoord.coord1.pix', 1.0)
    bp.clear('Chunk.position.axis.function.refCoord.coord1.val')
    bp.add_fits_attribute(
        'Chunk.position.axis.function.refCoord.coord1.val', 'RA_DEG'
    )
    bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 1.0)
    bp.clear('Chunk.position.axis.function.refCoord.coord2.val')
    bp.add_fits_attribute(
        'Chunk.position.axis.function.refCoord.coord2.val', 'DEC_DEG'
    )
    # CW
    # Fibre size is 1.6", i.e. 0.000444 deg
    bp.set('Chunk.position.axis.function.cd11', -0.000444)
    bp.set('Chunk.position.axis.function.cd12', 0.0)
    bp.set('Chunk.position.axis.function.cd21', 0.0)
    bp.set('Chunk.position.axis.function.cd22', 0.000444)

    bp.add_fits_attribute('Chunk.position.equinox', 'EQUINOX')

    bp.set(
        'Chunk.time.axis.function.delta',
        'get_espadons_time_refcoord_delta(params)',
    )
    bp.set(
        'Chunk.time.axis.function.refCoord.val',
        'get_espadons_time_refcoord_val(params)',
    )

    bp.set('Chunk.time.exposure', 'get_espadons_exptime(params)')
    bp.set('Chunk.time.resolution', 'get_espadons_exptime(params)')

    if cfht_name.suffix == 'p':
        bp.configure_polarization_axis(6)
        # caom2IngestEspadons.py, l209, lTODO
        bp.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
        bp.set('Chunk.polarization.axis.function.delta', 1)
        bp.set('Chunk.polarization.axis.function.naxis', 1)
        bp.set('Chunk.polarization.axis.function.refCoord.pix', 1)
        bp.set(
            'Chunk.polarization.axis.function.refCoord.val',
            'get_polarization_function_val(header)',
        )

    logging.debug('Done _accumulate_espadons_bp.')


def _accumulate_mega_bp(bp, uri, cfht_name):
    """Configure the MegaCam/MegaPrime-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_mega_bp.')

    bp.set(
        'Plane.provenance.lastExecuted',
        'get_mega_provenance_last_executed(header)',
    )
    bp.set_default('Plane.provenance.name', 'ELIXIR')
    bp.set_default(
        'Plane.provenance.reference',
        'http://www.cfht.hawaii.edu/Instruments/Elixir/',
    )

    logging.debug('Done _accumulate_mega_bp.')


def _accumulate_sitelle_bp(bp, uri, cfht_name):
    """Configure the Sitelle-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_sitelle_bp.')

    if cfht_name.suffix == 'v':
        bp.set('Observation.intent', ObservationIntentType.SCIENCE)
        bp.set('Observation.sequenceNumber', cfht_name.product_id[:-1])
        bp.set('Plane.dataRelease', 'get_sitelle_v_plane_data_release(header)')
        bp.clear('Plane.provenance.version')
        bp.add_fits_attribute('Plane.provenance.version', 'PROGRAM')
        bp.set('Artifact.productType', ProductType.SCIENCE)
    bp.set('Plane.dataProductType', 'get_sitelle_plane_data_product_type(uri)')
    bp.set_default('Plane.provenance.name', 'ORBS')
    bp.set_default('Plane.provenance.reference', 'http://ascl.net/1409.007')

    bp.set(
        'Chunk.energy.resolvingPower',
        'get_sitelle_energy_resolving_power(params)',
    )

    bp.set(
        'Chunk.time.axis.function.delta',
        'get_sitelle_time_refcoord_delta(params)',
    )
    bp.set('Chunk.time.axis.function.refCoord.val', '_get_mjd_start(header)')

    logging.debug('End _accumulate_sitelle_bp.')


def _accumulate_spirou_bp(bp, uri, cfht_name):
    """Configure the SPIRou-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    bp.set('Observation.target.targetID', '_get_gaia_target_id(header)')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'RADECSYS')

    if cfht_name.suffix == 'r':
        pass
    elif cfht_name.suffix in ['a', 'c', 'd', 'f', 'o', 'x']:
        bp.set('Plane.provenance.name', 'get_spirou_provenance_name(header)')
        bp.set(
            'Plane.provenance.reference',
            'http://www.cfht.hawaii.edu/Instruments/SPIRou/',
        )
        bp.set(
            'Plane.provenance.version',
            'get_spirou_provenance_version(header)',
        )
    else:
        bp.set('Plane.provenance.name', 'DRS')
        bp.set(
            'Plane.provenance.reference',
            (
                'https://www.cfht.hawaii.edu/Instruments/SPIRou/'
                'SPIRou_pipeline.php'
            ),
        )
        bp.clear('Plane.provenance.version')
        bp.add_fits_attribute('Plane.provenance.version', 'VERSION')

    bp.clear('Plane.provenance.lastExecuted')
    bp.add_fits_attribute('Plane.provenance.lastExecuted', 'DRSPDATE')

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
    bp.set(
        'Chunk.position.axis.function.refCoord.coord1.val',
        'get_ra_deg_from_0th_header(header)',
    )
    bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 1.0)
    bp.set(
        'Chunk.position.axis.function.refCoord.coord2.val',
        'get_dec_deg_from_0th_header(header)',
    )
    bp.set('Chunk.position.axis.function.cd11', -0.00035833)
    bp.set('Chunk.position.axis.function.cd12', 0.0)
    bp.set('Chunk.position.axis.function.cd21', 0.0)
    bp.set('Chunk.position.axis.function.cd22', 0.00035833)

    bp.set(
        'Chunk.position.coordsys',
        'get_position_coordsys_from_0th_header(header)',
    )
    bp.set(
        'Chunk.position.equinox',
        'get_position_equinox_from_0th_header(header)',
    )

    if cfht_name.suffix not in ['g', 'p']:
        bp.set(
            'Chunk.time.axis.function.delta',
            'get_spirou_time_refcoord_delta(params)',
        )
        bp.set(
            'Chunk.time.axis.function.naxis',
            'get_spirou_time_refcoord_naxis(params)',
        )
        bp.set('Chunk.time.exposure', 'get_spirou_exptime(params)')
        bp.set('Chunk.time.resolution', 'get_spirou_resolution(params)')

    if cfht_name.suffix == 'p':
        bp.configure_polarization_axis(7)
        bp.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
        bp.set('Chunk.polarization.axis.function.naxis', 1)
        bp.set('Chunk.polarization.axis.function.delta', 1.0)
        bp.set('Chunk.polarization.axis.function.refCoord.pix', 1.0)


def _accumulate_wircam_bp(bp, uri, cfht_name):
    """Configure the WIRCam-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_wircam_bp.')
    bp.set('Observation.type', 'get_wircam_obs_type(params)')

    bp.set('Plane.provenance.keywords', 'get_wircam_provenance_keywords(uri)')
    bp.set_default('Plane.provenance.name', 'IIWI')
    bp.set_default(
        'Plane.provenance.reference',
        'http://www.cfht.hawaii.edu/Instruments/Imaging/WIRCam',
    )

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

    headers = kwargs.get('headers')
    fqn = kwargs.get('fqn')
    uri = kwargs.get('uri')
    subject = kwargs.get('subject')

    ingesting_hdf5 = False

    # TODO - use the instrument name in the observation? might not work
    # for a redo
    if uri is None:
        instrument = md.Inst(observation.instrument.name)
    else:
        suffix = cn.cfht_names.get(uri).suffix
        if suffix == 'z':
            instrument = md.Inst.SITELLE
            ingesting_hdf5 = True
            logging.info(
                f'Ingesting the hdf5 plane for {observation.observation_id}'
            )
        else:
            instrument = cb.CFHTBuilder.get_instrument(headers, uri)

    if instrument is md.Inst.MEGACAM:
        # need the 'megacam' for the filter lookup at SVO, but there is only
        # a 'MegaPrime' instrument in the CAOM collection at CADC
        # see e.g. 2003A.frpts.z.36.00
        observation.instrument = cc.copy_instrument(
            observation.instrument, md.Inst.MEGAPRIME.value
        )

    # fqn is not always defined
    if uri is not None:
        cfht_name = cn.cfht_names.get(uri)
    elif fqn is not None:
        temp_uri = mc.build_uri(cn.ARCHIVE, os.path.basename(fqn))
        cfht_name = cn.cfht_names.get(temp_uri)
    else:
        raise mc.CadcException(
            f'Cannot define a CFHTName instance for '
            f'{observation.observation_id}'
        )

    if ingesting_hdf5:
        # avoid all the code that references undefined headers variable
        if not isinstance(observation, DerivedObservation):
            observation = cc.change_to_composite(
                observation, 'scan', cn.COLLECTION
            )
        _update_sitelle_plane(observation, uri)
        logging.debug('Done hdf5 update.')
        return observation

    is_derived, derived_type = _is_derived(
        headers, cfht_name, observation.observation_id
    )
    if is_derived and not isinstance(observation, DerivedObservation):
        logging.info(
            f'{observation.observation_id} will be changed to a Derived '
            f'Observation.'
        )
        algorithm_name = 'master_detrend'
        if observation.observation_id[-1] == 'p':
            if cfht_name.has_polarization:
                algorithm_name = 'polarization'
            else:
                algorithm_name = 'scan'
        observation = cc.change_to_composite(
            observation, algorithm_name, cn.COLLECTION
        )

    if instrument is md.Inst.SITELLE and cfht_name.suffix == 'v':
        idx = 0
    else:
        idx = _update_observation_metadata(
            observation, headers, cfht_name, fqn, uri, subject
        )
    x = instruments.instrument_factory(
        instrument,
        headers,
        idx,
        cfht_name,
        observation,
    )
    x.update_observation()
    for plane in observation.planes.values():
        if plane.product_id != cfht_name.product_id:
            # do only the work for the applicable plane
            continue

        x.plane = plane
        for artifact in plane.artifacts.values():
            if artifact.uri != cfht_name.file_uri:
                continue
            params = {
                'uri': artifact.uri,
                'header': headers[idx],
            }
            if get_calibration_level(params) == CalibrationLevel.RAW_STANDARD:
                time_delta = get_time_refcoord_delta_simple(params)
            else:
                time_delta = get_time_refcoord_delta_derived(headers[idx])

            for part in artifact.parts.values():
                if instrument is md.Inst.SPIROU and cfht_name.suffix == 's':
                    part.chunks = TypedList(
                        Chunk,
                    )

                for chunk in part.chunks:
                    cc.undo_astropy_cdfix_call(chunk, time_delta)
                    x.part = part
                    x.chunk = chunk
                    x.update_chunk()

        if isinstance(observation, DerivedObservation):
            if derived_type is ProvenanceType.IMCMB:
                cc.update_plane_provenance(
                    plane,
                    headers[1:],
                    derived_type.value,
                    cn.COLLECTION,
                    _repair_imcmb_provenance_value,
                    observation.observation_id,
                )
            elif derived_type is ProvenanceType.COMMENT:
                cc.update_plane_provenance_single(
                    plane,
                    headers,
                    derived_type.value,
                    cn.COLLECTION,
                    _repair_comment_provenance_value,
                    observation.observation_id,
                )
            else:
                cc.update_plane_provenance(
                    plane,
                    headers,
                    derived_type.value,
                    cn.COLLECTION,
                    _repair_filename_provenance_value,
                    observation.observation_id,
                )
        x.update_plane()
        # this is here, because the bits that are being copied have been
        # created/modified by the update_chunk call
        if instrument is md.Inst.SITELLE and cfht_name.suffix == 'p':
            _update_sitelle_plane(observation, uri)

    # relies on update_plane_provenance being called
    if isinstance(observation, DerivedObservation):
        cc.update_observation_members(observation)

    if cfht_name.suffix in ['p', 'y'] and instrument is md.Inst.WIRCAM:
        # complete the ingestion of the missing bits of a sky construct file
        _update_wircam_plane(observation, cfht_name)

    logging.debug('Done update.')
    return observation


def get_calibration_level(params):
    header, cfht_name, uri = _decompose_params(params)
    result = CalibrationLevel.CALIBRATED
    if (
        cfht_name.is_simple
        and not cfht_name.simple_by_suffix
        and not cfht_name.is_master_cal
    ):
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
    header, cfht_name, uri = _decompose_params(params)
    if _has_energy(header):
        if _is_espadons_energy(cfht_name):
            # caom2IngestEspadons.py l639
            result = 0.0031764
        elif _is_sitelle_energy(cfht_name):
            result = header.get('CDELT3')
            if result is None:
                # units in file are nm, units in blueprint are Angstroms
                result = 10.0 * mc.to_float(header.get('FILTERBW'))
        else:
            filter_name = header.get('FILTER')
            temp, ignore = _get_filter_md(
                cn.cfht_names.get(uri).instrument, filter_name, uri
            )
            result = ac.FilterMetadataCache.get_fwhm(temp)
    return result


def get_energy_function_naxis(params):
    result = 1.0
    header, cfht_name, uri = _decompose_params(params)
    if _is_espadons_energy(cfht_name):
        # caom2IngestEspadons.py l636
        result = 213542
    elif _is_sitelle_energy(cfht_name):
        result = header.get('NAXIS3')
        if result is None:
            result = 1.0
    return result


def get_energy_function_pix(params):
    result = None
    header, cfht_name, uri = _decompose_params(params)
    if _has_energy(header):
        result = 1.0
        if _is_espadons_energy(cfht_name):
            # caom2IngestEspadons.py l637
            result = 0.5
        elif _is_sitelle_energy(cfht_name):
            result = header.get('CRPIX3')
            if result is None:
                result = 0.5
    return result


def get_energy_function_val(params):
    result = None
    header, cfht_name, uri = _decompose_params(params)
    if _has_energy(header):
        if _is_espadons_energy(cfht_name):
            # caom2IngestEspadons.py l638
            result = 370.0
        elif _is_sitelle_energy(cfht_name):
            result = header.get('CRVAL3')
            if result is None:
                # units in file are nm, units in blueprint are Angstroms
                result = 10.0 * mc.to_float(header.get('FILTERLB'))
        else:
            filter_name = header.get('FILTER')
            temp, ignore = _get_filter_md(
                cn.cfht_names.get(uri).instrument, filter_name, uri
            )
            result = ac.FilterMetadataCache.get_central_wavelength(temp)
    return result


def get_energy_resolving_power(params):
    result = None
    header, cfht_name, uri = _decompose_params(params)
    if _has_energy(header):
        delta = get_energy_function_delta(params)
        val = get_energy_function_val(params)
        result = None
        if delta is not None and val is not None:
            result = val / delta
    return result


def get_environment_elevation(header):
    elevation = mc.to_float(header.get('TELALT'))
    if elevation is not None and not (0.0 <= elevation <= 90.0):
        logging.info(
            f'Setting elevation to None for {_get_filename(header)} '
            f'because the value is {elevation}.'
        )
        elevation = None
    return elevation


def get_exptime(params):
    header, cfht_name, uri = _decompose_params(params)
    exptime = mc.to_float(header.get('EXPTIME'))
    if cfht_name.instrument is md.Inst.SITELLE:
        if cfht_name.suffix == 'p':
            num_steps = header.get('STEPNB')
            exptime = exptime * num_steps
    # units are seconds
    if exptime is None:
        if cfht_name.is_simple:
            # caom2IngestMegacaomdetrend.py, l438
            exptime = 0.0
    return exptime


def get_espadons_energy_resolving_power(params):
    result = None
    header, cfht_name, uri = _decompose_params(params)
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
    header, cfht_name, uri = _decompose_params(params)
    exptime = mc.to_float(header.get('EXPTIME'))
    if cfht_name.suffix == 'p':
        # caom2IngestEspadons.py, l406
        exptime = 0.0
        polar_seq = mc.to_int(header.get('POLARSEQ'))
        for ii in range(1, polar_seq + 1):
            exptime += mc.to_float(header.get(f'EXPTIME{ii}'))
    # units are seconds
    if exptime is None:
        cfht_name = cn.cfht_names.get(uri)
        if cfht_name.is_simple:
            # caom2IngestMegacaomdetrend.py, l438
            exptime = 0.0
    return exptime


def get_espadons_provenance_keywords(params):
    header, cfht_name, ignore = _decompose_params(params)
    result = None
    if cfht_name.suffix in ['i', 'p']:
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


def get_espadons_provenance_reference(header):
    result = 'http://www.cfht.hawaii.edu/Instruments/Spectroscopy/Espadons/'
    temp = get_espadons_provenance_name(header)
    if temp == 'UPENA':
        result = 'http://www.cfht.hawaii.edu/Instruments/Upena/'
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
    header, cfht_name, uri = _decompose_params(params)
    if cfht_name.suffix == 'p':
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
            if (
                date_obs is None
                or time_obs is None
                or date_obs == '1970-01-01'
                or date_obs == '1970-00-01'
            ):
                hst_time = header.get('HSTTIME)')
                # fmt 'Mon Nov 27 15:58:17 HST 2006'
                mjd_obs = ac.get_datetime(hst_time)
            else:
                mjd_obs_str = f'{date_obs}T{time_obs}'
                mjd_obs = ac.get_datetime(mjd_obs_str)
            mjd_obs = mjd_obs.value
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
    if obs_type is None:
        # no 'OBSTYPE' keyword, so fits2caom2 will set the value to
        # science
        result = None
    elif obs_type == 'OBJECT':
        run_id = _get_run_id(header)
        if run_id[3].lower() != 'q':
            result = ObservationIntentType.SCIENCE
    return result


def get_obs_sequence_number(params):
    header, cfht_name, uri = _decompose_params(params)
    result = None
    # SF 09-01-20
    # *y files are produced from other files, I am guessing the sky
    # subtraction software at CFHT copies the header from one of the exposure
    # and does not update the EXPNUM.
    #
    # SGo - because of this, use the file name to find the sequence number,
    # not the 'EXPNUM' keyword as in the originating caom2Ingest*.py scripts.
    if (cfht_name.is_simple and not cfht_name.is_master_cal) or (
        cfht_name.instrument
        in [md.Inst.ESPADONS, md.Inst.SITELLE, md.Inst.SPIROU]
        and cfht_name.suffix == 'p'
    ):
        result = cfht_name.file_id[:-1]
    return result


def get_obs_type(header):
    result = _get_obstype(header)
    if result is not None:
        if result == 'FRPTS':
            result = 'FRINGE'
        elif result == 'scatter':
            result = 'FLAT'
    return result


def get_plane_data_product_type(params):
    header, cfht_name, uri = _decompose_params(params)
    # caom2wircam.default
    # caom2wircamdetrend.default
    result = DataProductType.IMAGE
    # caom2spirou.default
    if cfht_name.instrument in [md.Inst.ESPADONS, md.Inst.SPIROU]:
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
            elif (
                run_id[3].lower() == 'e' or run_id[3].lower() == 'q'
            ) and date_obs is not None:
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
                        f'REL_DATE not in header. Derive from RUNID {run_id}.'
                    )
                    semester = mc.to_int(run_id[0:2])
                    rel_year = 2000 + semester + 1
                    if run_id[2] == 'A':
                        result = f'{rel_year}-08-31T00:00:00'
                    else:
                        rel_year += 1
                        result = f'{rel_year}-02-28T00:00:00'
    return result


def get_sitelle_v_plane_data_release(header):
    # REL_DATE not in header, RUN_ID not in header, derive from DATE-OBS
    result = None
    rel_date = header.get('REL_DATE', header.get('DATE-OBS'))
    if rel_date is not None:
        rel_date_dt = mc.make_time(rel_date)
        # add approximately 13 months
        result = rel_date_dt + timedelta(days=13 * 30)
    return result


def get_polarization_function_val(header):
    lookup = {'I': 1, 'Q': 2, 'U': 3, 'V': 4, 'W': 5}
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
    header, cfht_name, ignore = _decompose_params(params)
    result = ProductType.SCIENCE
    obs_type = get_obs_intent(header)
    if obs_type == ObservationIntentType.CALIBRATION:
        result = ProductType.CALIBRATION
    if cfht_name.suffix in ['g', 'm', 'w', 'y']:
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
    result = None
    pi_name = header.get('PI_NAME')
    if pi_name is not None and 'CFHTLS' in pi_name:
        result = 'CFHTLS'
    else:
        run_id = header.get('RUNID')
        if run_id is not None:
            result = md.cache.get_program(run_id)
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
    header, cfht_name, uri = _decompose_params(params)
    if _has_energy(header):
        # from caom2IngestSitelle.py, l555+
        sitresol = header.get('SITRESOL')
        if sitresol is not None and sitresol > 0.0:
            result = sitresol
        if result is None:
            result = 1.0
            if cfht_name.suffix in ['a', 'c', 'f', 'o', 'x']:
                # from caom2IngestSitelle.py, l596
                crval3 = mc.to_float(header.get('FILTERLB'))
                cdelt3 = mc.to_float(header.get('FILTERBW'))
                if crval3 is not None and cdelt3 is not None:
                    result = crval3 / cdelt3
    return result


def get_sitelle_plane_data_product_type(uri):
    cfht_name = cn.cfht_names.get(uri)
    result = DataProductType.IMAGE
    if cfht_name.is_derived_sitelle:
        result = DataProductType.CUBE
    return result


def get_sitelle_time_refcoord_delta(params):
    header, cfht_name, uri = _decompose_params(params)
    if cfht_name.suffix == 'p':
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
    header, cfht_name, uri = _decompose_params(params)
    if cfht_name.suffix in ['a', 'c', 'd', 'f', 'o', 'r', 'x']:
        result = header.get('EXPTIME')
    elif cfht_name.suffix == 'p':
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
    header, cfht_name, uri = _decompose_params(params)
    result = get_spirou_time_refcoord_delta(params)
    if cfht_name.suffix == 'r':
        result = result * (24.0 * 3600.0)
    else:
        result = get_spirou_exptime(params)
    if result is None:
        logging.warning(f'No Time WCS resolution value for {uri}.')
    return result


def get_spirou_time_refcoord_delta(params):
    # caom2IngestSpirou.py, l530+
    header, cfht_name, uri = _decompose_params(params)
    result = None
    if cfht_name.suffix == 'r':
        temp = header.get('FRMTIME')
    elif cfht_name.suffix == 'p':
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
    header, cfht_name, uri = _decompose_params(params)
    result = 1.0
    if cfht_name.suffix == 'r':
        result = header.get('NREADS')
    if result is None:
        logging.warning(f'No Time WCS refcoord.naxis value for {uri}.')
    return result


def get_target_position_cval1(params):
    header, cfht_name, uri = _decompose_params(params)
    ra, ignore_dec = _get_ra_dec(header)
    if ra is None:
        if cfht_name.instrument is md.Inst.ESPADONS:
            ra = header.get('RA_DEG')
    return ra


def get_target_position_cval2(params):
    header, cfht_name, uri = _decompose_params(params)
    ignore_ra, dec = _get_ra_dec(header)
    if dec is None:
        if cfht_name.instrument is md.Inst.ESPADONS:
            dec = header.get('DEC_DEG')
    return dec


def get_target_standard(params):
    header, cfht_name, uri = _decompose_params(params)
    obs_type = _get_obstype(header)
    run_id = _get_run_id(header)
    result = None
    if run_id is not None:
        run_id_type = run_id[3].lower()
        if run_id_type == 'q' and obs_type == 'OBJECT':
            obj_name = header.get('OBJECT').lower()
            if cfht_name.instrument is md.Inst.SITELLE:
                if 'std' in obj_name:
                    result = True
                else:
                    result = False
            else:
                if (
                    'flat' in obj_name
                    or 'focus' in obj_name
                    or 'zenith' in obj_name
                ):
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
        mjd_end = mjd_end.value
        exp_time = mjd_end - mjd_obs
    return exp_time


def get_time_refcoord_delta_simple(params):
    # caom2IngestMegacam.py
    exp_time = get_exptime(params)
    if exp_time is None:
        header, cfht_name, uri = _decompose_params(params)
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
        logging.warning(
            f'Chunk.time.axis.function.refCoord.val is None for '
            f'{_get_filename(header)}'
        )
    else:
        mjd_obs = mjd_obs.value
    return mjd_obs


def get_time_refcoord_val_simple(header):
    result = _get_mjd_obs(header)
    if result is None:
        temp = header.get('DATE-OBS')
        if temp is None:
            # from caom2IngestMegacam.py, l549
            temp = header.get('DATE')
        result = ac.get_datetime(temp)
        result = result.value
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
    header, cfht_name, uri = _decompose_params(params)
    result = _get_obstype(header)
    # caom2IngestWircamdetrend.py, l369
    if 'weight' in uri:
        result = 'WEIGHT'
    elif 'badpix' in uri or 'hotpix' in uri or 'deadpix' in uri:
        result = 'BPM'
    elif cfht_name.suffix == 'g' and result is None:
        result = 'GUIDE'
    return result


def get_wircam_provenance_keywords(uri):
    suffix = cn.cfht_names.get(uri).suffix
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
    return header, cn.cfht_names.get(uri), uri


def _get_filename(header):
    return header.get('FILENAME')


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
            if time_obs is None:
                mjd_obs = date_obs
            else:
                mjd_obs = date_obs + time_obs
        mjd_obs = mjd_obs.value
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
        if obj_ra_dec == 'gappt' or obj_ra_dec == 'null':
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
            logging.warning(
                f'Setting RUNID to default 17BE for {header.get("FILENAME")}.'
            )
            run_id = '17BE'
        else:
            run_id = run_id.strip()
    return run_id


def _get_types(params):
    header, cfht_name, ignore = _decompose_params(params)
    dp_result = DataProductType.IMAGE
    pt_result = ProductType.SCIENCE
    obs_type = get_obs_intent(header)
    if obs_type == ObservationIntentType.CALIBRATION:
        pt_result = ProductType.CALIBRATION
    if cfht_name.suffix in ['m', 'w', 'y']:
        dp_result = DataProductType.AUXILIARY
        pt_result = ProductType.AUXILIARY
    return dp_result, pt_result


def _has_energy(header):
    obs_type = _get_obstype(header)
    # from conversation with CW, SF
    # also from caom2IngestEspadons.py, l393, despite an existing example
    # with energy information
    return obs_type not in ['BIAS', 'DARK']


def _is_espadons_energy(cfht_name):
    result = False
    if cfht_name.instrument is md.Inst.ESPADONS:
        if cfht_name.suffix in ['a', 'b', 'c', 'd', 'f', 'i', 'o', 'p', 'x']:
            result = True
    return result


def _is_sitelle_energy(cfht_name):
    result = False
    if cfht_name.instrument is md.Inst.SITELLE:
        if cfht_name.suffix in ['a', 'c', 'f', 'o', 'p', 'x']:
            result = True
    return result


def _get_filter_md(instrument, filter_name, entry):
    filter_md = md.filter_cache.get_svo_filter(instrument, filter_name)
    if not md.filter_cache.is_cached(instrument, filter_name):
        # want to stop ingestion if the filter name is not expected
        raise mc.CadcException(
            f'Could not find filter metadata for {filter_name} in {entry}.'
        )
    # CW - 15-05-20
    # some flats like this have filter names like ‘i’, instead of ‘i.MP9701’.
    # Even though there may only be a few of these, we want zero so they don’t
    # appear in the filters picklist and confuse users. If the header doesn’t
    # have the full filter name maybe you can hack it based on your knowledge
    # of which i filter was used during this era.
    #
    # SGo - hence the reverse lookup of FILTER_REPAIR CACHE
    updated_filter_name = mc.reverse_lookup(
        filter_name, md.cache.get_from(md.FILTER_REPAIR_CACHE)
    )
    if updated_filter_name is None:
        updated_filter_name = filter_name
    return filter_md, updated_filter_name


def _get_gaia_target_id(header):
    catalog_id = header.get('GAIAID')
    # catalog id looks like:
    # GAIAID  = 'Gaia DR2 470826482635704064'
    # should look like:
    # Gaia:DRX/SOURCE_ID
    # from JJK: 25-02-21
    # PD: 01-03-21
    # should be wary of upper-case letters in schemes, so use 'gaia'
    result = None
    if catalog_id is not None:
        if isinstance(catalog_id, int):
            catalog_dr = header.get('GAIADR')
            bits = catalog_dr.split()
            if len(bits) == 2:
                result = mc.build_uri(
                    scheme=bits[0].lower(),
                    archive=bits[1],
                    file_name=str(catalog_id),
                )
            else:
                logging.warning(f'Unexpected GAIADR value {catalog_dr}.')
        else:
            bits = catalog_id.split()
            if len(bits) == 3:
                result = mc.build_uri(
                    scheme=bits[0].lower(),
                    archive=bits[1],
                    file_name=bits[2],
                )
            else:
                logging.warning(f'Unexpected GAIAID value {catalog_id}.')
    return result


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
                f'Treating {obs_id} with filetype {file_type} as derived. '
            )
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


def _update_observation_metadata(obs, headers, cfht_name, fqn, uri, subject):
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
    logging.debug(
        f'Begin _update_observation_metadata for {cfht_name.file_name}'
    )
    idx = 0
    run_id = headers[0].get('RUNID')
    if run_id is None:
        run_id = headers[0].get('CRUNID')
        # xor
        if (
            run_id is None
            and not (
                cfht_name.instrument is md.Inst.SPIROU
                and cfht_name.suffix == 'g'
            )
        ) or (
            run_id is not None
            and (
                cfht_name.instrument is md.Inst.SPIROU
                and cfht_name.suffix == 'g'
            )
        ):
            if len(headers) > 1:
                idx = 1

                logging.warning(
                    f'Resetting the header/blueprint relationship for '
                    f'{cfht_name.file_name} in {obs.observation_id}'
                )
                if fqn is not None:
                    uri = cfht_name.file_uri
                    # this is the fits2caom2 implementation, which returns
                    # a list structure
                    unmodified_headers = get_cadc_headers(
                        f'file://{fqn}', subject=subject
                    )
                elif uri is not None:
                    # this is the fits2caom2 implementation, which returns
                    # a list structure
                    unmodified_headers = get_cadc_headers(uri, subject=subject)
                else:
                    raise mc.CadcException(
                        f'Cannot retrieve un-modified headers for {uri}'
                    )
                module = importlib.import_module(__name__)
                bp = ObsBlueprint(module=module)
                accumulate_bp(bp, uri)
                # TODO this is not a long-term implementation
                # re-read the headers from disk, because the first pass
                # through caom2gen will have modified the header content based
                # on the original blueprint the execution goes through and
                # sets the CDELT4 value, then from that sets the CD4_4 value.
                # Then the second time through, the CD4_4 value update is
                # expressly NOT done, because the CD4_4 value is already set
                # from the first pass-through - need to figure out a way to
                # fix this .... sigh
                tc.add_headers_to_obs_by_blueprint(
                    obs, unmodified_headers[1:], bp, uri, cfht_name.product_id
                )
            else:
                logging.debug(
                    f'Cannot reset the header/blueprint relationship for '
                    f'{cfht_name.file_name} in {obs.observation_id}'
                )

    logging.debug(f'End _update_observation_metadata.')
    return idx


def _update_sitelle_plane(observation, uri):
    logging.debug(
        f'Begin _update_sitelle_plane for {observation.observation_id}'
    )
    # if the 'p' plane exists, the observation id is the same as the plane id,
    # so copy the metadata to the 'z' plane
    z_plane_key = observation.observation_id.replace('p', 'z')
    temp_z_uri = uri.replace('p', 'z', 1)
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
            temp = uri.replace('.hdf5', '.fits.fz')
            if temp.count('z') == 1:
                # uri looks like: ad:CFHT/2384125p.fits.fz
                p_artifact_key = temp
            else:
                p_artifact_key = temp.replace('z', 'p', 1)
            if p_artifact_key not in p_plane.artifacts.keys():
                p_artifact_key = uri.replace('z', 'p', 1).replace(
                    '.hdf5', '.fits'
                )
                if p_artifact_key not in p_plane.artifacts.keys():
                    p_artifact_key = uri.replace('z', 'p', 1).replace(
                        '.hdf5', '.fits.gz'
                    )
                    if p_artifact_key not in p_plane.artifacts.keys():
                        p_artifact_key = uri.replace('z', 'p', 1).replace(
                            '.hdf5', '.fits.header'
                        )
                        if p_artifact_key not in p_plane.artifacts.keys():
                            raise mc.CadcException(
                                f'Unexpected extension name pattern for '
                                f'artifact URI {p_artifact_key} in '
                                f'{observation.observation_id}.'
                            )
            features = mc.Features()
            features.supports_latest_caom = True
            for part in p_plane.artifacts[p_artifact_key].parts.values():
                z_plane.artifacts[z_artifact_key].parts.add(cc.copy_part(part))
                for chunk in part.chunks:
                    z_plane.artifacts[z_artifact_key].parts[
                        part.name
                    ].chunks.append(cc.copy_chunk(chunk, features))
            z_plane.artifacts[
                z_artifact_key
            ].meta_producer = p_plane.artifacts[p_artifact_key].meta_producer
            z_plane.provenance = p_plane.provenance
            z_plane.calibration_level = p_plane.calibration_level
            z_plane.data_product_type = p_plane.data_product_type
            z_plane.data_release = p_plane.data_release
            z_plane.meta_producer = p_plane.meta_producer
            z_plane.meta_release = p_plane.meta_release

    logging.debug('End _update_sitelle_plane')


def _update_wircam_plane(observation, cfht_name):
    logging.debug(
        f'Begin _update_wircam_plane for {observation.observation_id}'
    )
    # for some 'y' files, that don't have enough metadata on their own,
    # if the 'p' plane exists, and the 'y' plane exists,
    # copy the metadata to the 'y' plane, because the 'y' file
    # will not have enough metadata to fill these things in alone
    copy_to_key = cfht_name.product_id
    copy_to_artifact_key = cfht_name.file_uri
    if cfht_name.suffix == 'p':
        copy_to_key = cfht_name.product_id.replace('p', 'y')
        copy_to_artifact_key = cfht_name.file_uri.replace('p', 'y', 1)

    copy_from_key = cfht_name.product_id
    copy_from_artifact_key = cfht_name.file_uri
    if cfht_name.suffix == 'y':
        copy_from_key = cfht_name.product_id.replace('y', 'p')
        copy_from_artifact_key = cfht_name.file_uri.replace('y', 'p', 1)

    if (
        copy_to_key in observation.planes.keys()
        and copy_from_key in observation.planes.keys()
    ):
        copy_to_plane = observation.planes[copy_to_key]
        copy_from_plane = observation.planes[copy_from_key]
        if (
            copy_from_artifact_key in copy_from_plane.artifacts.keys() and
            copy_to_artifact_key in copy_to_plane.artifacts.keys()
        ):
            copy_from_artifact = copy_from_plane.artifacts[
                copy_from_artifact_key
            ]
            copy_to_artifact = copy_to_plane.artifacts[copy_to_artifact_key]
            if copy_from_plane.provenance is not None:
                copy_to_plane.provenance = cc.copy_provenance(
                    copy_from_plane.provenance
                )
                # set to None, because caom2IngestWircam.py sets only for
                # 'p', 's' files: l1064, l1092
                while len(copy_to_plane.provenance.keywords) > 0:
                    copy_to_plane.provenance.keywords.pop()
            _semi_deep_copy_plane(
                copy_from_plane,
                copy_to_plane,
                copy_from_artifact,
                copy_to_artifact,
            )
    logging.debug('End _update_wircam_plane')


def _semi_deep_copy_plane(from_plane, to_plane, from_artifact, to_artifact):
    to_plane.calibration_level = from_plane.calibration_level
    to_plane.data_product_type = from_plane.data_product_type
    to_plane.data_release = from_plane.data_release
    to_plane.meta_producer = from_plane.meta_producer
    to_plane.meta_release = from_plane.meta_release
    for part in from_artifact.parts.values():
        to_artifact.parts.add(cc.copy_part(part))
        for chunk in part.chunks:
            to_artifact.parts[part.name].chunks.append(cc.copy_chunk(chunk))


def _identify_instrument(uri, cfht_builder):
    logging.debug(f'Begin _identify_instrument for uri {uri}.')
    storage_name = cfht_builder.build(mc.CaomName(uri).file_name)
    logging.debug(
        f'End _identify_instrument for uri {uri} instrument is '
        f'{storage_name.instrument}.'
    )
    return storage_name.instrument


def _build_blueprints(storage_names):
    """This application relies on the caom2utils fits2caom2 ObsBlueprint
    definition for mapping FITS file values to CAOM model element
    attributes. This method builds the DRAO-ST blueprint for a single
    artifact.

    The blueprint handles the mapping of values with cardinality of 1:1
    between the blueprint entries and the model attributes.

    :param storage_names list of CFHTName instances for the files to be
        processed."""
    module = importlib.import_module(__name__)
    blueprints = {}
    for storage_name in storage_names:
        blueprint = ObsBlueprint(module=module)
        if not mc.StorageName.is_preview(
            storage_name.file_uri
        ) and not mc.StorageName.is_hdf5(storage_name.file_uri):
            accumulate_bp(blueprint, storage_name.file_uri)
        blueprints[storage_name.file_uri] = blueprint
    return blueprints


def _get_storage_names(args, cfht_builder):
    result = []
    if args.local:
        for ii in args.local:
            storage_name = cfht_builder.build(ii)
            result.append(storage_name)
    elif args.lineage:
        for ii in args.lineage:
            ignore_product_id, uri = mc.decompose_lineage(ii)
            ignore_scheme, ignore_path, file_name = mc.decompose_uri(uri)
            storage_name = cfht_builder.build(file_name)
            result.append(storage_name)
    else:
        raise mc.CadcException(f'Could not define uri from these args {args}')
    return result


def _get_keyword(lookup, headers):
    result = headers[0].get(lookup)
    if result is None:
        result = headers[1].get(lookup)
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
    """
    This function is called by pipeline execution. It must have this name.
    """
    args = _cfht_args_parser()
    config = mc.Config()
    config.get_executors()
    clients = client_composable.ClientCollection(config)
    cfht_builder = cb.CFHTBuilder(
        clients.data_client, config.archive, config.use_local_files
    )
    storage_names = _get_storage_names(args, cfht_builder)
    cn.cfht_names = {}
    for storage_name in storage_names:
        cn.cfht_names[storage_name.file_uri] = storage_name
    blueprints = _build_blueprints(storage_names)
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
        logging.debug(tb)
        sys.exit(-1)
