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

import logging

from caom2 import ObservationIntentType, ProductType
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc

from cfht2caom2 import metadata as md


__all__ = ['APPLICATION']


APPLICATION = 'cfht2caom2'

# All comments denoted 'CW' are copied from
# ssh://goliaths@gimli3/srv/cadc/git/wcaom2archive.git
# from /cfht2caom2/scripts


def accumulate_bp(bp, cfht_name):
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
    if not (
        cfht_name.suffix == 'p' and cfht_name.instrument == md.Inst.SPIROU
    ):
        # TODO if instrument == SITELLE and suffix == 'p', energy_axis == 3
        # it will remove all the if instrument == sitelle and suffix = 'p'
        # in the get_energy_* functions
        bp.configure_time_axis(3)
    if cfht_name.has_energy:
        bp.configure_energy_axis(4)
    bp.configure_observable_axis(6)

    bp.set('Observation.intent', 'get_obs_intent()')

    meta_producer = mc.get_version(APPLICATION)
    bp.set('Observation.metaProducer', meta_producer)
    bp.set('Observation.metaRelease', 'get_meta_release()')

    bp.set('Observation.sequenceNumber', 'get_obs_sequence_number()')
    bp.set('Observation.type', 'get_obs_type()')

    bp.set_default('Observation.algorithm.name', None)

    bp.set(
        'Observation.environment.elevation',
        'get_environment_elevation()',
    )
    bp.set(
        'Observation.environment.humidity',
        'get_obs_environment_humidity()',
    )

    # title is select title from runid_title where proposal_id = 'runid'
    # this obtained from cache.yml now
    bp.clear('Observation.proposal.id')
    bp.add_fits_attribute('Observation.proposal.id', 'RUNID')
    bp.clear('Observation.proposal.pi')
    bp.add_fits_attribute('Observation.proposal.pi', 'PI_NAME')
    bp.set('Observation.proposal.project', 'get_proposal_project()')
    bp.set('Observation.proposal.title', 'get_proposal_title()')

    bp.set('Observation.instrument.name', cfht_name.instrument.value)
    bp.set(
        'Observation.instrument.keywords',
        'get_instrument_keywords()',
    )

    bp.set('Observation.target.standard', 'get_target_standard()')

    bp.clear('Observation.target_position.coordsys')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'OBJRADEC')
    bp.clear('Observation.target_position.equinox')
    bp.add_fits_attribute('Observation.target_position.equinox', 'OBJEQN')
    bp.add_fits_attribute('Observation.target_position.equinox', 'OBJEQUIN')
    bp.set(
        'Observation.target_position.point.cval1',
        'get_target_position_cval1()',
    )
    bp.set(
        'Observation.target_position.point.cval2',
        'get_target_position_cval2()',
    )

    bp.set('Observation.telescope.name', 'CFHT 3.6m')
    x, y, z = ac.get_geocentric_location('cfht')
    bp.set('Observation.telescope.geoLocationX', x)
    bp.set('Observation.telescope.geoLocationY', y)
    bp.set('Observation.telescope.geoLocationZ', z)

    bp.set('Plane.dataProductType', 'get_plane_data_product_type()')
    bp.set('Plane.calibrationLevel', 'get_calibration_level()')
    bp.set('Plane.dataRelease', 'get_plane_data_release()')
    bp.set('Plane.metaRelease', 'get_meta_release()')

    bp.set(
        'Plane.provenance.lastExecuted',
        'get_provenance_last_executed()',
    )
    bp.set('Plane.metaProducer', meta_producer)
    bp.set_default('Plane.provenance.producer', 'CFHT')
    bp.set('Plane.provenance.project', 'STANDARD PIPELINE')
    bp.clear('Plane.provenance.runID')
    bp.add_fits_attribute('Plane.provenance.runID', 'CRUNID')
    bp.set('Plane.provenance.version', 'get_provenance_version()')

    bp.set('Artifact.metaProducer', meta_producer)
    bp.set('Artifact.productType', 'get_product_type()')
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
        bp.set('Chunk.energy.axis.axis.ctype', 'get_energy_ctype()')
        bp.set('Chunk.energy.axis.axis.cunit', 'get_energy_cunit()')
        bp.set('Chunk.energy.axis.error.rnder', 1.0)
        bp.set('Chunk.energy.axis.error.syser', 1.0)
        bp.set(
            'Chunk.energy.axis.function.delta',
            'get_energy_function_delta()',
        )
        bp.set(
            'Chunk.energy.axis.function.naxis',
            'get_energy_function_naxis()',
        )
        bp.set(
            'Chunk.energy.axis.function.refCoord.pix',
            'get_energy_function_pix()',
        )
        bp.set(
            'Chunk.energy.axis.function.refCoord.val',
            'get_energy_function_val()',
        )
        bp.clear('Chunk.energy.bandpassName')
        bp.add_fits_attribute('Chunk.energy.bandpassName', 'FILTER')
        bp.set(
            'Chunk.energy.resolvingPower',
            'get_energy_resolving_power()',
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
        bp.set('Chunk.time.exposure', 'get_exptime()')
        bp.set('Chunk.time.resolution', 'get_exptime()')
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
            'get_time_refcoord_delta_simple()',
        )
        bp.set(
            'Chunk.time.axis.function.refCoord.val',
            'get_time_refcoord_val_simple()',
        )
    else:
        bp.set(
            'Chunk.time.axis.function.delta',
            'get_time_refcoord_delta_derived()',
        )
        bp.set(
            'Chunk.time.axis.function.refCoord.val',
            'get_time_refcoord_val_derived()',
        )
    bp.set('Chunk.time.axis.function.refCoord.pix', 0.5)

    if cfht_name.instrument is md.Inst.ESPADONS:
        _accumulate_espadons_bp(bp, cfht_name)
    elif cfht_name.instrument in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
        _accumulate_mega_bp(bp, cfht_name)
    elif cfht_name.instrument is md.Inst.SITELLE:
        _accumulate_sitelle_bp(bp, cfht_name)
    elif cfht_name.instrument is md.Inst.SPIROU:
        _accumulate_spirou_bp(bp, cfht_name)
    elif cfht_name.instrument is md.Inst.WIRCAM:
        _accumulate_wircam_bp(bp, cfht_name)

    logging.debug('Done accumulate_bp.')


def _accumulate_espadons_bp(bp, cfht_name):
    """Configure the ESPaDOnS-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_espadons_bp.')

    # bp.set('Observation.target.targetID', '_get_gaia_target_id()')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'RADECSYS')

    bp.set(
        'Plane.provenance.keywords',
        'get_espadons_provenance_keywords()',
    )
    bp.set(
        'Plane.provenance.lastExecuted',
        'get_espadons_provenance_last_executed()',
    )
    bp.set('Plane.provenance.name', 'get_espadons_provenance_name()')
    bp.set(
        'Plane.provenance.project',
        'get_espadons_provenance_project()',
    )
    bp.set(
        'Plane.provenance.reference',
        'get_espadons_provenance_reference()',
    )
    bp.set(
        'Plane.provenance.version',
        'get_espadons_provenance_version()',
    )

    bp.set(
        'Chunk.energy.resolvingPower',
        'get_espadons_energy_resolving_power()',
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
        'get_espadons_time_refcoord_delta()',
    )
    bp.set(
        'Chunk.time.axis.function.refCoord.val',
        'get_espadons_time_refcoord_val()',
    )

    bp.set('Chunk.time.exposure', 'get_espadons_exptime()')
    bp.set('Chunk.time.resolution', 'get_espadons_exptime()')

    if cfht_name.suffix == 'p':
        bp.configure_polarization_axis(6)
        # caom2IngestEspadons.py, l209, lTODO
        bp.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
        bp.set('Chunk.polarization.axis.function.delta', 1)
        bp.set('Chunk.polarization.axis.function.naxis', 1)
        bp.set('Chunk.polarization.axis.function.refCoord.pix', 1)
        bp.set(
            'Chunk.polarization.axis.function.refCoord.val',
            'get_polarization_function_val()',
        )

    logging.debug('Done _accumulate_espadons_bp.')


def _accumulate_mega_bp(bp, cfht_name):
    """Configure the MegaCam/MegaPrime-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_mega_bp.')

    bp.set(
        'Plane.provenance.lastExecuted',
        'get_mega_provenance_last_executed()',
    )
    bp.set_default('Plane.provenance.name', 'ELIXIR')
    bp.set_default(
        'Plane.provenance.reference',
        'http://www.cfht.hawaii.edu/Instruments/Elixir/',
    )

    logging.debug('Done _accumulate_mega_bp.')


def _accumulate_sitelle_bp(bp, cfht_name):
    """Configure the Sitelle-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_sitelle_bp.')

    if cfht_name.suffix == 'v':
        bp.set('Observation.intent', ObservationIntentType.SCIENCE)
        bp.set('Observation.sequenceNumber', cfht_name.product_id[:-1])
        bp.set('Plane.dataRelease', 'get_sitelle_v_plane_data_release()')
        bp.clear('Plane.provenance.version')
        bp.add_fits_attribute('Plane.provenance.version', 'PROGRAM')
        bp.set('Artifact.productType', ProductType.SCIENCE)
    bp.set('Plane.dataProductType', 'get_sitelle_plane_data_product_type()')
    bp.set_default('Plane.provenance.name', 'ORBS')
    bp.set_default('Plane.provenance.reference', 'http://ascl.net/1409.007')

    bp.set(
        'Chunk.energy.resolvingPower',
        'get_sitelle_energy_resolving_power()',
    )

    bp.set(
        'Chunk.time.axis.function.delta',
        'get_sitelle_time_refcoord_delta()',
    )
    bp.set('Chunk.time.axis.function.refCoord.val', '_get_mjd_start()')

    logging.debug('End _accumulate_sitelle_bp.')


def _accumulate_spirou_bp(bp, cfht_name):
    """Configure the SPIRou-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    bp.set('Observation.target.targetID', '_get_gaia_target_id()')
    bp.add_fits_attribute('Observation.target_position.coordsys', 'RADECSYS')

    if cfht_name.suffix == 'r':
        pass
    elif cfht_name.suffix in ['a', 'c', 'd', 'f', 'o', 'x']:
        bp.set('Plane.provenance.name', 'get_spirou_provenance_name()')
        bp.set(
            'Plane.provenance.reference',
            'http://www.cfht.hawaii.edu/Instruments/SPIRou/',
        )
        bp.set(
            'Plane.provenance.version',
            'get_spirou_provenance_version()',
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
        'get_ra_deg_from_0th_header()',
    )
    bp.set('Chunk.position.axis.function.refCoord.coord2.pix', 1.0)
    bp.set(
        'Chunk.position.axis.function.refCoord.coord2.val',
        'get_dec_deg_from_0th_header()',
    )
    bp.set('Chunk.position.axis.function.cd11', -0.00035833)
    bp.set('Chunk.position.axis.function.cd12', 0.0)
    bp.set('Chunk.position.axis.function.cd21', 0.0)
    bp.set('Chunk.position.axis.function.cd22', 0.00035833)

    bp.set(
        'Chunk.position.coordsys',
        'get_position_coordsys_from_0th_header()',
    )
    bp.set(
        'Chunk.position.equinox',
        'get_position_equinox_from_0th_header()',
    )

    if cfht_name.suffix not in ['g', 'p']:
        bp.set(
            'Chunk.time.axis.function.delta',
            'get_spirou_time_refcoord_delta()',
        )
        bp.set(
            'Chunk.time.axis.function.naxis',
            'get_spirou_time_refcoord_naxis()',
        )
        bp.set('Chunk.time.exposure', 'get_spirou_exptime()')
        bp.set('Chunk.time.resolution', 'get_spirou_resolution()')

    if cfht_name.suffix == 'p':
        bp.configure_polarization_axis(7)
        bp.set('Chunk.polarization.axis.axis.ctype', 'STOKES')
        bp.set('Chunk.polarization.axis.function.naxis', 1)
        bp.set('Chunk.polarization.axis.function.delta', 1.0)
        bp.set('Chunk.polarization.axis.function.refCoord.pix', 1.0)


def _accumulate_wircam_bp(bp, cfht_name):
    """Configure the WIRCam-specific ObsBlueprint at the CAOM model
    Observation level.
    """
    logging.debug('Begin _accumulate_wircam_bp.')
    bp.set('Observation.type', 'get_wircam_obs_type()')

    bp.set('Plane.provenance.keywords', 'get_wircam_provenance_keywords()')
    bp.set_default('Plane.provenance.name', 'IIWI')
    bp.set_default(
        'Plane.provenance.reference',
        'http://www.cfht.hawaii.edu/Instruments/Imaging/WIRCam',
    )

    bp.set('Chunk.energy.bandpassName', 'get_wircam_bandpass_name(header)')

    logging.debug('Done _accumulate_wircam_bp.')


def _get_keyword(lookup, headers):
    result = headers[0].get(lookup)
    if result is None:
        result = headers[1].get(lookup)
    return result
