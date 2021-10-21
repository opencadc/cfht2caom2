# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
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
#  : 4 $
#
# ***********************************************************************
#

import copy
import logging
import math

from caom2 import Axis, Slice, ObservableAxis, Chunk, DataProductType
from caom2 import CoordAxis2D, CoordRange2D, RefCoord, SpatialWCS, Coord2D
from caom2 import TemporalWCS, CoordAxis1D, CoordFunction1D, CoordError
from caom2utils.fits2caom2 import WcsParser
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import metadata as md

__all__ = ['instrument_factory', 'InstrumentType']


class InstrumentType:
    def __init__(
        self,
        name,
        cfht_name,
        observation,
        headers,
        extension,
    ):
        self._name = name
        self._cfht_name = cfht_name
        self._chunk = None
        self._obs_id = observation.observation_id
        self._observation = observation
        self._observation_intent = observation.intent
        self._observation_type = observation.type
        self._part = None
        self._plane = None
        self._headers = headers
        self._extension = extension
        self._logger = logging.getLogger(name.value)

    @property
    def chunk(self):
        return self._chunk

    @chunk.setter
    def chunk(self, value):
        self._chunk = value

    @property
    def part(self):
        return self._part

    @part.setter
    def part(self, value):
        self._part = value

    @property
    def plane(self):
        return self._plane

    @plane.setter
    def plane(self, value):
        self._plane = value

    def clean_up_energy(self):
        if self._chunk.energy is not None:
            if self._chunk.energy.bandpass_name in ['NONE', 'Open']:
                # CW
                # If no or "open" filter then set filter name to
                # null
                self._chunk.energy.bandpass_name = None

            if (
                self._chunk.energy.axis is not None
                and self._chunk.energy.axis.axis is not None
                and self._chunk.energy.axis.axis.ctype is not None
            ):
                # PD 08-04-20
                # the correct way to express "inverse meter" is
                # either  m**-1 or m^-1
                #
                # we support both exponentiations but convert ^
                # to ** so I guess at that time we thought ** was
                # the more common style.
                if self._chunk.energy.axis.axis.cunit == '1 / m':
                    self._chunk.energy.axis.axis.cunit = 'm**-1'

    def get_filter_md(self, filter_name):
        filter_md = md.filter_cache.get_svo_filter(self._name, filter_name)
        if not md.filter_cache.is_cached(self._name, filter_name):
            # want to stop ingestion if the filter name is not expected
            raise mc.CadcException(
                f'Could not find filter metadata for {filter_name} in '
                f'{self._cfht_name.file_uri}.'
            )
        # CW - 15-05-20
        # some flats like this have filter names like ‘i’, instead of
        # ‘i.MP9701’. Even though there may only be a few of these, we want
        # zero so they don’t appear in the filters picklist and confuse users.
        # If the header doesn’t have the full filter name maybe you can hack
        # it based on your knowledge of which i filter was used during this
        # era.
        #
        # SGo - hence the reverse lookup of FILTER_REPAIR CACHE
        updated_filter_name = mc.reverse_lookup(
            filter_name, md.cache.get_from(md.FILTER_REPAIR_CACHE)
        )
        if updated_filter_name is None:
            updated_filter_name = filter_name
        return filter_md, updated_filter_name

    def make_axes_consistent(self):
        pass

    def reset_energy(self):
        pass

    def reset_position(self):
        pass

    def update_chunk(self):
        self.clean_up_energy()
        self.update_observable()
        self.update_polarization()
        self.update_time()
        self.update_position()
        self.update_energy()
        self.reset_energy()
        self.reset_position()
        self.make_axes_consistent()

    def update_energy(self):
        pass

    def update_observable(self):
        pass

    def update_observation(self):
        if (
            self._observation.proposal is not None
            and self._observation.proposal.pi_name is None
        ):
            self._observation.proposal.pi_name = 'CFHT'

    def update_plane(self):
        if self._observation.algorithm.name == 'scan':
            self.plane.data_product_type = DataProductType.CUBE
            if self.plane.provenance is not None:
                self.plane.provenance.last_executed = mc.make_time(
                    self._headers[self._extension].get('DATE'),
                )

    def _update_plane_provenance(self):
        self._logger.debug(
            f'Begin _update_plane_provenance for {self._obs_id}'
        )
        obs_uri_ignore, plane_uri = cc.make_plane_uri(
            self._obs_id, f'{self._obs_id}o', cn.COLLECTION
        )
        # obs_member_str = mc.CaomName.make_obs_uri_from_obs_id(
        #     cn.COLLECTION, obs_id
        # )
        # obs_member = ObservationURI(obs_member_str)
        # plane_uri = PlaneURI.get_plane_uri(
        #     obs_member, f'{self._obs_id}{self._cfht_name.suffix}')
        self.plane.provenance.inputs.add(plane_uri)
        self._logger.debug(f'End _update_plane_provenance for {self._obs_id}')

    def update_polarization(self):
        pass

    def update_position(self):
        pass

    def update_time(self):
        pass


class Espadons(InstrumentType):
    def __init__(
        self,
        headers,
        extension,
        cfht_name,
        observation,
    ):
        super(Espadons, self).__init__(
            md.Inst.ESPADONS,
            cfht_name,
            observation,
            headers,
            extension,
        )

    def _get_espadons_energy_resolving_power(self):
        instmode = self._headers[self._extension].get('INSTMODE')
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

    def make_axes_consistent(self):
        if not (
            self._chunk.naxis is not None and self._chunk.position is None
        ):
            if self._chunk.energy is not None:
                self._chunk.energy_axis = None
            if self._chunk.time is not None:
                self._chunk.time_axis = None
        if self._chunk.naxis is None:
            if self._chunk.observable_axis is not None:
                self._chunk.observable_axis = None
            if self._chunk.polarization_axis is not None:
                self._chunk.polarization_axis = None

    def reset_position(self):
        if self._chunk.position is not None:
            # conform to stricter WCS validation
            self._chunk.position_axis_1 = None
            self._chunk.position_axis_2 = None
            self._chunk.naxis = None
            # CW - Ignore position wcs if a calibration file
            # suffix list from caom2IngestEspadons.py, l389
            # 'b', 'd', 'c', 'f', 'x'
            # with missing spatial indicator keywords
            radecsys = self._headers[self._extension].get('RADECSYS')
            if not (
                self._cfht_name.suffix in ['a', 'i', 'o', 'p']
                and (radecsys is None or radecsys.lower() != 'null')
                and self._headers[self._extension].get('RA_DEG') is not None
                and self._headers[self._extension].get('DEC_DEG') is not None
            ):
                cc.reset_position(self._chunk)

    def update_energy(self):
        self._logger.debug(f'Begin _update_energy_espadons for {self._obs_id}')
        if (
            self._cfht_name.suffix
            in [
                'b',
                'd',
                'i',
                'p',
            ]
            or self._observation_type in ['BIAS', 'DARK']
        ):
            # i, p, are done in the espadons energy data visitor, and b, d are
            # not done at all
            # CW caom2IngestEspadons.py, l393
            # Ignore energy wcs if some type of calibration file
            return

        if self._cfht_name.suffix in ['a', 'c', 'f', 'o', 'x']:
            from caom2 import Axis, RefCoord, CoordRange1D, CoordAxis1D
            from caom2 import SpectralWCS

            # caom2IngestEspadons.py, l818
            axis = Axis('WAVE', 'nm')
            # caom2IngestEspadons.py l636
            naxis1 = 213542
            # caom2IngestEspadons.py l639
            cdelt1 = 0.0031764
            # caom2IngestEspadons.py l638
            crval1 = 370.0
            ref_coord_1 = RefCoord(0.5, crval1)
            ref_coord_2 = RefCoord(1.5, crval1 + float(naxis1) * cdelt1)
            coord_range = CoordRange1D(ref_coord_1, ref_coord_2)
            coord_axis = CoordAxis1D(axis=axis, range=coord_range)
            resolving_power = self._get_espadons_energy_resolving_power()
            self._chunk.energy = SpectralWCS(
                coord_axis,
                specsys='TOPOCENT',
                ssyssrc='TOPOCENT',
                resolving_power=resolving_power,
            )
            self._chunk.energy_axis = 1
        self._logger.debug(f'End _update_energy_espadons for {self._obs_id}')

    def update_observable(self):
        self._logger.debug(f'Begin _update_observable for {self._obs_id}')
        if self._cfht_name.suffix in ['i', 'p']:
            # caom2IngestEspadons.py, l828
            # CW Set up observable axes, row 1 is wavelength, row 2 is
            # normalized flux, row 3 ('p' only) is Stokes spectrum

            # this check is here, because it's quite difficult to find the
            # 'right' chunk, and blind updating generally causes both chunks
            # to have the same metadata values.
            if self._chunk.observable is None:
                independent_axis = Axis('WAVE', 'nm')
                independent = Slice(independent_axis, 1)
                dependent_axis = Axis('flux', 'counts')
                dependent = Slice(dependent_axis, 2)
                self._chunk.observable = ObservableAxis(dependent, independent)
                self._chunk.observable_axis = 2

                if (
                    self._cfht_name.suffix == 'p'
                    and len(self.part.chunks) == 1
                ):
                    # caom2IngestEspadons.py, l863
                    dependent_axis = Axis('polarized flux', 'percent')
                    dependent = Slice(dependent_axis, 3)
                    new_chunk = copy.deepcopy(self._chunk)
                    new_chunk.observable = ObservableAxis(
                        dependent, independent
                    )
                    new_chunk._id = Chunk._gen_id()
                    self.part.chunks.append(new_chunk)
        self._logger.debug(f'End _update_observable for {self._obs_id}')

    def update_observation(self):
        super(Espadons, self).update_observation()
        if self._observation.target_position is not None:
            if self._observation.target_position.equinox is not None:
                self._observation.target_position.equinox = (
                    md.cache.get_repair(
                        'Observation.target_position.equinox',
                        self._observation.target_position.equinox,
                    )
                )
            if self._observation.target_position.coordsys is not None:
                self._observation.target_position.coordsys = (
                    md.cache.get_repair(
                        'Observation.target_position.coordsys',
                        self._observation.target_position.coordsys,
                    )
                )

    def update_plane(self):
        super(Espadons, self).update_plane()
        # caom2IngestEspadons.py, l714
        if self._cfht_name.suffix == 'i':
            self._update_plane_provenance()

    def update_position(self):
        if self._chunk.position is not None:
            if self._chunk.position.equinox is not None:
                self._chunk.position.equinox = md.cache.get_repair(
                    'Chunk.position.equinox',
                    self._chunk.position.equinox,
                )
            if self._chunk.position.coordsys is not None:
                self._chunk.position.coordsys = md.cache.get_repair(
                    'Chunk.position.coordsys',
                    self._chunk.position.coordsys,
                )


class Mega(InstrumentType):
    def __init__(
        self,
        name,
        headers,
        extension,
        cfht_name,
        observation,
    ):
        super(Mega, self).__init__(
            name,
            cfht_name,
            observation,
            headers,
            extension,
        )
        filter_name = self._headers[self._extension].get('FILTER')
        if filter_name is None and len(self._headers) > self._extension + 1:
            filter_name = self._headers[self._extension + 1].get('FILTER')
        self._filter_name = filter_name

    def make_axes_consistent(self):
        # PD - in general, do not assign, unless the wcs
        # metadata is in the FITS header
        if self._chunk.time_axis is not None:
            self._chunk.time_axis = None

    def reset_energy(self):
        # CW
        # Ignore energy wcs if some type of calibration file
        # or filter='None' or 'Open' or there is no filter
        # match
        filter_md, updated_filter_name = self.get_filter_md(self._filter_name)
        if (
            self._filter_name is None
            or self._filter_name in ['Open', 'NONE']
            or ac.FilterMetadataCache.get_fwhm(filter_md) is None
            or self._cfht_name.suffix in ['b', 'l', 'd']
            or self._observation_type in ['DARK']
        ):
            cc.reset_energy(self._chunk)

    def reset_position(self):
        # CW
        # Ignore position wcs if a calibration file (except 'x'
        # type calibration) and/or position info not in header
        # or binned 8x8
        ccdbin = self._headers[self._extension].get('CCBIN1')
        radecsys = self._headers[self._extension].get('RADECSYS')
        ctype1 = self._headers[self._extension].get('CTYPE1')
        if (
            self._cfht_name.suffix in ['b', 'l', 'd', 'f']
            or ccdbin == 8
            or radecsys is None
            or ctype1 is None
            # TODO - figure out if this should be called
            # or self.observation_intent is ObservationIntentType.CALIBRATION
        ):
            cc.reset_position(self._chunk)
            self._chunk.naxis = None

    def update_energy(self):
        # SGo - use range for energy with filter information
        filter_md, updated_filter_name = self.get_filter_md(self._filter_name)
        if not (
            self._filter_name is None
            or self._filter_name in ['Open', 'NONE']
            or ac.FilterMetadataCache.get_fwhm(filter_md) is None
            or self._cfht_name.suffix in ['b', 'l', 'd']
            or self._observation_type in ['DARK']
        ):
            cc.build_chunk_energy_range(
                self._chunk, updated_filter_name, filter_md
            )
            if self._chunk.energy is not None:
                from caom2 import CoordError

                self._chunk.energy.ssysobs = 'TOPOCENT'
                self._chunk.energy.ssyssrc = 'TOPOCENT'
                # values from caom2megacam.default, caom2megacamdetrend.default
                self._chunk.energy.axis.error = CoordError(1.0, 1.0)

    def update_plane(self):
        super(Mega, self).update_plane()
        # caom2IngestMegacam.py, l142
        if self._cfht_name.suffix == 'p':
            self._update_plane_provenance()


class Sitelle(InstrumentType):
    def __init__(
        self,
        headers,
        extension,
        cfht_name,
        observation,
    ):
        super(Sitelle, self).__init__(
            md.Inst.SITELLE,
            cfht_name,
            observation,
            headers,
            extension,
        )

    def make_axes_consistent(self):
        self._chunk.time_axis = None
        if (
            self._cfht_name.suffix in ['a', 'o', 'x']
            and self._chunk.position is None
            and self._chunk.naxis is not None
            and self._chunk.naxis == 3
            and self._chunk.energy is not None
        ):
            self._chunk.energy_axis = 3

        # because SITELLE has the information from two
        # detectors amalgamated into one set of HDUs
        if self._chunk.naxis is not None and self._chunk.naxis <= 2:
            if self._chunk.position_axis_1 is None:
                self._chunk.naxis = None
            self._chunk.energy_axis = None
        else:
            self._chunk.energy_axis = 3

    def reset_energy(self):
        if self._chunk.energy_axis is not None:
            # CW, SF 18-12-19 - SITELLE biases and darks have
            # no energy, all other observation types do
            if self._observation_type in ['BIAS', 'DARK']:
                cc.reset_energy(self._chunk)
            else:
                self._chunk.energy_axis = 3

        if self._cfht_name.suffix == 'v':
            cc.reset_energy(self._chunk)

    def update_position(self):
        if (
            self._cfht_name.suffix in ['a', 'o', 'x']
            and self._chunk.position is None
        ):
            self._update_range_position()

        if self._cfht_name.suffix == 'p':
            if self._chunk.position.axis.function is None:
                self._update_function_position()

    def _update_function_position(self):
        self._logger.debug(
            f'Begin update_position_function for {self._obs_id}'
        )
        header = self._headers[self._extension]
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

        wcs_parser = WcsParser(header, self._obs_id, self._extension)
        if self._chunk is None:
            self._chunk = Chunk()
        wcs_parser.augment_position(self._chunk)
        self._logger.debug(f'End update_function_position for {self._obs_id}')

    def _update_range_position(self):
        self._logger.debug(f'Begin _update_position for {self._obs_id}')
        # from caom2IngestSitelle.py l894
        header = self._headers[self._extension]
        obs_ra = header.get('RA_DEG')
        obs_dec = header.get('DEC_DEG')
        if obs_ra is None or obs_dec is None:
            logging.error(
                'RA_DEG {obs_ra} DEC_DEG {obs_dec} for {obs_id} are not set.'
            )
            return

        pix_scale2 = mc.to_float(header.get('PIXSCAL2'))
        delta_dec = (1024.0 * pix_scale2) / 3600.0
        obs_dec_bl = obs_dec - delta_dec
        obs_dec_tr = obs_dec + delta_dec

        pix_scale1 = mc.to_float(header.get('PIXSCAL1'))
        # obs_dec_tr == obs_dec_tl
        delta_ra_top = (1024.0 * pix_scale1) / (
            3600.0 * (math.cos(obs_dec_tr * 3.14159 / 180.0))
        )
        delta_ra_bot = (1024.0 * pix_scale1) / (
            3600.0 * (math.cos(obs_dec_bl * 3.14159 / 180.0))
        )
        obs_ra_bl = obs_ra + delta_ra_bot
        obs_ra_tr = obs_ra - delta_ra_top

        axis = CoordAxis2D(Axis('RA', 'deg'), Axis('DEC', 'deg'))
        axis.range = CoordRange2D(
            Coord2D(RefCoord(0.5, obs_ra_bl), RefCoord(0.5, obs_dec_bl)),
            Coord2D(RefCoord(2048.5, obs_ra_tr), RefCoord(2048.5, obs_dec_tr)),
        )
        position = SpatialWCS(axis)
        position.coordsys = 'FK5'
        position.equinox = 2000.0
        self._chunk.position = position
        self._chunk.position_axis_1 = 1
        self._chunk.position_axis_2 = 2
        self._logger.debug(f'End _update_position for {self._obs_id}')


class Spirou(InstrumentType):
    def __init__(
        self,
        headers,
        extension,
        cfht_name,
        observation,
    ):
        super(Spirou, self).__init__(
            md.Inst.SPIROU,
            cfht_name,
            observation,
            headers,
            extension,
        )
        self._header = headers[extension]

    def make_axes_consistent(self):
        # stricter WCS validation
        self._chunk.naxis = None
        if self._chunk.energy is not None:
            self._chunk.energy_axis = None
            if self._observation_type == 'DARK':
                # caom2IngestSpirou.py, l514
                self._chunk.energy = None
        if self._chunk.time is not None:
            self._chunk.time_axis = None
        if self._chunk.position is not None:
            self._chunk.position_axis_1 = None
            self._chunk.position_axis_2 = None

    def reset_position(self):
        self._logger.debug(f'Begin reset_position for {self._obs_id}')
        # from caom2IngestSpirou.py, l499+
        # CW - ignore position wcs if a calibration file
        ra_deg = self._header.get('RA_DEG')
        dec_deg = self._header.get('DEC_DEG')
        ra_dec_sys = self._header.get('RADECSYS')
        if self._observation_type not in ['OBJECT', 'ALIGN'] or (
            ra_deg is None
            and dec_deg is None
            and (ra_dec_sys is None or ra_dec_sys.lower() == 'null')
        ):
            cc.reset_position(self._chunk)
        self._logger.debug(f'End reset_position for {self._obs_id}')

    def update_observation(self):
        super(Spirou, self).update_observation()
        if self._cfht_name.suffix != 'r':
            cc.rename_parts(self._observation, self._headers)

    def update_plane(self):
        super(Spirou, self).update_plane()
        # caom2IngestSpirou.py, l584
        if self._cfht_name.suffix in ['e', 's', 't']:
            self._update_plane_provenance()
        elif self._cfht_name.suffix == 'g':
            self.plane.data_product_type = DataProductType.IMAGE

    def update_polarization(self):
        self._logger.debug(f'End update_polarization for {self._obs_id}')
        if self._cfht_name.suffix == 'p':
            header = None
            for h in self._headers:
                if h.get('EXTNAME') == self.part.name:
                    header = h
                    break

            stokes_param = header.get('STOKES')
            if stokes_param is None:
                self._logger.warning(
                    f'No STOKES value for HDU {self.part.name} in '
                    f'{self._obs_id}. No polarization.'
                )
                self._chunk.polarization = None
                self._chunk.polarization_axis = None
            else:
                lookup = {
                    'I': 1.0,
                    'Q': 2.0,
                    'U': 3.0,
                    'V': 4.0,
                    'W': 5.0,
                }
                crval = lookup.get(stokes_param, 0.0)
                if crval == 0.0:
                    self._logger.warning(
                        f'STOKES value is {crval}. No polarization.'
                    )
                    self._chunk.polarization = None
                    self._chunk.polarization_axis = None
                else:
                    if (
                        self._chunk.polarization is not None
                        and self._chunk.polarization.axis is not None
                        and self._chunk.polarization.axis.function is not None
                    ):
                        self._chunk.polarization.axis.function.ref_coord.val = (
                            crval
                        )
            # check with Dustin on what a polarization cut-out
            # looks like before deciding this is semi-ok
            self._chunk.naxis = None
            self._chunk.position_axis_1 = None
            self._chunk.position_axis_2 = None
            self._chunk.energy_axis = None
            self._chunk.time_axis = None
            self._chunk.polarization_axis = None
        self._logger.debug(f'End update_polarization for {self._obs_id}')

    def update_time(self):
        if self._cfht_name.suffix == 'g':
            self._update_time_g()
        elif self._cfht_name.suffix == 'p':
            self._update_time_p()

    def _update_time_g(self):
        self._logger.debug(f'Begin update_time for {self._obs_id}')
        # construct TemporalWCS for 'g' files from the CAOM2 pieces
        # because the structure of 'g' files is so varied, it's not
        # possible to hand over even part of the construction to the
        # blueprint.
        # SF - 22-09-20 - use ETIME

        ref_coord_val = mc.get_keyword(self._headers, 'DATE')
        ref_coord_mjd = ac.get_datetime(ref_coord_val).value

        if self._chunk.time is None:
            self._chunk.time = TemporalWCS(
                CoordAxis1D(Axis('TIME', 'd')), timesys='UTC'
            )
        if self._chunk.time.axis is None:
            self._chunk.time.axis = CoordAxis1D(
                axis=Axis('TIME', 'd'),
                error=None,
                range=None,
                bounds=None,
                function=None,
            )

        time_index = self._headers[0].get('ZNAXIS')
        if time_index is None or time_index == 0:
            time_index = self._headers[1].get('ZNAXIS')
            if time_index is None or time_index == 0:
                time_index = self._headers[0].get('NAXIS')
                if time_index is None or time_index == 0:
                    time_index = self._headers[1].get('NAXIS')
        naxis_key = f'ZNAXIS{time_index}'
        time_naxis = mc.get_keyword(self._headers, naxis_key)
        if time_naxis is None:
            naxis_key = f'NAXIS{time_index}'
            time_naxis = mc.get_keyword(self._headers, naxis_key)

        ref_coord = RefCoord(pix=0.5, val=mc.to_float(ref_coord_mjd))

        e_time = mc.get_keyword(self._headers, 'ETIME')
        if e_time is None:
            self._logger.warning(
                f'No exposure found for {self._cfht_name.file_name}. '
                f'No Temporal WCS.'
            )
        else:
            self._chunk.time.exposure = e_time / 1000.0
            self._chunk.time.resolution = self._chunk.time.exposure
            time_delta = self._chunk.time.exposure / 86400.0
            self._chunk.time.axis.function = CoordFunction1D(
                naxis=time_naxis, delta=time_delta, ref_coord=ref_coord
            )
        self._logger.debug(f'End _update_time_g for {self._obs_id}')

    def _update_time_p(self):
        self._logger.debug(f'Begin _update_time_p for {self._obs_id}')
        # TOTETIME is not in all the HDUs, so copy it from the HDUs that have
        # it, and use it everywhere - this matches existing CFHT SPIRou 'p'
        # behaviour

        def _find_keywords_in_header(header, lookup):
            values = []
            for keyword in header:
                if keyword.startswith(lookup) and len(keyword) > len(lookup):
                    values.append(header.get(keyword))
            return values

        header = None
        for h in self._headers:
            if h.get('EXTNAME') == self.part.name:
                header = h
                break

        tot_e_time = header.get('TOTETIME')

        lower = _find_keywords_in_header(header, 'MJDATE')
        upper = _find_keywords_in_header(header, 'MJDEND')

        if len(lower) > 0:
            for ii, entry in enumerate(lower):
                self._chunk.time = cc.build_temporal_wcs_append_sample(
                    self._chunk.time, entry, upper[ii]
                )

        if self._chunk.time is None:
            self._logger.warning(
                f'No chunk time metadata for {self._obs_id}, part '
                f'{self.part.name}.'
            )
        elif tot_e_time is None:
            self._logger.warning(
                f'Cannot find time metadata for {self._obs_id}, part '
                f'{self.part.name}.'
            )
        else:
            self._chunk.time.exposure = tot_e_time
            self._chunk.time.resolution = tot_e_time
        self._logger.debug(f'End _update_time_p for {self._obs_id}')


class Wircam(InstrumentType):
    def __init__(
        self,
        headers,
        extension,
        cfht_name,
        observation,
    ):
        super(Wircam, self).__init__(
            md.Inst.WIRCAM,
            cfht_name,
            observation,
            headers,
            extension,
        )

    def make_axes_consistent(self):
        # position axis check is to determine if naxis should
        # be set
        if (
            self._cfht_name.suffix in ['d', 'f', 'g', 'o']
            and self._chunk.position is None
        ):
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
            self._chunk.naxis = None
            self._chunk.position_axis_1 = None
            self._chunk.position_axis_2 = None
            self._chunk.energy_axis = None
            self._chunk.time_axis = None

        if self._chunk.naxis == 2:
            self._chunk.time_axis = None
            self._chunk.energy_axis = None

    def reset_energy(self):
        if self._cfht_name.suffix == 'g':
            temp_bandpass_name = self._headers[self._extension].get('FILTER')
            if temp_bandpass_name == 'FakeBlank':
                cc.reset_energy(self._chunk)

    def reset_position(self):
        if (
            self._chunk.position is not None
            and self._chunk.position.coordsys.lower() == 'null'
        ):
            cc.reset_position(self._chunk)
            self._chunk.naxis = None
            self._chunk.energy_axis = None
            self._chunk.time_axis = None

        if self._cfht_name.suffix in ['f'] or self._observation_type in [
            'BPM',
            'DARK',
            'FLAT',
            'WEIGHT',
        ]:
            cc.reset_position(self._chunk)
            self._chunk.naxis = None

    def update_energy(self):
        filter_name = self._headers[self._extension].get('FILTER')
        if filter_name is None and len(self._headers) > self._extension + 1:
            filter_name = self._headers[self._extension + 1].get('FILTER')
        filter_md, updated_filter_name = self.get_filter_md(filter_name)
        cc.build_chunk_energy_range(
            self._chunk, updated_filter_name, filter_md
        )
        if self._chunk.energy is not None:
            from caom2 import CoordError

            self._chunk.energy.ssysobs = 'TOPOCENT'
            self._chunk.energy.ssyssrc = 'TOPOCENT'
            # values from caom2megacam.default, caom2megacamdetrend.default
            self._chunk.energy.axis.error = CoordError(1.0, 1.0)

    # def update_observation(self):
    #     if self._cfht_name.suffix in ['p', 'y']:

    def update_plane(self):
        super(Wircam, self).update_plane()
        # caom2IngestWircam.py, l193
        # CW 09-01-20
        # Only the 'o' is input
        if self._cfht_name.suffix in ['p', 's']:
            self._update_plane_provenance()

    def update_position(self):
        if self._cfht_name.suffix == 'g':
            self._update_position_g()
        elif self._cfht_name.suffix == 'o':
            self._update_position_o()

    def _update_position_g(self):
        """'g' file position handling."""
        self._logger.debug(f'Begin _update_position_g for {self._obs_id}')
        header = self._headers[self._extension]
        obj_name = header.get('OBJNAME')
        part_index = mc.to_int(self.part.name)
        if obj_name == 'zenith' or part_index >= 5:
            # SGo - the second clause is here, because there are only four
            # sets of position information in headers (RA/DEC of guide start
            # on arrays 1 2 3 4), and that's the only thing that is calculated
            # in the original code. Values for HDU 5+ are not written as part
            # of the override file, and thus default to 0, which fails
            # ingestion to the service.
            self._logger.warning(
                f'obj_name is {obj_name}. part_num is {part_index} No '
                f'position for {self._obs_id}'
            )
            cc.reset_position(self._chunk)
            return

        header = self._headers[part_index]
        cd1_1 = None
        cd2_2 = None
        if (
            header.get('CRVAL2') is not None
            or self._headers[0].get('CRVAL2') is not None
        ):
            cd1_1 = mc.to_float(header.get('PIXSCAL1')) / 3600.0
            cd2_2 = mc.to_float(header.get('PIXSCAL2')) / 3600.0

        if cd1_1 is None or cd2_2 is None:
            self._logger.warning(
                f'cd1_1 is {cd1_1}, cd2_2 is {cd2_2}, part_index is '
                f'{part_index}. No position for this part for {self._obs_id}.'
            )
            cc.reset_position(self._chunk)
            return

        wcgd_ra = header.get(f'WCGDRA{self._part.name}')
        if wcgd_ra is None:
            wcgd_ra = self._headers[0].get(f'WCGDRA{self._part.name}')
        wcgd_dec = header.get(f'WCGDDEC{self._part.name}')
        if wcgd_dec is None:
            wcgd_dec = self._headers[0].get(f'WCGDDEC{self._part.name}')
        if wcgd_ra is None or wcgd_dec is None:
            self._logger.warning(
                f'WCGDRA{self._part.name} and WCGDDEC{self._part.name} are '
                f'undefined. No position.'
            )
            cc.reset_position(self._chunk)
            return
        cr_val1, cr_val2 = ac.build_ra_dec_as_deg(
            wcgd_ra, wcgd_dec, frame='fk5'
        )
        if math.isclose(cr_val1, 0.0) and math.isclose(cr_val2, 0.0):
            self._logger.warning(
                f'WCGDRA{self._part.name} and WCGDDEC{self._part.name} are '
                f'close to 0. No position.'
            )
            cc.reset_position(self._chunk)
            return

        naxis_1 = header.get('NAXIS1')
        if naxis_1 is None:
            naxis_1 = header.get('ZNAXIS1')
        naxis_2 = header.get('NAXIS2')
        if naxis_2 is None:
            naxis_2 = header.get('ZNAXIS2')

        if mc.to_float(self._obs_id) < 980000:
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

        wcs_parser = WcsParser(header, self._obs_id, self._extension)
        if self._chunk is None:
            self._chunk = Chunk()
        wcs_parser.augment_position(self._chunk)
        if self._chunk.position is not None:
            self._chunk.naxis = 2
            self._chunk.position_axis_1 = 1
            self._chunk.position_axis_2 = 2
        self._logger.debug(f'End _update_position_g for {self._obs_id}')

    def _update_position_o(self):
        self._logger.debug(f'Begin _update_position_o for {self._obs_id}')
        part_index = mc.to_int(self.part.name)
        header = self._headers[part_index]
        ra_deg = header.get('RA_DEG')
        dec_deg = header.get('DEC_DEG')
        if (
            self._chunk.position is None
            and ra_deg is not None
            and dec_deg is not None
        ):
            self._logger.info(
                f'Adding position information for {self._obs_id}'
            )
            header['CTYPE1'] = 'RA---TAN'
            header['CTYPE2'] = 'DEC--TAN'
            header['CUNIT1'] = 'deg'
            header['CUNIT2'] = 'deg'
            header['CRVAL1'] = ra_deg
            header['CRVAL2'] = dec_deg
            wcs_parser = WcsParser(header, self._obs_id, self._extension)
            if self._chunk is None:
                self._chunk = Chunk()
            wcs_parser.augment_position(self._chunk)
            if self._chunk.position is not None:
                self._chunk.position_axis_1 = 1
                self._chunk.position_axis_2 = 2
        self._logger.debug(f'End _update_position_o')

    def update_time(self):
        self._logger.debug(f'Begin _update_time for {self._obs_id}')
        if self._cfht_name.suffix == 'g':
            # construct TemporalWCS for 'g' files from the CAOM2 pieces
            # because the structure of 'g' files is so varied, it's not
            # possible to hand over even part of the construction to the
            # blueprint.

            # SF - 07-05-20
            # so NAXIS here (where 'here' is a 'g' file) is ZNAXIS=3: time
            # sequence of images of the guiding camera => this means try
            # ZNAXIS* keyword values before trying NAXIS*, hence the header
            # lookup code

            ref_coord_val = mc.get_keyword(self._headers, 'MJD-OBS')
            # ref_coord_val = self._headers[0].get('MJD-OBS')
            # if ref_coord_val is None:
            #     ref_coord_val = self._headers[1].get('MJD-OBS')
            part_index = mc.to_int(self.part.name)
            part_header = self._headers[part_index]

            if self._chunk.time is None:
                self._chunk.time = TemporalWCS(
                    CoordAxis1D(Axis('TIME', 'd')), timesys='UTC'
                )

            if self._chunk.time.axis is None:
                self._chunk.time.axis = CoordAxis1D(
                    axis=Axis('TIME', 'd'),
                    error=None,
                    range=None,
                    bounds=None,
                    function=None,
                )

            if self._chunk.time.axis.error is None:
                self._chunk.time.axis.error = CoordError(
                    rnder=0.0000001, syser=0.0000001
                )

            if self._chunk.time.axis.function is None:
                ref_coord = RefCoord(pix=0.5, val=mc.to_float(ref_coord_val))

                time_index = part_header.get('ZNAXIS')
                if time_index is None:
                    time_index = part_header.get('NAXIS')
                naxis_key = f'ZNAXIS{time_index}'
                time_naxis = part_header.get(naxis_key)
                if time_naxis is None:
                    naxis_key = f'NAXIS{time_index}'
                    time_naxis = part_header.get(naxis_key)

                if (
                    time_naxis is not None
                    and time_index is not None
                    and time_index == 3
                ):
                    # caom2.4 wcs validation conformance
                    self._chunk.time_axis = 3

                # CW
                # Define time samples for guidecube data
                # Guiding time doesn't seem to match up very well, so just
                # say that use MJD-OBS gnaxis3 and WCPERIOD
                #
                # code from caom2IngestWircam.py, l876+
                wcgdra1_0 = self._headers[0].get('WCGDRA1')
                wc_period_0 = self._headers[0].get('WCPERIOD')
                wc_period = None
                if wcgdra1_0 is not None and wc_period_0 is not None:
                    wc_period = wc_period_0
                else:
                    wcgdra1_1 = self._headers[1].get('WCGDRA1')
                    wc_period_1 = self._headers[1].get('WCPERIOD')
                    if wcgdra1_1 is not None and wc_period_1 is not None:
                        wc_period = wc_period_1

                time_delta = None
                if wc_period is not None:
                    if wc_period < 0.0:
                        wc_period = 100.0
                    # caom2IngestWircam.py, l375
                    self._chunk.time.exposure = wc_period / 1000.0
                    self._chunk.time.resolution = self._chunk.time.exposure
                    time_delta = self._chunk.time.exposure / 86400.0
                self._chunk.time.axis.function = CoordFunction1D(
                    naxis=time_naxis, delta=time_delta, ref_coord=ref_coord
                )

        # fits2caom2 prefers ZNAXIS to NAXIS, but the originating scripts
        # seem to prefer NAXIS, so odd as this placement seems, do not rely
        # on function execution, because it affects NAXIS, not ZNAXIS - sigh
        if (
            self._observation_type not in ['BPM', 'DARK', 'FLAT', 'WEIGHT']
            and self._cfht_name.suffix != 'g'
        ):
            n_exp = self._headers[self._extension].get('NEXP')
            if n_exp is not None:
                # caom2IngestWircam.py, l843
                self._chunk.time.axis.function.naxis = mc.to_int(n_exp)

        self._logger.debug(f'End _update_time for {self._obs_id}')


def instrument_factory(
    name,
    headers,
    extension,
    cfht_name,
    observation,
):
    if name is md.Inst.ESPADONS:
        temp = Espadons(headers, extension, cfht_name, observation)
    elif name in [md.Inst.MEGAPRIME, md.Inst.MEGACAM]:
        temp = Mega(name, headers, extension, cfht_name, observation)
    elif name is md.Inst.SITELLE:
        temp = Sitelle(headers, extension, cfht_name, observation)
    elif name is md.Inst.SPIROU:
        temp = Spirou(headers, extension, cfht_name, observation)
    elif name is md.Inst.WIRCAM:
        temp = Wircam(headers, extension, cfht_name, observation)
    else:
        temp = InstrumentType(name, cfht_name, observation, headers, extension)
    return temp
