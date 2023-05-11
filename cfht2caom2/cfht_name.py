# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2020.                            (c) 2020.
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

from logging import getLogger
from os.path import basename
from re import match
from urllib.parse import urlparse

from caom2pipe.manage_composable import build_uri, CadcException, StorageName
from cfht2caom2.metadata import Inst


__all__ = ['CFHTName']


# TODO - is this still required?
# minimize the number of times the code tries to figure out which instrument
# generated the file being ingested
#
# declare the global here, so that it survives the importlib.import_module
# done by fits2caom2
# key - Artifact uri
# value - CFHTName instance
cfht_names = {}


class CFHTName(StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support fz files in storage,
    - decompress all .gz files, recompress if it's lossless
    - product id == file id
    - the file_name attribute has ALL the extensions, including compression
      type.

    Decisions to make based on the name:
    - Simple vs Derived
    - Calibration Level (RAW vs CALIBRATED)
    - Data Product Type (image, cube, eventlist, spectrum, timeseries, visibility, measurements, catalog, event, sed)
    - Product Type (SCIENCE, CALIBRATION, AUXILIARY, PREVIEW, THUMBNAIL)
    - Algorithm Name (exposure, scan, master_detrend, polarization)
    - WCS Axes to configure
    - Time WCS Axis calculation
    - Espadons/Mega Energy participation
    - Espadons Polarization participation
    - Espadons/Mega/WIRCam Position resetting
    - Sitelle Derived
    - Sitelle energy definition
    - Sitelle data_release value
    - SPIRou 'r' file handling
    - WIRCam OBSTYPE
    - WIRCam ingestion order
    - WIRCam provenance
    """

    CFHT_NAME_PATTERN = '*'

    def __init__(
        self,
        obs_id=None,
        file_name=None,
        instrument=None,
        source_names=[],
        bitpix=None,
    ):
        self._instrument = Inst(instrument)
        self._file_id = None
        self._file_name = None
        self._obs_id = None
        self._suffix = None
        # make recompression decisions based on bitpix
        self._bitpix = bitpix
        super().__init__(
            obs_id=obs_id,
            file_name=file_name.replace('.header', ''),
            source_names=source_names,
        )
        self._logger = getLogger(self.__class__.__name__)
        self._logger.debug(self)

    def __str__(self):
        return (
            f'\n'
            f'      instrument: {self.instrument}\n'
            f'          bitpix: {self.bitpix}\n'
            f'          obs_id: {self.obs_id}\n'
            f'         file_id: {self.file_id}\n'
            f'       file_name: {self.file_name}\n'
            f'    source_names: {self.source_names}\n'
            f'destination_uris: {self.destination_uris}\n'
            f'        file_uri: {self.file_uri}\n'
            f'      product_id: {self.product_id}\n'
            f'          suffix: {self.suffix}\n'
        )

    def _get_uri(self, file_name, scheme):
        """
        SF - 20-05-22
        it looks like compression will be:
        if bitpix==(-32|-64):
             gunzip file.fits.gz
        else:
             imcopy file.fits.gz file.fits.fz[compress]

        """
        if file_name.endswith('.gz'):
            if (
                self._bitpix is None
                or self._bitpix in [-32, -64]
                # fpack failures for unsupported instruments
                or self._instrument is Inst.UNSUPPORTED
            ):
                # use fpack only if bitpix is known
                result = build_uri(
                    scheme=StorageName.scheme,
                    archive=StorageName.collection,
                    file_name=file_name,
                ).replace('.gz', '')
            else:
                result = build_uri(
                    scheme=StorageName.scheme,
                    archive=StorageName.collection,
                    file_name=file_name,
                ).replace('.gz', '.fz')
                self._logger.debug(
                    f'{self._file_name} will be recompressed with fpack.'
                )
        else:
            result = build_uri(scheme=scheme, archive=StorageName.collection, file_name=file_name)
        # .header is mostly for test execution
        return result.replace('.header', '')

    def is_valid(self):
        return True

    @property
    def bitpix(self):
        return self._bitpix

    @bitpix.setter
    def bitpix(self, value):
        self._bitpix = value

    @property
    def file_uri(self):
        """
        The CADC Storage URI for the file.
        """
        return self._get_uri(self._file_name, StorageName.scheme)

    @property
    def instrument(self):
        return self._instrument

    @instrument.setter
    def instrument(self, value):
        self._instrument = value

    @property
    def prev(self):
        """The preview file name for the file."""
        return '{}_preview_1024.jpg'.format(self.product_id)

    @property
    def sequence_number(self):
        result = None
        # SF 09-01-20
        # *y files are produced from other files, I am guessing the sky
        # subtraction software at CFHT copies the header from one of the
        # exposure and does not update the EXPNUM.
        #
        # SGo - because of this, as a secondary measure, try the file name for the sequence number
        temp = match('^[0-9]{5,7}', self._file_name)
        if temp:
            result = self._file_name[:temp.end()]
        return result

    @property
    def thumb(self):
        """The thumbnail file name for the file."""
        return '{}_preview_256.jpg'.format(self.product_id)

    @property
    def zoom(self):
        """The zoom preview file name for the file."""
        return '{}_preview_zoom_1024.jpg'.format(self.product_id)

    @property
    def zoom_uri(self):
        """The zoom preview URI."""
        return self._get_uri(self.zoom, StorageName.scheme)

    @property
    def has_different_destination_name(self):
        if len(self._source_names) == 0:
            result = (
                basename(self._file_name) !=
                basename(self._destination_uris[0])
            )
        else:
            for index, entry in enumerate(self._destination_uris):
                temp = basename(entry)
                result = temp != basename(self._source_names[index])
        return result

    @property
    def is_feasible(self):
        """
        Executing parts of the pipeline is not feasible for hdf5 files at this
        time - make it possible to know whether or not to try for each
        StorageName instance.
        :return:
        """
        return not StorageName.is_hdf5(self._file_name)

    @property
    def simple(self):
        """
        :return: True if the file should be represented as a SimpleObservation, False otherwise
        """
        mega_suffix_list = ['b', 'd', 'f', 'l', 'o', 'x']
        s = {
            Inst.ESPADONS: ['a', 'b', 'c', 'd', 'f', 'o', 'x'],
            Inst.MEGACAM: mega_suffix_list,
            Inst.MEGAPRIME: mega_suffix_list,
            Inst.SITELLE: ['a', 'b', 'c', 'd', 'f', 'o', 'x'],
            Inst.SPIROU: ['a', 'c', 'd', 'f', 'g', 'o', 'r', 'x'],
            Inst.WIRCAM: ['a', 'd', 'f', 'g', 'm', 'o', 'x', 'w', 'v'],
        }
        if self._suffix is not None and self._suffix in s.get(self._instrument):
            result = True
        else:
            # _flag has inputs/members metadata
            result = '_' in self._file_id and '_flag' not in self._file_id
        return result

    @property
    def derived(self):
        """
        Prefer !simple, as this method is here for testing completeness.

        :return: True if the file should be represented in a DerivedObservation, False otherwise
        """
        d = {
            Inst.ESPADONS: ['i', 'p'],
            # caom2IngestMegacam.py, l142 provenance.inputs = caom:CFHT/%s/%so
            Inst.MEGAPRIME: ['p'],
            Inst.SITELLE: ['p', 'v', 'z'],
            Inst.SPIROU: ['e', 'p', 's', 't', 'v'],
            Inst.WIRCAM: ['p', 's', 'y'],
        }
        # diag is not Derived because it's treated as an Auxiliary file, and has no inputs/members metadata
        if self._suffix is not None and self._suffix in d.get(self._instrument) and '_diag' not in self._file_id:
            result = True
        else:
            result = '.' in self._file_id
        return result

    @property
    def raw_time(self):
        """
        :return: True for those processed file naming patterns with Temporal WCS defined according to raw keywords.
        """
        d = {
            Inst.ESPADONS: ['i'],
            Inst.MEGACAM: ['p', 's'],
            Inst.MEGAPRIME: ['p', 's'],
            Inst.SPIROU: ['e', 's', 't', 'v'],
            Inst.WIRCAM: ['p', 's', 'y'],
        }
        if self.simple:
            result = '_' not in self._file_id
        else:
            result = self._instrument in d and self._suffix is not None and self._suffix in d.get(self._instrument)
        return result

    @property
    def suffix(self):
        return self._suffix

    def set_destination_uris(self):

        def _set_extension(for_entry):
            temp = basename(urlparse(for_entry).path)
            if self._bitpix is None or self._bitpix in [-32, -64]:
                result = self._get_uri(temp, StorageName.scheme).replace('.header', '').replace('.gz', '')
            else:
                result = self._get_uri(temp, StorageName.scheme).replace('.gz', '.fz')
            return result

        if len(self._source_names) == 0:
            self._destination_uris.append(self.file_uri)
        else:
            for entry in self._source_names:
                temp = _set_extension(entry)
                if temp.endswith('.fz'):
                    temp2 = temp.replace('.fz', '')
                    if temp2 in self.destination_uris:
                        self._destination_uris.remove(temp2)
                if temp not in self.destination_uris:
                    self._destination_uris.append(temp)

    def set_file_id(self):
        self._file_id = CFHTName.remove_extensions(self._file_name)
        self._suffix = None
        if self.sequence_number is not None:
            # for file names that have _flag or _diag in them
            self._suffix = self._file_id.split('_')[0][-1]

    def set_obs_id(self):
        self._obs_id = CFHTName.get_obs_id(self._file_id)
        if self._obs_id == self._file_id and self.sequence_number is not None:
            self._obs_id = self.sequence_number

    def set_product_id(self, **kwargs):
        if '_diag' in self._file_name:
            # SF 16-03-23
            # artifact of the *p ones
            self._product_id = self._file_id.split('_')[0]
        else:
            super().set_product_id(**kwargs)

    @staticmethod
    def get_obs_id(file_id):
        """
        This function is used a lot for provenance obs id finding, so make a more-efficient version that avoids
        all the source/destination handling.

        # SF - 14-09-22
        # observation ID should be run ID
        # Laurie Rosseau-Nepton - 11-08-22 - remove the 'p'
        #
        # which over-rides this:
        #
        # SF - slack - 02-04-20
        # - MegaCam - the logic should be probably be 2 planes: p
        # and o for science. - all cfht exposures are sorted by EXPNUM
        # if i understand their data acquisition. b,f,d,x should be 1
        # plane observations. - my assumption is that the b,f,d,x have
        # no reason to have a processed equivalent.
        """
        try:
            ignore = int(file_id[:-1])
            obs_id = file_id[:-1]
        except ValueError as e:
            obs_id = file_id
        return obs_id

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        # ESPaDOnS files have a .gz extension
        # SITELLE has hdf5 files
        return (
            name.replace('.fits', '')
            .replace('.fz', '')
            .replace('.header', '')
            .replace('.gz', '')
            .replace('.hdf5', '')
        )
