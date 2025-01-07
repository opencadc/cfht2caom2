# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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

import logging

from os.path import basename, join
from re import match
from urllib.parse import urlparse

from caom2utils.data_util import get_local_file_info, get_local_file_headers
from caom2pipe.execute_composable import CaomExecuteRunnerMeta
from caom2pipe.execute_composable import MetaVisitRunnerMeta, NoFheadStoreVisitRunnerMeta, OrganizeExecutesRunnerMeta
from caom2pipe.execute_composable import NoFheadScrapeRunnerMeta, NoFheadVisitRunnerMeta
from caom2pipe.manage_composable import build_uri, CadcException, get_keyword, StorageName, TaskType
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
        bitpix=None,
        instrument=None,
        source_names=[],
    ):
        self._instrument = Inst(instrument)
        self._suffix = None
        # make recompression decisions based on bitpix
        self._bitpix = bitpix
        self._descriptors = {}
        super().__init__(source_names=source_names)

    def __str__(self):
        f_info_keys = '\n                  '.join(ii for ii in self.file_info.keys())
        metadata_keys = '\n                  '.join(ii for ii in self.metadata.keys())
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
            f'  file_info keys: {f_info_keys}\n'
            f'   metadata keys: {metadata_keys}'
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

    def descriptor(self, key):
        return self._descriptors.get(key)

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

    def set_metadata(self, **kwargs):
        self._instrument = get_instrument(self._metadata.get(self.file_uri), self._file_name)
        if not self.hdf5:
            self._bitpix = get_keyword(self._metadata.get(self.file_uri), 'BITPIX')

    def set_obs_id(self, **kwargs):
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


def get_instrument(headers, entry):
    """
    SF - 15-04-20 - slack - what if there's no INSTRUME or DETECTOR
    keyword?
    something like: if CFHT: if has(INSTRUME) then
    MEGA else if NEXTEND>30 then MEGA else...

    On MegaPrime vs MegaCam values from SVO:
    megaprime == full instrument (optics+camera) ,
    megacam == camera

    SVO has mistakes (they say megacam=CCD+mirror+optics), but i think the
    Megaprime entry has a lot more filters. So to be consistent, we should
    always use the Megaprime entry

    SGo - make the 'always use Megaprime' happen here by setting the
    instrument to MegaPrime.

    :param headers: astropy fits headers
    :param entry: string for error logging
    :return: md.Inst instance
    """
    if StorageName.is_hdf5(entry):
        inst = Inst.SITELLE
    else:
        nextend = None
        detector = None
        instrument = headers[0].get('INSTRUME')
        if instrument is None:
            if len(headers) > 1:
                instrument = headers[1].get('INSTRUME')
            if instrument is None:
                instrument = headers[0].get('DETECTOR')
                if instrument is None:
                    if len(headers) > 1:
                        instrument = headers[1].get('DETECTOR')
                    if instrument is None:
                        nextend = headers[0].get('NEXTEND')
        elif instrument == 'Unknown':
            detector = headers[0].get('DETECTOR')
        if (instrument is None and nextend is not None and nextend > 30) or '_diag' in entry:
            inst = Inst.MEGAPRIME
        else:
            msg = (
                f'Unknown value for instrument {instrument}, detector '
                f'{detector} and nextend {nextend} for {entry}.'
            )

            try:
                inst = Inst(instrument)
            except ValueError:
                if (
                    instrument == 'CFHT MegaPrime'
                    or instrument == 'megacam'
                ):
                    inst = Inst.MEGAPRIME
                elif instrument == 'Unknown' and detector is not None:
                    if detector == 'OLAPA':
                        # SF 24-06-20
                        # ok to hack ESPADONS name for OLAPA
                        inst = Inst.ESPADONS
                    else:
                        try:
                            inst = Inst(detector)
                        except ValueError:
                            # SF 24-06-20
                            # nasty hacks: all the 48 failed espadons
                            # files have a PATHNAME with espadon in it
                            #
                            # so you may add PATHNAME check if everything
                            # else fails
                            pathname = headers[0].get('PATHNAME')
                            if 'espadons' in pathname:
                                inst = Inst.ESPADONS
                            else:
                                raise CadcException(msg)
                else:
                    raise CadcException(msg)
    logging.debug(f'Instrument is {inst}')
    return inst


class CFHTMetaVisitRunnerMeta(MetaVisitRunnerMeta):
    """
    Defines the pipeline step for Collection creation or augmentation by a visitor of metadata into CAOM.
    """

    def __init__(
        self,
        clients,
        config,
        meta_visitors,
        reporter,
    ):
        super().__init__(
            clients=clients,
            config=config,
            meta_visitors=meta_visitors,
            reporter=reporter,
        )

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        for index, source_name in enumerate(self._storage_name.source_names):
            uri = self._storage_name.destination_uris[index]
            if uri not in self._storage_name.file_info:
                self._storage_name.file_info[uri] = self._clients.data_client.info(uri)
            if uri not in self._storage_name.metadata:
                self._storage_name.metadata[uri] = []
                if '.fits' in source_name:
                    self._storage_name._metadata[uri] = self._clients.data_client.get_head(uri)
                elif self._storage_name.hdf5:
                    if uri not in self._storage_name._descriptors:
                        # local import to limit exposure in Docker builds
                        import h5py
                        f_in = h5py.File(source_name)
                        self._storage_name._descriptors[uri] = f_in
                        # Laurie Rosseau-Nepton - 26-04-23
                        # The standard_spectrum is related to flux calibration used on the data. The other one is
                        # for the science data and is the one that should be used.
                        if len(f_in.attrs) >= 20 and 'NAXIS' in f_in.attrs:
                            # 20 is the not-so-random boundary I picked to differentiate between the hdf5 files
                            # that have sufficient metadata to describe WCS.
                            self._logger.debug(f'Found attrs for {source_name}')
                            self._storage_name._metadata[uri] = [f_in.attrs]
                        else:
                            self._logger.warning(f'No attrs for {source_name}.')
                            self._storage_name._metadata[uri] = []

            self._storage_name._instrument = get_instrument(
                self._storage_name.metadata.get(uri), self._storage_name._file_name
            )
            if not self._storage_name.hdf5:
                self._storage_name._bitpix = get_keyword(self._storage_name.metadata.get(uri), 'BITPIX')

        # ensure the destination uris have the correct extensions based on BITPIX values
        self._storage_name.set_destination_uris()
        self._logger.debug('End _set_preconditions')


class CFHTNoFheadStoreVisitRunnerMeta(NoFheadStoreVisitRunnerMeta):

    def __init__(self, clients, config, data_visitors, meta_visitors, reporter, store_transferrer):
        super().__init__(config, clients, store_transferrer, meta_visitors, data_visitors, reporter)

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        #
        # because of decompression, the destination URIs are unknown until after access to keywords in the headers
        # is available, so don't rely on the destination URIs
        #
        # source names do not change, so start with them
        #
        for source_name in self._storage_name.source_names:
            if source_name not in self._storage_name.file_info:
                self._storage_name.file_info[source_name] = get_local_file_info(source_name)
            if source_name not in self._storage_name.metadata:
                self._storage_name.metadata[source_name] = []
                if '.fits' in source_name:
                    self._storage_name._metadata[source_name] = get_local_file_headers(source_name)
                elif self._storage_name.hdf5:
                    if source_name not in self._storage_name._descriptors:
                        # local import to limit exposure in Docker builds
                        import h5py
                        f_in = h5py.File(source_name)
                        self._storage_name._descriptors[source_name] = f_in
                        # Laurie Rosseau-Nepton - 26-04-23
                        # The standard_spectrum is related to flux calibration used on the data. The other one is
                        # for the science data and is the one that should be used.
                        if len(f_in.attrs) >= 20 and 'NAXIS' in f_in.attrs:
                            # 20 is the not-so-random boundary I picked to differentiate between the hdf5 files
                            # that have sufficient metadata to describe WCS.
                            self._logger.debug(f'Found attrs for {source_name}')
                            self._storage_name._metadata[source_name] = [f_in.attrs]
                        else:
                            self._logger.warning(f'No attrs for {source_name}.')
                            self._storage_name._metadata[source_name] = []

            self._storage_name._instrument = get_instrument(
                self._storage_name.metadata.get(source_name), self._storage_name._file_name
            )
            if not self._storage_name.hdf5:
                self._storage_name._bitpix = get_keyword(self._storage_name.metadata.get(source_name), 'BITPIX')

        # ensure the destination uris have the correct extensions based on BITPIX values
        self._storage_name.set_destination_uris()

        # now reset the storage name instance to use the destination URIs, because common code expects URIs
        file_info = self._storage_name.file_info
        metadata = self._storage_name.metadata
        descriptors = self._storage_name._descriptors
        self._storage_name.file_info = {}
        self._storage_name.metadata = {}
        self._storage_name._descriptors = {}
        for index, source_name in enumerate(self._storage_name.source_names):
            uri = self._storage_name.destination_uris[index]
            self._storage_name.file_info[uri] = file_info.get(source_name)
            self._storage_name.metadata[uri] = metadata.get(source_name)
            self._storage_name._descriptors[uri] = descriptors.get(source_name)
        self._logger.debug('End _set_preconditions')


class CFHTNoFheadScrapeRunnerMeta(NoFheadScrapeRunnerMeta):
    """Defines a pipeline step for all the operations that require access to the file on disk for metdata and data
    operations. This is to support HDF5 operations, since at the time of writing, there is no --fhead metadata
    retrieval option for HDF5 files.

    """

    def __init__(self, config, data_visitors, meta_visitors, reporter):
        super().__init__(config, data_visitors, meta_visitors, reporter)

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        #
        # because of decompression, the destination URIs are unknown until after access to keywords in the headers
        # is available, so don't rely on the destination URIs
        #
        # source names do not change, so start with them
        #
        for source_name in self._storage_name.source_names:
            if source_name not in self._storage_name.file_info:
                self._storage_name.file_info[source_name] = get_local_file_info(source_name)
            if source_name not in self._storage_name.metadata:
                self._storage_name.metadata[source_name] = []
                if '.fits' in source_name:
                    self._storage_name._metadata[source_name] = get_local_file_headers(source_name)
                elif self._storage_name.hdf5:
                    if source_name not in self._storage_name._descriptors:
                        # local import to limit exposure in Docker builds
                        import h5py
                        f_in = h5py.File(source_name)
                        self._storage_name._descriptors[source_name] = f_in
                        # Laurie Rosseau-Nepton - 26-04-23
                        # The standard_spectrum is related to flux calibration used on the data. The other one is
                        # for the science data and is the one that should be used.
                        if len(f_in.attrs) >= 20 and 'NAXIS' in f_in.attrs:
                            # 20 is the not-so-random boundary I picked to differentiate between the hdf5 files
                            # that have sufficient metadata to describe WCS.
                            self._logger.debug(f'Found attrs for {source_name}')
                            self._storage_name._metadata[source_name] = [f_in.attrs]
                        else:
                            self._logger.warning(f'No attrs for {source_name}.')
                            self._storage_name._metadata[source_name] = []

            self._storage_name._instrument = get_instrument(
                self._storage_name.metadata.get(source_name), self._storage_name._file_name
            )
            if not self._storage_name.hdf5:
                self._storage_name._bitpix = get_keyword(self._storage_name.metadata.get(source_name), 'BITPIX')

        # ensure the destination uris have the correct extensions based on BITPIX values
        self._storage_name.set_destination_uris()

        # now reset the storage name instance to use the destination URIs, because common code expects URIs
        file_info = self._storage_name.file_info
        metadata = self._storage_name.metadata
        descriptors = self._storage_name._descriptors
        self._storage_name.file_info = {}
        self._storage_name.metadata = {}
        self._storage_name._descriptors = {}
        for index, source_name in enumerate(self._storage_name.source_names):
            uri = self._storage_name.destination_uris[index]
            self._storage_name.file_info[uri] = file_info.get(source_name)
            self._storage_name.metadata[uri] = metadata.get(source_name)
            self._storage_name._descriptors[uri] = descriptors.get(source_name)

        self._logger.debug('End _set_preconditions')


class CFHTNoFheadVisitRunnerMeta(NoFheadVisitRunnerMeta):
    """Defines a pipeline step for all the operations that require access to the file on disk for metdata and data
    operations. This is to support HDF5 operations, since at the time of writing, there is no --fhead metadata
    retrieval option for HDF5 files.

    """

    def __init__(
        self,
        clients,
        config,
        data_visitors,
        meta_visitors,
        modify_transferrer,
        reporter,
    ):
        super().__init__(clients, config, data_visitors, meta_visitors, modify_transferrer, reporter)

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        #
        # for the INGEST + MODIFY combination, it is ok to rely on the destination URIs, since they are all known
        # now, and will not be further affected by decompression work. This URI will work with SI retrieval.
        #
        for index, source_name in enumerate(self._storage_name.source_names):
            uri = self._storage_name.destination_uris[index]
            local_fqn = join(self._working_dir, basename(uri))
            if uri not in self._storage_name.file_info:
                self._storage_name.file_info[uri] = get_local_file_info(local_fqn)
            if uri not in self._storage_name.metadata:
                self._storage_name.metadata[uri] = []
                if '.fits' in source_name:
                    self._storage_name._metadata[uri] = get_local_file_headers(local_fqn)
                elif self._storage_name.hdf5:
                    if uri not in self._storage_name._descriptors:
                        # local import to limit exposure in Docker builds
                        import h5py
                        f_in = h5py.File(local_fqn)
                        self._storage_name._descriptors[uri] = f_in
                        # Laurie Rosseau-Nepton - 26-04-23
                        # The standard_spectrum is related to flux calibration used on the data. The other one is
                        # for the science data and is the one that should be used.
                        if len(f_in.attrs) >= 20 and 'NAXIS' in f_in.attrs:
                            # 20 is the not-so-random boundary I picked to differentiate between the hdf5 files
                            # that have sufficient metadata to describe WCS.
                            self._logger.debug(f'Found attrs for {source_name}')
                            self._storage_name._metadata[uri] = [f_in.attrs]
                        else:
                            self._logger.warning(f'No attrs for {source_name}.')
                            self._storage_name._metadata[uri] = []

        self._storage_name.set_metadata()
        self._logger.debug('End _set_preconditions')


class CFHTNoFheadLocalVisitRunnerMeta(CaomExecuteRunnerMeta):

    def __init__(self, clients, config, data_visitors, meta_visitors, reporter):
        super().__init__(clients, config, meta_visitors, reporter)
        self._data_visitors = data_visitors

    def _set_preconditions(self):
        """This is probably not the best approach, but I want to think about where the optimal location for the
        retrieve_file_info and retrieve_headers methods will be long-term. So, for the moment, use them here."""
        self._logger.debug(f'Begin _set_preconditions for {self._storage_name.file_name}')
        #
        # for the INGEST + MODIFY combination, it is ok to rely on the destination URIs, since they are all known
        # now, and will not be further affected by decompression work. This URI will work with SI retrieval.
        #
        for index, source_name in enumerate(self._storage_name.source_names):
            uri = self._storage_name.destination_uris[index]
            if uri not in self._storage_name.file_info:
                self._storage_name.file_info[uri] = get_local_file_info(source_name)
            if uri not in self._storage_name.metadata:
                self._storage_name.metadata[uri] = []
                if '.fits' in source_name:
                    self._storage_name._metadata[uri] = get_local_file_headers(source_name)
                elif self._storage_name.hdf5:
                    if uri not in self._storage_name._descriptors:
                        # local import to limit exposure in Docker builds
                        import h5py
                        f_in = h5py.File(source_name)
                        self._storage_name._descriptors[uri] = f_in
                        # Laurie Rosseau-Nepton - 26-04-23
                        # The standard_spectrum is related to flux calibration used on the data. The other one is
                        # for the science data and is the one that should be used.
                        if len(f_in.attrs) >= 20 and 'NAXIS' in f_in.attrs:
                            # 20 is the not-so-random boundary I picked to differentiate between the hdf5 files
                            # that have sufficient metadata to describe WCS.
                            self._logger.debug(f'Found attrs for {source_name}')
                            self._storage_name._metadata[uri] = [f_in.attrs]
                        else:
                            self._logger.warning(f'No attrs for {source_name}.')
                            self._storage_name._metadata[uri] = []

        self._storage_name.set_metadata()
        self._logger.debug('End _set_preconditions')

    def execute(self, context):
        self._logger.debug('begin execute with the steps:')
        self.storage_name = context.get('storage_name')

        self._logger.debug('set the preconditions')
        self._set_preconditions()

        self._logger.debug('get the observation for the existing model')
        self._caom2_read()

        self._logger.debug('execute the meta visitors')
        self._visit_meta()

        self._logger.debug('execute the data visitors')
        self._visit_data()

        self._logger.debug('write the observation to disk for debugging')
        self._write_model()

        self._logger.debug('store the updated xml')
        self._caom2_store()

        self._logger.debug('End execute.')


class CFHTStoreIngestRunnerMeta(CFHTNoFheadStoreVisitRunnerMeta):

    def __init__(self, clients, config, meta_visitors, reporter, store_transferrer):
        super().__init__(clients, config, None, meta_visitors, reporter, store_transferrer)

    def execute(self, context):
        self._logger.debug('begin execute with the steps:')
        self.storage_name = context.get('storage_name')

        self._logger.debug('set the preconditions')
        self._set_preconditions()

        self._logger.debug('store the input files')
        self._store_data()

        self._logger.debug('get the observation for the existing model')
        self._caom2_read()

        self._logger.debug('execute the meta visitors')
        self._visit_meta()

        self._logger.debug('write the observation to disk for debugging')
        self._write_model()

        self._logger.debug('store the updated xml')
        self._caom2_store()

        self._logger.debug('End execute.')


class CFHTOrganizeExecutesRunnerMeta(OrganizeExecutesRunnerMeta):
    """A class that extends OrganizeExecutes to handle the choosing of the correct executors based on the config.yml.
    Attributes:
        _needs_delete (bool): if True, the CAOM repo action is delete/create instead of update.
        _reporter: An instance responsible for reporting the execution status.
    Methods:
        _choose():
            Determines which descendants of CaomExecute to instantiate based on the content of the config.yml
            file for an application.
    """

    def __init__(
            self,
            config,
            meta_visitors,
            data_visitors,
            needs_delete=False,
            store_transfer=None,
            modify_transfer=None,
            clients=None,
            reporter=None,
    ):
        super().__init__(
            config,
            meta_visitors,
            data_visitors,
            store_transfer=store_transfer,
            modify_transfer=modify_transfer,
            clients=clients,
            reporter=reporter,
            needs_delete=needs_delete,
        )

    def _choose(self):
        """The logic that decides which descendants of CaomExecute to instantiate. This is based on the content of
        the config.yml file for an application.
        """
        super()._choose()
        if self._needs_delete:
            raise CadcException('No need identified for this yet.')

        if self.can_use_single_visit():
            if TaskType.SCRAPE in self.task_types:
                self._logger.debug(
                    f'Over-riding with executor CFHTNoFheadScrapeRunnerMeta for tasks {self.task_types}.'
                )
                self._executors = []
                self._executors.append(
                    CFHTNoFheadScrapeRunnerMeta(
                        self.config,
                        self._meta_visitors,
                        self._data_visitors,
                        self._reporter,
                    )
                )
            elif TaskType.STORE in self.task_types:
                if TaskType.MODIFY in self.task_types:
                    self._logger.debug(
                        f'Over-riding with executor CFHTNoFheadStoreVisitRunnerMeta for {self.task_types}.'
                    )
                    self._executors = []
                    self._executors.append(
                        CFHTNoFheadStoreVisitRunnerMeta(
                            self._clients,
                            self.config,
                            self._data_visitors,
                            self._meta_visitors,
                            self._reporter,
                            self._store_transfer,
                        )
                    )
                else:
                    self._logger.debug(
                        f'Over-riding with executor CFHTStoreIngestRunnerMeta for {self.task_types}.'
                    )
                    self._executors = []
                    self._executors.append(
                        CFHTStoreIngestRunnerMeta(
                            self._clients,
                            self.config,
                            self._meta_visitors,
                            self._reporter,
                            self._store_transfer,
                        )
                    )
            else:
                if self.config.use_local_files:
                    self._logger.debug(
                        f'Over-riding with executor CFHTNoFheadLocalVisitRunnerMeta for tasks {self.task_types}.'
                    )
                    self._executors = []
                    self._executors.append(
                        CFHTNoFheadLocalVisitRunnerMeta(
                            self._clients,
                            self.config,
                            self._data_visitors,
                            self._meta_visitors,
                            self._reporter,
                        )
                    )
                else:
                    self._logger.debug(
                        f'Over-riding with executor CFHTNoFheadVisitRunnerMeta for tasks {self.task_types}.'
                    )
                    self._executors = []
                    self._executors.append(
                        CFHTNoFheadVisitRunnerMeta(
                            self._clients,
                            self.config,
                            self._data_visitors,
                            self._meta_visitors,
                            self._modify_transfer,
                            self._reporter,
                        )
                    )
        else:
            for task_type in self.task_types:
                if task_type == TaskType.INGEST:
                    self._logger.debug(f'Over-riding with executor CFHTMetaVisitRunnerMeta for {task_type}.')
                    self._executors = []
                    self._executors.append(
                        CFHTMetaVisitRunnerMeta(self._clients, self.config, self._meta_visitors, self._reporter)
                    )
