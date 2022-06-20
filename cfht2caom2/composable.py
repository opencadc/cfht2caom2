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

import logging
import sys
import traceback

from os.path import basename

from caom2pipe import client_composable as clc
from caom2pipe import data_source_composable as dsc
from caom2pipe.manage_composable import Config, get_keyword, StorageName
from caom2pipe.reader_composable import (
    FileMetadataReader,
    StorageClientReader,
)
from caom2pipe import run_composable as rc
from cfht2caom2 import cleanup_augmentation
from cfht2caom2 import espadons_energy_augmentation, preview_augmentation
from cfht2caom2 import fits2caom2_augmentation
from cfht2caom2.cfht_builder import CFHTBuilder
from cfht2caom2.cfht_name import CFHTName
from cfht2caom2.metadata import Inst


META_VISITORS = [fits2caom2_augmentation]
DATA_VISITORS = [
    espadons_energy_augmentation,
    preview_augmentation,
    cleanup_augmentation,
]

CFHT_BOOKMARK = 'cfht_timestamp'


def _common_init():
    config = Config()
    config.get_executors()
    StorageName.collection = config.collection
    StorageName.scheme = (
        'cadc' if config.features.supports_latest_client else 'ad'
    )
    clients = clc.ClientCollection(config)
    if config.use_local_files:
        reader = FileMetadataReader()
    else:
        reader = StorageClientReader(clients.data_client)
    builder = CFHTBuilder(
        config.archive,
        config.use_local_files,
        reader,
        config.features.supports_latest_client,
    )
    source = None
    if config.use_local_files:
        source = dsc.LocalFilesDataSource(
            config,
            clients.data_client,
            reader,
            recursive=config.recurse_data_sources,
        )
    return config, clients, reader, builder, source


def _run_state():
    config, clients, reader, builder, source = _common_init()
    return rc.run_by_state(
        config=config,
        name_builder=builder,
        bookmark_name=CFHT_BOOKMARK,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        clients=clients,
        source=source,
        metadata_reader=reader,
    )


def run_state():
    """Wraps _run_state in exception handling."""
    try:
        _run_state()
        sys.exit(0)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_by_builder():
    """Run the processing for observations using a todo file to identify the
    work to be done, but with the support of a Builder, so that StorageName
    instances can be provided. This is important here, because the
    instrument name needs to be provided to the StorageName constructor.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    config, clients, reader, builder, source = _common_init()
    return rc.run_by_todo(
        config,
        builder,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        clients=clients,
        source=source,
        metadata_reader=reader,
    )


def run_by_builder():
    try:
        result = _run_by_builder()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run_decompress():
    """
    Run decompression for CFHT. This is interim code, hence the
    location of the DecompressReader class definition.
    """

    class DecompressReader(StorageClientReader):
        """
        Using the source_names for the set_* implementations means
        the information can be obtained for the compressed file name, and
        since the header information is required to decide whether or not
        to recompress a file, this is important distinction during
        execution. In the case, the chicken in the compressed file name,
        and the egg is the "does it need to be recompressed" decision.
        """

        def __init__(self, client):
            super().__init__(client)

        def set_file_info(self, storage_name):
            """Retrieves FileInfo information to memory."""
            for entry in storage_name.source_names:
                if entry not in self._file_info.keys():
                    self._file_info[entry] = self._client.info(entry)

        def set_headers(self, storage_name):
            """Retrieves the Header information to memory."""
            self._logger.debug(
                f'Begin set_headers for {storage_name.file_name}'
            )
            for entry in storage_name.source_names:
                if entry not in self._headers.keys():
                    if '.fits' in entry:
                        self._headers[entry] = self._client.get_head(entry)
                    else:
                        self._headers[entry] = []
            self._logger.debug('End set_headers')

    class DecompressBuilder(CFHTBuilder):
        def __init__(
            self,
            archive,
            use_local_files,
            metadata_reader,
            supports_latest_client,
        ):
            super().__init__(
                archive,
                use_local_files,
                metadata_reader,
                supports_latest_client,
            )

        def build(self, entry):
            """
            :param entry an entry is a file name, complete with the appropriate
            compression extension, that is sufficient to retrieve file header
            information from CADC's storage system.
            """

            # retrieve the header information, extract the instrument name
            self._logger.debug(f'Build a StorageName instance for {entry}.')
            bitpix = None
            if StorageName.is_hdf5(entry):
                instrument = Inst.SITELLE
            else:
                file_name = basename(entry).replace('.header', '')
                # the separate construction of file name for the uri supports
                # unit testing
                storage_name = StorageName(
                    file_name=file_name, source_names=[entry]
                )
                self._metadata_reader.set(storage_name)
                headers = self._metadata_reader.headers.get(entry)
                instrument = CFHTBuilder.get_instrument(headers, entry)
                bitpix = get_keyword(headers, 'BITPIX')
            result = CFHTName(
                file_name=basename(entry),
                source_names=[entry],
                instrument=instrument,
                bitpix=bitpix,
            )
            self._logger.debug('End build.')
            return result

    config, clients, reader_ignore, builder_ignore, source = _common_init()
    # the headers have to come from AD, because SI doesn't do that atm
    original_feature = config.features.supports_latest_client
    config.features.supports_latest_client = False
    ad_reading = clc.ClientCollection(config)
    decompress_reader = DecompressReader(ad_reading.data_client)
    config.features.supports_latest_client = original_feature
    builder = DecompressBuilder(
        config.archive,
        config.use_local_files,
        decompress_reader,
        config.features.supports_latest_client,
    )

    return rc.run_by_todo(
        config,
        builder,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        clients=clients,
        source=source,
        metadata_reader=decompress_reader,
    )


def run_decompress():
    try:
        result = _run_decompress()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
