# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2022.                            (c) 2022.
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

from traceback import format_exc

from caom2pipe.data_source_composable import LocalFilesDataSource
from caom2pipe.manage_composable import get_now


class CFHTLocalFilesDataSource(LocalFilesDataSource):
    """

    With this implementation, all files in a directory listing will be added
    to the set of work to be done. Downstream, the "don't replace files that
    are the same size" checks and the replace functionality of the Storage
    Inventory client, will stop a duplicate STORE execution, but there are a
    lot of checks before that happens.

Test cases:

Source Name   | Interim Name       | Name at               | Move
at CFHT:      | for checking       | CADC:                 | Name at
              | md5sum:            |                       | CFHT
              |                    |                       |
/a/ab.fits    | /a/ab.fits         | cadc:CFHT/ab.fits     | /a/s/ab.fits
/a/cd.fits.gz | /tmp/cd/cd.fits    | cadc:CFHT/cd.fits     | /a/s/cd.fits.gz
/a/ef.fits.gz | /tmp/ef/ef.fits.fz | cadc:CFHT/ef.fits.fz  | /a/s/ef.fits.gz
/a/gh.fits.fz | /a/gh.fits.fz      | cadc:CFHT/gh.fits.fz  | /a/s/gh.fits.fz

    """

    def __init__(
        self,
        config,
        cadc_client,
        metadata_reader,
        recursive,
        builder,
    ):
        super().__init__(config, cadc_client, metadata_reader, recursive)
        # use a builder so there's no need to guess CFHT naming patterns
        self._builder = builder
        self._end_dt = get_now()

    def _post_store_check_md5sum(self, entry_path):
        """
        May decompress/recompress the files as part of the check against
        CADC.

        :return: boolean False if the metadata is the same locally as at
            CADC, True otherwise
        """
        # get the metadata locally
        result = True
        if self._is_connected:
            # get the CADC FileInfo
            storage_name = self._builder.build(entry_path)
            try:
                cadc_meta = self._cadc_client.info(storage_name.file_uri)
            except Exception as e:
                self._logger.error(
                    f'info call failed for {storage_name.destination_uris[0]} '
                    f'with {e}'
                )
                self._logger.debug(format_exc())
                cadc_meta = None

            if cadc_meta is None:
                result = True
            else:
                local_meta = self._metadata_reader.file_info.get(storage_name.file_uri).md5sum.replace('md5:', '')
                self._logger.debug(
                    f'Finding post-store metadata {local_meta} from '
                    f'{storage_name.file_uri}.'
                )
                if local_meta == cadc_meta.md5sum.replace('md5:', ''):
                    result = False
        else:
            self._logger.debug(
                f'SCRAPE\'ing data - no md5sum checking with CADC for '
                f'{entry_path}.'
            )
        temp_text = 'different' if result else 'same'
        self._logger.debug(
            f'Done _post_store_check_md5sum for {entry_path} result is '
            f'{temp_text} at CADC.'
        )
        return result

    def clean_up(self, entry, execution_result, current_count=0):
        self._logger.debug(f'Begin clean_up with {entry}')
        if self._cleanup_when_storing and (
            (not self._retry_failures)
            or (self._retry_failures and current_count >= self._retry_count)
            or (self._retry_failures and execution_result == 0)
        ):
            if isinstance(entry, str):
                fqn = entry
            else:
                fqn = entry.entry_name
            self._logger.debug(f'Clean up {fqn}')
            if self._post_store_check_md5sum(fqn):
                # the transfer itself failed, so track as a failure
                self._move_action(fqn, self._cleanup_failure_directory)
            else:
                self._move_action(fqn, self._cleanup_success_directory)
        self._logger.debug('End clean_up.')
