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
import logging
import shutil
import traceback

from os.path import basename
from caom2pipe import astro_composable as ac
from caom2pipe import client_composable as clc
from caom2pipe import data_source_composable as dsc
from caom2pipe import manage_composable as mc

__all__ = ['CFHTTimeBoxDataSource']


class CFHTTimeBoxDataSource(dsc.ListDirTimeBoxDataSource):
    def __init__(self, config, cadc_client):
        super(CFHTTimeBoxDataSource, self).__init__(config)
        self._cadc_client = cadc_client
        self._cleanup_when_storing = config.cleanup_files_when_storing
        self._cleanup_failure_directory = config.cleanup_failure_destination
        self._cleanup_success_directory = config.cleanup_success_destination
        self._store_modified_files_only = config.store_modified_files_only
        self._supports_latest_client = config.features.supports_latest_client
        self._archive = config.archive
        self._collection = config.collection
        self._logger = logging.getLogger(self.__class__.__name__)
        self._logger.error('hello sailor')

    def clean_up(self):
        if self._cleanup_when_storing:
            for entry in self._work:
                self._move_action(
                    entry.entry_name, self._cleanup_success_directory
                )

    def default_filter(self, entry):
        copy_file = True
        if super(CFHTTimeBoxDataSource, self).default_filter(entry):
            if entry.name.startswith('.'):
                # skip dot files
                copy_file = False
            elif ac.check_fits(entry.path):
                # only transfer files that pass the FITS verification
                if self._store_modified_files_only:
                    # only transfer files with a different MD5 checksum
                    copy_file = self._check_md5sum(entry.path)
                    if not copy_file and self._cleanup_when_storing:
                        # KW - 23-06-21
                        # if the file already exists, with the same
                        # checksum, at CADC, Kanoa says move it to the
                        # 'succeeded' directory.
                        self._move_action(
                            entry.path, self._cleanup_success_directory
                        )
            else:
                self._logger.warning(
                    f'Moving {entry.path} to '
                    f'{self._cleanup_failure_directory}'
                )
                self._move_action(entry.path, self._cleanup_failure_directory)
                copy_file = False
        return copy_file

    def _check_md5sum(self, entry_path):
        # get the metadata locally
        result = True
        local_meta = mc.get_file_meta(entry_path)

        # get the metadata at CADC
        f_name = basename(entry_path)
        if self._supports_latest_client:
            destination_name = mc.build_uri(self._collection, f_name, 'cadc')
            cadc_meta = clc.get_cadc_meta_client_v(
                destination_name, self._cadc_client
            )
            cadc_md5sum = cadc_meta.md5sum
        else:
            cadc_meta = clc.get_cadc_meta_client(
                self._cadc_client, self._archive, f_name,
            )
            cadc_md5sum = cadc_meta.md5sum

        if local_meta.get('md5sum') == cadc_md5sum:
            self._logger.warning(
                f'{entry_path} has the same md5sum at CADC. Not transferring.'
            )
            result = False
        return result

    def _move_action(self, fqn, destination):
        # if move when storing is enabled, move to the failure location
        if self._cleanup_when_storing:
            # shutil.move is atomic if it's within a file system, which I
            # believe is the description Kanoa gave. It also supports
            # the same behaviour as
            # https://www.gnu.org/software/coreutils/manual/html_node/
            # mv-invocation.html#mv-invocation when copying between
            # file systems for moving a single file.
            try:
                shutil.move(fqn, destination)
            except Exception as e:
                self._logger.debug(traceback.format_exc())
                self._logger.error(
                    f'Failed to move {fqn} to {destination}'
                )
                raise mc.CadcException(e)