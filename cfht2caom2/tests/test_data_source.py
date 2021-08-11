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

import pytest
import shutil

from datetime import datetime, timedelta, timezone
from mock import Mock, patch
from pathlib import Path

from cadcdata import FileInfo
from caom2pipe import manage_composable as mc
from cfht2caom2 import data_source


def test_cfht_transfer_check_fits_verify():
    # how things should probably work at CFHT
    delta = timedelta(minutes=30)
    # half an hour ago
    test_start_time = datetime.now(tz=timezone.utc) - delta
    # half an hour from now
    test_end_time = test_start_time + delta + delta
    test_source_directory = Path('/cfht_source')
    test_failure_directory = Path('/cfht_transfer_failure')
    test_success_directory = Path('/cfht_transfer_success')
    test_empty_file = Path('/cfht_source/empty_file.fits.fz')
    test_broken_file = Path('/cfht_source/broken.fits')
    test_broken_source = Path('/test_files/broken.fits')
    test_correct_file = Path('/cfht_source/correct.fits.gz')
    test_correct_source = Path('/test_files/correct.fits.gz')
    test_dot_file = Path('/cfht_source/.dot_file.fits')
    test_same_source = Path('/test_files/same_file.fits')
    test_same_file = Path('/cfht_source/same_file.fits')
    test_already_successful = Path('/cfht_source/already_successful.fits')
    test_already_successful_source = Path(
        '/test_files/already_successful.fits'
    )

    def mock_info(uri):
        return FileInfo(
            id=uri,
            size=12,
            md5sum='e4e153121805745792991935e04de322',
        )

    def _at_cfht(test_start_ts, test_end_ts):
        test_config.cleanup_files_when_storing = True
        cadc_client_mock = Mock(autospec=True)
        cadc_client_mock.info.side_effect = mock_info
        test_subject = data_source.CFHTUseLocalFilesDataSource(
            test_config, cadc_client_mock
        )

        assert test_subject is not None, 'expect construction to work'
        test_result = test_subject.get_time_box_work(test_start_ts, test_end_ts)
        assert len(test_result) == 1, 'wrong number of results returned'
        assert (
                test_result[0].entry_name == '/cfht_source/correct.fits.gz'
        ), 'wrong result'

        for entry in [test_empty_file, test_broken_file]:
            assert not entry.exists(), 'file at source'
            moved = Path(test_failure_directory, entry.name)
            assert moved.exists(), 'file at destination'

        # same file has been moved
        assert not test_same_file.exists(), 'same file at source'
        moved = Path(test_success_directory, test_same_file.name)
        assert moved.exists(), 'same file at destination'

        # correct file - the file stays where it is until it's transferred
        assert test_correct_file.exists(), 'correct file at source'
        moved = Path(test_success_directory, test_correct_file.name)
        assert not moved.exists(), 'correct file at destination'

        # and after the transfer
        test_subject.clean_up()
        assert not test_correct_file.exists(), 'correct file at source'
        assert moved.exists(), 'correct file at destination'

        # dot file - which should be ignored
        assert test_dot_file.exists(), 'correct file at source'
        moved = Path(test_success_directory, test_dot_file.name)
        assert not moved.exists(), 'dot file at destination'

    def _at_cadc(test_start_ts, test_end_ts):
        test_config.cleanup_files_when_storing = False
        test_config.features.supports_latest_client = True
        cadc_client_mock = Mock(autospec=True)
        cadc_client_mock.info.side_effect = mock_info
        test_subject = data_source.CFHTUseLocalFilesDataSource(
            test_config, cadc_client_mock
        )
        assert test_subject is not None, 'expect construction to work'
        test_result = test_subject.get_time_box_work(
            test_start_ts, test_end_ts
        )
        assert len(test_result) == 1, 'wrong number of results returned'
        assert (
            test_result[0].entry_name == '/cfht_source/correct.fits.gz'
        ), 'wrong result'
        for f in [
            test_empty_file,
            test_broken_file,
            test_correct_file,
            test_dot_file,
            test_same_file,
        ]:
            assert f.exists(), 'file at source'
            moved = Path(test_failure_directory, f.name)
            assert not moved.exists(), 'file at destination'
        # clean up should do nothing
        test_subject.clean_up()
        for f in [
            test_empty_file,
            test_broken_file,
            test_correct_file,
            test_dot_file,
            test_same_file,
        ]:
            assert f.exists(), 'file at source'
            moved = Path(test_failure_directory, f.name)
            assert not moved.exists(), 'file at destination'

    def _move_failure(test_start_ts, test_end_ts):
        def _move_mock():
            raise mc.CadcException('move mock')

        move_orig = shutil.move
        shutil.move = Mock(side_effect=_move_mock)
        try:
            test_config.cleanup_files_when_storing = True
            cadc_client_mock = Mock(autospec=True)
            cadc_client_mock.info.side_effect = mock_info
            test_subject = data_source.CFHTUseLocalFilesDataSource(
                test_config, cadc_client_mock
            )
            with pytest.raises(mc.CadcException):
                test_result = test_subject.get_time_box_work(
                    test_start_ts, test_end_ts
                )
        finally:
            shutil.move = move_orig

    for test in [_at_cfht, _at_cadc, _move_failure]:
        for entry in [
            test_failure_directory,
            test_success_directory,
            test_source_directory
        ]:
            if not entry.exists():
                entry.mkdir()
            for child in entry.iterdir():
                child.unlink()
        for entry in [test_empty_file, test_dot_file]:
            entry.touch()
        for source in [
            test_broken_source,
            test_correct_source,
            test_same_source,
            test_already_successful_source,
        ]:
            shutil.copy(source, test_source_directory)

        # CFHT test case - try to move a file that would have the effect
        # of replacing a file already in the destination directory
        shutil.copy(test_already_successful_source, test_success_directory)

        test_config = mc.Config()
        test_config.use_local_files = True
        test_config.data_sources = [test_source_directory.as_posix()]
        test_config.data_source_extensions = ['.fits', '.fits.gz', '.fits.fz']
        test_config.cleanup_success_destination = \
            test_success_directory.as_posix()
        test_config.cleanup_failure_destination = \
            test_failure_directory.as_posix()
        test_config.store_modified_files_only = True

        test(
            test_start_time.timestamp(),
            test_end_time.timestamp(),
        )
