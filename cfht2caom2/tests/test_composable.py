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

import os
import shutil
import warnings

from astropy.io import fits
from astropy.utils.exceptions import AstropyUserWarning
from collections import deque
from datetime import datetime, timedelta, timezone
from hashlib import md5
from logging import error
from tempfile import TemporaryDirectory
from traceback import format_exc

from mock import ANY, call, Mock, patch

from cadcdata import FileInfo
from caom2utils import data_util
from caom2 import SimpleObservation, Algorithm, Instrument
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe.data_source_composable import StateRunnerMeta
from caom2pipe import manage_composable as mc
from caom2pipe.manage_composable import exec_cmd_array
from caom2pipe.run_composable import run_by_state
from cfht2caom2 import APPLICATION
from cfht2caom2 import composable, CFHT_BOOKMARK, cfht_name, metadata
from cfht2caom2.cfht_data_source import CFHTLocalFilesDataSource
import test_fits2caom2_augmentation

TEST_DIR = f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/composable_test'


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.check_fits')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
def test_run_by_builder(
    repo_mock, access_mock, util_headers_mock, fits_mock, cache_mock
):
    util_headers_mock.side_effect = ac.make_headers_from_file
    # files are valid FITS
    fits_mock.return_value = True

    test_fqn = os.path.join(TEST_DIR, 'test_files')
    if not os.path.exists(test_fqn):
        os.mkdir(test_fqn)
    try:
        access_mock.return_value = 'https://localhost'
        # should attempt to run MetaVisit
        repo_mock.return_value.read.side_effect = _mock_repo_read
        repo_mock.return_value.create.side_effect = Mock()
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=TEST_DIR)
        try:
            # execution
            test_result = composable._run_by_builder()
            assert test_result == 0, 'wrong result'
        finally:
            os.getcwd = getcwd_orig

        assert repo_mock.return_value.read.called, 'repo read not called'
        assert repo_mock.return_value.create.called, 'repo create not called'
    finally:
        _cleanup(TEST_DIR)


@patch('caom2pipe.execute_composable.FitsForCADCCompressor.fix_compression')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.check_fits')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
@patch('caom2pipe.client_composable.CadcTapClient')
def test_run_store(
    tap_mock,
    data_client_mock,
    repo_client_mock,
    check_fits_mock,
    cache_mock,
    headers_mock,
    compression_mock,
):
    compression_mock.side_effect = _mock_fix_compression
    test_dir_fqn = os.path.join(
        test_fits2caom2_augmentation.TEST_DATA_DIR, 'store_test'
    )
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=test_dir_fqn)
    repo_client_mock.return_value.read.side_effect = _mock_repo_read_not_none
    data_client_mock.return_value.info.side_effect = _mock_get_file_info
    headers_mock.side_effect = ac.make_headers_from_file
    check_fits_mock.return_value = True
    try:
        # execution
        test_result = composable._run_by_builder()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig

    assert data_client_mock.return_value.put.called, 'expect a file put'
    data_client_mock.return_value.put.assert_called_with(
        test_dir_fqn, 'cadc:CFHT/1000003f.fits.fz', None
    ), 'wrong put_file args'


@patch('caom2pipe.execute_composable.FitsForCADCCompressor.fix_compression')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.check_fits')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
@patch('caom2pipe.client_composable.CadcTapClient')
def test_run_store_retry(
    tap_mock,
    data_client_mock,
    repo_client_mock,
    check_fits_mock,
    cache_mock,
    compression_mock,
):
    compression_mock.side_effect = _mock_fix_compression
    test_dir_fqn = os.path.join(
        test_fits2caom2_augmentation.TEST_DATA_DIR, 'store_retry_test'
    )
    test_failure_dir = os.path.join(
        test_fits2caom2_augmentation.TEST_DATA_DIR,
        'store_retry_test/failure',
    )
    test_success_dir = os.path.join(
        test_fits2caom2_augmentation.TEST_DATA_DIR, 'store_retry_test/success'
    )
    for ii in [test_failure_dir, test_success_dir]:
        if not os.path.exists(ii):
            os.mkdir(ii)
    test_failure_fqn = os.path.join(test_failure_dir, '1000003f.fits.fz')
    test_success_fqn = os.path.join(test_success_dir, '1000003f.fits.fz')

    try:
        getcwd_orig = os.getcwd
        get_local_orig = data_util.get_local_headers_from_fits
        os.getcwd = Mock(return_value=test_dir_fqn)
        repo_client_mock.return_value.read.side_effect = (
            _mock_repo_read_not_none
        )
        data_client_mock.return_value.info.side_effect = _mock_get_file_info
        data_client_mock.return_value.put.side_effect = OSError
        data_util.get_local_headers_from_fits = Mock(
            side_effect=_mock_header_read
        )
        check_fits_mock.return_value = True
        try:
            # execution
            test_result = composable._run_by_builder()
            assert test_result == -1, 'all the puts should fail'
        finally:
            os.getcwd = getcwd_orig
            data_util.get_local_headers_from_fits = get_local_orig

        assert data_client_mock.return_value.put.called, 'expect a file put'
        data_client_mock.return_value.put.assert_called_with(
            f'{test_dir_fqn}/new', 'cadc:CFHT/1000003f.fits.fz', None
        ), 'wrong put_file args'
        assert os.path.exists(test_failure_fqn), 'expect failure move'
        success_content = os.listdir(test_success_dir)
        assert len(success_content) == 0, 'should be no success files'
    finally:
        if os.path.exists(test_failure_fqn):
            test_new_fqn = os.path.join(test_dir_fqn, 'new/1000003f.fits.fz')
            shutil.move(test_failure_fqn, test_new_fqn)
        _cleanup(test_dir_fqn)


@patch('caom2pipe.client_composable.ClientCollection')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch(
    'caom2pipe.data_source_composable.ListDirTimeBoxDataSource.'
    'get_time_box_work',
    autospec=True,
)
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run_state(
    run_mock,
    get_work_mock,
    util_headers_mock,
    cache_mock,
    clients_mock,
):
    try:
        util_headers_mock.side_effect = ac.make_headers_from_file
        run_mock.return_value = 0
        get_work_mock.side_effect = _mock_dir_listing
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=TEST_DIR)

        # it's a WIRCAM file
        test_obs_id = '2281792'
        test_f_name = f'{test_obs_id}p.fits.fz'
        try:
            # execution
            test_result = composable._run_state()
            assert test_result == 0, 'mocking correct execution'
        finally:
            os.getcwd = getcwd_orig

        assert run_mock.called, 'should have been called'
        args, kwargs = run_mock.call_args
        test_storage = args[0]
        assert isinstance(test_storage, cfht_name.CFHTName), type(
            test_storage
        )
        assert (
            test_storage.obs_id == test_obs_id
        ), f'wrong obs id {test_storage.obs_id}'
        assert test_storage.file_name == test_f_name, 'wrong file name'
        assert (
            test_storage.file_uri == f'cadc:CFHT/{test_f_name}'
        ), 'wrong uri'
    finally:
        _cleanup(TEST_DIR)


# common definitions for the next two tests
class LocalFilesDataSourceCleanupTest(CFHTLocalFilesDataSource):

        def __init__(self, config, cadc_client, metadata_reader, builder):
            super().__init__(
                config, cadc_client, metadata_reader, False, builder
            )

        def _move_action(self, source, destination):
            uris = {
                '781920i.fits.gz': [
                    '/test_files/781920i.fits.gz',
                    '/test_files/success',
                ],
                '1681594g.fits.gz': [
                    '/test_files/1681594g.fits.gz',
                    '/test_files/success',
                ],
                '1028439o.fits': [
                    '/test_files/1028439o.fits',
                    '/test_files/success',
                ],
                '2359320o.fits.fz': [
                    '/test_files/2359320o.fits.fz',
                    '/test_files/success',
                ],
            }
            f_name = os.path.basename(source)
            lookup = uris.get(f_name)
            assert lookup is not None, f'unexpected f_name {f_name}'
            assert lookup[0] == source, f'source {source}'
            assert lookup[1] == destination, f'destination {destination}'


def _mock_dir_list(
    arg1, output_file='', data_only=True, response_format='arg4'
):
    result = deque()
    result.append(
        StateRunnerMeta(
            '/test_files/781920i.fits.gz',
            '2019-10-23T16:27:19.000',
        ),  # BITPIX -32, no recompression
    )
    result.append(
        StateRunnerMeta(
            '/test_files/1681594g.fits.gz',
            '2019-10-23T16:27:20.000',
        ),  # BITPIX 16, recompression
    )
    result.append(
        StateRunnerMeta(
            '/test_files/1028439o.fits',
            '2019-10-23T16:27:21.000',
        ),  # already uncompressed, no decompression or recompression
    )
    result.append(
        StateRunnerMeta(
            '/test_files/2359320o.fits.fz',
            '2019-10-23T16:27:22.000',
        ),  # already compressed, no decompression or recompression
    )
    return result

uris = {
    '781920': FileInfo(
        'cadc:CFHT/781920i.fits',
        size=10301760,  # not the compressed size of 7630485
        file_type='application/fits',
        md5sum='md5:24cf5c193a312d9aa76d94a5e2cf39c3',
    ),
    '1681594': FileInfo(
        'cadc:CFHT/1681594g.fits.fz',
        size=967680,  # not the .gz compressed size of 197442
        file_type='application/fits',
        md5sum='md5:b1c65d8b1cf5282dcc4444f9c23b7281',
    ),
    '1028439': FileInfo(
        'cadc:CFHT/1028439o.fits',
        size=67317120,  # original size
        file_type='application/fits',
        md5sum='md5:70243ab5f189209d1d74e08edd4a85ae',
    ),
    '2359320': FileInfo(
        'cadc:CFHT/2359320o.fits.fz',
        size=8585280,  # original size
        file_type='application/fits',
        md5sum='md5:1061786cb4da268512e89e252ea26882',
    ),
}

info_returns = [
    uris.get('781920'),
    uris.get('1681594'),
    uris.get('1028439'),
    uris.get('2359320'),
]

info_calls = [
    call('cadc:CFHT/781920i.fits'),
    call('cadc:CFHT/1681594g.fits.fz'),
    call('cadc:CFHT/1028439o.fits'),
    call('cadc:CFHT/2359320o.fits.fz'),
]


@patch('caom2pipe.client_composable.ClientCollection')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch(
    'caom2pipe.data_source_composable.ListDirTimeBoxDataSource.'
    'get_time_box_work',
    autospec=True,
)
def test_run_state_compression_cleanup(
    get_work_mock,
    cache_mock,
    clients_mock,
):
    # this test works with FITS files, not header-only versions of FITS
    # files, because it's testing the decompression/recompression cycle

    get_work_mock.side_effect = _mock_dir_list

    clients_mock.return_value.data_client.info.side_effect = info_returns

    def _check_uris(obs):
        file_info = uris.get(obs.observation_id)
        assert file_info is not None, 'wrong observation id'
        for plane in obs.planes.values():
            for artifact in plane.artifacts.values():
                if artifact.uri == file_info.id:
                    assert (
                        artifact.content_type == file_info.file_type
                    ), 'type'
                    assert artifact.content_length == file_info.size, 'size'
                    assert (
                        artifact.content_checksum.uri == file_info.md5sum
                    ), 'md5'
                    return

        assert False, f'observation id not found {obs.observation_id}'

    cwd = os.getcwd()
    with TemporaryDirectory() as tmp_dir_name:
        os.chdir(tmp_dir_name)

        def _mock_read(p1, p2):
            fqn = f'{tmp_dir_name}/logs/{p2}.xml'
            if os.path.exists(fqn):
                # mock the modify task
                return mc.read_obs_from_file(fqn)
            else:
                # mock the ingest task
                return None

        clients_mock.return_value.metadata_client.read.side_effect = _mock_read

        test_state_fqn = f'{tmp_dir_name}/state.yml'
        start_time = datetime.now(tz=timezone.utc) - timedelta(minutes=5)
        start_file_content = (
            f'bookmarks:\n  cfht_timestamp:\n    last_record: {start_time}\n'
        )
        with open(test_state_fqn, 'w') as f:
            f.write(start_file_content)

        test_config = mc.Config()
        test_config.working_directory = tmp_dir_name
        test_config.task_types = [
            mc.TaskType.STORE, mc.TaskType.INGEST, mc.TaskType.MODIFY
        ]
        test_config.logging_level = 'DEBUG'
        test_config.log_to_file = True
        test_config.collection = 'CFHT'
        test_config.proxy_file_name = 'cadcproxy.pem'
        test_config.proxy_fqn = f'{tmp_dir_name}/cadcproxy.pem'
        test_config.features.supports_latest_client = True
        test_config.features.supports_decompression = True
        test_config.use_local_files = True
        test_config.log_file_directory = f'{tmp_dir_name}/logs'
        test_config.data_sources = '/test_files'
        test_config.state_file_name = 'state.yml'
        test_config.retry_failures = False
        test_config.cleanup_files_when_storing = True
        test_config.cleanup_success_destination = '/test_files/success'
        test_config.cleanup_failure_destination = '/test_files/failure'
        mc.Config.write_to_file(test_config)
        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=tmp_dir_name)
        try:
            # execution
            try:
                (
                    test_config,
                    test_clients,
                    test_reader,
                    test_build,
                    test_source_ignore,
                ) = composable._common_init()
                test_source = LocalFilesDataSourceCleanupTest(
                    test_config,
                    test_clients.data_client,
                    test_reader,
                    test_build,
                )
                test_result = run_by_state(
                    config=test_config,
                    name_builder=test_build,
                    bookmark_name=CFHT_BOOKMARK,
                    meta_visitors=composable.META_VISITORS,
                    data_visitors=composable.DATA_VISITORS,
                    clients=test_clients,
                    source=test_source,
                    metadata_reader=test_reader,
                    application=APPLICATION,
                )
                assert test_result == 0, 'expecting correct execution'
            except Exception as e:
                error(e)
                error(format_exc())
                raise e

            clients_mock.return_value.data_client.put.assert_called(), 'put'
            assert (
                clients_mock.return_value.data_client.put.call_count == 15
            ), 'put call count, including the previews'
            put_calls = [
                call(
                    f'{tmp_dir_name}/781920', 'cadc:CFHT/781920i.fits', None
                ),
                call(
                    f'{tmp_dir_name}/781920',
                    'cadc:CFHT/781920i_preview_1024.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/781920',
                    'cadc:CFHT/781920i_preview_256.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/1681594',
                    'cadc:CFHT/1681594g.fits.fz',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/1681594',
                    'cadc:CFHT/1681594g_preview_256.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/1681594',
                    'cadc:CFHT/1681594g_preview_1024.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/1681594',
                    'cadc:CFHT/1681594g_preview_zoom_1024.jpg',
                    None,
                ),
                call('/test_files', 'cadc:CFHT/1028439o.fits', None),
                call(
                    f'{tmp_dir_name}/1028439',
                    'cadc:CFHT/1028439o_preview_256.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/1028439',
                    'cadc:CFHT/1028439o_preview_1024.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/1028439',
                    'cadc:CFHT/1028439o_preview_zoom_1024.jpg',
                    None,
                ),
                call('/test_files', 'cadc:CFHT/2359320o.fits.fz', None),
                call(
                    f'{tmp_dir_name}/2359320',
                    'cadc:CFHT/2359320o_preview_256.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/2359320',
                    'cadc:CFHT/2359320o_preview_1024.jpg',
                    None,
                ),
                call(
                    f'{tmp_dir_name}/2359320',
                    'cadc:CFHT/2359320o_preview_zoom_1024.jpg',
                    None,
                ),
            ]
            clients_mock.return_value.data_client.put.assert_has_calls(
                put_calls
            ), 'wrong put args'

            clients_mock.return_value.data_client.info.assert_called(), 'info'
            assert (
                clients_mock.return_value.data_client.info.call_count == 4
            ), 'info call count, only checking _post_store_check_md5sum'
            clients_mock.return_value.data_client.info.assert_has_calls(
                info_calls
            ), 'wrong info args'

            clients_mock.return_value.data_client.get_head.assert_not_called(
            ), 'LocalStore, get_head should not be called'
            clients_mock.return_value.data_client.get.assert_not_called(
            ), 'LocalStore, get should not be called'

            clients_mock.return_value.metadata_client.read.assert_called(
            ), 'read'
            assert (
                clients_mock.return_value.metadata_client.read.call_count == 8
            ), 'meta read call count, ingest + modify'
            read_calls = [
                call('CFHT', '781920'),
                call('CFHT', '781920'),
                call('CFHT', '1681594'),
                call('CFHT', '1681594'),
                call('CFHT', '1028439'),
                call('CFHT', '1028439'),
                call('CFHT', '2359320'),
                call('CFHT', '2359320'),
            ]
            clients_mock.return_value.metadata_client.read.assert_has_calls(
                read_calls,
            ), 'wrong read args'

            clients_mock.return_value.metadata_client.create.assert_called(
            ), 'create'
            assert (
                clients_mock.return_value.metadata_client.create.call_count
                == 4
            ), 'meta create call count'
            create_calls = [
                call(ANY),
            ]
            clients_mock.return_value.metadata_client.create.assert_has_calls(
                create_calls,
            ), 'wrong create args'

            for obs_id in ['781920', '1681594', '1028439', '2359320']:
                test_obs = mc.read_obs_from_file(
                    f'{test_config.working_directory}/logs/{obs_id}.xml'
                )
                _check_uris(test_obs)
        finally:
            os.chdir(cwd)
            os.getcwd = getcwd_orig


@patch('caom2pipe.manage_composable.exec_cmd_array')
@patch('caom2pipe.client_composable.ClientCollection')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch(
    'caom2pipe.data_source_composable.ListDirTimeBoxDataSource.'
    'get_time_box_work',
    autospec=True,
)
def test_run_state_compression_commands(
    get_work_mock,
    cache_mock,
    clients_mock,
    exec_mock,
):
    # this test works with FITS files, not header-only versions of FITS
    # files, because it's testing the decompression/recompression cycle
    # but it's checking that the commands to the exec_cmd_array call are
    # correct

    get_work_mock.side_effect = _mock_dir_list
    clients_mock.return_value.data_client.info.side_effect = info_returns

    def _mock_exec_cmd_array(arg1, arg2):
        exec_cmd_array(arg1, arg2)

    exec_mock.side_effect = _mock_exec_cmd_array

    cwd = os.getcwd()
    with TemporaryDirectory() as tmp_dir_name:
        os.chdir(tmp_dir_name)
        test_state_fqn = f'{tmp_dir_name}/state.yml'
        start_time = datetime.now(tz=timezone.utc) - timedelta(minutes=5)
        start_file_content = (
            f'bookmarks:\n  cfht_timestamp:\n    last_record: {start_time}\n'
        )
        with open(test_state_fqn, 'w') as f:
            f.write(start_file_content)

        test_config = mc.Config()
        test_config.working_directory = tmp_dir_name
        test_config.task_types = [mc.TaskType.STORE]
        test_config.logging_level = 'INFO'
        test_config.collection = 'CFHT'
        test_config.proxy_file_name = 'cadcproxy.pem'
        test_config.proxy_fqn = f'{tmp_dir_name}/cadcproxy.pem'
        test_config.features.supports_latest_client = True
        test_config.features.supports_decompression = True
        test_config.use_local_files = True
        test_config.data_sources = '/test_files'
        test_config.state_file_name = 'state.yml'
        test_config.retry_failures = False
        test_config.cleanup_files_when_storing = True
        test_config.cleanup_success_destination = '/test_files/success'
        test_config.cleanup_failure_destination = '/test_files/failure'
        mc.Config.write_to_file(test_config)
        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=tmp_dir_name)
        try:
            # execution
            try:
                (
                    test_config,
                    test_clients,
                    test_reader,
                    test_build,
                    test_source_ignore,
                ) = composable._common_init()
                test_source = LocalFilesDataSourceCleanupTest(
                    test_config,
                    test_clients.data_client,
                    test_reader,
                    test_build,
                )
                test_result = run_by_state(
                    config=test_config,
                    name_builder=test_build,
                    bookmark_name=CFHT_BOOKMARK,
                    meta_visitors=composable.META_VISITORS,
                    data_visitors=composable.DATA_VISITORS,
                    clients=test_clients,
                    source=test_source,
                    metadata_reader=test_reader,
                    application=APPLICATION,
                )
                assert test_result == 0, 'expecting correct execution'
            except Exception as e:
                error(e)
                error(format_exc())
                raise e

            clients_mock.return_value.data_client.put.assert_called(), 'put'
            assert (
                clients_mock.return_value.data_client.put.call_count == 4
            ), 'put call count, including the previews'
            put_calls = [
                call(
                    f'{tmp_dir_name}/781920', 'cadc:CFHT/781920i.fits', None
                ),
                call(
                    f'{tmp_dir_name}/1681594',
                    'cadc:CFHT/1681594g.fits.fz',
                    None,
                ),
                call('/test_files', 'cadc:CFHT/1028439o.fits', None),
                call('/test_files', 'cadc:CFHT/2359320o.fits.fz', None),
            ]
            clients_mock.return_value.data_client.put.assert_has_calls(
                put_calls
            ), 'wrong put args'

            exec_mock.assert_called(), 'exec_cmd_array'
            assert exec_mock.call_count == 1, 'exec_cmd_array call count'
            exec_mock.assert_called_with(
                [
                    '/bin/bash',
                    '-c',
                    f"imcopy /test_files/1681594g.fits.gz "
                    f"'{tmp_dir_name}/1681594/1681594g.fits.fz[compress]'",
                ],
                ANY,
            ), 'exec_cmd_array args'

            clients_mock.return_value.data_client.info.assert_called(), 'info'
            assert (
                clients_mock.return_value.data_client.info.call_count == 4
            ), 'info call count, only checking _post_store_check_md5sum'
            clients_mock.return_value.data_client.info.assert_has_calls(
                info_calls
            ), 'wrong info args'

            clients_mock.return_value.data_client.get_head.assert_not_called(
            ), 'LocalStore, get_head should not be called'
            clients_mock.return_value.data_client.get.assert_not_called(
            ), 'LocalStore, get should not be called'
            clients_mock.return_value.metadata_client.read.assert_not_called(
            ), 'read'
        finally:
            os.chdir(cwd)
            os.getcwd = getcwd_orig


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.CadcTapClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
def test_run_by_builder_hdf5_first(
    data_mock, tap_mock, repo_mock, cache_mock
):
    # create a new observation with an hdf5 file, just using scrape
    # to make sure the observation is writable to an ams service
    #
    # also make sure the SCRAPE task works, so ensure
    # there's no need for credentials, or CADC library clients

    test_obs_id = '2384125p'
    test_dir = f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/hdf5_test'
    fits_fqn = f'{test_dir}/{test_obs_id}.fits.header'
    hdf5_fqn = f'{test_dir}/2384125z.hdf5'
    actual_fqn = f'{test_dir}/logs/{test_obs_id}.xml'
    expected_hdf5_only_fqn = f'{test_dir}/hdf5_only.expected.xml'

    # clean up existing observation
    if os.path.exists(actual_fqn):
        os.unlink(actual_fqn)
    if os.path.exists(fits_fqn):
        os.unlink(fits_fqn)

    # make sure expected files are present
    if not os.path.exists(hdf5_fqn):
        shutil.copy(
            f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/'
            f'multi_plane/2384125z.hdf5',
            hdf5_fqn,
        )

    _common_execution(test_dir, actual_fqn, expected_hdf5_only_fqn)


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('caom2pipe.astro_composable.check_fits')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.CadcTapClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
def test_run_by_builder_hdf5_added_to_existing(
    data_mock, tap_mock, repo_mock, fits_check_mock, header_mock, cache_mock
):
    test_dir = f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/hdf5_test'
    try:
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        fits_check_mock.return_value = True
        header_mock.side_effect = ac.make_headers_from_file

        # add to an existing observation with an hdf5 file, just using scrape
        # to make sure the observation is writable to an ams service, and the
        # 'p' metadata gets duplicated correctly

        test_obs_id = '2384125p'
        hdf5_fqn = f'{test_dir}/2384125z.hdf5'
        fits_fqn = f'{test_dir}/{test_obs_id}.fits.header'
        actual_fqn = f'{test_dir}/logs/{test_obs_id}.xml'
        expected_fqn = f'{test_dir}/all.expected.xml'
        expected_hdf5_only_fqn = f'{test_dir}/hdf5_only.expected.xml'

        # make sure expected files are present
        if not os.path.exists(actual_fqn):
            if not os.path.exists(os.path.join(test_dir, 'logs')):
                os.mkdir(os.path.join(test_dir, 'logs'))
            shutil.copy(expected_hdf5_only_fqn, actual_fqn)
        if not os.path.exists(fits_fqn):
            shutil.copy(
                f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/multi_plane/'
                f'{test_obs_id}.fits.header',
                fits_fqn,
            )

        # clean up unexpected files
        if os.path.exists(hdf5_fqn):
            os.unlink(hdf5_fqn)

        _common_execution(test_dir, actual_fqn, expected_fqn)
    finally:
        _cleanup(test_dir)


@patch('cadcutils.net.ws.WsCapabilities.get_access_url')
@patch('caom2pipe.execute_composable.CaomExecute._caom2_store')
@patch('caom2pipe.execute_composable.CaomExecute._visit_meta')
@patch('caom2pipe.data_source_composable.TodoFileDataSource.get_work')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
def test_run_ingest(
    data_client_mock,
    repo_client_mock,
    data_source_mock,
    meta_visit_mock,
    caom2_store_mock,
    access_url_mock,
):
    access_url_mock.return_value = 'https://localhost:8080'
    temp_deque = deque()
    test_f_name = '1319558w.fits.fz'
    temp_deque.append(test_f_name)
    data_source_mock.return_value = temp_deque
    repo_client_mock.return_value.read.return_value = None
    data_client_mock.return_value.get_head.return_value = [
        {'INSTRUME': 'WIRCam'},
    ]

    data_client_mock.return_value.info.return_value = FileInfo(
        id=test_f_name,
        file_type='application/fits',
        md5sum='abcdef',
    )

    cwd = os.getcwd()
    with TemporaryDirectory() as tmp_dir_name:
        os.chdir(tmp_dir_name)
        test_config = mc.Config()
        test_config.working_directory = tmp_dir_name
        test_config.task_types = [mc.TaskType.INGEST]
        test_config.logging_level = 'INFO'
        test_config.collection = 'CFHT'
        test_config.proxy_file_name = 'cadcproxy.pem'
        test_config.proxy_fqn = f'{tmp_dir_name}/cadcproxy.pem'
        test_config.features.supports_latest_client = True
        test_config.features.supports_decompression = True
        test_config.use_local_files = False
        mc.Config.write_to_file(test_config)
        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=tmp_dir_name)
        try:
            test_result = composable._run_by_builder()
            assert test_result is not None, 'expect result'
            assert test_result == 0, 'expect success'
            assert repo_client_mock.return_value.read.called, 'read called'
            assert data_client_mock.return_value.info.called, 'info'
            assert (
                data_client_mock.return_value.info.call_count == 1
            ), 'wrong number of info calls'
            data_client_mock.return_value.info.assert_called_with(
                f'cadc:CFHT/{test_f_name}',
            )
            assert (
                data_client_mock.return_value.get_head.called
            ), 'get_head should be called'
            assert (
                data_client_mock.return_value.get_head.call_count == 1
            ), 'wrong number of get_heads'
            data_client_mock.return_value.get_head.assert_called_with(
                f'cadc:CFHT/{test_f_name}',
            )
            assert meta_visit_mock.called, '_visit_meta call'
            assert meta_visit_mock.call_count == 1, '_visit_meta call count'
            assert caom2_store_mock.called, '_caom2_store call'
            assert caom2_store_mock.call_count == 1, '_caom2_store call count'
        finally:
            os.getcwd = getcwd_orig
            os.chdir(cwd)


def _cleanup(test_dir_fqn):
    for ii in [
        f'{test_dir_fqn}/logs',
        f'{test_dir_fqn}/logs_0',
        f'{test_dir_fqn}/metrics',
        f'{test_dir_fqn}/rejected',
    ]:
        if os.path.exists(ii):
            for entry in os.scandir(ii):
                os.unlink(entry.path)
            os.rmdir(ii)


def _common_execution(test_dir, actual_fqn, expected_fqn):
    # set up mocks
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=test_dir)
    try:
        # execution
        test_result = composable._run_by_builder()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig
    assert os.path.exists(actual_fqn), f'expect {actual_fqn}'
    compare_result = mc.compare_observations(actual_fqn, expected_fqn)
    if compare_result is not None:
        raise AssertionError(compare_result)


def _mock_repo_create(arg1):
    # arg1 is an Observation instance
    act_fqn = f'{TEST_DIR}/{arg1.observation_id}.xml'
    ex_fqn = f'{TEST_DIR}/{arg1.observation_id}.expected.xml'
    mc.write_obs_to_file(arg1, act_fqn)
    result = cc.compare(ex_fqn, act_fqn, arg1.observation_id)
    if result is not None:
        assert False, result
    pass


def _mock_repo_read(arg1, arg2):
    return None


def _mock_repo_read_not_none(arg1, arg2):
    return SimpleObservation(
        observation_id='TEST_OBS_ID',
        collection='TEST',
        algorithm=Algorithm(name='exposure'),
        instrument=Instrument(name=metadata.Inst.MEGAPRIME.value),
    )


def _mock_repo_update(ignore1):
    return None


def _mock_get_file_info(file_id):
    if '_prev' in file_id:
        return FileInfo(
            id=file_id,
            size=10290,
            md5sum='md5:{}'.format(md5('-37'.encode()).hexdigest()),
            file_type='image/jpeg',
            lastmod=datetime(year=2019, month=3, day=4, hour=19, minute=5),
        )
    else:
        return FileInfo(
            id=file_id,
            size=665345,
            md5sum='md5:a347f2754ff2fd4b6209e7566637efad',
            file_type='application/fits',
            lastmod=datetime(year=2019, month=3, day=4, hour=19, minute=5),
        )


def _mock_dir_listing(
    arg1, output_file='', data_only=True, response_format='arg4'
):
    return [
        StateRunnerMeta(
            os.path.join(
                os.path.join(TEST_DIR, 'test_files'), '2281792p.fits.fz'
            ),
            '2019-10-23T16:27:19.000',
        ),
    ]


def _mock_header_read(file_name):
    x = """SIMPLE  =                    T / 
BITPIX  =                  -32 / Bits per pixel
NAXIS   =                    2 / Number of dimensions
NAXIS1  =                 2048 /
NAXIS2  =                 2048 /
INSTRUME= 'WIRCam  '           /
END
"""
    # TODO return data_util.make_headers_from_string(x)
    delim = '\nEND'
    extensions = [e + delim for e in x.split(delim) if e.strip()]
    headers = [fits.Header.fromstring(e, sep='\n') for e in extensions]
    return headers


def _mock_fix_compression(fqn):
    return fqn
