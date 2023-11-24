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

import h5py
import os
import shutil
import warnings

from astropy.io import fits
from astropy.utils.exceptions import AstropyUserWarning
from collections import deque
from datetime import datetime, timedelta
from hashlib import md5
from logging import error
from traceback import format_exc

from mock import ANY, call, Mock, patch, PropertyMock

from cadcdata import FileInfo
from caom2utils import data_util
from caom2 import SimpleObservation, Algorithm, Instrument
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe.data_source_composable import StateRunnerMeta
from caom2pipe import manage_composable as mc
from caom2pipe.manage_composable import exec_cmd_array
from caom2pipe.run_composable import run_by_state
from cfht2caom2 import composable, cfht_name, metadata
from cfht2caom2.cfht_data_source import CFHTLocalFilesDataSource
import test_caom_gen_visit

TEST_DIR = f'{test_caom_gen_visit.TEST_DATA_DIR}/composable_test'


@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.check_fitsverify')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('caom2pipe.client_composable.ClientCollection')
def test_run_by_builder(
    clients_mock, util_headers_mock, fits_mock, cache_mock, vo_mock, preview_mock
):
    # use_local_files = True, INGEST + MODIFY
    util_headers_mock.side_effect = ac.make_headers_from_file
    # files are valid FITS
    fits_mock.return_value = True
    # have a file with only header values, so preview generation will fail, make it look like it succeeds
    preview_mock.side_effect = _visit_mock

    test_fqn = os.path.join(TEST_DIR, 'test_files')
    if not os.path.exists(test_fqn):
        os.mkdir(test_fqn)
    try:
        # should attempt to run NoFheadVisit
        clients_mock.return_value.metadata_client.read.side_effect = _mock_repo_read
        clients_mock.return_value.metadata_client.create.side_effect = Mock()
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=TEST_DIR)
        try:
            # execution
            test_result = composable._run_by_builder()
            assert test_result == 0, 'wrong result'
        finally:
            os.getcwd = getcwd_orig

        assert clients_mock.return_value.metadata_client.read.called, 'repo read not called'
        assert clients_mock.return_value.metadata_client.create.called, 'repo create not called'
    finally:
        _cleanup(TEST_DIR, 'test_obs_id')


@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.manage_composable.ExecutionReporter')
@patch('caom2pipe.execute_composable.FitsForCADCCompressor.fix_compression')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.check_fitsverify')
@patch('caom2pipe.client_composable.ClientCollection', autospec=True)
def test_run_store(
    clients_mock,
    check_fits_mock,
    cache_mock,
    vo_mock,
    headers_mock,
    compression_mock,
    reporter_mock,
    preview_mock,
    test_data_dir,
):
    compression_mock.side_effect = _mock_fix_compression
    test_dir_fqn = os.path.join(test_data_dir, 'store_test')
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=test_dir_fqn)
    clients_mock.return_value.metadata_client.read.side_effect = _mock_repo_read_not_none
    clients_mock.return_value.data_client.info.side_effect = _mock_get_file_info
    headers_mock.side_effect = ac.make_headers_from_file
    check_fits_mock.return_value = True
    # have a file with only header values, so preview generation will fail, make it look like it succeeds
    preview_mock.side_effect = _visit_mock

    # make sure the rejected.yml file doesn't contain the file under test
    rejected_fqn = f'{test_dir_fqn}/rejected/rejected.yml'
    if os.path.exists(rejected_fqn):
        os.unlink(rejected_fqn)

    try:
        # execution
        test_result = composable._run_by_builder()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig

    assert clients_mock.return_value.data_client.put.called, 'expect a file put'
    assert clients_mock.return_value.data_client.put.call_count == 1, 'file put count'
    clients_mock.return_value.data_client.put.assert_called_with(test_dir_fqn, 'cadc:CFHT/1000003f.fits.fz'), 'put'
    assert reporter_mock.return_value.capture_success.called, 'capture_success'
    assert reporter_mock.return_value.capture_success.call_count == 1, 'capture_success call count'
    assert reporter_mock.return_value.capture_todo.called, 'capture_todo'
    reporter_mock.return_value.capture_todo.assert_called_with(1, 0, 0), 'wrong capture_todo args'
    assert clients_mock.return_value.metadata_client.update.called, 'caom2 update'
    assert clients_mock.return_value.metadata_client.read.called, 'caom2 read'
    clients_mock.return_value.metadata_client.read.assert_called_with('CFHT', '1000003'), 'caom2 read'


@patch('caom2pipe.astro_composable.get_vo_table')
@patch('caom2pipe.execute_composable.FitsForCADCCompressor.fix_compression')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.check_fitsverify')
@patch('caom2pipe.client_composable.ClientCollection')
def test_run_store_retry(clients_mock, check_fits_mock, cache_mock, compression_mock, vo_mock, test_data_dir):
    compression_mock.side_effect = _mock_fix_compression
    test_dir_fqn = os.path.join(test_data_dir, 'store_retry_test')
    test_failure_dir = os.path.join(test_data_dir, 'store_retry_test/failure')
    test_success_dir = os.path.join(test_data_dir, 'store_retry_test/success')
    test_rejected_fqn = os.path.join(test_data_dir, 'store_retry_test/rejected/rejected.yml')
    for ii in [test_failure_dir, test_success_dir]:
        if not os.path.exists(ii):
            os.mkdir(ii)
    test_failure_fqn = os.path.join(test_failure_dir, '1000003f.fits.fz')

    try:
        getcwd_orig = os.getcwd
        get_local_orig = data_util.get_local_headers_from_fits
        os.getcwd = Mock(return_value=test_dir_fqn)
        clients_mock.return_value.metadata_client.read.side_effect = _mock_repo_read_not_none
        clients_mock.return_value.data_client.info.side_effect = _mock_get_file_info
        clients_mock.return_value.data_client.put.side_effect = OSError
        data_util.get_local_headers_from_fits = Mock(side_effect=_mock_header_read)
        check_fits_mock.return_value = True
        try:
            # execution
            test_result = composable._run_by_builder()
            assert test_result == -1, 'all the puts should fail'
        finally:
            os.getcwd = getcwd_orig
            data_util.get_local_headers_from_fits = get_local_orig

        assert clients_mock.return_value.data_client.put.called, 'expect a file put'
        clients_mock.return_value.data_client.put.assert_called_with(
            f'{test_dir_fqn}/new', 'cadc:CFHT/1000003f.fits.fz'
        ), 'wrong put_file args'
        assert os.path.exists(test_failure_fqn), 'expect failure move'
        success_content = os.listdir(test_success_dir)
        assert len(success_content) == 0, 'should be no success files'
        assert os.path.exists(test_rejected_fqn), 'should be an empty rejected file'
        test_rejected = mc.Rejected(test_rejected_fqn)
        for reason, f_name in test_rejected.content.items():
            assert f_name != '1000003f.fits.fz', f'should be no {reason} failure'
    finally:
        if os.path.exists(test_failure_fqn):
            test_new_fqn = os.path.join(test_dir_fqn, 'new/1000003f.fits.fz')
            shutil.move(test_failure_fqn, test_new_fqn)
        _cleanup(test_dir_fqn, '1000003')


@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.client_composable.ClientCollection')
def test_run_store_retry_rejected_entry(clients_mock, cache_mock, vo_mock, test_config, tmp_path):
    # this was originally a test to determine why a file was ending up in the rejected.yml listing
    # when it was in ok shape. It shook out an additional error in the blueprint handling,
    # so now it's just here for regression

    test_failure_dir = os.path.join(tmp_path, 'failure')
    test_success_dir = os.path.join(tmp_path, 'success')

    test_config.change_working_directory(tmp_path)
    test_config.data_sources = [tmp_path.as_posix()]
    test_config.data_source_extensions = ['.fits.gz']
    test_config.use_local_files = True
    test_config.store_modified_files_only = True
    test_config.log_to_file = True
    test_config.retry_failures = True
    test_config.retry_decay = 0
    test_config.cleanup_files_when_storing = True
    test_config.cleanup_failure_destination = test_failure_dir
    test_config.cleanup_success_destination = test_success_dir
    test_config.task_types = [mc.TaskType.STORE, mc.TaskType.INGEST, mc.TaskType.MODIFY]
    test_rejected_fqn = os.path.join(tmp_path, 'rejected/rejected.yml')
    for ii in [test_failure_dir, test_success_dir]:
        os.mkdir(ii)
    test_f_name = '2615124g.fits.gz'
    shutil.copyfile(f'/test_files/{test_f_name}', f'{tmp_path}/{test_f_name}')
    test_failure_fqn = os.path.join(test_failure_dir, test_f_name)
    test_success_fqn = os.path.join(test_success_dir, test_f_name)

    try:
        getcwd_orig = os.getcwd()
        os.chdir(tmp_path)
        mc.Config.write_to_file(test_config)
        clients_mock.return_value.metadata_client.read.side_effect = _mock_repo_read_not_none
        clients_mock.return_value.data_client.info.side_effect = _mock_get_file_info
        clients_mock.return_value.data_client.put.side_effect = [
            None,
            mc.CadcException(
                'Failed to store data with Read timeout on '
                'https://ws-uv.canfar.net/minoc/files/cadc:CFHT/2615124g_preview_1024.jpg'
            ),
            None,  # file retry
            None,  # thumbnail
            None,  # preview
            None,  # zoom
        ]
        data_util.get_local_headers_from_fits = Mock(side_effect=_mock_header_read)

        # execution
        test_result = composable._run_by_builder()
        assert test_result == -1, 'the first preview put should fail'

        assert clients_mock.return_value.data_client.put.called, 'expect a file put'
        assert clients_mock.return_value.data_client.put.call_count == 6, 'file put call count'
        clients_mock.return_value.data_client.put.assert_has_calls(
            [
                call(f'{tmp_path}/2615124', f'cadc:CFHT/{test_f_name.replace(".gz", "")}'),
                call(f'{tmp_path}/2615124', f'cadc:CFHT/2615124g_preview_256.jpg'),
                call(f'{tmp_path}/2615124', f'cadc:CFHT/{test_f_name.replace(".gz", "")}'),
                call(f'{tmp_path}/2615124', f'cadc:CFHT/2615124g_preview_256.jpg'),
                call(f'{tmp_path}/2615124', f'cadc:CFHT/2615124g_preview_1024.jpg'),
                call(f'{tmp_path}/2615124', f'cadc:CFHT/2615124g_preview_zoom_1024.jpg'),
            ],
        ), 'wrong put_file args'
        assert not os.path.exists(test_failure_fqn), 'no failure move'
        assert os.path.exists(test_success_fqn), 'expect success move'
        success_content = os.listdir(test_success_dir)
        assert len(success_content) == 1, 'expect success files'
        assert os.path.exists(test_rejected_fqn), 'should be an empty rejected file'
        test_rejected = mc.Rejected(test_rejected_fqn)
        for reason, f_name in test_rejected.content.items():
            assert f_name != test_f_name, f'should be no {reason} failure'
    finally:
        os.chdir(getcwd_orig)


@patch(
    'cfht2caom2.cfht_data_source.CFHTLocalFilesDataSource.end_dt',
    new_callable=PropertyMock(return_value=datetime(year=2019, month=3, day=7, hour=19, minute=5))
)
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
    end_time_mock,
    test_config,
):
    # it's a WIRCAM file
    test_obs_id = '2281792'
    test_f_name = f'{test_obs_id}p.fits.fz'

    try:
        test_state_fqn = f'{TEST_DIR}/state.yml'
        start_time = datetime(year=2019, month=3, day=3, hour=19, minute=5)
        mc.State.write_bookmark(test_state_fqn, test_config.bookmark, start_time)
        util_headers_mock.side_effect = ac.make_headers_from_file
        run_mock.return_value = 0
        get_work_mock.side_effect = _mock_dir_listing
        getcwd_orig = os.getcwd
        os.getcwd = Mock(return_value=TEST_DIR)

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
        _cleanup(TEST_DIR, test_obs_id)


# common definitions for the test_run_state_compression* tests
class LocalFilesDataSourceCleanupTest(CFHTLocalFilesDataSource):
    def __init__(self, config, cadc_client, metadata_reader, builder):
        super().__init__(config, cadc_client, metadata_reader, False, builder)

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
    # BITPIX -32, no recompression
    result.append(StateRunnerMeta('/test_files/781920i.fits.gz', datetime(2019, 10, 23, 16, 27, 19)))
    # BITPIX 16, recompression
    result.append(StateRunnerMeta('/test_files/1681594g.fits.gz', datetime(2019, 10, 23, 16, 27, 20)))
    # already uncompressed, no decompression or recompression
    result.append(StateRunnerMeta('/test_files/1028439o.fits', datetime(2019, 10, 23, 16, 27, 21)))
    # already compressed, no decompression or recompression
    result.append(StateRunnerMeta('/test_files/2359320o.fits.fz', datetime(2019, 10, 23, 16, 27, 22)))
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


@patch('caom2pipe.manage_composable.ExecutionReporter')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('caom2pipe.client_composable.ClientCollection')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.data_source_composable.ListDirTimeBoxDataSource.get_time_box_work', autospec=True)
def test_run_state_compression_cleanup(
    get_work_mock,
    cache_mock,
    clients_mock,
    vo_table_mock,
    reporter_mock,
    test_config,
    tmp_path,
):
    # this test works with FITS files, not header-only versions of FITS
    # files, because it's testing the decompression/recompression cycle

    get_work_mock.side_effect = _mock_dir_list

    clients_mock.return_value.data_client.info.side_effect = info_returns
    vo_table_mock.side_effect = test_caom_gen_visit._vo_mock

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

    def _mock_read(p1, p2):
        fqn = f'{tmp_path}/logs/{p2}.xml'
        if os.path.exists(fqn):
            # mock the modify task
            return mc.read_obs_from_file(fqn)
        else:
            # mock the ingest task
            return None

    clients_mock.return_value.metadata_client.read.side_effect = _mock_read

    test_config.change_working_directory(tmp_path)
    start_time = datetime.now() - timedelta(minutes=5)
    mc.State.write_bookmark(test_config.state_fqn, test_config.bookmark, start_time)
    test_config.task_types = [mc.TaskType.STORE, mc.TaskType.INGEST, mc.TaskType.MODIFY]
    test_config.log_to_file = True
    test_config.proxy_fqn = f'{tmp_path}/cadcproxy.pem'
    test_config.use_local_files = True
    test_config.data_sources = ['/test_files']
    test_config.retry_failures = False
    test_config.cleanup_files_when_storing = True
    test_config.cleanup_success_destination = '/test_files/success'
    test_config.cleanup_failure_destination = '/test_files/failure'
    test_config.logging_level = 'INFO'

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        mc.Config.write_to_file(test_config)
        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')
        # execution
        try:
            test_config, test_clients, test_reader, test_build, test_source_ignore = composable._common_init()
            test_source = LocalFilesDataSourceCleanupTest(
                test_config, test_clients.data_client, test_reader, test_build
            )
            test_result = run_by_state(
                config=test_config,
                name_builder=test_build,
                meta_visitors=composable.META_VISITORS,
                data_visitors=composable.DATA_VISITORS,
                clients=test_clients,
                sources=[test_source],
                metadata_reader=test_reader,
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
            call(f'{tmp_path}/781920', 'cadc:CFHT/781920i.fits'),
            call(f'{tmp_path}/781920', 'cadc:CFHT/781920i_preview_1024.jpg'),
            call(f'{tmp_path}/781920', 'cadc:CFHT/781920i_preview_256.jpg'),
            call(f'{tmp_path}/1681594', 'cadc:CFHT/1681594g.fits.fz'),
            call(f'{tmp_path}/1681594', 'cadc:CFHT/1681594g_preview_256.jpg'),
            call(f'{tmp_path}/1681594', 'cadc:CFHT/1681594g_preview_1024.jpg'),
            call(f'{tmp_path}/1681594', 'cadc:CFHT/1681594g_preview_zoom_1024.jpg'),
            call('/test_files', 'cadc:CFHT/1028439o.fits'),
            call(f'{tmp_path}/1028439', 'cadc:CFHT/1028439o_preview_256.jpg'),
            call(f'{tmp_path}/1028439', 'cadc:CFHT/1028439o_preview_1024.jpg'),
            call(f'{tmp_path}/1028439', 'cadc:CFHT/1028439o_preview_zoom_1024.jpg'),
            call('/test_files', 'cadc:CFHT/2359320o.fits.fz'),
            call(f'{tmp_path}/2359320', 'cadc:CFHT/2359320o_preview_256.jpg'),
            call(f'{tmp_path}/2359320', 'cadc:CFHT/2359320o_preview_1024.jpg'),
            call(f'{tmp_path}/2359320', 'cadc:CFHT/2359320o_preview_zoom_1024.jpg'),
        ]
        clients_mock.return_value.data_client.put.assert_has_calls(put_calls), 'wrong put args'
        clients_mock.return_value.data_client.info.assert_called(), 'info'
        assert (
            clients_mock.return_value.data_client.info.call_count == 4
        ), 'info call count, only checking _post_store_check_md5sum'
        clients_mock.return_value.data_client.info.assert_has_calls(info_calls), 'wrong info args'

        # LocalStore, get_head should not be called
        clients_mock.return_value.data_client.get_head.assert_not_called()
        # LocalStore, get should not be called
        clients_mock.return_value.data_client.get.assert_not_called()

        assert clients_mock.return_value.metadata_client.read.called, 'read'
        assert clients_mock.return_value.metadata_client.read.call_count == 4, 'meta read call count, ingest + modify'
        read_calls = [
            call('CFHT', '781920'),
            call('CFHT', '1681594'),
            call('CFHT', '1028439'),
            call('CFHT', '2359320'),
        ]
        clients_mock.return_value.metadata_client.read.assert_has_calls(read_calls), 'wrong read args'

        assert clients_mock.return_value.metadata_client.create.called, 'create'
        assert clients_mock.return_value.metadata_client.create.call_count == 4, 'meta create call count'
        create_calls = [call(ANY)]
        clients_mock.return_value.metadata_client.create.assert_has_calls(create_calls), 'wrong create args'

        for obs_id in ['781920', '1681594', '1028439', '2359320']:
            test_obs = mc.read_obs_from_file(f'{test_config.working_directory}/logs/{obs_id}.xml')
            _check_uris(test_obs)
        # capture_todo is not called because the data source is mocked
        assert reporter_mock.return_value.capture_success.called, 'capture_success'
        assert reporter_mock.return_value.capture_success.call_count == 4, 'capture_success call count'
    finally:
        os.chdir(cwd)


@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('caom2pipe.manage_composable.exec_cmd_array')
@patch('caom2pipe.client_composable.ClientCollection')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.data_source_composable.ListDirTimeBoxDataSource.get_time_box_work', autospec=True)
def test_run_state_compression_commands(
    get_work_mock,
    cache_mock,
    clients_mock,
    exec_mock,
    vo_table_mock,
    visit_mock,
    test_config,
    tmp_path,
):
    # this test works with FITS files, not header-only versions of FITS
    # files, because it's testing the decompression/recompression cycle
    # but it's checking that the commands to the exec_cmd_array call are
    # correct

    get_work_mock.side_effect = _mock_dir_list
    clients_mock.return_value.data_client.info.side_effect = info_returns
    clients_mock.return_value.metadata_client.read.return_value = None
    vo_table_mock.side_effect = test_caom_gen_visit._vo_mock
    visit_mock.side_effect = _visit_mock

    def _mock_exec_cmd_array(arg1, arg2):
        exec_cmd_array(arg1, arg2)

    exec_mock.side_effect = _mock_exec_cmd_array

    cwd = os.getcwd()
    test_config.change_working_directory(tmp_path)
    test_config.task_types = [mc.TaskType.STORE, mc.TaskType.INGEST]
    test_config.proxy_file_name = 'cadcproxy.pem'
    test_config.proxy_fqn = f'{tmp_path}/cadcproxy.pem'
    test_config.use_local_files = True
    test_config.data_sources = ['/test_files']
    test_config.state_file_name = 'state.yml'
    test_config.retry_failures = False
    test_config.cleanup_files_when_storing = True
    test_config.cleanup_success_destination = '/test_files/success'
    test_config.cleanup_failure_destination = '/test_files/failure'

    try:
        os.chdir(tmp_path)
        mc.Config.write_to_file(test_config)
        start_time = datetime.now() - timedelta(minutes=5)
        mc.State.write_bookmark(test_config.state_fqn, test_config.bookmark, start_time)
        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')

        # execution
        try:
            test_config, test_clients, test_reader, test_build, test_source_ignore = composable._common_init()
            test_source = LocalFilesDataSourceCleanupTest(
                test_config, test_clients.data_client, test_reader, test_build
            )
            test_result = run_by_state(
                config=test_config,
                name_builder=test_build,
                meta_visitors=composable.META_VISITORS,
                data_visitors=composable.DATA_VISITORS,
                clients=test_clients,
                sources=[test_source],
                metadata_reader=test_reader,
            )
            assert test_result == 0, 'expecting correct execution'
        except Exception as e:
            error(e)
            error(format_exc())
            raise e

        clients_mock.return_value.data_client.put.assert_called(), 'put'
        assert clients_mock.return_value.data_client.put.call_count == 4, 'put call count, including the previews'
        put_calls = [
            call(f'{tmp_path}/781920', 'cadc:CFHT/781920i.fits'),
            call(f'{tmp_path}/1681594', 'cadc:CFHT/1681594g.fits.fz'),
            call('/test_files', 'cadc:CFHT/1028439o.fits'),
            call('/test_files', 'cadc:CFHT/2359320o.fits.fz'),
        ]
        clients_mock.return_value.data_client.put.assert_has_calls(put_calls), 'wrong put args'

        exec_mock.assert_called(), 'exec_cmd_array'
        assert exec_mock.call_count == 1, 'exec_cmd_array call count'
        exec_mock.assert_called_with(
            [
                '/bin/bash',
                '-c',
                f"imcopy /test_files/1681594g.fits.gz "
                f"'{tmp_path}/1681594/1681594g.fits.fz[compress]'",
            ],
            ANY,
        ), 'exec_cmd_array args'

        clients_mock.return_value.data_client.info.assert_called(), 'info'
        assert (
            clients_mock.return_value.data_client.info.call_count == 4
        ), 'info call count, only checking _post_store_check_md5sum'
        clients_mock.return_value.data_client.info.assert_has_calls(info_calls), 'wrong info args'

        # NoFheadStoreVisit, get_head should not be called
        clients_mock.return_value.data_client.get_head.assert_not_called()
        # NoFheadStoreVisit, get should be called
        clients_mock.return_value.data_client.get.assert_not_called()
        assert clients_mock.return_value.metadata_client.read.called, 'read'
        test_report_result = mc.ExecutionSummary.read_report_file(test_config.report_fqn)
        assert test_report_result is not None, 'expect content'
        # 0 is ok here, because the capture_todo call is in the mocked `get_time_box_work` call
        assert test_report_result.entries == 0, f'entries {test_report_result}'
        assert test_report_result.success == 4, f'success {test_report_result}'
        assert test_report_result._errors_sum == 0, f'failure {test_report_result}'
        assert test_report_result._rejected_sum == 0, f'rejected {test_report_result}'
        assert test_report_result._retry_sum == 0, f'retry {test_report_result}'
        assert test_report_result._skipped_sum == 0, f'skipped {test_report_result}'
        assert test_report_result._timeouts_sum == 0, f'timeouts {test_report_result}'
    finally:
        os.chdir(cwd)


@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.client_composable.ClientCollection')
def test_run_by_builder_hdf5_first(clients_mock, cache_mock, vo_mock, preview_mock, test_data_dir):
    # create a new observation with an hdf5 file, just using scrape
    # to make sure the observation is writable to an ams service
    #
    # also make sure the SCRAPE task works, so ensure
    # there's no need for credentials, or CADC library clients

    test_obs_id = '2384125'
    test_dir = f'{test_data_dir}/hdf5_test'
    fits_fqn = f'{test_dir}/{test_obs_id}p.fits.header'
    hdf5_fqn = f'{test_dir}/{test_obs_id}z.hdf5'
    actual_fqn = f'{test_dir}/logs/{test_obs_id}.xml'
    reject_fqn = f'{test_dir}/rejected/rejected.yml'
    expected_hdf5_only_fqn = f'{test_dir}/hdf5_only.expected.xml'
    preview_mock.side_effect = _visit_mock

    # clean up existing observation
    if os.path.exists(actual_fqn):
        os.unlink(actual_fqn)
    if os.path.exists(fits_fqn):
        os.unlink(fits_fqn)
    if os.path.exists(reject_fqn):
        os.unlink(reject_fqn)

    # make sure expected files are present
    # this file does not have the necessary attrs key/value pairs to set the metadata, but is has a few, so it's
    # a test case all on its own
    if not os.path.exists(hdf5_fqn):
        f_in = h5py.File('/test_files/2384125z.hdf5')
        f_out = h5py.File(hdf5_fqn, 'w')
        for k, v in f_in.attrs.items():
            f_out.attrs.create(k, v)
        f_in.close()
        f_out.close()

    try:
        _common_execution(test_dir, actual_fqn, expected_hdf5_only_fqn)
    finally:
        _cleanup(test_dir, test_obs_id)


@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2utils.data_util.get_local_file_info')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('caom2pipe.astro_composable.check_fitsverify')
@patch('caom2pipe.client_composable.ClientCollection')
def test_run_by_builder_hdf5_added_to_existing(
    clients_mock, fits_check_mock, header_mock, file_info_mock, cache_mock, vo_mock, visit_mock, test_data_dir
):
    test_dir = f'{test_data_dir}/hdf5_test'
    visit_mock.side_effect = _visit_mock

    def _make_headers(fqn):
        if not fqn.endswith('.header)'):
            fqn = f'{fqn}.header'
        return ac.make_headers_from_file(fqn)

    try:
        warnings.simplefilter('ignore', category=AstropyUserWarning)
        fits_check_mock.return_value = True
        header_mock.side_effect = _make_headers

        # add to an existing observation with an hdf5 file, just using scrape
        # to make sure the observation is writable to an ams service, and the
        # 'p' metadata gets duplicated correctly

        # can't use the real 'p' file and can't do modify, because it runs out of memory

        test_obs_id = '2384125'
        hdf5_fqn = f'{test_dir}/2384125z.hdf5'
        fits_fqn = f'{test_dir}/{test_obs_id}p.fits.header'
        actual_fqn = f'{test_dir}/logs/{test_obs_id}.xml'
        expected_fqn = f'{test_dir}/all.expected.xml'
        expected_hdf5_only_fqn = f'{test_dir}/hdf5_only.expected.xml'

        file_info_mock.return_value = FileInfo(
            id=f'cadc:CFHT/{test_obs_id}p.fits',
            size=56537,
            file_type='application/fits',
            md5sum='abd123',
        )

        # make sure expected files are present
        if not os.path.exists(actual_fqn):
            if not os.path.exists(os.path.join(test_dir, 'logs')):
                os.mkdir(os.path.join(test_dir, 'logs'))
            shutil.copy(expected_hdf5_only_fqn, actual_fqn)
        if not os.path.exists(fits_fqn):
            shutil.copy(f'{test_data_dir}/multi_plane/sitelle/{test_obs_id}p.fits.header', fits_fqn)

        # clean up unexpected files
        if os.path.exists(hdf5_fqn):
            os.unlink(hdf5_fqn)

        _common_execution(test_dir, actual_fqn, expected_fqn)
    finally:
        _cleanup(test_dir, test_obs_id)


@patch('cfht2caom2.cleanup_augmentation.visit')
@patch('cfht2caom2.espadons_energy_augmentation.visit')
@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.execute_composable.CaomExecute._caom2_store')
@patch('caom2pipe.execute_composable.CaomExecute._visit_meta')
@patch('caom2pipe.data_source_composable.TodoFileDataSource.get_work')
@patch('caom2pipe.client_composable.ClientCollection')
def test_run_ingest(
    clients_mock,
    data_source_mock,
    meta_visit_mock,
    caom2_store_mock,
    cache_mock,
    vo_mock,
    visit_mock,
    espadons_mock,
    cleanup_mock,
    test_config,
    tmp_path,
):
    # access_url_mock.return_value = 'https://localhost:8080'
    temp_deque = deque()
    test_f_name = '1319558w.fits.fz'
    temp_deque.append(test_f_name)
    visit_mock.side_effect = _visit_mock
    espadons_mock.side_effect = _visit_mock
    cleanup_mock.side_effect = _visit_mock
    data_source_mock.return_value = temp_deque
    clients_mock.return_value.metadata_client.read.return_value = None
    clients_mock.return_value.data_client.get_head.return_value = [{'INSTRUME': 'WIRCam'}]

    def _mock_get(arg1, arg2):
        obs_id_dir = f'{tmp_path}/1319558'
        shutil.copy(f'/test_files/{test_f_name}', obs_id_dir)
    clients_mock.return_value.data_client.get.side_effect = _mock_get

    clients_mock.return_value.data_client.info.return_value = FileInfo(
        id=test_f_name,
        file_type='application/fits',
        md5sum='abcdef',
    )

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        test_config.change_working_directory(tmp_path)
        test_config.task_types = [mc.TaskType.INGEST, mc.TaskType.MODIFY]
        test_config.proxy_file_name = 'cadcproxy.pem'
        test_config.proxy_fqn = f'{tmp_path}/cadcproxy.pem'
        test_config.use_local_files = False
        test_config.logging_level = 'DEBUG'
        mc.Config.write_to_file(test_config)
        with open(test_config.proxy_fqn, 'w') as f:
            f.write('test content')
        test_result = composable._run_by_builder()
        assert test_result is not None, 'expect result'
        assert test_result == 0, 'expect success'
        assert clients_mock.return_value.metadata_client.read.called, 'read called'
        assert not clients_mock.return_value.data_client.info.called, 'no info'
        assert not clients_mock.return_value.data_client.get_head.called, 'no get_head'
        assert meta_visit_mock.called, '_visit_meta call'
        assert meta_visit_mock.call_count == 1, '_visit_meta call count'
        assert caom2_store_mock.called, '_caom2_store call'
        assert caom2_store_mock.call_count == 1, '_caom2_store call count'
        assert visit_mock.called, 'preview'
        assert espadons_mock.called, 'espadons'
        assert cleanup_mock.called, 'cleanup'
    finally:
        os.chdir(cwd)


@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.client_composable.ClientCollection')
def test_store_ingest_hdf5(clients_mock, cache_mock, vo_mock, preview_mock, test_config, test_data_dir, tmp_path):
    # create a new observation with an hdf5 file that has attributes
    test_input_dir = f'{test_data_dir}/hdf5_with_metadata_test'
    test_config.change_working_directory(tmp_path)
    test_config.use_local_files = True
    test_config.task_types = [mc.TaskType.STORE, mc.TaskType.INGEST]
    test_config.data_sources = [test_input_dir]
    test_config.data_source_extensions = ['.hdf5']
    test_config.logging_level = 'INFO'
    # test_config.logging_level = 'DEBUG'

    test_hdf5_obs = mc.read_obs_from_file(f'{test_input_dir}/hdf5.xml')
    clients_mock.return_value.metadata_client.read.side_effect = [None, test_hdf5_obs]
    # preview generation will fail because there's no data in the test hdf5 file
    preview_mock.side_effect = _visit_mock

    orig_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        mc.Config.write_to_file(test_config)

        # execution
        test_result = composable._run_by_builder()
        assert test_result == 0, 'wrong result'
        assert clients_mock.return_value.data_client.put.called, 'STORE'
        clients_mock.return_value.data_client.put.assert_called_with(test_input_dir, 'cadc:CFHT/testz.hdf5'), 'STORE'
        assert clients_mock.return_value.metadata_client.create.called, 'STORE metadata'
        test_report_result = mc.ExecutionSummary.read_report_file(test_config.report_fqn)
        assert test_report_result is not None, 'expect content'
        assert test_report_result.entries == 1, f'entries {test_report_result}'
        assert test_report_result.success == 1, f'success {test_report_result}'
        assert test_report_result._errors_sum == 0, f'failure {test_report_result}'
        assert test_report_result._rejected_sum == 0, f'rejected {test_report_result}'
        assert test_report_result._retry_sum == 0, f'retry {test_report_result}'
        assert test_report_result._skipped_sum == 0, f'skipped {test_report_result}'
        assert test_report_result._timeouts_sum == 0, f'timeouts {test_report_result}'
    finally:
        os.chdir(orig_cwd)


@patch('cfht2caom2.preview_augmentation.visit')
@patch('caom2pipe.astro_composable.get_vo_table')
@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.client_composable.ClientCollection')
def test_ingest_modify_hdf5(clients_mock, cache_mock, vo_mock, preview_mock, test_config, test_data_dir, tmp_path):
    # update an observation with an hdf5 file that has attributes
    test_input_dir = f'{test_data_dir}/hdf5_with_metadata_test'
    test_config.change_working_directory(tmp_path)
    test_config.task_types = [mc.TaskType.INGEST, mc.TaskType.MODIFY]
    test_config.logging_level = 'INFO'
    test_config.logging_level = 'DEBUG'

    test_hdf5_obs = mc.read_obs_from_file(f'{test_input_dir}/hdf5.xml')
    clients_mock.return_value.metadata_client.read.return_value = test_hdf5_obs

    def _get_mock(working_directory, ign1):
        import logging
        logging.error('when is this called')
        shutil.copy(f'{test_input_dir}/testz.hdf5', working_directory)
    clients_mock.return_value.data_client.get.side_effect = _get_mock

    preview_mock.side_effect = _visit_mock

    orig_cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        mc.Config.write_to_file(test_config)
        with open(test_config.work_fqn, 'w') as f_out:
            f_out.write('testz.hdf5')

        # execution
        test_result = composable._run_by_builder()
        assert test_result == 0, 'wrong result'
        assert clients_mock.return_value.data_client.get.called, 'MODIFY'
        clients_mock.return_value.data_client.get.assert_called_with(
            f'{tmp_path}/testz', 'cadc:CFHT/testz.hdf5',
        ), 'MODIFY args'
        assert clients_mock.return_value.metadata_client.read.called, 'MODIFY metadata'
        clients_mock.return_value.metadata_client.read.assert_called_with('CFHT', 'testz'), 'MODIFY metadata'
        assert clients_mock.return_value.metadata_client.update.called, 'MODIFY metadata'

    finally:
        os.chdir(orig_cwd)


def _cleanup(test_dir_fqn, obs_id):
    for ii in [
        f'{test_dir_fqn}/logs',
        f'{test_dir_fqn}/logs_0',
        f'{test_dir_fqn}/metrics',
        f'{test_dir_fqn}/rejected',
        f'{test_dir_fqn}/{obs_id}',
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


def _visit_mock(arg1, **kwargs):
    # return the unchanged observation
    return arg1


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
    elif '2615124g' in file_id:
        return FileInfo(
        id=f'cadc:CFHT/2615124g.fits.fz',
        size=1,
        md5sum='e2f542abdc2712bc3bfac91137cba1d6',
        file_type='application/fits',
        encoding=None,
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
            os.path.join(os.path.join(TEST_DIR, 'test_files'), '2281792p.fits.fz'),
            datetime(year=2019, month=10, day=23, hour=16, minute=27, second=19),
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
