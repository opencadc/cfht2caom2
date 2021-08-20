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

from astropy.io import fits
from datetime import datetime
from hashlib import md5

from mock import Mock, patch

from cadcdata import FileInfo
from caom2utils import data_util
from caom2 import SimpleObservation, Algorithm, Instrument
from caom2pipe import caom_composable as cc
from caom2pipe import data_source_composable as dsc
from caom2pipe import manage_composable as mc
from cfht2caom2 import composable, cfht_name, metadata
import test_main_app, cfht_mocks

TEST_DIR = f'{test_main_app.TEST_DATA_DIR}/composable_test'


@patch('caom2pipe.execute_composable.CaomExecute._fits2caom2_cmd_local')
@patch('caom2pipe.client_composable.CAOM2RepoClient')
def test_run_by_builder(repo_mock, exec_mock):
    # should attempt to run LocalMetaCreate
    repo_mock.return_value.read.side_effect = _mock_repo_read
    repo_mock.return_value.create.side_effect = Mock()
    repo_mock.return_value.update.side_effect = _mock_repo_update
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
    assert exec_mock.called, 'expect to be called'
    exec_mock.assert_called_with(), 'wrong exec args'


@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
@patch('caom2pipe.client_composable.CadcTapClient')
def test_run_store(
    tap_mock, data_client_mock, repo_client_mock
):
    test_dir_fqn = os.path.join(test_main_app.TEST_DATA_DIR, 'store_test')
    getcwd_orig = os.getcwd
    get_local_orig = data_util.get_local_file_headers
    os.getcwd = Mock(return_value=test_dir_fqn)
    repo_client_mock.return_value.read.side_effect = _mock_repo_read_not_none
    data_client_mock.return_value.info.side_effect = (
        _mock_get_file_info
    )
    data_util.get_local_file_headers = Mock(
        side_effect=_mock_header_read
    )
    try:
        # execution
        test_result = composable._run_by_builder()
        assert test_result == 0, 'wrong result'
    finally:
        os.getcwd = getcwd_orig
        data_util.get_local_file_headers = get_local_orig

    assert data_client_mock.return_value.put.called, 'expect a file put'
    data_client_mock.return_value.put.assert_called_with(
        test_dir_fqn, 'ad:CFHT/1000003f.fits.fz', 'default'
    ), 'wrong put_file args'


@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.CadcTapClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
@patch(
    'caom2pipe.data_source_composable.ListDirTimeBoxDataSource.'
    'get_time_box_work',
    autospec=True,
)
@patch('caom2pipe.execute_composable.OrganizeExecutes.do_one')
def test_run_state(
    run_mock, tap_mock, data_client_mock, tap_client_mock, repo_mock
):
    run_mock.return_value = 0
    tap_mock.side_effect = _mock_dir_listing
    data_client_mock.return_value.get_head.side_effect = cfht_mocks._mock_get_head
    data_client_mock.return_value.info.side_effect = (
        _mock_get_file_info
    )
    getcwd_orig = os.getcwd
    os.getcwd = Mock(return_value=TEST_DIR)

    test_obs_id = '2359320p'
    test_f_name = f'{test_obs_id}.fits'
    try:
        # execution
        test_result = composable._run_state()
        assert test_result == 0, 'mocking correct execution'
    finally:
        os.getcwd = getcwd_orig

    assert run_mock.called, 'should have been called'
    args, kwargs = run_mock.call_args
    test_storage = args[0]
    assert isinstance(test_storage, cfht_name.CFHTName), type(test_storage)
    assert (
        test_storage.obs_id == test_obs_id
    ), f'wrong obs id {test_storage.obs_id}'
    assert test_storage.file_name == test_f_name, 'wrong file name'
    assert test_storage.fname_on_disk == test_f_name, 'wrong fname on disk'
    assert test_storage.url is None, 'wrong url'
    assert (
        test_storage.lineage == f'{test_obs_id}/ad:CFHT/{test_f_name}'
    ), 'wrong lineage'
    assert test_storage.external_urls is None, 'wrong external urls'


@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.CadcTapClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
def test_run_by_builder_hdf5_first(data_mock, tap_mock, repo_mock):
    # create a new observation with an hdf5 file, just using scrape
    # to make sure the observation is writable to an ams service
    #
    # also make sure the SCRAPE task works, so ensure
    # there's no need for credentials, or CADC library clients

    test_obs_id = '2384125p'
    test_dir = f'{test_main_app.TEST_DATA_DIR}/hdf5_test'
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
            f'{test_main_app.TEST_DATA_DIR}/multi_plane/2384125z.hdf5',
            hdf5_fqn,
        )

    _common_execution(test_dir, actual_fqn, expected_hdf5_only_fqn)


@patch('caom2pipe.client_composable.CAOM2RepoClient')
@patch('caom2pipe.client_composable.CadcTapClient')
@patch('caom2pipe.client_composable.StorageClientWrapper')
def test_run_by_builder_hdf5_added_to_existing(data_mock, tap_mock, repo_mock):
    # add to an existing observation with an hdf5 file, just using scrape
    # to make sure the observation is writable to an ams service, and the
    # 'p' metadata gets duplicated correctly

    test_obs_id = '2384125p'
    test_dir = f'{test_main_app.TEST_DATA_DIR}/hdf5_test'
    hdf5_fqn = f'{test_dir}/2384125z.hdf5'
    fits_fqn = f'{test_dir}/{test_obs_id}.fits.header'
    actual_fqn = f'{test_dir}/logs/{test_obs_id}.xml'
    expected_fqn = f'{test_dir}/all.expected.xml'
    expected_hdf5_only_fqn = f'{test_dir}/hdf5_only.expected.xml'

    # make sure expected files are present
    shutil.copy(expected_hdf5_only_fqn, actual_fqn)
    if not os.path.exists(fits_fqn):
        shutil.copy(
            f'{test_main_app.TEST_DATA_DIR}/multi_plane/'
            f'{test_obs_id}.fits.header',
            fits_fqn,
        )

    # clean up unexpected files
    if os.path.exists(hdf5_fqn):
        os.unlink(hdf5_fqn)

    _common_execution(test_dir, actual_fqn, expected_fqn)


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
            lastmod=datetime(
                year=2019, month=3, day=4, hour=19, minute=5
            ),
        )
    else:
        return FileInfo(
            id=file_id,
            size=665345,
            md5sum='md5:a347f2754ff2fd4b6209e7566637efad',
            file_type='application/fits',
            lastmod=datetime(
                year=2019, month=3, day=4, hour=19, minute=5
            ),
        )


def _mock_dir_listing(
    arg1, output_file='', data_only=True, response_format='arg4'
):
    return [
        dsc.StateRunnerMeta(
            os.path.join(
                os.path.join(TEST_DIR, 'test_files'), '2359320p.fits'
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
    delim = '\nEND'
    extensions = [e + delim for e in x.split(delim) if e.strip()]
    headers = [fits.Header.fromstring(e, sep='\n') for e in extensions]
    return headers
