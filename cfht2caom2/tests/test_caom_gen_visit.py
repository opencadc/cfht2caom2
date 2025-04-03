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
#  : 4 $
#
# ***********************************************************************
#

import glob
import logging
import traceback
import warnings

from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import FITSFixedWarning
from mock import Mock, patch
from os import unlink
from os.path import basename, dirname, exists, join, realpath

from astropy.io.votable import parse_single_table
from caom2.diff import get_differences
from caom2pipe.manage_composable import (
    CadcException,
    ExecutionReporter2,
    Observable,
    read_obs_from_file,
    TaskType,
    write_obs_to_file,
)
from caom2utils.data_util import get_local_file_headers, get_local_file_info
from cfht2caom2.cfht_name import CFHTName, CFHTMetaVisitRunnerMeta
from cfht2caom2 import file2caom2_augmentation
from cfht2caom2 import metadata as md


THIS_DIR = dirname(realpath(__file__))
TEST_DATA_DIR = join(THIS_DIR, 'data')
SINGLE_PLANE_DIR = join(TEST_DATA_DIR, 'single_plane')


def pytest_generate_tests(metafunc):
    obs_id_list = glob.glob(f'{SINGLE_PLANE_DIR}/**/*.fits.header')
    metafunc.parametrize('test_name', obs_id_list)


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.get_vo_table')
def test_visitor(vo_mock, cache_mock, test_name, test_config, tmp_path, change_test_dir):
    # cache_mock there so there are no update cache calls - so the tests work without a network connection
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    vo_mock.side_effect = _vo_mock
    # during cfht2caom2 operation, want to use astropy on FITS files
    # but during testing want to use headers and built-in Python file
    # operations
    instr = dirname(basename(test_name))
    instrument = {
        'espadons': md.Inst.ESPADONS,
        'mega': md.Inst.MEGAPRIME,
        'sitelle': md.Inst.SITELLE,
        'spirou': md.Inst.SPIROU,
        'wircam': md.Inst.WIRCAM,
    }.get(instr)
    test_config.change_working_directory(tmp_path.as_posix())
    test_config.proxy_file_name = 'test_proxy.pem'
    test_config.task_types = [TaskType.SCRAPE]
    test_config.use_local_files = True
    test_config.log_to_file = True
    test_config.data_sources = [tmp_path.as_posix()]
    test_reporter = ExecutionReporter2(test_config)
    expected_fqn = f'{test_name}/{basename(test_name)}.expected.xml'
    in_fqn = expected_fqn.replace('.expected', '.in')
    actual_fqn = expected_fqn.replace('expected', 'actual')
    if exists(actual_fqn):
        unlink(actual_fqn)

    observation = None
    if exists(in_fqn):
        observation = read_obs_from_file(in_fqn)

    with open(test_config.proxy_fqn, 'w') as f:
        f.write('test content')

    clients_mock = Mock()
    test_subject = CFHTMetaVisitRunnerMeta(clients_mock, test_config, [file2caom2_augmentation], test_reporter)
    test_subject._observation = observation

    test_set = [test_name]
    for entry in test_set:
        def _mock_repo_read(collection, obs_id):
            return test_subject._observation
        clients_mock.metadata_client.read.side_effect = _mock_repo_read

        def _read_header_mock(ignore1):
            return get_local_file_headers(entry)
        clients_mock.data_client.get_head.side_effect = _read_header_mock

        def _info_mock(uri):
            temp = get_local_file_info(entry)
            if storage_name.hdf5:
                temp.file_type = 'application/x-hdf5'
            else:
                temp.file_type = 'application/fits'
            temp.size = None  # because the previous tests did not set file size
            return temp
        clients_mock.data_client.info.side_effect = _info_mock

        storage_name = CFHTName(source_names=[entry], instrument=instrument)
        context = {'storage_name': storage_name}
        try:
            test_subject.execute(context)
        except CadcException as e:
            logging.error(traceback.format_exc())
            assert False

    _compare(test_name, test_subject._observation, storage_name.obs_id)
    if storage_name.file_name == '19BMfr.fringe.gri.40.00.fits.header':
        assert len(test_reporter._observable.rejected.get_bad_metadata()) == 1, f'invalid TV_STOP'
    else:
        assert (
            len(test_reporter._observable.rejected.get_bad_metadata()) == 0
        ), f'should be no rejections {test_reporter._observable.rejected.get_bad_metadata()}'

    # assert False


def _compare(test_name, observation, obs_id):
    if observation is None:
        assert False, f'No observation for {obs_id}'
    else:
        expected_fqn = f'{dirname(test_name)}/{obs_id}.expected.xml'
        actual_fqn = expected_fqn.replace('expected', 'actual')
        if exists(expected_fqn):
            expected = read_obs_from_file(expected_fqn)
            compare_result = get_differences(expected, observation)
            if compare_result is None:
                if exists(actual_fqn):
                    unlink(actual_fqn)
            else:
                write_obs_to_file(observation, actual_fqn)
                compare_text = '\n'.join([r for r in compare_result])
                if observation.instrument is None:
                    msg = f'Differences found in {expected.observation_id}\n{compare_text}'
                else:
                    msg = f'Differences found in {expected.observation_id} {observation.instrument.name}\n{compare_text}'
                raise AssertionError(msg)
        else:
            write_obs_to_file(observation, actual_fqn)
            assert False, f'No expected file {expected_fqn} for {obs_id}'


def _vo_mock(url):
    try:
        x = url.split('/')
        filter_name = x[-1].replace('&VERB=0', '')
        votable = parse_single_table(
            f'{TEST_DATA_DIR}/votable/{filter_name}.xml'
        )
        return votable, None
    except Exception as e:
        logging.error(f'get_vo_table failure for url {url}')
