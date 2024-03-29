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

import glob
import logging
import warnings

from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import FITSFixedWarning
from mock import patch
from os import unlink
from os.path import basename, dirname, exists, join, realpath

from astropy.io.votable import parse_single_table
from caom2.diff import get_differences
from cadcdata import FileInfo
from caom2pipe import astro_composable as ac
from caom2pipe.manage_composable import get_keyword, Observable, read_obs_from_file, write_obs_to_file
from caom2pipe import reader_composable as rdc
from caom2utils import data_util
from cfht2caom2 import CFHTName, fits2caom2_augmentation
from cfht2caom2 import metadata as md


THIS_DIR = dirname(realpath(__file__))
TEST_DATA_DIR = join(THIS_DIR, 'data')
SINGLE_PLANE_DIR = join(TEST_DATA_DIR, 'single_plane')


def pytest_generate_tests(metafunc):
    obs_id_list = glob.glob(f'{SINGLE_PLANE_DIR}/**/*.fits.header')
    metafunc.parametrize('test_name', obs_id_list)


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('cfht2caom2.instruments.get_local_headers_from_fits')
@patch('caom2pipe.astro_composable.get_vo_table')
def test_visitor(vo_mock, local_headers_mock, cache_mock, test_name, test_config, tmp_path):
    # cache_mock there so there are no update cache calls - so the tests work without a network connection
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    vo_mock.side_effect = _vo_mock
    # during cfht2caom2 operation, want to use astropy on FITS files
    # but during testing want to use headers and built-in Python file
    # operations
    local_headers_mock.side_effect = _local_headers
    instr = dirname(basename(test_name))
    instrument = {
        'espadons': md.Inst.ESPADONS,
        'mega': md.Inst.MEGAPRIME,
        'sitelle': md.Inst.SITELLE,
        'spirou': md.Inst.SPIROU,
        'wircam': md.Inst.WIRCAM,
    }.get(instr)
    storage_name = CFHTName(
        file_name=basename(test_name).replace('.header', ''),
        instrument=instrument,
        source_names=[test_name],
    )
    file_info = FileInfo(
        id=storage_name.file_uri, file_type='application/fits'
    )
    headers = ac.make_headers_from_file(test_name)
    metadata_reader = rdc.Hdf5FileMetadataReader()
    metadata_reader._headers = {storage_name.file_uri: headers}
    metadata_reader._file_info = {storage_name.file_uri: file_info}
    test_config.rejected_file_name = 'rejected.yml'
    test_config.rejected_directory = tmp_path.as_posix()
    test_observable = Observable(test_config)
    kwargs = {
        'storage_name': storage_name,
        'metadata_reader': metadata_reader,
        'observable': test_observable,
        'config': test_config,
    }
    storage_name._bitpix = get_keyword(headers, 'BITPIX')
    observation = None
    observation = fits2caom2_augmentation.visit(observation, **kwargs)

    _compare(test_name, observation, storage_name.obs_id)
    if storage_name.file_name == '19BMfr.fringe.gri.40.00.fits':
        assert len(test_observable.rejected.get_bad_metadata()) == 1, f'invalid TV_STOP'
    else:
        assert (
            len(test_observable.rejected.get_bad_metadata()) == 0
        ), f'should be no rejections {test_observable.rejected.get_bad_metadata()}'

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


def _local_headers(fqn):
    from urllib.parse import urlparse
    from astropy.io import fits

    file_uri = urlparse(fqn)
    try:
        fits_header = open(file_uri.path).read()
        headers = data_util.make_headers_from_string(fits_header)
    except UnicodeDecodeError:
        hdulist = fits.open(fqn, memmap=True, lazy_load_hdus=True)
        hdulist.verify('fix')
        hdulist.close()
        headers = [h.header for h in hdulist]
    return headers


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
