# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2023.                            (c) 2023.
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

from cadcdata import FileInfo
from caom2.diff import get_differences
from caom2pipe.astro_composable import make_headers_from_file
from caom2pipe.manage_composable import Observable, read_obs_from_file, write_obs_to_file
from caom2pipe.reader_composable import Hdf5FileMetadataReader
from cfht2caom2.metadata import Inst
from cfht2caom2.cfht_name import CFHTName
from cfht2caom2 import fits2caom2_augmentation

import pytest
from mock import Mock, patch
from test_fits2caom2_augmentation import _vo_mock


@pytest.fixture()
def test_kwargs(test_data_dir, test_config):    
    test_f_name = '1013501p_flag.fits.fz'
    test_fqn = f'{test_data_dir}/{test_f_name}.txt'
    test_storage_name = CFHTName(
        file_name = test_f_name,
        instrument = Inst.MEGAPRIME,
        source_names = [f'{test_config.scheme}:{test_config.collection}/{test_f_name}'],
    )
    file_info = FileInfo(id=test_storage_name.file_uri, md5sum='abc')
    headers = make_headers_from_file(test_fqn)
    reader = Hdf5FileMetadataReader()
    reader._headers = {test_storage_name.file_uri: headers}
    reader._file_info = {test_storage_name.file_uri: file_info}
    test_clients = Mock(autospec=True)

    def _get_head_mock(ignore):
        return make_headers_from_file(test_fqn)

    test_clients.data_client.get_head.side_effect = _get_head_mock
    test_config.rejected_file_name = 'rejected.yml'
    test_config.rejected_directory = '/tmp'
    kwargs = {
        'clients': test_clients,
        'metadata_reader': reader,
        'observable': Observable(test_config),
        'storage_name': test_storage_name,
        'config': test_config,
    }
    return kwargs


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.get_vo_table')
def test_new_flag_observation(vo_mock, cache_mock, test_kwargs, test_data_dir):
    vo_mock.side_effect = _vo_mock
    test_obs = fits2caom2_augmentation.visit(None, **test_kwargs)
    write_obs_to_file(test_obs, f'{test_data_dir}/x.xml')
    expected_obs = read_obs_from_file(f'{test_data_dir}/flag_no_obs.xml')
    test_result = get_differences(test_obs, expected_obs)
    msg = None
    if test_result is not None:
        msg = '\n'.join(ii for ii in test_result)
    assert test_result is None, f'wrong {msg}'


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.get_vo_table')
def test_add_flag_to_existing_observation(vo_mock, cache_mock, test_kwargs, test_data_dir):
    vo_mock.side_effect = _vo_mock
    test_obs = read_obs_from_file(f'{test_data_dir}/flag_existing_obs_start.xml')
    test_obs = fits2caom2_augmentation.visit(test_obs, **test_kwargs)
    expected_obs = read_obs_from_file(f'{test_data_dir}/flag_existing_obs_end.xml')
    test_result = get_differences(test_obs, expected_obs)
    msg = None
    if test_result is not None:
        msg = '\n'.join(ii for ii in test_result)
    assert test_result is None, f'wrong {msg}'
