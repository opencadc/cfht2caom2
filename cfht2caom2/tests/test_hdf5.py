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
import h5py
import warnings

from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import FITSFixedWarning
from mock import Mock, patch
from os.path import dirname,join, realpath

from cadcdata import FileInfo
from caom2pipe.manage_composable import ExecutionReporter2
from cfht2caom2 import CFHTName
from cfht2caom2 import file2caom2_augmentation, metadata
import test_caom_gen_visit


THIS_DIR = dirname(realpath(__file__))
TEST_DATA_DIR = join(THIS_DIR, 'data')
SITELLE_DIR = join(TEST_DATA_DIR, 'sitelle')


def pytest_generate_tests(metafunc):
    obs_id_list = glob.glob(f'{SITELLE_DIR}/*.hdf5')
    metafunc.parametrize('test_name', obs_id_list)


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.get_vo_table')
def test_visitor(vo_mock, cache_mock, test_name, test_config, tmp_path, change_test_dir):
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    vo_mock.side_effect = test_caom_gen_visit._vo_mock
    # cache_mock there so there are no update cache calls - so the tests
    # work without a network connection
    test_config.change_working_directory(tmp_path.as_posix())
    storage_name = CFHTName(
        instrument=metadata.Inst.SITELLE,
        source_names=[test_name],
    )
    file_info = FileInfo(id=storage_name.file_uri, file_type='application/x-hdf5')
    storage_name.file_info[storage_name.file_uri] = file_info
    f_in = h5py.File(test_name)
    storage_name.metadata[storage_name.file_uri] = [f_in.attrs]
    storage_name._descriptors[storage_name.file_uri] = f_in
    test_reporter = ExecutionReporter2(test_config)
    kwargs = {
        'working_directory': tmp_path,
        'config': test_config,
        'clients': Mock(),
        'storage_name': storage_name,
        'reporter': test_reporter,
    }
    storage_name._bitpix = -32
    observation = None
    observation = file2caom2_augmentation.visit(observation, **kwargs)

    test_caom_gen_visit._compare(test_name, observation, storage_name.obs_id)
    # assert False
