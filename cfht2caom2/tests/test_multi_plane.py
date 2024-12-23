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

import logging
import traceback
import warnings

from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import FITSFixedWarning
from caom2utils.data_util import get_local_file_headers, get_local_file_info
from caom2pipe.manage_composable import CadcException, ExecutionReporter, Observable
from cfht2caom2 import cfht_name, file2caom2_augmentation, metadata
from glob import glob
from os.path import basename, dirname

from mock import Mock, patch
import test_caom_gen_visit


# structured by observation id, list of file ids that make up a multi-plane
# observation
DIR_NAME = 'multi_plane'
LOOKUP = {
    '979339': ['979339i.fits', '979339o.fits'],
    '2281792': [
        '2281792s.fits',
        '2281792p.fits',
        '2281792o.fits',
        '2281792g.fits',
    ],
    '1151210': [
        '1151210g.fits',
        '1151210m.fits',
        '1151210w.fits',
    ],
    '979412': ['979412o.fits', '979412p.fits'],
    '1257365': ['1257365o.fits', '1257365p.fits'],
    '840066': ['840066o.fits', '840066g.fits'],
    '1927963': ['1927963f.fits', '1927963o.fits', '1927963p.fits'],
    '2384125': ['2384125p.fits', '2384125z.hdf5'],
    '2401734': [
        '2401734o.fits',
        '2401734e.fits',
        '2401734r.fits',
        '2401734s.fits',
        '2401734t.fits',
        '2401734v.fits',
    ],
    '1979958': ['1979958p.fits', '1979958y.fits'],
    '2460606': ['2460606i.fits', '2460606o.fits']
}


def pytest_generate_tests(metafunc):
    dir_listing = glob(f'{test_caom_gen_visit.TEST_DATA_DIR}/multi_plane/**/*.expected.xml')
    metafunc.parametrize('expected_fqn', dir_listing)


@patch('cfht2caom2.metadata.CFHTCache._try_to_append_to_cache')
@patch('caom2pipe.astro_composable.get_vo_table')
def test_visitor(vo_mock, cache_mock, expected_fqn, test_data_dir, test_config, tmp_path, change_test_dir):
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    vo_mock.side_effect = test_caom_gen_visit._vo_mock
    # cache_mock there so there are no update cache calls - so the tests
    # work without a network connection
    instr = basename(dirname(expected_fqn))
    instrument = {
        'espadons': metadata.Inst.ESPADONS,
        'mega': metadata.Inst.MEGAPRIME,
        'sitelle': metadata.Inst.SITELLE,
        'spirou': metadata.Inst.SPIROU,
        'wircam': metadata.Inst.WIRCAM,
    }.get(instr)
    test_name = basename(expected_fqn).replace('.expected.xml', '')

    test_config.change_working_directory(tmp_path.as_posix())
    test_config.use_local_files = True
    test_config.log_to_file = True
    test_config.data_sources = [tmp_path.as_posix()]

    test_reporter = ExecutionReporter(test_config, Observable(test_config))
    clients_mock = Mock()
    test_subject = cfht_name.CFHTMetaVisitRunnerMeta(
        clients_mock, test_config, [file2caom2_augmentation], test_reporter
    )
    test_subject._observation = None

    for f_name in LOOKUP[test_name]:
        if 'hdf5' in f_name:
            source_names = [f'/test_files/{f_name}']
        else:
            source_names = [f'{test_caom_gen_visit.TEST_DATA_DIR}/multi_plane/{instr}/{f_name}.header' ]

        def _mock_repo_read(collection, obs_id):
            return test_subject._observation
        clients_mock.metadata_client.read.side_effect = _mock_repo_read

        def _read_header_mock(ignore1):
            return get_local_file_headers(source_names[0])
        clients_mock.data_client.get_head.side_effect = _read_header_mock

        def _info_mock(uri):
            temp = get_local_file_info(source_names[0])
            if storage_name.hdf5:
                temp.file_type = 'application/x-hdf5'
            else:
                temp.file_type = 'application/fits'
            return temp
        clients_mock.data_client.info.side_effect = _info_mock

        storage_name = cfht_name.CFHTName(source_names=source_names, instrument=instrument)
        context = {'storage_name': storage_name}
        try:
            test_subject.execute(context)
        except CadcException as e:
            logging.error(traceback.format_exc())
            assert False

    test_caom_gen_visit._compare(expected_fqn, test_subject._observation, test_name)
