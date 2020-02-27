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
import sys

from caom2.diff import get_differences
from caom2pipe import manage_composable as mc
from cfht2caom2 import cfht_name, main_app

from mock import patch
import test_main_app

# structured by observation id, list of file ids that make up a multi-plane
# observation
DIR_NAME = 'multi_plane'
LOOKUP = {'979339': ['979339i.fits.gz', '979339o.fits.gz'],
          '2281792': ['2281792s.fits.fz', '2281792p.fits.fz',
                      '2281792o.fits.fz', '2281792g.fits.gz'],
          '1151210': ['1151210g.fits.fz', '1151210m.fits.fz',
                      '1151210w.fits.gz'],
          # '2384125': ['2384125p.fits.fz', '2384125v.fits.fz', '2384125z.hdf5']
          '2384125p': ['2384125p.fits.fz', '2384125z.hdf5']
          }


def pytest_generate_tests(metafunc):
    obs_id_list = []
    for ii in LOOKUP:
        obs_id_list.append(ii)
    metafunc.parametrize('test_name', obs_id_list)


@patch('cfht2caom2.main_app._identify_instrument')
def test_multi_plane(inst_mock, test_name):
    inst_mock.side_effect = test_main_app._identify_inst_mock
    obs_id = test_name
    lineage = _get_lineage(obs_id)
    actual_fqn = '{}/{}/{}.actual.xml'.format(
        test_main_app.TEST_DATA_DIR, DIR_NAME, obs_id)

    local = _get_local(test_name)
    plugin = test_main_app.PLUGIN

    with patch('caom2utils.fits2caom2.CadcDataClient') as data_client_mock, \
            patch('caom2pipe.astro_composable.get_vo_table') as svofps_mock:
        def get_file_info(archive, file_id):
            return {'type': 'application/fits'}

        data_client_mock.return_value.get_file_info.side_effect = get_file_info
        svofps_mock.side_effect = test_main_app._vo_mock

        if os.path.exists(actual_fqn):
            os.remove(actual_fqn)

        sys.argv = \
            (f'{main_app.APPLICATION} --quiet --no_validate --observation '
             f'{cfht_name.COLLECTION} {test_name} --local {local} '
             f'--plugin {plugin} --module {plugin} --out {actual_fqn} '
             f'--lineage {lineage}').split()
        print(sys.argv)
        main_app.to_caom2()
        expected_fqn = '{}/{}/{}.expected.xml'.format(
            test_main_app.TEST_DATA_DIR, DIR_NAME, obs_id)
        actual = mc.read_obs_from_file(actual_fqn)
        expected = mc.read_obs_from_file(expected_fqn)
        result = get_differences(expected, actual, 'Observation')
        if result:
            msg = 'Differences found in observation {} test name {}\n{}'. \
                format(expected.observation_id, test_name, '\n'.join(
                [r for r in result]))
            raise AssertionError(msg)
        # assert False  # cause I want to see logging messages


def _get_lineage(obs_id):
    result = ''
    for ii in LOOKUP[obs_id]:
        fits = mc.get_lineage(cfht_name.ARCHIVE,
                              cfht_name.CFHTName.remove_extensions(ii), ii)
        result = f'{result } {fits}'
    return result


def _get_local(obs_id):
    result = ''
    root = f'{test_main_app.TEST_DATA_DIR}/{DIR_NAME}'
    if '979339' in obs_id:
        result = f'{root}/979339i.fits ' \
                 f'{root}/979339o.fits.header'
    elif '2384125p' in obs_id:
        # result = f'{root}/2384125p.fits.header ' \
        #          f'{root}/2384125v.fits.header' \
        #          f'{root}/2384125z.hdf5'
        # result = f'{root}/2384125z.hdf5'
        result = f'{root}/2384125p.fits.header ' \
                 f'{root}/2384125z.hdf5'
    else:
        for ii in LOOKUP[obs_id]:
            result = f'{result} {root}/' \
                     f'{cfht_name.CFHTName.remove_extensions(ii)}.fits.header'
    return result