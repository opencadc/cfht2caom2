# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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
import glob
import logging

from astropy.io.votable import parse_single_table

from mock import patch

from cfht2caom2 import main_app, APPLICATION, COLLECTION, CFHTName
from cfht2caom2 import ARCHIVE
from cfht2caom2 import metadata as md
from caom2.diff import get_differences
from caom2pipe import manage_composable as mc

import os
import sys

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
PLUGIN = os.path.join(os.path.dirname(THIS_DIR), 'main_app.py')

LOOKUP = {'key': ['fileid1', 'fileid2']}

WIRCAM_TEST_EXPANSE = {
    'a': '2462928a',
    'd': '787355d',
    'f': '1706150f',
    'g': '1007126g',
    'm': '981337m',
    'o': '2157095o',
    'p': '2281792p',
    's': '2281792s',
    'w': '1151210w',
    'x': '1758254x',
    'y': '1681594y'
}


def pytest_generate_tests(metafunc):
    # obs_id_list = []
    # for ii in LOOKUP:
    #     obs_id_list.append(ii)
    obs_id_list = glob.glob(f'{TEST_DATA_DIR}/*.fits.header')
    metafunc.parametrize('test_name', obs_id_list)


@ patch('cfht2caom2.main_app._identify_instrument')
def test_main_app(inst_mock, test_name):
    inst_mock.side_effect = _identify_inst_mock
    basename = os.path.basename(test_name)
    instrument = _identify_inst_mock(test_name)
    extension = '.fz'
    if instrument is md.Inst.ESPADONS:
        extension = '.gz'
    file_name = basename.replace('.header', extension)
    cfht_name = CFHTName(file_name=file_name,
                         instrument=instrument)
    if '979339' in test_name:
        obs_path = f'{TEST_DATA_DIR}/979339.expected.xml'
        output_file = f'{TEST_DATA_DIR}/979339.actual.xml'
    else:
        obs_path = f'{TEST_DATA_DIR}/{cfht_name.obs_id}.expected.xml'
        output_file = f'{TEST_DATA_DIR}/{basename}.actual.xml'
    expected = mc.read_obs_from_file(obs_path)

    local = _get_local(basename)

    with patch('caom2utils.fits2caom2.CadcDataClient') as data_client_mock, \
        patch('caom2pipe.astro_composable.get_vo_table') as vo_mock:
        def get_file_info(archive, file_id):
            return {'type': 'application/fits'}

        data_client_mock.return_value.get_file_info.side_effect = get_file_info
        vo_mock.side_effect = _vo_mock

        sys.argv = \
            ('{} --no_validate --local {} --observation {} {} -o {} '
             '--plugin {} --module {} --lineage {}'.
             format(APPLICATION, local, COLLECTION,
                    cfht_name.obs_id, output_file, PLUGIN, PLUGIN,
                    _get_lineage(cfht_name))).split()
        print(sys.argv)
        try:
            main_app.to_caom2()
        except Exception as e:
            import logging
            import traceback
            logging.error(traceback.format_exc())

    actual = mc.read_obs_from_file(output_file)
    result = get_differences(expected, actual, 'Observation')
    if result:
        msg = 'Differences found in observation {} test name {}\n{}'. \
            format(expected.observation_id, test_name, '\n'.join(
            [r for r in result]))
        raise AssertionError(msg)
    # assert False  # cause I want to see logging messages


def _get_lineage(cfht_name):
    if '979339' in cfht_name.product_id:
        temp_i = mc.get_lineage(ARCHIVE, '979339i', '979339i.fits.gz')
        temp_o = mc.get_lineage(ARCHIVE, '979339o', '979339o.fits.gz')
        result = f'{temp_i} {temp_o}'
    elif '2281792' in cfht_name.product_id:
        temp_s = mc.get_lineage(ARCHIVE, '2281792s', '2281792s.fits.fz')
        temp_p = mc.get_lineage(ARCHIVE, '2281792p', '2281792p.fits.fz')
        temp_o = mc.get_lineage(ARCHIVE, '2281792o', '2281792o.fits.fz')
        temp_g = mc.get_lineage(ARCHIVE, '2281792g', '2281792g.fits.gz')
        result = f'{temp_s} {temp_p} {temp_o} {temp_g}'
    elif '1151210' in cfht_name.product_id:
        temp_g = mc.get_lineage(ARCHIVE, '1151210g', '1151210g.fits.fz')
        temp_m = mc.get_lineage(ARCHIVE, '1151210m', '1151210m.fits.fz')
        temp_w = mc.get_lineage(ARCHIVE, '1151210w', '1151210w.fits.gz')
        result = f'{temp_g} {temp_m} {temp_w}'
    else:
        result = mc.get_lineage(ARCHIVE, cfht_name.product_id,
                                f'{cfht_name.file_name}')
    return result


def _get_local(test_name):
    if '979339' in test_name:
        result = f'{TEST_DATA_DIR}/979339i.fits ' \
                 f'{TEST_DATA_DIR}/979339o.fits.header'
    elif '1151210' in test_name:
        result = f'{TEST_DATA_DIR}/1151210g.fits.header ' \
                 f'{TEST_DATA_DIR}/1151210m.fits.header ' \
                 f'{TEST_DATA_DIR}/1151210w.fits.header'
    elif '2460503' in test_name:
        result = f'{TEST_DATA_DIR}/2460503p.fits'
    elif '2281792' in test_name:
        result = f'{TEST_DATA_DIR}/2281792s.fits.header ' \
                 f'{TEST_DATA_DIR}/2281792p.fits.header ' \
                 f'{TEST_DATA_DIR}/2281792o.fits.header ' \
                 f'{TEST_DATA_DIR}/2281792g.fits.header'
    else:
        result = f'{TEST_DATA_DIR}/{test_name}'
    return result


def _vo_mock(url):
    try:
        x = url.split('/')
        filter_name = x[-1].replace('&VERB=0', '')
        votable = parse_single_table(
            f'{TEST_DATA_DIR}/votable/{filter_name}.xml')
        return votable, None
    except Exception as e:
        logging.error(f'get_vo_table failure for url {url}')


def _identify_inst_mock(uri):
    result = md.Inst.SITELLE  # 1944968p, 2445397p
    if '2463796o' in uri:
        result = md.Inst.MEGACAM
    elif ('2281792p' in uri or '2157095o' in uri or 'weight' in uri or
          '2281792' in uri):
        result = md.Inst.WIRCAM
    elif ('1001063b' in uri or '1001836x' in uri or '1003681' in uri or
            '1219059' in uri or '1883829c' in uri or
            '2460602a' in uri or
            '760296f' in uri or '881162d' in uri or '979339' in uri or
            '2460503p' in uri):
        result = md.Inst.ESPADONS
    return result
