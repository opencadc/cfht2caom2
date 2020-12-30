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

from mock import patch
from cfht2caom2 import metadata as md

import test_main_app


@patch('caom2pipe.manage_composable.query_endpoint')
def test_project_titles_cache(query_mock):
    query_mock.side_effect = _mock_query

    # already cached value
    test_subject = md.cache
    assert test_subject is not None, 'expect a cache'
    test_run_id = '09BC26'
    test_result = test_subject.get_title(test_run_id)
    assert test_result is not None, 'expect a result'
    assert test_result == 'The Next Generation Virgo Survey -- ' \
                          'Infrared: K_s-Band Observations of the Central 4 ' \
                          'deg^2', 'wrong result'

    # empty string does nothing - support case with no RUN_ID in header
    test_run_id = ''
    test_result = test_subject.get_title(test_run_id)
    assert test_result is None, 'expect a result'

    # value not found in cache
    test_run_id = '20AS19'
    test_result = test_subject.get_title(test_run_id)
    assert test_result is not None, f'expect a result {test_run_id}'
    assert test_result == 'Revealing the origin of rprocess elements from ' \
                          'the renhanced stars', 'wrong result'

    test_run_id = '19BS10'
    test_result = test_subject.get_title(test_run_id)
    assert test_result is not None, f'expect a result {test_run_id}'
    assert test_result == 'RVxTESS: Photometric and Spectropolarimetric ' \
                          'studies of M dwarfs with simultaneous TESS and ' \
                          'CFHT/SPIRou Observations', 'wrong result'

    test_run_id = '14BC11'
    test_result = test_subject.get_title(test_run_id)
    assert test_result is not None, f'expect a result {test_run_id}'
    assert test_result == 'Characterizing the magnetic fields of two ' \
                          'recently discovered rare Sigma Ori E type ' \
                          'stars', 'wrong result'


def test_get_repair():
    test_subject = md.cache
    assert test_subject is not None, 'expect a cache'
    # successful repair test case
    test_key = 'Observation.target_position.coordsys'
    test_value = 'FKS'
    test_result = test_subject.get_repair(test_key, test_value)
    assert test_result is not None, 'expect a result'
    assert test_result == 'FK5', 'wrong result'
    test_key = 'Chunk.position.equinox'
    test_value = 200.0
    test_result = test_subject.get_repair(test_key, test_value)
    assert test_result is not None, 'expect a result'
    assert test_result == 2000.0, 'wrong result'

    # no need for repair test case
    test_key = 'Chunk.position.equinox'
    test_value = 201.0
    test_result = test_subject.get_repair(test_key, test_value)
    assert test_result is not None, 'expect a result'
    assert test_result == 201.0, 'wrong result'

    # not found test case
    test_key = 'Not.Found'
    test_result = test_subject.get_repair(test_key, test_value)
    assert test_result is not None, 'expect a result'
    assert test_result == 201.0, 'wrong result'


def _mock_query(url):
    class Object(object):
        def __init__(self):
            pass

        def close(self):
            pass

    response = Object()
    if url == 'http://www.cfht.hawaii.edu/en/science/QSO/2020A/':
        with open(f'{test_main_app.TEST_DATA_DIR}/programs/'
                  f'2020a_index.htm') as f:
            response.text = f.read()
    elif url == 'http://www.cfht.hawaii.edu/en/science/QSO/2020A/' \
                'qso_prog_ESP_2020A.html':
        with open(f'{test_main_app.TEST_DATA_DIR}/programs/'
                  f'espadons_2020a.htm') as f:
            response.text = f.read()
    else:
        response.text = ''
    return response
