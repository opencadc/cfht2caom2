# -*- coding: utf-8 -*-
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

import os

from caom2pipe import manage_composable as mc
from caom2pipe import run_composable as rc
from cfht2caom2 import CFHTBuilder, metadata

from mock import patch
import cfht_mocks


@patch('caom2pipe.run_composable.RunnerClients', autospec=True)
def test_cfht_builder(client_init_mock):
    test_config = mc.Config()
    test_config.use_local_files = True
    test_run_clients = rc.RunnerClients(test_config)
    test_subject = CFHTBuilder(
        test_run_clients,
        test_config.archive,
        test_config.use_local_files,
    )
    assert test_subject is not None, 'ctor failure'

    test_result = test_subject.build('123p.hdf5')
    assert test_result is not None, 'expect a result'
    assert test_result.file_name == '123p.hdf5', 'wrong local hdf5 name'
    assert (
        test_result.instrument == metadata.Inst.SITELLE
    ), 'wrong hdf5 instrument'

    test_fqn = os.path.join(
        cfht_mocks.TEST_DATA_DIR,
        'composable_test/test_files/2281792p.fits.fz',
    )

    test_result = test_subject.build(test_fqn)
    assert test_result is not None, 'local fits file failed'
    assert (
        test_result.file_name == os.path.basename(test_fqn)
    ), 'wrong local file name'
    assert test_result.instrument == metadata.Inst.WIRCAM
    assert not (
        client_init_mock.return_value.data_client.called
    ), 'should not use client'

    client_init_mock.return_value.data_client.get_file.side_effect = \
        cfht_mocks._mock_get_file
    test_config.use_local_files = False
    test_subject = CFHTBuilder(
        test_run_clients,
        test_config.archive,
        test_config.use_local_files,
    )
    assert test_subject is not None, 'ctor failure 2'
    test_result = test_subject.build(test_fqn)
    assert test_result is not None, 'remote fits file failed'
    assert (
            test_result.file_name == os.path.basename(test_fqn)
    ), 'wrong remote file name'
    assert test_result.instrument == metadata.Inst.WIRCAM
    assert (
        client_init_mock.return_value.data_client.get_file.called
    ), 'should use client'
