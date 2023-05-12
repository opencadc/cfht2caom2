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

from os.path import basename, join

from glob import glob
from caom2pipe import astro_composable as ac
from cfht2caom2.cfht_builder import CFHTBuilder, CFHTLocalBuilder, set_storage_name_values
from cfht2caom2 import metadata

from mock import Mock
import pytest


@pytest.fixture()
def fqn(test_data_dir):
    return join(test_data_dir, 'composable_test/test_files/2281792p.fits.fz')


def test_cfht_local_builder(fqn, test_config):
    headers_mock = Mock(autospec=True)
    headers_mock.headers.get.side_effect = lambda ignore: ac.make_headers_from_file(fqn)
    test_config.use_local_files = True
    test_subject = CFHTLocalBuilder(test_config.collection, test_config.use_local_files, headers_mock)
    assert test_subject is not None, 'ctor failure'

    test_result = test_subject.build('123p.hdf5')
    assert test_result is not None, 'expect a result'
    assert test_result.file_name == '123p.hdf5', 'wrong local hdf5 name'
    assert (test_result.instrument == metadata.Inst.SITELLE), 'wrong hdf5 instrument'

    test_result = test_subject.build(fqn)
    assert test_result is not None, 'local fits file failed'
    assert test_result.file_name == basename(fqn), 'wrong local file name'
    assert test_result.instrument == metadata.Inst.WIRCAM


def test_cfht_builder(fqn, test_config):
    test_config.use_local_files = False
    test_subject = CFHTBuilder(test_config.collection)
    assert test_subject is not None, 'ctor failure 2'

    test_result = test_subject.build(fqn)
    assert test_result is not None, 'remote fits file failed'
    assert test_result.file_name == basename(fqn), 'wrong remote file name'
    headers = ac.make_headers_from_file(fqn)
    set_storage_name_values(test_result, headers)
    assert test_result.instrument == metadata.Inst.WIRCAM
    test_uri = 'cadc:CFHT/2281792p.fits.fz'
    assert test_result.source_names[0] == fqn, 'wrong file name'
    assert test_result.destination_uris[0] == test_uri, 'wrong uri'


def test_diag(test_data_dir, test_config):
    test_subject = CFHTBuilder(test_config.collection)
    test_storage_name = test_subject.build('695816p_diag.fits')
    headers = ac.make_headers_from_file(f'{test_data_dir}/single_plane/mega/695816p_diag.fits.header')
    set_storage_name_values(test_storage_name, headers)
    assert test_storage_name.instrument == metadata.Inst.MEGAPRIME


def test_suffixes(test_data_dir, test_config):
    # ensure every test file can be identified as Simple or Derived
    for plane_name in ['single_plane', 'multi_plane']:
        for instrument in ['espadons', 'mega', 'sitelle', 'spirou' ,'wircam']:
            headers_mock = Mock(autospec=True)

            def _mock_get(uri):
                fqn = f'{test_data_dir}/{plane_name}/{instrument}/{basename(uri)}.header'
                return ac.make_headers_from_file(fqn)

            headers_mock.headers.get.side_effect = _mock_get
            test_config.use_local_files = True
            test_builder = CFHTLocalBuilder(test_config.collection, test_config.use_local_files, headers_mock)
            assert test_builder is not None, 'ctor failure'

            plane_list = glob(f'{test_data_dir}/{plane_name}/{instrument}/*.header')
            for entry in plane_list:
                test_subject = test_builder.build(entry.replace('.header', ''))
                assert test_subject is not None, 'ctor'
                found_one = False
                if test_subject.simple:
                    assert not test_subject.derived, f'not derived {test_subject}'
                    found_one = True
                if test_subject.derived:
                    assert not test_subject.simple, f'not simple {test_subject}'
                    found_one = True

                assert found_one, f'{entry} neither derived nor simple {test_subject.is_master_cal}'
