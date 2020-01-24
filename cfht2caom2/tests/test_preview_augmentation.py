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
import pytest

from datetime import datetime
from mock import patch, Mock

from caom2 import ChecksumURI, Dimension2D, Artifact, ReleaseType, ProductType
from gem2caom2 import preview_augmentation, pull_augmentation, SCHEME
from gem2caom2 import ARCHIVE
from caom2pipe import manage_composable as mc

pytest.main(args=['-s', os.path.abspath(__file__)])
THIS_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_DIR = os.path.join(THIS_DIR, 'data')
TEST_OBS = 'GN-2013B-Q-28-150-002'
TEST_FILE = 'N20131203S0006.jpg'
TEST_FP_OBS = 'GN-2015A-C-1-20-001'
TEST_FP_FILE = 'N20150404S0726.fits'
TEST_OBS_FILE = '{}/{}'.format(TEST_DATA_DIR, 'visit_obs_start.xml')
TEST_PRODUCT_ID = 'GN2001BQ013-04'
REJECTED_FILE = os.path.join(TEST_DATA_DIR, 'rejected.yml')


def test_preview_augment():
    # this should result in two new artifacts being added to the plane
    # one for a thumbnail and one for a preview

    obs = mc.read_obs_from_file(TEST_OBS_FILE)
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.utcnow()
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    test_rejected = mc.Rejected(REJECTED_FILE)
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))
    cadc_client_mock = Mock()
    kwargs = {'working_directory': TEST_DATA_DIR,
              'cadc_client': cadc_client_mock,
              'stream': 'stream',
              'observable': test_observable}

    test_prev = '{}/{}.jpg'.format(TEST_DATA_DIR, TEST_PRODUCT_ID)
    if os.path.exists(test_prev):
        os.unlink(test_prev)

    try:
        with patch('caom2pipe.manage_composable.http_get') as http_mock, \
                patch('caom2pipe.manage_composable.data_put') as ad_put_mock, \
                patch('caom2pipe.manage_composable.get_artifact_metadata') as \
                art_mock, \
                patch('caom2pipe.manage_composable.exec_cmd') as exec_mock:
            def _art_mock(fqn, product_type, release_type, uri, temp):
                if '_th' in uri:
                    return Artifact(uri,
                                    ProductType.PREVIEW,
                                    ReleaseType.DATA,
                                    'image/jpeg',
                                    13,
                                    ChecksumURI('md5:13'))
                else:
                    return Artifact(uri,
                                    ProductType.PREVIEW,
                                    ReleaseType.DATA,
                                    'image/jpeg',
                                    12,
                                    ChecksumURI('md5:12'))

            def _get_mock(url_ignore, fqn_ignore):
                with open(test_prev, 'w') as f:
                    f.write('abc')

            cadc_client_mock.return_value.data_get.return_value = \
                mc.CadcException('test')
            art_mock.side_effect = _art_mock
            http_mock.side_effect = _get_mock
            result = preview_augmentation.visit(obs, **kwargs)
            test_url = '{}{}.fits'.format(preview_augmentation.PREVIEW_URL,
                                          TEST_PRODUCT_ID)
            assert http_mock.called, 'http mock should be called'
            http_mock.assert_called_with(test_url, test_prev),  \
            'mock not called'
            assert ad_put_mock.called, 'ad put mock not called'
            assert art_mock.called, 'art mock not called'
            assert exec_mock.called, 'exec mock not called'
            assert result is not None, 'expect a result'
            assert result['artifacts'] == 2, 'artifacts should be added'
            assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 3, \
                'two new artifacts'
            prev_uri = mc.build_uri(
                ARCHIVE, '{}.jpg'.format(TEST_PRODUCT_ID), SCHEME)
            thumb_uri = mc.build_uri(ARCHIVE, '{}_th.jpg'.format(
                TEST_PRODUCT_ID))
            assert prev_uri in obs.planes[TEST_PRODUCT_ID].artifacts.keys(), \
                'no preview'
            assert thumb_uri in obs.planes[TEST_PRODUCT_ID].artifacts, \
                'no thumbnail'
    finally:
        if os.path.exists(test_prev):
            os.unlink(test_prev)


def test_preview_augment_known_no_preview():
    # rejected file exists that says there's a preview known to not
    # exist, so trying to generate a thumbnail will result in no
    # change to the plane/artifact structure

    try:
        obs = mc.read_obs_from_file(TEST_OBS_FILE)
        obs.planes[TEST_PRODUCT_ID].data_release = datetime.utcnow()
        assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, \
            'initial condition'

        if os.path.exists(REJECTED_FILE):
            os.unlink(REJECTED_FILE)
        test_rejected = mc.Rejected(REJECTED_FILE)
        test_rejected.record(
            mc.Rejected.NO_PREVIEW, '{}.jpg'.format(TEST_PRODUCT_ID))
        test_config = mc.Config()
        test_observable = mc.Observable(
            test_rejected, mc.Metrics(test_config))

        cadc_client_mock = Mock()
        kwargs = {'working_directory': TEST_DATA_DIR,
                  'cadc_client': cadc_client_mock,
                  'stream': 'stream',
                  'observable': test_observable}

        with patch('caom2pipe.manage_composable.http_get') as http_mock, \
                patch('caom2pipe.manage_composable.data_put') as ad_put_mock, \
                patch('caom2pipe.manage_composable.get_artifact_metadata') as \
                art_mock, \
                patch('caom2pipe.manage_composable.exec_cmd') as exec_mock:
            cadc_client_mock.return_value.data_get.return_value = mc.CadcException(
                'test')
            result = preview_augmentation.visit(obs, **kwargs)
            assert not http_mock.called, 'http mock should not be called'
            assert not ad_put_mock.called, 'ad put mock should not be called'
            assert not art_mock.called, 'art mock should not be called'
            assert not exec_mock.called, 'exec mock should not be called'
            assert result is not None, 'expect a result'
            assert result['artifacts'] == 0, 'no artifacts should be updated'
            assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, \
                'no new artifacts'

        test_rejected.persist_state()
        assert os.path.exists(REJECTED_FILE)
    finally:
        if os.path.exists(REJECTED_FILE):
            os.unlink(REJECTED_FILE)


def test_preview_augment_unknown_no_preview():
    # what happens when it's not known that there's no preview
    obs = mc.read_obs_from_file(TEST_OBS_FILE)
    obs.planes[TEST_PRODUCT_ID].data_release = datetime.utcnow()
    assert len(obs.planes[TEST_PRODUCT_ID].artifacts) == 1, 'initial condition'

    # make sure the rejected file is empty
    if os.path.exists(REJECTED_FILE):
        os.unlink(REJECTED_FILE)
    test_rejected = mc.Rejected(REJECTED_FILE)
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))

    cadc_client_mock = Mock()
    kwargs = {'working_directory': TEST_DATA_DIR,
              'cadc_client': cadc_client_mock,
              'stream': 'stream',
              'observable': test_observable}

    with patch('caom2pipe.manage_composable.http_get',
               side_effect=mc.CadcException(
                   '404 Client Error: Not Found for url')) as http_mock, \
            patch('caom2pipe.manage_composable.data_put') as ad_put_mock, \
            patch('caom2pipe.manage_composable.get_artifact_metadata') as \
                art_mock, \
            patch('caom2pipe.manage_composable.exec_cmd') as exec_mock:
        cadc_client_mock.return_value.data_get.return_value = mc.CadcException(
            'test')
        result = preview_augmentation.visit(obs, **kwargs)
        assert result is not None, 'expect result'
        # 0 because the preview artifact doesn't already exist
        assert result['artifacts'] == 0, 'wrong result'
        test_url = '{}{}.fits'.format(preview_augmentation.PREVIEW_URL,
                                      TEST_PRODUCT_ID)
        test_prev = '{}/{}.jpg'.format(TEST_DATA_DIR, TEST_PRODUCT_ID)
        http_mock.assert_called_with(test_url, test_prev),  'mock not called'
        assert not ad_put_mock.called, 'ad put mock should not be called'
        assert not art_mock.called, 'art mock should not be called'
        assert not exec_mock.called, 'exec mock should not be called'


def test_preview_augment_delete_preview():
    # plane starts with a preview artifact, but it represents a non-existent
    # file, so remove the artifact from the CAOM observation
    test_product_id = 'S20080610S0045'
    fqn = os.path.join(TEST_DATA_DIR, 'GS-2008A-C-5-35-002.fits.xml')
    obs = mc.read_obs_from_file(fqn)
    assert len(obs.planes[test_product_id].artifacts) == 2, 'initial condition'
    test_rejected = mc.Rejected('/tmp/nonexistent')
    test_rejected.content = {
        'bad_metadata': [],
        'no_preview':
            ['S20080610S0043.jpg',
             'S20080610S0041.jpg',
             'S20080610S0044.jpg',
             'S20080610S0045.jpg']}
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))
    kwargs = {'working_directory': TEST_DATA_DIR,
              'cadc_client': None,
              'stream': 'stream',
              'observable': test_observable}
    result = preview_augmentation.visit(obs, **kwargs)
    assert result is not None, 'expect a result'
    assert result['artifacts'] == 1, 'wrong result'
    assert len(obs.planes[test_product_id].artifacts) == 1, 'post condition'
