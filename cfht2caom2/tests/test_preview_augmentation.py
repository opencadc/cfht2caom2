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

import glob
import os
import pytest

from datetime import datetime
from mock import patch, Mock

from caom2 import ChecksumURI, Artifact, ReleaseType, ProductType
from cfht2caom2 import preview_augmentation, ARCHIVE
from cfht2caom2 import metadata as md
from cfht2caom2 import cfht_name
from caom2pipe import manage_composable as mc

import test_main_app

TEST_FILES_DIR = '/test_files'
REJECTED_FILE = os.path.join(test_main_app.TEST_DATA_DIR, 'rejected.yml')


@patch('caom2pipe.manage_composable.data_put')
def test_preview_augment(ad_put_mock):

    # this should result in three new artifacts being added to every plane:
    # one for a thumbnail and two for previews (one zoom)

    test_rejected = mc.Rejected(REJECTED_FILE)
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))
    cadc_client_mock = Mock()

    test_files = {
        'visit_obs_start_wircam.xml':
            ['1151210o.fits.fz', '1151210s.fits.fz', '1151210m.fits.fz',
             '1151210p.fits.fz', '1151210w.fits.fz', '1151210y.fits.fz',
             '1151210g.fits.fz'],
        'visit_obs_start_megacam_sci.xml':
            ['1927963f.fits.fz', '1927963o.fits.fz', '1927963p.fits.fz'],
        'visit_obs_start_megacam_cal.xml':
            ['979412b.fits.fz', '979412o.fits.fz', '979412p.fits.fz'],
        'visit_obs_start_sitelle_calibrated_cube.xml':
            ['2359320p.fits'],
        'visit_obs_start_sitelle.xml':
            ['2359320o.fits.fz'],
        'visit_obs_start_espadons.xml':
            ['2460606i.fits.gz', '2460606o.fits.gz'],
        'visit_obs_start_espadons_cal.xml':
            ['1001063b.fits.gz'],
        'visit_obs_start_spirou.xml': ['2401734o.fits', '2401734e.fits',
                                       '2401734r.fits', '2401734s.fits',
                                       '2401734t.fits']
                                       # '2401734t.fits', '2401734v.fits']
    }

    test_checksums = {
        'ad:CFHT/2460606i_preview_1024.jpg':
            'md5:931f8219dd742c05a02aaeed9bc80bae',
        'ad:CFHT/2460606i_preview_256.jpg':
            'md5:bd829d7f5b0b12924894c1cd97c1ccb7',
        'ad:CFHT/2460606o_preview_1024.jpg':
            'md5:eeeb436f0e6b1cb48aca5a9ba5d1ab3a',
        'ad:CFHT/2460606o_preview_256.jpg':
            'md5:3879c153776ee0569ee33bc8a35b5e75',
        'ad:CFHT/2460606o_preview_zoom_1024.jpg':
            'md5:dfe23ed11655ed3597ffb3557dbd2f63',
        'ad:CFHT/1001063b_preview_1024.jpg':
            'md5:d457fa8ecdd397790de65f8bd3a80d06',
        'ad:CFHT/1001063b_preview_256.jpg':
            'md5:3ae90e1e120fb97216ddb2d965e9423e',
        'ad:CFHT/1001063b_preview_zoom_1024.jpg':
            'md5:be1a546f858800746348419e5e4d1ec0',
        'ad:CFHT/979412b_preview_1024.jpg':
            'md5:8b3e4ff2c210b3861748423339a71cc5',
        'ad:CFHT/979412b_preview_256.jpg':
            'md5:893a9cd18f66fab907573bf1e664148d',
        'ad:CFHT/979412b_preview_zoom_1024.jpg':
            'md5:194750cb8d718c2572ae8c41c0a01d89',
        'ad:CFHT/979412o_preview_1024.jpg':
            'md5:379a7b21022a0e758fbdea6130c07233',
        'ad:CFHT/979412o_preview_256.jpg':
            'md5:59c2ecedda2a2560cafdca543533e190',
        'ad:CFHT/979412o_preview_zoom_1024.jpg':
            'md5:706e5b9e175ff4998c4fce34e14a255f',
        'ad:CFHT/979412p_preview_1024.jpg':
            'md5:d449d26987a4d9e220d249f8573025ee',
        'ad:CFHT/979412p_preview_256.jpg':
            'md5:102be79a6ee26a7e5c2ee7e40623d452',
        'ad:CFHT/979412p_preview_zoom_1024.jpg':
            'md5:a26455916377b90aa94a4650a05cd69e',
        'ad:CFHT/1927963f_preview_1024.jpg':
            'md5:9afd99f2d60de3bc168832ce37200508',
        'ad:CFHT/1927963f_preview_256.jpg':
            'md5:5de596f4b5ab1ae8bfab1de5550ec56d',
        'ad:CFHT/1927963f_preview_zoom_1024.jpg':
            'md5:87ab244b836390c2bfa761349fb33301',
        'ad:CFHT/1927963o_preview_1024.jpg':
            'md5:d8a9c9c6dde50c59f75a1b384300d41a',
        'ad:CFHT/1927963o_preview_256.jpg':
            'md5:84cfc76e094678e61eaa9f9421e30dd9',
        'ad:CFHT/1927963o_preview_zoom_1024.jpg':
            'md5:6fa913f58470c6da13a02aa808bbd278',
        'ad:CFHT/1927963p_preview_1024.jpg':
            'md5:61b903863106438cfa2771731ed8690d',
        'ad:CFHT/1927963p_preview_256.jpg':
            'md5:01e4de796ff3f20949d02fe776b5f699',
        'ad:CFHT/1927963p_preview_zoom_1024.jpg':
            'md5:fed58d46c3ead044d047f98818c3452e',
        'ad:CFHT/2359320o_preview_1024.jpg':
            'md5:45b164799c67751c5a90255c11ccd9e5',
        'ad:CFHT/2359320o_preview_256.jpg':
            'md5:dce23df89c3424bc67d733a7784a13d3',
        'ad:CFHT/2359320o_preview_zoom_1024.jpg':
            'md5:6379ffcd7369a3cb20249223cde7642b',
        # 'ad:CFHT/2359320p_preview_1024.jpg':  ####
        #     'md5:252c08c95acb5bacf973685b22ad2d86',
        # 'ad:CFHT/2359320p_preview_256.jpg':
        #     'md5:b0a999abfca69cf40b157ee1f65e0a55',
        # 'ad:CFHT/2359320p_preview_zoom_1024.jpg':
        #     'md5:0d027b4d2feee6073344490f65e573df',
        'ad:CFHT/2401734e_preview_1024.jpg':
            'md5:11d15d1bbe9ccf618559e0221646b346',
        'ad:CFHT/2401734e_preview_256.jpg':
            'md5:3d9e4c1f59c8aa76f9fe6613ed3c54b2',
        'ad:CFHT/2401734o_preview_1024.jpg':
            'md5:1b84806bcd9ed45e850afa802499e239',
        'ad:CFHT/2401734o_preview_256.jpg':
            'md5:b1e8aa9b1f0275505660435d3ee6443e',
        'ad:CFHT/2401734o_preview_zoom_1024.jpg':
            'md5:397ecc64a367226d73970d7a6e992f47',
        'ad:CFHT/2401734r_preview_1024.jpg':
            'md5:3c0dc874898cefe2ca46195e2e19b763',
        'ad:CFHT/2401734r_preview_256.jpg':
            'md5:3b06668db0e8d9819cc1c50100bd85f4',
        'ad:CFHT/2401734r_preview_zoom_1024.jpg':
            'md5:ab7dbc11c15b63420a61da45e81d86dd',
        'ad:CFHT/2401734s_preview_1024.jpg':
            'md5:22f0abc0b9b1fce9511fd7ab990925fc',
        'ad:CFHT/2401734s_preview_256.jpg':
            'md5:b07abd8a3f8d50bfd36b30b7f1289883',
        'ad:CFHT/2401734t_preview_1024.jpg':
            'md5:6872763c63f0dfa9f027941319714e1b',
        'ad:CFHT/2401734t_preview_256.jpg':
            'md5:d9ee166527bc17b70ee4a66dd8b73d05',
        'ad:CFHT/1151210g_preview_1024.jpg':
            'md5:cfe4f124a401e81256c040f4b8345cd2',
        'ad:CFHT/1151210g_preview_256.jpg':
            'md5:02c9421541356b8852734ba36870686a',
        'ad:CFHT/1151210g_preview_zoom_1024.jpg':
            'md5:fc76e8305754b61d8f1e94d26834a2ab',
        'ad:CFHT/1151210m_preview_1024.jpg':
            'md5:02abd143e36acd593cab8eeea458eaf2',
        'ad:CFHT/1151210m_preview_256.jpg':
            'md5:99d5b739c4ab1ef69aa0b34a0696f9cf',
        'ad:CFHT/1151210m_preview_zoom_1024.jpg':
            'md5:63c8499bd8caec17fca13395ca7feb77',
        'ad:CFHT/1151210w_preview_1024.jpg':
            'md5:542b11ae4e21edc19380e0b37dec129b',
        'ad:CFHT/1151210w_preview_256.jpg':
            'md5:a10b7c6b8f101f52ab16d0da0649451a',
        'ad:CFHT/1151210w_preview_zoom_1024.jpg':
            'md5:bf5c787cd7590ee28f7242de3dd7de50',
        'ad:CFHT/1151210o_preview_1024.jpg':
            'md5:302ce34edd7c69da23119d001e85d39a',
        'ad:CFHT/1151210o_preview_256.jpg':
            'md5:137fa8a405cdf03f827000908831ece9',
        'ad:CFHT/1151210o_preview_zoom_1024.jpg':
            'md5:135ec5d59f6f21ad15f68a3268bb91af',
        'ad:CFHT/1151210p_preview_1024.jpg':
            'md5:f05cb0ad85ac610ef9bc3cf1ea00c853',
        'ad:CFHT/1151210p_preview_256.jpg':
            'md5:73aa50e0d59a69452ae8be06de9c6317',
        'ad:CFHT/1151210p_preview_zoom_1024.jpg':
            'md5:4706d8b8147dfb51c4211f94b9b1f298',
        'ad:CFHT/1151210s_preview_1024.jpg':
            'md5:7d7019efd828ff748fece2b778c5fecb',
        'ad:CFHT/1151210s_preview_256.jpg':
            'md5:43051edaa8dbb6436d4dcb32fc5a2bd0',
        'ad:CFHT/1151210s_preview_zoom_1024.jpg':
            'md5:36a993e92d1120e0c8444bc4fa1fd0b6',
        'ad:CFHT/1151210y_preview_1024.jpg':
            'md5:931dd9cabffe9e470a556327c08c6b96',
        'ad:CFHT/1151210y_preview_256.jpg':
            'md5:45df1e0a5dd4b4b058b7e20cdc55b9f7',
        'ad:CFHT/1151210y_preview_zoom_1024.jpg':
            'md5:009f7aa763d9296341358af372fc167f'
    }
    kwargs = {'working_directory': TEST_FILES_DIR,
              'cadc_client': cadc_client_mock,
              'stream': 'stream',
              'observable': test_observable}

    for entry in glob.glob(f'{TEST_FILES_DIR}/*.jpg'):
        os.unlink(entry)

    for key, value in test_files.items():
        obs = mc.read_obs_from_file(f'{test_main_app.TEST_DATA_DIR}/{key}')
        if 'wircam' in key:
            instrument = md.Inst.WIRCAM
        elif 'mega' in key:
            instrument = md.Inst.MEGACAM
        elif 'sitelle' in key:
            instrument = md.Inst.SITELLE
        elif 'espadons' in key:
            instrument = md.Inst.ESPADONS
        elif 'spirou' in key:
            instrument = md.Inst.SPIROU
        else:
            assert False, 'do not understand instrument'
        for f_name in value:
            kwargs['science_file'] = f_name

            test_name = cfht_name.CFHTName(file_name=f_name,
                                           instrument=instrument)
            check_number = 1
            if test_name.suffix == 'p' and instrument is md.Inst.SITELLE:
                check_number = 3
            assert len(obs.planes[test_name.product_id].artifacts) == \
                check_number, f'initial condition {f_name}'

            try:
                test_result = preview_augmentation.visit(obs, **kwargs)
            except Exception as e:
                import logging
                logging.error(e)
                import traceback
                logging.error(traceback.format_exc())
                assert False

            assert test_result is not None, f'expect a result {f_name}'

            check_number = 3
            end_artifact_count = 4
            expected_call_count = 3
            f_name_list = [test_name.prev_uri, test_name.thumb_uri,
                           test_name.zoom_uri]
            if ((instrument is md.Inst.ESPADONS and test_name.suffix == 'i') or
                    (instrument is md.Inst.SPIROU and
                     test_name.suffix in ['e', 'p', 's', 't'])):
                check_number = 2
                expected_call_count = 2
                end_artifact_count = 3
                f_name_list = [test_name.prev_uri, test_name.thumb_uri]

            assert test_result['artifacts'] == check_number, \
                f'artifacts should be added {f_name}'

            mc.write_obs_to_file(obs, f'/test_files/act.{f_name}.xml')
            if test_name.suffix == 'p' and instrument is md.Inst.SITELLE:
                end_artifact_count = 6
            assert len(obs.planes[test_name.product_id].artifacts) == \
                end_artifact_count, f'new artifacts {f_name}'

            for p in f_name_list:
                assert p in \
                       obs.planes[test_name.product_id].artifacts.keys(), \
                       f'no preview {p}'
                artifact = obs.planes[test_name.product_id].artifacts[p]
                # because 2359320p_preview_1024 keeps changing ....
                # if artifact.uri in test_checksums:
                #     assert artifact.content_checksum.uri == test_checksums[p], \
                #         f'wrong checksum {p} {artifact.content_checksum} ' \
                #         f'{test_checksums[p]}'

            assert ad_put_mock.called, f'ad put mock not called {f_name}'
            assert ad_put_mock.call_count == expected_call_count, \
                f'ad put called wrong number of times {f_name}'
            ad_put_mock.reset_mock()
    # assert False
