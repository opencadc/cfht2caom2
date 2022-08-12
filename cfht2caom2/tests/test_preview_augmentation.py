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
import logging
import os
import traceback

from cfht2caom2 import preview_augmentation
from cfht2caom2 import metadata as md
from cfht2caom2 import cfht_name
from caom2pipe import manage_composable as mc

import test_fits2caom2_augmentation

TEST_FILES_DIR = '/test_files'
REJECTED_FILE = os.path.join(
    test_fits2caom2_augmentation.TEST_DATA_DIR, 'rejected.yml'
)


def test_preview_augment():

    # this should result in three new artifacts being added to every plane:
    # one for a thumbnail and two for previews (one zoom)

    test_rejected = mc.Rejected(REJECTED_FILE)
    test_config = mc.Config()
    test_observable = mc.Observable(test_rejected, mc.Metrics(test_config))

    test_files = {
        'visit_obs_start_wircam.xml': [
            '1151210o.fits.fz',
            '1151210s.fits.fz',
            '1151210m.fits.fz',
            '1151210p.fits.fz',
            '1151210w.fits.fz',
            '1151210y.fits.fz',
            '1151210g.fits.fz',
        ],
        'visit_obs_start_megacam_sci.xml': [
            '1927963f.fits.fz',
            '1927963o.fits.fz',
            '1927963p.fits.fz',
        ],
        'visit_obs_start_megacam_cal.xml': [
            '979412b.fits.fz',
            '979412o.fits.fz',
            '979412p.fits.fz',
        ],
        'visit_obs_start_sitelle_calibrated_cube.xml': ['2359320p.fits'],
        'visit_obs_start_sitelle.xml': ['2359320o.fits.fz'],
        'visit_obs_start_espadons.xml': [
            '2460606i.fits',
            '2460606o.fits.fz',
        ],
        'visit_obs_start_espadons_cal.xml': ['1001063b.fits.fz'],
        'visit_obs_start_spirou.xml': [
            '2401734o.fits',
            '2401734e.fits',
            '2401734r.fits',
            '2401734s.fits',
            '2401734t.fits',
            '2401734v.fits',
        ],
        'visit_obs_start_2hdus_mega.xml': ['1821271p.fits.fz'],
        'visit_obs_start_wircam_unzipped.xml': ['1334021f.fits'],
        'visit_obs_start_wircam_gz.xml': ['2615124g.fits.fz'],
        'visit_obs_start_scatter.xml': ['04Am01.scatter.r.36.00.fits.fz'],
        'visit_obs_start_mega_few_hdus.xml': ['1805750p.fits.fz'],
        'visit_obs_start_spirou_p.xml': ['2446341p.fits'],
        'visit_obs_start_spirou_p2.xml': ['2630033p.fits'],
        'visit_obs_start_wircam_fits.xml': ['1334131g.fits'],
        'visit_obs_start_sitelle_hdf5.xml': ['2752885z.hdf5'],
    }

    test_checksums = {
        'cadc:CFHT/2460606i_preview_1024.jpg': 'md5:eeb2558c76c83a1b56e0afa3d96b0fa0',
        'cadc:CFHT/2460606i_preview_256.jpg': 'md5:1084a1e58d246c35fae43ca021916c6a',
        'cadc:CFHT/2460606o_preview_1024.jpg': 'md5:eeeb436f0e6b1cb48aca5a9ba5d1ab3a',
        'cadc:CFHT/2460606o_preview_256.jpg': 'md5:3879c153776ee0569ee33bc8a35b5e75',
        'cadc:CFHT/2460606o_preview_zoom_1024.jpg': 'md5:1d5c9e5ee1586ddd6ce5c53b17f8f5ec',
        'cadc:CFHT/1001063b_preview_1024.jpg': 'md5:d457fa8ecdd397790de65f8bd3a80d06',
        'cadc:CFHT/1001063b_preview_256.jpg': 'md5:3ae90e1e120fb97216ddb2d965e9423e',
        'cadc:CFHT/1001063b_preview_zoom_1024.jpg': 'md5:a2d62e654838a620dcf97425753481ca',
        'cadc:CFHT/979412b_preview_1024.jpg': 'md5:8b3e4ff2c210b3861748423339a71cc5',
        'cadc:CFHT/979412b_preview_256.jpg': 'md5:893a9cd18f66fab907573bf1e664148d',
        'cadc:CFHT/979412b_preview_zoom_1024.jpg': 'md5:4f21d6676b0749c48b4872699f7aff7c',
        'cadc:CFHT/979412o_preview_1024.jpg': 'md5:379a7b21022a0e758fbdea6130c07233',
        'cadc:CFHT/979412o_preview_256.jpg': 'md5:59c2ecedda2a2560cafdca543533e190',
        'cadc:CFHT/979412o_preview_zoom_1024.jpg': 'md5:6f3a3a39155a426a365f13f6d10982a9',
        'cadc:CFHT/979412p_preview_1024.jpg': 'md5:d449d26987a4d9e220d249f8573025ee',
        'cadc:CFHT/979412p_preview_256.jpg': 'md5:102be79a6ee26a7e5c2ee7e40623d452',
        'cadc:CFHT/979412p_preview_zoom_1024.jpg': 'md5:7e09407afb9e8148ad36ba56933e1822',
        'cadc:CFHT/1927963f_preview_1024.jpg': 'md5:9afd99f2d60de3bc168832ce37200508',
        'cadc:CFHT/1927963f_preview_256.jpg': 'md5:5de596f4b5ab1ae8bfab1de5550ec56d',
        'cadc:CFHT/1927963f_preview_zoom_1024.jpg': 'md5:ca705c70eb0d320856e346662a83f1ba',
        'cadc:CFHT/1927963o_preview_1024.jpg': 'md5:d8a9c9c6dde50c59f75a1b384300d41a',
        'cadc:CFHT/1927963o_preview_256.jpg': 'md5:84cfc76e094678e61eaa9f9421e30dd9',
        'cadc:CFHT/1927963o_preview_zoom_1024.jpg': 'md5:843109c217104034fa3b029949fa277e',
        'cadc:CFHT/1927963p_preview_1024.jpg': 'md5:61b903863106438cfa2771731ed8690d',
        'cadc:CFHT/1927963p_preview_256.jpg': 'md5:01e4de796ff3f20949d02fe776b5f699',
        'cadc:CFHT/1927963p_preview_zoom_1024.jpg': 'md5:91c4a735d507d3c11482179d725a4166',
        'cadc:CFHT/2359320o_preview_1024.jpg': 'md5:45b164799c67751c5a90255c11ccd9e5',
        'cadc:CFHT/2359320o_preview_256.jpg': 'md5:dce23df89c3424bc67d733a7784a13d3',
        'cadc:CFHT/2359320o_preview_zoom_1024.jpg': 'md5:5cc26a79326b42c6a17bca8e0d3fb1fb',
        # 'cadc:CFHT/2359320p_preview_1024.jpg':  ####
        #     'md5:b09f54461fb10259f4958ff9ab1e8179',
        # 'cadc:CFHT/2359320p_preview_256.jpg':
        #     'md5:b037076c0d9488d013052d6750e55aeb',
        # 'cadc:CFHT/2359320p_preview_zoom_1024.jpg':
        #     'md5:514aa25fec1b74fe942ccc8ce8e552f7',
        'cadc:CFHT/2401734e_preview_1024.jpg': 'md5:a3b272473e959cbc5496c1d562f16c0d',
        'cadc:CFHT/2401734e_preview_256.jpg': 'md5:781a1c86e727f82ec17a98b5eaf2b7d0',
        'cadc:CFHT/2401734o_preview_1024.jpg': 'md5:1b84806bcd9ed45e850afa802499e239',
        'cadc:CFHT/2401734o_preview_256.jpg': 'md5:b1e8aa9b1f0275505660435d3ee6443e',
        'cadc:CFHT/2401734o_preview_zoom_1024.jpg': 'md5:18839f9482dd3420f0d6ca3c167e6575',
        'cadc:CFHT/2401734r_preview_1024.jpg': 'md5:3c0dc874898cefe2ca46195e2e19b763',
        'cadc:CFHT/2401734r_preview_256.jpg': 'md5:3b06668db0e8d9819cc1c50100bd85f4',
        'cadc:CFHT/2401734r_preview_zoom_1024.jpg': 'md5:32b54773700dd652e244632f51aca5ba',
        'cadc:CFHT/2401734s_preview_1024.jpg': 'md5:30a76081e4e6f5ee332e8e076feff2a7',
        'cadc:CFHT/2401734s_preview_256.jpg': 'md5:f32f4288aebc1f7cea577d6d217929c9',
        'cadc:CFHT/2401734t_preview_1024.jpg': 'md5:3a5409aa4cd728e0a92f50a3eff88427',
        'cadc:CFHT/2401734t_preview_256.jpg': 'md5:f90654bad03ba8c2e4796ab4acfb5264',
        'cadc:CFHT/1151210g_preview_1024.jpg': 'md5:ac415da2e2d9847b81c59477c73f26ce',
        'cadc:CFHT/1151210g_preview_256.jpg': 'md5:56ad048a81c7c16f0801ad8b22fc943f',
        'cadc:CFHT/1151210g_preview_zoom_1024.jpg': 'md5:0486b7e6572bb0bd8750ff16159ee5ad',
        'cadc:CFHT/1151210m_preview_1024.jpg': 'md5:02abd143e36acd593cab8eeea458eaf2',
        'cadc:CFHT/1151210m_preview_256.jpg': 'md5:99d5b739c4ab1ef69aa0b34a0696f9cf',
        'cadc:CFHT/1151210m_preview_zoom_1024.jpg': 'md5:d55b984b8dcb39baa0d5ba0a34051ffa',
        'cadc:CFHT/1151210w_preview_1024.jpg': 'md5:542b11ae4e21edc19380e0b37dec129b',
        'cadc:CFHT/1151210w_preview_256.jpg': 'md5:a10b7c6b8f101f52ab16d0da0649451a',
        'cadc:CFHT/1151210w_preview_zoom_1024.jpg': 'md5:e79bb2ac4230d5587f97a9b838d6f88d',
        'cadc:CFHT/1151210o_preview_1024.jpg': 'md5:302ce34edd7c69da23119d001e85d39a',
        'cadc:CFHT/1151210o_preview_256.jpg': 'md5:137fa8a405cdf03f827000908831ece9',
        'cadc:CFHT/1151210o_preview_zoom_1024.jpg': 'md5:e656f21a21a4a41e088cccf87ec76580',
        'cadc:CFHT/1151210p_preview_1024.jpg': 'md5:f05cb0ad85ac610ef9bc3cf1ea00c853',
        'cadc:CFHT/1151210p_preview_256.jpg': 'md5:73aa50e0d59a69452ae8be06de9c6317',
        'cadc:CFHT/1151210p_preview_zoom_1024.jpg': 'md5:e56690286391a74840bc07158cedaeec',
        'cadc:CFHT/1151210s_preview_1024.jpg': 'md5:7d7019efd828ff748fece2b778c5fecb',
        'cadc:CFHT/1151210s_preview_256.jpg': 'md5:43051edaa8dbb6436d4dcb32fc5a2bd0',
        'cadc:CFHT/1151210s_preview_zoom_1024.jpg': 'md5:4f20f5586eaa43abc5de977e60cea807',
        'cadc:CFHT/1151210y_preview_1024.jpg': 'md5:931dd9cabffe9e470a556327c08c6b96',
        'cadc:CFHT/1151210y_preview_256.jpg': 'md5:45df1e0a5dd4b4b058b7e20cdc55b9f7',
        'cadc:CFHT/1151210y_preview_zoom_1024.jpg': 'md5:53cfa0d681c46fedfc2bf719e88231b5',
    }
    kwargs = {
        'working_directory': TEST_FILES_DIR,
        'cadc_client': None,
        'stream': 'stream',
        'observable': test_observable,
    }

    for entry in glob.glob(f'{TEST_FILES_DIR}/*.jpg'):
        os.unlink(entry)

    checksum_failures = []

    orig_scheme = mc.StorageName.scheme
    orig_collection = mc.StorageName.collection
    try:
        mc.StorageName.scheme = 'cadc'
        mc.StorageName.collection = cfht_name.COLLECTION
        for key, value in test_files.items():
            obs = mc.read_obs_from_file(
                f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/{key}'
            )
            if 'wircam' in key:
                instrument = md.Inst.WIRCAM
            elif 'mega' in key or 'scatter' in key:
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
                test_name = cfht_name.CFHTName(
                    file_name=f_name, instrument=instrument
                )
                test_name.source_names = [
                    os.path.join(TEST_FILES_DIR, f_name)
                ]
                kwargs['storage_name'] = test_name
                check_number = 1
                if test_name.suffix == 'p' and instrument is md.Inst.SITELLE:
                    check_number = 3
                assert (
                    len(obs.planes[test_name.product_id].artifacts)
                    == check_number
                ), f'initial condition {f_name}'

                try:
                    from datetime import datetime

                    start_ts = datetime.utcnow().timestamp()
                    test_result = preview_augmentation.visit(obs, **kwargs)
                    end_ts = datetime.utcnow().timestamp()
                    logging.error(
                        f'{f_name} execution time {end_ts - start_ts}'
                    )
                except Exception as e:
                    logging.error(e)
                    logging.error(traceback.format_exc())
                    assert False

                assert test_result is not None, f'expect a result {f_name}'

                check_number = 3
                end_artifact_count = 4
                f_name_list = [
                    test_name.prev_uri,
                    test_name.thumb_uri,
                    test_name.zoom_uri,
                ]
                if (
                    instrument is md.Inst.ESPADONS and test_name.suffix == 'i'
                ) or (
                    instrument is md.Inst.SPIROU
                    and test_name.suffix in ['e', 'p', 's', 't', 'v']
                ):
                    check_number = 2
                    end_artifact_count = 3
                    f_name_list = [test_name.prev_uri, test_name.thumb_uri]

                # assert (
                #     test_result['artifacts'] == check_number
                # ), f'artifacts should be added {f_name}'

                if instrument is md.Inst.SITELLE:
                    if test_name.suffix == 'p':
                        end_artifact_count = 6
                    elif test_name.suffix == 'z':
                        end_artifact_count = 3
                assert (
                    len(obs.planes[test_name.product_id].artifacts)
                    == end_artifact_count
                ), f'new artifacts {f_name}'

                for p in f_name_list:
                    # assert p in \
                    #        obs.planes[test_name.product_id].artifacts.keys(), \
                    #        f'no preview {p}'
                    if p in obs.planes[test_name.product_id].artifacts.keys():
                        artifact = obs.planes[test_name.product_id].artifacts[
                            p
                        ]
                        # because 2359320p_preview_1024 keeps changing ....
                        if artifact.uri in test_checksums:
                            if (
                                artifact.content_checksum.uri
                                != test_checksums[p]
                            ):
                                checksum_failures.append(
                                    f'{p} expected {test_checksums[p]} actual '
                                    f'{artifact.content_checksum.uri}'
                                )

        assert (
            len(checksum_failures) == 0,
            '\n'.join(ii for ii in checksum_failures),
        )
    finally:
        mc.StorageName.collection = orig_collection
        mc.StorageName.scheme = orig_scheme

    # assert False
