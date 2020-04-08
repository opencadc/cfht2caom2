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
            ['2359320o.fits.fz', '2359320v.fits.fz']
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
        else:
            assert False, 'do not understand instrument'
        for f_name in value:
            kwargs['science_file'] = f_name

            test_name = cfht_name.CFHTName(file_name=f_name,
                                           instrument=instrument)
            check_number = 1
            if test_name.suffix == 'g':
                check_number = 4
            elif test_name.suffix == 'p' and instrument is md.Inst.SITELLE:
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
            if test_name.suffix == 'p' and instrument is md.Inst.WIRCAM:
                check_number = 6
            elif test_name.suffix == 'g':
                check_number = 0
            assert test_result['artifacts'] == check_number, \
                f'artifacts should be added {f_name}'
            assert len(obs.planes[test_name.product_id].artifacts) == 4, \
                f'new artifacts {f_name}'

            if test_name.suffix == 'g':
                assert not ad_put_mock.called, f'ad put mock called for g'
            else:
                for p in [test_name.prev_uri, test_name.thumb_uri,
                          test_name.zoom_uri]:
                    assert p in \
                           obs.planes[test_name.product_id].artifacts.keys(), \
                           f'no preview {p}'

                assert ad_put_mock.called, f'ad put mock not called {f_name}'
                assert ad_put_mock.call_count == 3, \
                    f'ad put called wrong number of times {f_name}'
                ad_put_mock.reset_mock()
