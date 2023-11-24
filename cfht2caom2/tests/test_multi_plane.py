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

import warnings

from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import FITSFixedWarning
from caom2pipe.manage_composable import Observable
from caom2pipe import reader_composable as rdc
from cfht2caom2 import cfht_name, fits2caom2_augmentation, metadata
from glob import glob
from os.path import basename, dirname

from mock import patch
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
@patch('cfht2caom2.instruments.get_local_headers_from_fits')
@patch('caom2utils.data_util.get_local_headers_from_fits')
@patch('caom2pipe.astro_composable.get_vo_table')
def test_visitor(
    vo_mock, local_headers_mock1, local_headers_mock2, cache_mock, expected_fqn, test_data_dir, test_config
):
    warnings.simplefilter('ignore', category=AstropyUserWarning)
    warnings.simplefilter('ignore', category=FITSFixedWarning)
    vo_mock.side_effect = test_caom_gen_visit._vo_mock
    # during cfht2caom2 operation, want to use astropy on FITS files
    # but during testing want to use headers and built-in Python file
    # operations
    local_headers_mock1.side_effect = test_caom_gen_visit._local_headers
    local_headers_mock2.side_effect = test_caom_gen_visit._local_headers
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
    observation = None
    for f_name in LOOKUP[test_name]:
        if 'hdf5' in f_name:
            source_names = [f'/test_files/{f_name}']
        else:
            source_names = [
                f'{test_caom_gen_visit.TEST_DATA_DIR}/multi_plane/{instr}/{f_name}.header'
            ]
        storage_name = cfht_name.CFHTName(
            file_name=f_name,
            instrument=instrument,
            source_names=source_names,
        )
        metadata_reader = rdc.FileMetadataReader()
        metadata_reader.set(storage_name)
        if 'hdf5' in f_name:
            metadata_reader.file_info[
                storage_name.file_uri
            ].file_type = 'application/x-hdf5'
        else:
            metadata_reader.file_info[
                storage_name.file_uri
            ].file_type = 'application/fits'
        test_config.rejected_file_name = 'rejected.yml'
        test_config.rejected_directory = test_data_dir
        test_observable = Observable(test_config)
        kwargs = {
            'storage_name': storage_name,
            'metadata_reader': metadata_reader,
            'observable': test_observable,
            'config': test_config,
        }

        observation = fits2caom2_augmentation.visit(observation, **kwargs)

    test_caom_gen_visit._compare(expected_fqn, observation, test_name)
