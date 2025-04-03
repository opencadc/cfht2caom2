# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
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

from glob import glob

from caom2utils.data_util import get_local_file_headers
from caom2pipe.manage_composable import StorageName
from cfht2caom2 import CFHTName


def test_is_valid(test_config):
    assert CFHTName(source_names=['1944968p.fits.fz'], instrument='SITELLE').is_valid()

    test_subject = CFHTName(source_names=['2463796o.fits.fz'], instrument='MegaCam')
    assert test_subject.obs_id == '2463796', 'wrong obs id'
    assert test_subject.file_id == '2463796o', 'wrong file id'
    assert (test_subject.file_uri == 'cadc:CFHT/2463796o.fits.fz'), 'wrong uri'
    assert test_subject.source_names == ['2463796o.fits.fz'], 'not local'
    assert test_subject.simple, f'should be simple {test_subject}'
    assert test_subject.suffix == 'o', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    StorageName.scheme = 'cadc'
    test_subject = CFHTName(source_names=['1944968p.fits.fz'], instrument='SITELLE')
    assert test_subject.obs_id == '1944968', 'wrong obs id'
    assert test_subject.file_id == '1944968p', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/1944968p.fits.fz', 'uri'
    assert not test_subject.simple, 'should be composite'
    assert test_subject.suffix == 'p', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    StorageName.scheme = 'cadc'
    test_subject = CFHTName(source_names=['2460503p.fits.gz'], instrument='ESPaDOnS')
    assert test_subject.obs_id == '2460503', 'wrong obs id'
    assert test_subject.file_id == '2460503p', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/2460503p.fits', 'wrong uri'
    assert not test_subject.simple, 'should be composite'
    assert test_subject.file_name == '2460503p.fits.gz', 'no decomp'
    assert test_subject.suffix == 'p', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['2452990p.fits.fz'], instrument='MegaPrime')
    assert test_subject.obs_id == '2452990', 'wrong obs id'
    assert test_subject.file_id == '2452990p', 'wrong file id'
    assert (test_subject.file_uri == 'cadc:CFHT/2452990p.fits.fz'), 'wrong uri'
    assert not test_subject.simple, 'should be derived'
    assert test_subject.suffix == 'p', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['2384125z.hdf5'], instrument='SITELLE')
    assert test_subject.obs_id == '2384125', 'wrong obs id'
    assert test_subject.file_id == '2384125z', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/2384125z.hdf5', 'wrong uri'
    assert not test_subject.simple, 'should be derived'
    assert test_subject.suffix == 'z', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['2384125p.fits.fz'], instrument='SITELLE')
    assert test_subject.obs_id == '2384125', 'wrong obs id'
    assert test_subject.file_id == '2384125p', 'wrong file id'
    assert (test_subject.file_uri == 'cadc:CFHT/2384125p.fits.fz'), 'wrong uri'
    assert not test_subject.simple, 'should be derived'
    assert test_subject.suffix == 'p', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['979412p.fits.fz'], instrument='MegaPrime')
    assert test_subject.obs_id == '979412', 'wrong obs id'
    assert test_subject.file_id == '979412p', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/979412p.fits.fz', 'wrong uri'
    assert not test_subject.simple, 'should be derived'
    assert test_subject.suffix == 'p', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['979412b.fits.fz'], instrument='MegaPrime')
    assert test_subject.obs_id == '979412', 'wrong obs id'
    assert test_subject.file_id == '979412b', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/979412b.fits.fz', 'wrong uri'
    assert test_subject.simple, 'should be simple'
    assert test_subject.suffix == 'b', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['2003A.frpts.z.36.00.fits.fz'], instrument='MegaPrime')
    assert test_subject.obs_id == '2003A.frpts.z.36.00', 'wrong obs id'
    assert test_subject.file_id == '2003A.frpts.z.36.00', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/2003A.frpts.z.36.00.fits.fz', 'wrong uri'
    assert not test_subject.simple, 'should not be simple'
    assert test_subject.derived, 'should be derived'
    assert test_subject.suffix is None, f'suffix {test_subject}'
    assert test_subject.sequence_number is None, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['02AE10.bias.0.36.00.fits'], instrument='MegaPrime')
    assert test_subject.obs_id == '02AE10.bias.0.36.00', 'wrong obs id'
    assert test_subject.file_id == '02AE10.bias.0.36.00', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/02AE10.bias.0.36.00.fits', 'wrong uri'
    assert not test_subject.simple, 'should not be simple'
    assert test_subject.derived, 'should be derived'
    assert test_subject.suffix is None, f'suffix {test_subject}'
    assert test_subject.sequence_number is None, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['2455409p.fits'], instrument='SPIRou')
    assert test_subject.obs_id == '2455409', 'wrong obs id'
    assert test_subject.file_id == '2455409p', 'wrong file id'
    assert test_subject.file_uri == 'cadc:CFHT/2455409p.fits', 'wrong uri'
    assert not test_subject.simple, 'should be derived'
    assert test_subject.suffix == 'p', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['2238502i.fits.fz'], instrument='ESPaDOnS')
    assert test_subject.obs_id == '2238502', 'wrong obs id'
    assert test_subject.suffix == 'i', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert not test_subject.simple, 'should not be simple'
    assert test_subject.derived, 'should be derived'
    assert test_subject.raw_time, 'raw time'

    StorageName.scheme = 'cadc'
    test_subject = CFHTName(
        source_names=['2602045r.fits.fz'],
        instrument='SPIRou',
        bitpix=-32,
    )
    assert test_subject.obs_id == '2602045', 'wrong obs id'
    assert test_subject.product_id == '2602045r', 'wrong product id'
    assert test_subject.file_uri == 'cadc:CFHT/2602045r.fits.fz', 'wrong file uri'
    assert test_subject.destination_uris[0] == 'cadc:CFHT/2602045r.fits.fz', 'wrong destination uri'
    assert test_subject.thumb_uri == 'cadc:CFHT/2602045r_preview_256.jpg', 'wrong thumb uri'
    assert test_subject.prev_uri == 'cadc:CFHT/2602045r_preview_1024.jpg', 'wrong preview uri'
    assert test_subject.zoom_uri == 'cadc:CFHT/2602045r_preview_zoom_1024.jpg', 'wrong zoom uri'
    assert not test_subject.has_different_destination_name, f'de/re'
    assert test_subject.suffix == 'r', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.simple, 'should be simple'
    assert not test_subject.derived, 'should not be derived'
    assert test_subject.raw_time, 'raw time'

    # decompression, no recompression
    StorageName.scheme = 'cadc'
    test_subject = CFHTName(
        source_names=['cadc:CFHT/2602045r.fits.gz'],
        instrument='SPIRou',
        bitpix=-32,
    )
    assert test_subject.file_uri == 'cadc:CFHT/2602045r.fits', 'wrong file uri'
    assert test_subject.destination_uris[0] == 'cadc:CFHT/2602045r.fits', 'wrong destination uri'
    assert test_subject.thumb_uri == 'cadc:CFHT/2602045r_preview_256.jpg', 'wrong thumb uri'
    assert test_subject.prev_uri == 'cadc:CFHT/2602045r_preview_1024.jpg', 'wrong preview uri'
    assert test_subject.zoom_uri == 'cadc:CFHT/2602045r_preview_zoom_1024.jpg', 'wrong zoom uri'
    assert test_subject.has_different_destination_name, 'de/re'
    assert test_subject.suffix == 'r', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    # decompression plus recompression
    test_subject = CFHTName(
        source_names=['2602045r.fits.gz'],
        instrument='SPIRou',
        bitpix=32,
    )
    assert test_subject.file_uri == 'cadc:CFHT/2602045r.fits.fz', 'wrong file uri'
    assert test_subject.destination_uris[0] == 'cadc:CFHT/2602045r.fits.fz', 'wrong destination uri'
    assert test_subject.thumb_uri == 'cadc:CFHT/2602045r_preview_256.jpg', 'wrong thumb uri'
    assert test_subject.prev_uri == 'cadc:CFHT/2602045r_preview_1024.jpg', 'wrong preview uri'
    assert test_subject.zoom_uri == 'cadc:CFHT/2602045r_preview_zoom_1024.jpg', 'wrong zoom uri'
    assert test_subject.has_different_destination_name, f'de/re'
    assert test_subject.suffix == 'r', f'suffix {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    # flag/diag
    test_subject = CFHTName(
        source_names=['1013552p_flag.fits.fz'], instrument='MegaPrime'
    )
    assert test_subject.file_uri == 'cadc:CFHT/1013552p_flag.fits.fz', 'wrong file uri'
    assert test_subject.destination_uris[0] == 'cadc:CFHT/1013552p_flag.fits.fz', 'wrong destination uri'
    assert test_subject.thumb_uri == 'cadc:CFHT/1013552p_flag_preview_256.jpg', 'wrong thumb uri'
    assert test_subject.suffix == 'p', 'wrong suffix'
    assert test_subject.obs_id == '1013552', f'wrong obs id {test_subject}'
    assert not test_subject.simple, f'should not be simple {test_subject}'
    assert test_subject.derived, f'should be derived {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    # flag/diag
    test_subject = CFHTName(source_names=['695816p_diag.fits'], instrument='MegaPrime')
    assert test_subject.file_uri == 'cadc:CFHT/695816p_diag.fits', 'wrong file uri'
    assert test_subject.destination_uris[0] == 'cadc:CFHT/695816p_diag.fits', 'wrong destination uri'
    assert test_subject.thumb_uri == 'cadc:CFHT/695816p_preview_256.jpg', 'wrong thumb uri'
    assert test_subject.suffix == 'p', 'wrong suffix'
    assert test_subject.obs_id == '695816', 'wrong obs id'
    assert test_subject.product_id == '695816p', 'wrong product id'
    assert test_subject.simple, f'should be simple {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert not test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['2513728g.fits'], instrument='SPIRou')
    assert test_subject.file_uri == 'cadc:CFHT/2513728g.fits', 'wrong file uri'
    assert test_subject.destination_uris[0] == 'cadc:CFHT/2513728g.fits', 'wrong destination uri'
    assert test_subject.thumb_uri == 'cadc:CFHT/2513728g_preview_256.jpg', 'wrong thumb uri'
    assert test_subject.suffix == 'g', 'wrong suffix'
    assert test_subject.obs_id == '2513728', 'wrong obs id'
    assert test_subject.product_id == '2513728g', 'wrong product id'
    assert test_subject.simple, f'should be simple {test_subject}'
    assert test_subject.sequence_number == test_subject.obs_id, f'sequence number {test_subject}'
    assert test_subject.raw_time, 'raw time'

    test_subject = CFHTName(source_names=['695816p_diag.fits'], instrument='MegaPrime')
    assert test_subject.obs_id == '695816', 'wrong obs id'
    assert test_subject.product_id == '695816p', 'wrong product id'
    assert test_subject.suffix == 'p', 'wrong suffix'


def test_suffixes(test_data_dir, test_config):
    # ensure every test file can be identified as Simple or Derived
    for plane_name in ['single_plane', 'multi_plane']:
        for instrument in ['espadons', 'mega', 'sitelle', 'spirou' ,'wircam']:
            test_config.use_local_files = True

            plane_list = glob(f'{test_data_dir}/{plane_name}/{instrument}/*.header')
            for entry in plane_list:
                test_subject = CFHTName(source_names=[entry])
                assert test_subject is not None, 'ctor'
                test_subject.metadata = {test_subject.file_uri: get_local_file_headers(entry)}
                test_subject.set_metadata()
                found_one = False
                if test_subject.simple:
                    assert not test_subject.derived, f'not derived {test_subject}'
                    found_one = True
                if test_subject.derived:
                    assert not test_subject.simple, f'not simple {test_subject}'
                    found_one = True
                assert found_one, f'{entry} neither derived nor simple {test_subject}'
