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

from os.path import join

from caom2pipe.manage_composable import read_obs_from_file, StorageName
from cfht2caom2 import espadons_energy_augmentation
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import metadata as md

import test_fits2caom2_augmentation

TEST_FILES_DIR = '/test_files'


def test_visit(test_config):
    product_id = '2460606i'
    f_name = f'{product_id}.fits.gz'
    obs_fqn = f'{test_fits2caom2_augmentation.TEST_DATA_DIR}/multi_plane/espadons/2460606.expected.xml'
    obs = read_obs_from_file(obs_fqn)

    # pre-conditions
    uri = f'cadc:CFHT/{product_id}.fits'
    assert (
        obs.planes[product_id].artifacts[uri].parts['0'].chunks[0].energy
        is None
    ), 'expect to assign'
    test_storage_name = cn.CFHTName(
        file_name=f_name, instrument=md.Inst.ESPADONS
    )
    test_storage_name.source_names = [
        join(TEST_FILES_DIR, f_name).replace('.gz', ''),
    ]
    kwargs = {
        'storage_name': test_storage_name,
        'working_directory': TEST_FILES_DIR,
    }
    cn.cfht_names[uri] = test_storage_name

    test_obs = espadons_energy_augmentation.visit(obs, **kwargs)
    assert test_obs is not None, 'expect a result'
    # assert test_result.get('chunks') == 1, 'expect 1 updated chunk'
    test_reference = (
        obs.planes[product_id].artifacts[uri].parts['0'].chunks[0]
    )
    assert test_reference is not None, 'expect to assign'
    assert test_reference.energy is not None, 'expect to assign energy'
    assert test_reference.naxis == 2, 'wrong naxis'
    assert test_reference.energy_axis == 1, 'wrong energy axis'
    assert test_reference.observable_axis == 2, 'wrong observable axis'
    assert test_reference.position_axis_1 is None, 'wrong position 1 axis'
    assert test_reference.position_axis_2 is None, 'wrong position 2 axis'
    assert test_reference.time_axis is None, 'wrong time axis'
    assert test_reference.custom_axis is None, 'wrong custom axis'
    assert test_reference.polarization_axis is None, 'wrong pol axis'
