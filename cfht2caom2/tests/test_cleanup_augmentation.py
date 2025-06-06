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
#  : 4 $
#
# ***********************************************************************
#

from os import path
from caom2pipe.caom_composable import get_all_artifact_keys
from caom2pipe import manage_composable as mc
from cfht2caom2 import cleanup_augmentation, cfht_name, Inst


def test_cleanup_augmentation(test_data_dir):
    test_obs = mc.read_obs_from_file(path.join(test_data_dir, 'visit/visit_obs_start_cleanup.xml'))
    assert len(test_obs.planes) == 3, 'initial conditions failed'
    storage_name = cfht_name.CFHTName(source_names=['abc.fits.fz'])
    kwargs = {'storage_name': storage_name}
    test_obs = cleanup_augmentation.visit(test_obs, **kwargs)
    assert len(test_obs.planes) == 2, 'post-test conditions failed'
    assert '1927963p' in test_obs.planes.keys(), 'wrong plane deleted'
    assert '1927963f' in test_obs.planes.keys(), 'wrong plane deleted'
    assert '1927963og' not in test_obs.planes.keys(), 'wrong plane deleted'


def test_cleanup_augmentation_png(test_data_dir):
    test_obs = mc.read_obs_from_file(path.join(test_data_dir, 'visit/visit_obs_start_png_cleanup.xml'))
    assert len(test_obs.planes) == 5, 'initial conditions failed'
    storage_name = cfht_name.CFHTName(source_names=['2368534s.fits'])
    storage_name._instrument = Inst.SPIROU
    kwargs = {'storage_name': storage_name}
    test_obs = cleanup_augmentation.visit(test_obs, **kwargs)
    assert len(test_obs.planes) == 5, 'post-test conditions failed - should be no plane clean-up'
    artifact_keys = get_all_artifact_keys(test_obs)
    # remove the artifacts from the plane that's been visited, since otherwise there might not be preview images
    assert 'cadc:CFHT/2368534e_preview_1024.png' in artifact_keys, 's 1024'
    assert 'cadc:CFHT/2368534e_preview_256.png' in artifact_keys, 'e 256'
    assert 'cadc:CFHT/2368534s_preview_256.png' not in artifact_keys, 's 256'
    assert 'cadc:CFHT/2368534s_preview_1024.png' not in artifact_keys, 's 1024'


def test_cleanup_augmentation_bad_artifact_uris(test_data_dir):
    test_obs = mc.read_obs_from_file(path.join(test_data_dir, 'visit/visit_obs_duplicate_uris_cleanup.xml'))
    assert len(test_obs.planes) == 2, 'initial conditions failed'
    f1 = '2255229o.fits.fz'
    f2 = '2255229p.fits.fz'
    test_plane_1 = test_obs.planes[f1.split('.')[0]]
    test_plane_2 = test_obs.planes[f2.split('.')[0]]
    assert (
        len(test_plane_1.artifacts) == 7
    ), 'plane 1 initial conditions failed'
    assert (
        len(test_plane_2.artifacts) == 7
    ), 'plane 2 initial conditions failed'

    for f_name in [f1, f2]:
        storage_name = cfht_name.CFHTName(source_names=[f_name])
        kwargs = {'storage_name': storage_name}
        test_obs = cleanup_augmentation.visit(test_obs, **kwargs)
    assert len(test_obs.planes) == 2, 'post-test conditions failed'
    assert len(test_plane_1.artifacts) == 4, 'plane 1 post conditions failed'
    assert len(test_plane_2.artifacts) == 4, 'plane 2 post conditions failed'
    assert 'cadc:CFHT/2255229p_preview_256.jpg' in test_plane_2.artifacts.keys()
    assert 'cadc:CFHT/2255229o_preview_256.jpg' in test_plane_1.artifacts.keys()


def test_cleanup_edge_case(test_data_dir):
    test_obs = mc.read_obs_from_file(path.join(test_data_dir, 'visit/visit_obs_edge_case_cleanup.xml'))
    storage_name = cfht_name.CFHTName(source_names=['2014318p.fits.fz'], instrument=Inst.WIRCAM)
    kwargs = {'storage_name': storage_name}
    test_plane = test_obs.planes[storage_name.product_id]
    assert len(test_plane.artifacts) == 7, 'initial condition'
    test_obs = cleanup_augmentation.visit(test_obs, **kwargs)
    assert len(test_plane.artifacts) == 4, 'post condition'
