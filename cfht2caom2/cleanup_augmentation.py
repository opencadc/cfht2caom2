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

import logging

from collections import defaultdict
from caom2 import Observation
from caom2pipe import manage_composable as mc
from cfht2caom2 import metadata as md


def visit(observation, **kwargs):
    mc.check_param(observation, Observation)

    storage_name = kwargs.get('storage_name')

    count = 0
    artifact_count = 0
    delete_list = []
    artifact_delete_list = defaultdict(list)

    if observation.instrument is None:
        # a SITELLE observation with only an HDF5 artifact, no cleanup required
        return observation

    for plane in observation.planes.values():
        if plane.product_id.endswith('og'):
            delete_list.append(plane.product_id)

        if storage_name.product_id != plane.product_id:
            continue

        # prior to cfht2caom2, artifact URIs were named after the observation
        # ID, with cfht2caom2 they're named after the product ID - clean up the
        # obsolete ones
        if (
            (
                observation.instrument.name is md.Inst.ESPADONS
                and storage_name.suffix == 'i'
                and len(plane.artifacts) > 3
            )
            or (
                observation.instrument.name is md.Inst.SPIROU
                and storage_name.suffix in ['e', 'p', 's', 't', 'v']
                and len(plane.artifacts) > 3
            )
            or (len(plane.artifacts) > 4)
        ):
            for artifact in plane.artifacts.values():
                if '.jpg' not in artifact.uri:
                    continue
                if f':CFHT/{storage_name.product_id}_' not in artifact.uri:
                    artifact_delete_list[plane.product_id].append(
                        artifact.uri
                    )

        for artifact in plane.artifacts.values():
            if '.png' in artifact.uri:
                # there are duplicate .png and .jpg previews, prefer the smaller jpgs
                artifact_delete_list[plane.product_id].append(artifact.uri)

    for entry in delete_list:
        logging.info(
            f'Removing plane {entry} from {observation.observation_id}'
        )
        count += 1
        observation.planes.pop(entry)

    for product_id, entries in artifact_delete_list.items():
        for plane in observation.planes.values():
            if plane.product_id == product_id:
                for entry in entries:
                    artifact_count += 1
                    logging.info(f'Removing artifact {entry} from plane {product_id} in {observation.observation_id}')
                    plane.artifacts.pop(entry)

    logging.info(
        f'Completed cleanup augmentation for {observation.observation_id}. '
        f'Remove {count} planes and {artifact_count} artifacts from the '
        f'observation.'
    )
    return observation
