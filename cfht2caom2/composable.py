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

import logging
import sys
import traceback

from caom2pipe import client_composable as clc
from caom2pipe.data_source_composable import LocalFilesDataSourceRunnerMeta
from caom2pipe.manage_composable import Config, StorageName, TaskType
from caom2pipe import run_composable as rc
from cfht2caom2 import cleanup_augmentation
from cfht2caom2 import espadons_energy_augmentation, preview_augmentation
from cfht2caom2 import file2caom2_augmentation
from cfht2caom2.cfht_name import CFHTName


META_VISITORS = [file2caom2_augmentation]
DATA_VISITORS = [
    espadons_energy_augmentation,
    preview_augmentation,
    cleanup_augmentation,
]

def can_use_single_visit(task_types):
    return (
        len(task_types) > 1
        and (
            (TaskType.STORE in task_types and TaskType.INGEST in task_types)
            or (TaskType.INGEST in task_types and TaskType.MODIFY in task_types)
            or (TaskType.SCRAPE in task_types and TaskType.MODIFY in task_types)
        )
    )


def _common_init():
    config = Config()
    config.get_executors()
    StorageName.collection = config.collection
    StorageName.scheme = config.scheme
    StorageName.preview_scheme = config.preview_scheme
    StorageName.data_source_extensions = config.data_source_extensions
    clients = clc.ClientCollection(config)
    sources = []
    if config.use_local_files:
        source = LocalFilesDataSourceRunnerMeta(config, clients.data_client, storage_name_ctor=CFHTName)
        sources.append(source)
    return config, clients, sources


def _run_state():
    config, clients, sources = _common_init()
    return rc.run_by_state_runner_meta(
        config=config,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        sources=sources,
        clients=clients,
        organizer_module_name='cfht2caom2.cfht_name',
        organizer_class_name='CFHTOrganizeExecutesRunnerMeta',
        storage_name_ctor=CFHTName,
    )


def run_state():
    """Wraps _run_state in exception handling."""
    try:
        _run_state()
        sys.exit(0)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)


def _run():
    """Run the processing for observations using a todo file to identify the work to be done. StorageName
    construction is incomplete with a todo file, because the instrument name and BITPIX are required.

    :return 0 if successful, -1 if there's any sort of failure. Return status
        is used by airflow for task instance management and reporting.
    """
    config, clients, sources = _common_init()
    return rc.run_by_todo_runner_meta(
        config,
        sources=sources,
        meta_visitors=META_VISITORS,
        data_visitors=DATA_VISITORS,
        clients=clients,
        organizer_module_name='cfht2caom2.cfht_name',
        organizer_class_name='CFHTOrganizeExecutesRunnerMeta',
        storage_name_ctor=CFHTName,
    )


def run():
    try:
        result = _run()
        sys.exit(result)
    except Exception as e:
        logging.error(e)
        tb = traceback.format_exc()
        logging.debug(tb)
        sys.exit(-1)
