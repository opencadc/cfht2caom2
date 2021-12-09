# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2021.                            (c) 2021.
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

import logging

from importlib import import_module

from caom2 import SimpleObservation, DerivedObservation, Algorithm
from caom2utils import ObsBlueprint, GenericParser, FitsParser
from cfht2caom2 import instruments, main_app


class Fits2caom2Visitor:
    def __init__(self, observation, **kwargs):
        self._observation = observation
        self._storage_name = kwargs.get('storage_name')
        self._file_metadata = kwargs.get('file_metadata')
        self._cadc_client = kwargs.get('cadc_client')
        self._dump_config = False
        self._logger = logging.getLogger(self.__class__.__name__)

    def visit(self):
        for uri, file_metadata in self._file_metadata.items():
            instrument_data = instruments.instrument_blueprint_factory(
                file_metadata.headers, self._storage_name
            )
            blueprint = ObsBlueprint(instantiated_class=instrument_data)
            main_app.accumulate_bp(blueprint, self._storage_name)

            if len(file_metadata.headers) == 0:
                parser = GenericParser(blueprint, uri)
            else:
                parser = FitsParser(file_metadata.headers, blueprint, uri)

            if self._dump_config:
                print(f'Blueprint for {uri}: {blueprint}')

            if self._observation is None:
                if blueprint._get('DerivedObservation.members') is None:
                    self._logger.debug('Build a SimpleObservation')
                    self._observation = SimpleObservation(
                        collection=self._storage_name.collection,
                        observation_id=self._storage_name.obs_id,
                        algorithm=Algorithm('exposure'),
                    )
                else:
                    self._logger.debug('Build a DerivedObservation')
                    self._observation = DerivedObservation(
                        collection=self._storage_name.collection,
                        observation_id=self._storage_name.obs_id,
                        algorithm=Algorithm('composite'),
                    )

            parser.augment_observation(
                observation=self._observation,
                artifact_uri=uri,
                product_id=self._storage_name.product_id,
            )

            self._observation = instrument_data.update(
                self._observation,
                file_metadata.file_info,
                self._cadc_client,
            )
        return self._observation


def visit(observation, **kwargs):
    s = Fits2caom2Visitor(observation, **kwargs)
    return s.visit()
