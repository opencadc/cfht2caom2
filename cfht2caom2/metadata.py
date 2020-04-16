# -*- coding: utf-8 -*-
# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2019.                            (c) 2019.
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
import re

from bs4 import BeautifulSoup
from enum import Enum

from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc


class Inst(Enum):
    ESPADONS = 'ESPaDOnS'
    MEGACAM = 'MegaCam'
    MEGAPRIME = 'MegaPrime'
    SITELLE = 'SITELLE'
    WIRCAM = 'WIRCam'
    NONE = None


# key is CFHT metadata
# value is SVO url piece
INSTRUMENT_REPAIR_LOOKUP = {'WIRCam': 'Wircam'}

# CW/SF 17-12-19 - content from conversation
# CW
# If no or "open" filter or pinhole mask, put in 20% on-off cuts from
# MegaCam_QE_data.txt
#
# SF 28-01-20 - from conversation - if there's no energy information,
# e.g. 1007126g.fits.gz with filter name == 'FakeBlank', leave the
# energy information empty

# keys in the cache:
FILTER_REPAIR_CACHE = 'filter_repair_lookup'
ENERGY_DEFAULTS_CACHE = 'energy_defaults'
PROJECT_TITLES_CACHE = 'project_titles'


class CFHTCache(mc.Cache):
    def __init__(self):
        super(CFHTCache, self).__init__()
        self._project_titles = self.get_from(PROJECT_TITLES_CACHE)
        self._cached_semesters = self._fill_cached_semesters()
        self._logger = logging.getLogger(__name__)

    def _fill_cached_semesters(self):
        result = []
        for key in self._project_titles.keys():
            # format of a run id is 20AS19, or 09BC99:
            result.append(CFHTCache.semester(key))
        return list(set(result))

    def _semester_cached(self, run_id):
        return CFHTCache.semester(run_id) in self._cached_semesters

    def _try_to_append_to_cache(self, run_id):
        if run_id is None or run_id == '':
            # the case of no value for run_id
            return

        sem = CFHTCache.semester(run_id)
        sem_int = mc.to_int(sem[:-1])
        updated_content = False
        if len(sem) < 3 or not sem[0] in ['0', '1', '2'] or sem_int < 9:
            # the URL here only works from 2009B on
            return
        base_url = 'http://www.cfht.hawaii.edu/en/science/QSO/'
        semester_url = f'{base_url}20{sem}/'
        self._logger.info(
            f'Checking for semester information at {semester_url}')
        response = mc.query_endpoint(semester_url)
        soup = BeautifulSoup(response.text, features='lxml')
        response.close()

        html_table = soup.find('table')
        rows = html_table.find_all('a', string=re.compile('\\.html'))
        for row in rows:
            inst_url = row.get('href')
            if not isinstance(inst_url, str) and 'qso_prog_' not in inst_url:
                continue
            if not inst_url.startswith('http'):
                inst_url = f'{semester_url}/{inst_url}'

            self._logger.info(
                f'Querying {inst_url} for new project information.')
            inst_response = mc.query_endpoint(inst_url)
            inst_soup = BeautifulSoup(inst_response.text, features='lxml')
            inst_response.close()
            table_rows = inst_soup.find_all('tr')
            for table_row in table_rows:
                tds = table_row.find_all('td')
                count = 0
                for td in tds:
                    if count == 0:
                        program_id = td.text
                    if count == 5:
                        title = td.text
                        updated_content = True
                        self._project_titles[program_id] = title
                        break
                    count += 1
        if updated_content:
            self.save()
        self._cached_semesters.append(CFHTCache.semester(run_id))

    def get_title(self, run_id):
        result = self._project_titles.get(run_id)
        if result is None:
            self._logger.warning(
                f'Could not find project information for run id {run_id}.')
            if (not self._semester_cached(run_id) and
                    run_id not in ['SMEARING', '00']):
                self._try_to_append_to_cache(run_id)
                # in case the cache was updated
                result = self._project_titles.get(run_id)
        return result

    @staticmethod
    def semester(run_id):
        return run_id[:3]


cache = CFHTCache()

filter_cache = ac.FilterMetadataCache(cache.get_from(FILTER_REPAIR_CACHE),
                                      INSTRUMENT_REPAIR_LOOKUP,
                                      'CFHT',
                                      cache.get_from(ENERGY_DEFAULTS_CACHE),
                                      'NONE')
