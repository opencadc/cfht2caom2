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

import logging
import os

from cadcdata import CadcDataClient
from caom2repo import CAOM2RepoClient
from caom2utils import fits2caom2
from caom2pipe import name_builder_composable as nbc
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import metadata as md


__all__ = ['CFHTBuilder']


class CFHTBuilder(nbc.StorageNameBuilder):

    def __init__(self, config):
        super(CFHTBuilder, self).__init__()
        self._config = config
        self._data_client = None
        self._repo_client = None
        self._metrics = mc.Metrics(self._config)
        if not self._config.use_local_files:
            subject = mc.define_subject(self._config)
            self._data_client = CadcDataClient(subject)
            self._repo_client = CAOM2RepoClient(
                subject, resource_id=self._config.resource_id)
        self._logger = logging.getLogger(__name__)

    def build(self, entry):
        """
        :param entry an entry is a file name, complete with the appropriate
        compression extension, that is sufficient to retrieve file header
        information from CADC's storage system.
        """

        # retrieve the header information, extract the instrument name
        self._logger.debug(f'Build a StorageName instance for {entry}.')
        if (mc.TaskType.INGEST_OBS in self._config.task_types and
                ('.fits' not in entry and not mc.StorageName.is_hdf5(entry))):
            obs = mc.repo_get(self._repo_client, self._config.collection,
                              entry, self._metrics)
            instrument = md.Inst(obs.instrument.name)
            result = cn.CFHTName(obs_id=entry, instrument=instrument)
        else:
            if mc.StorageName.is_hdf5(entry):
                headers = []
            else:
                if self._config.use_local_files:
                    cwd = os.getcwd()
                    headers = fits2caom2.get_cadc_headers(
                        f'file://{cwd}/{entry}')
                else:
                    headers_str = mc.get_cadc_headers_client(
                        self._config.archive, entry, self._data_client)
                    headers = ac.make_headers_from_string(headers_str)

            instrument = CFHTBuilder.get_instrument(headers, entry)
            result = cn.CFHTName(file_name=entry, instrument=instrument,
                                 entry=entry)
        return result

    @staticmethod
    def get_instrument(headers, entry):
        """
        SF - 15-04-20 - slack - what if there's no INSTRUME or DETECTOR
        keyword?
        something like: if CFHT: if has(INSTRUME) then
        MEGA else if NEXTEND>30 then MEGA else...

        On MegaPrime vs MegaCam values from SVO:
        megaprime == full instrument (optics+camera) ,
        megacam == camera

        SVO has mistakes (they say megacam=CCD+mirror+optics), but i think the
        Megaprime entry has a lot more filters. So to be consistent, we should
        always use the Megaprime entry

        SGo - make the 'always use Megaprime' happen here by setting the
        instrument to MegaPrime.

        :param headers: astropy fits headers
        :param entry: string for error logging
        :return: md.Inst instance
        """
        if mc.StorageName.is_hdf5(entry):
            inst = md.Inst.SITELLE
        else:
            nextend = None
            detector = None
            instrument = headers[0].get('INSTRUME')
            if instrument is None:
                instrument = headers[1].get('INSTRUME')
                if instrument is None:
                    instrument = headers[0].get('DETECTOR')
                    if instrument is None:
                        instrument = headers[1].get('DETECTOR')
                        if instrument is None:
                            nextend = headers[0].get('NEXTEND')
                            if nextend is None:
                                raise mc.CadcException(f'Could not identify '
                                                       f'instrument for '
                                                       f'{entry}.')
            elif instrument == 'Unknown':
                detector = headers[0].get('DETECTOR')
            if instrument is None and nextend is not None and nextend > 30:
                inst = md.Inst.MEGAPRIME
            else:
                msg = f'Unknown value for instrument {instrument}, detector ' \
                      f'{detector} and nextend {nextend} for {entry}.'

                try:
                    inst = md.Inst(instrument)
                except ValueError:
                    if (instrument == 'CFHT MegaPrime' or
                            instrument == 'megacam'):
                        inst = md.Inst.MEGAPRIME
                    elif instrument == 'Unknown' and detector is not None:
                        if detector == 'OLAPA':
                            # SF 24-06-20
                            # ok to hack ESPADONS name for OLAPA
                            inst = md.Inst.ESPADONS
                        else:
                            try:
                                inst = md.Inst(detector)
                            except ValueError:
                                # SF 24-06-20
                                # nasty hacks: all the 48 failed espadons
                                # files have a PATHNAME with espadon in it
                                #
                                # so you may add PATHNAME check if everything
                                # else fails
                                pathname = headers[0].get('PATHNAME')
                                if 'espadons' in pathname:
                                    inst = md.Inst.ESPADONS
                                else:
                                    raise mc.CadcException(msg)
                    else:
                        raise mc.CadcException(msg)
        return inst
