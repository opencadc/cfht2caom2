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

from caom2pipe import manage_composable as mc
from cfht2caom2 import metadata as md


__all__ = ['CFHTName', 'COLLECTION', 'ARCHIVE']


COLLECTION = 'CFHT'
ARCHIVE = 'CFHT'


class CFHTName(mc.StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support fz and gz files in storage
    - product id == file id
    - the file_name attribute has ALL the extensions, including compression
      type.
    """

    CFHT_NAME_PATTERN = '*'

    def __init__(
        self,
        obs_id=None,
        file_name=None,
        instrument=None,
        ad_uri=None,
        entry=None,
        source_names=[],
    ):
        # set compression to an empty string so the file uri method still
        # works, since the file_name element will have all extensions,
        # including the .fz | .gz | '' to indicate compression type
        if obs_id is None:
            super(CFHTName, self).__init__(
                None,
                COLLECTION,
                CFHTName.CFHT_NAME_PATTERN,
                file_name,
                compression='',
                entry=entry,
                source_names=source_names,
            )
            self._instrument = md.Inst(instrument)
            if ad_uri is not None and file_name is None:
                file_name = mc.CaomName(ad_uri).file_name
            self._file_name = file_name
            self._file_id = CFHTName.remove_extensions(file_name)
            self._suffix = self._file_id[-1]
            if self._instrument in [md.Inst.MEGAPRIME, md.Inst.MEGACAM]:
                # SF - slack - 02-04-20
                # - MegaCam - the logic should be probably be 2 planes: p
                # and o for science. - all cfht exposures are sorted by EXPNUM
                # if i understand their data acquisition. b,f,d,x should be 1
                # plane observations. - my assumption is that the b,f,d,x have
                # no reason to have a processed equivalent.
                if (
                    self._suffix in ['b', 'd', 'f', 'x']
                    or self._suffix.isnumeric()
                ):
                    self._obs_id = self._file_id
                else:
                    self._obs_id = self._file_id[:-1]
            else:
                if self.is_simple and not self.is_master_cal:
                    self.obs_id = self._file_id[:-1]
                else:
                    self.obs_id = self._file_id
                    if self.is_derived_sitelle:
                        self.obs_id = self.obs_id.replace(self._suffix, 'p')
        else:
            super(CFHTName, self).__init__(
                obs_id,
                COLLECTION,
                CFHTName.CFHT_NAME_PATTERN,
                compression='',
            )
            self.obs_id = obs_id
            self._instrument = None
            self._file_id = None
            self._file_name = None
            self._file_id = None
            self._suffix = None

    def __str__(self):
        return (
            f'\n'
            f'  instrument: {self.instrument}\n'
            f'      obs_id: {self.obs_id}\n'
            f'     file_id: {self.file_id}\n'
            f'   file_name: {self.file_name}\n'
            f'source_names: {self.source_names}\n'
            f'    file_uri: {self.file_uri}\n'
            f'     lineage: {self.lineage}\n'
        )

    def is_valid(self):
        return True

    @property
    def file_id(self):
        """The file id - no extensions."""
        return self._file_id

    @property
    def file_name(self):
        """The file name."""
        return self._file_name

    @file_name.setter
    def file_name(self, value):
        """The file name."""
        self._file_name = value

    @property
    def instrument(self):
        return self._instrument

    @property
    def product_id(self):
        result = self._file_id
        return result

    @property
    def prev(self):
        """The preview file name for the file."""
        return '{}_preview_1024.jpg'.format(self.product_id)

    @property
    def thumb(self):
        """The thumbnail file name for the file."""
        return '{}_preview_256.jpg'.format(self.product_id)

    @property
    def zoom(self):
        """The zoom preview file name for the file."""
        return '{}_preview_zoom_1024.jpg'.format(self.product_id)

    @property
    def zoom_uri(self):
        """The zoom preview URI."""
        return self._get_uri(self.zoom)

    @property
    def is_master_cal(self):
        return (
            'weight' in self._file_id
            or 'master' in self._file_id
            or 'hotpix' in self._file_id
            or 'badpix' in self._file_id
            or 'deadpix' in self._file_id
            or 'dark' in self._file_id
            or 'scatter' in self._file_id
        )

    @property
    def has_energy(self):
        return not (
            self._instrument is md.Inst.ESPADONS and self._suffix in ['i', 'p']
        )

    @property
    def has_polarization(self):
        return self._suffix in ['p'] and self._instrument in [
            md.Inst.ESPADONS,
            md.Inst.SPIROU,
        ]

    @property
    def is_derived_sitelle(self):
        return self._instrument == md.Inst.SITELLE and self._suffix in [
            'p',
            'v',
            'z',
        ]

    @property
    def is_feasible(self):
        """
        Executing parts of the pipeline is not feasible for hdf5 files at this
        time - make it possible to know whether or not to try for each
        StorageName instance.
        :return:
        """
        return not mc.StorageName.is_hdf5(self._file_name)

    @property
    def is_simple(self):
        result = False
        if (
            self._suffix
            in [
                'a',
                'b',
                'c',
                'd',
                'f',
                'g',
                'l',
                'm',
                'o',
                's',
                'w',
                'x',
                'y',
            ]
            or self.simple_by_suffix
            or self.is_master_cal
            or (
                self._instrument is md.Inst.SPIROU
                and self._suffix in ['a', 'c', 'd', 'e', 'f', 'o', 'r', 'v']
            )
        ):
            result = True
        return result

    @property
    def simple_by_suffix(self):
        return (
            (
                self._suffix in ['p', 's']
                and self._instrument
                in [md.Inst.MEGACAM, md.Inst.MEGAPRIME, md.Inst.WIRCAM]
            )
            or (self._suffix == 'i' and self._instrument is md.Inst.ESPADONS)
            or (
                self._suffix in ['e', 's', 't', 'v']
                and self._instrument is md.Inst.SPIROU
            )
        )

    @property
    def suffix(self):
        return self._suffix

    @staticmethod
    def remove_extensions(name):
        """How to get the file_id from a file_name."""
        # ESPaDOnS files have a .gz extension ;)
        # SITELLE has hdf5 files
        return (
            name.replace('.fits', '')
            .replace('.fz', '')
            .replace('.header', '')
            .replace('.gz', '')
            .replace('.hdf5', '')
        )
