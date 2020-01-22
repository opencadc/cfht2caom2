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

from enum import Enum

from caom2pipe import astro_composable as ac


class Inst(Enum):
    ESPADONS = 'ESPaDOnS'
    MEGACAM = 'MegaCam'
    MEGAPRIME = 'MetaPrime'
    SITELLE = 'SITELLE'
    WIRCAM = 'WIRCam'
    NONE = None


# key is CFHT metadata
# value is SVO url piece
FILTER_REPAIR_LOOKUP = {
    # MegaPrime
    'u.MP9302': 'u_sdss',
    'u.MP9301': 'u',
    'CaHK.MP9303': 'CaHK',
    'g.MP9401': 'g',
    'g.MP9402': 'g_sdss',
    'OIII.MP9501': 'OIII',
    'OIII.MP9502': 'OIII_off',
    'gri.MP9605': 'gri',
    'Halpha.MP7604': 'Halpha_off',
    'Halpha.MP7605': 'Halpha_on',
    'Halpha.MP9603': 'Halpha',
    'Halpha.MP9604': 'Halpha_Off_2',
    'r.MP9601': 'r',
    'r.MP9602': 'r_sdss',
    'i.MP9701': 'i1',
    'i.MP9702': 'i',
    'i.MP9703': 'i_sdss',
    'TiO.MP7701': 'TiO',
    'CN.MP7803': 'CN',
    'NB920.MP7902': 'NB920',
    'z.MP9801': 'z',
    'z.MP9901': 'z_sdss',
    # WIRCam
    'BrG': 'Brackett_gamma',
    'Ks.WC8302': 'Ks'}

# key is CFHT metadata
# value is SVO url piece
INSTRUMENT_REPAIR_LOOKUP = {Inst.WIRCAM: 'Wircam'}

# CW/SF 17-12-19 - content from conversation
# CW
# If no or "open" filter or pinhole mask, put in 20% on-off cuts from
# MegaCam_QE_data.txt
#
CFHT_CACHE = {
    # pinhole mask, so use full energy range
    'PHG.MP9999': {'cw': 6200., 'fwhm': 6000.},
    # test filter - just use full energy range
    'N393.MP1111': {'cw': 6200., 'fwhm': 6000.},
    # CFH12K filter
    'OIII:MP7504': {'cw': 0.0, 'fwhm': 0.0},
    # 'NONE': {'cw': 6200., 'fwhm': 6000.},
    # 'OPEN': {'cw': 6200., 'fwhm': 6000.},
    'MegaPrime.None': {'cw': 6200., 'fwhm': 6000.},
    # SITELLE - not available at SVO
    'C1': {'cw': None, 'fwhm': None},
    'C2': {'cw': None, 'fwhm': None},
    'C3': {'cw': None, 'fwhm': None},
    'C4': {'cw': None, 'fwhm': None},
    'SN1': {'cw': None, 'fwhm': None},
    'SN2': {'cw': None, 'fwhm': None},
    'SN3': {'cw': None, 'fwhm': None},
    # WIRCam
    'WIRCam.None': {'cw': 17000.0, 'fwhm': 14000.0},
    'FakeBlank': {'cw': 17000.0, 'fwhm': 14000.0}
}

filter_cache = ac.FilterMetadataCache(FILTER_REPAIR_LOOKUP,
                                      INSTRUMENT_REPAIR_LOOKUP,
                                      'CFHT', CFHT_CACHE, 'NONE')
