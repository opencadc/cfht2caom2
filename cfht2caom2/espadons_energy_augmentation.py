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

import logging

from caom2 import Observation, CoordAxis1D, SpectralWCS, Axis
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import main_app
from cfht2caom2 import metadata as md


def visit(observation, **kwargs):
    mc.check_param(observation, Observation)

    working_dir = kwargs.get('working_directory', './')
    science_file = kwargs.get('science_file')
    if science_file is None:
        raise mc.CadcException('Visitor needs a science_file parameter.')

    cfht_name = cn.CFHTName(file_name=science_file,
                            instrument=observation.instrument.name)
    count = 0
    if (cfht_name.instrument is md.Inst.ESPADONS and
            cfht_name.suffix in ['i', 'p']):
        for plane in observation.planes.values():
            for artifact in plane.artifacts.values():
                logging.error(cfht_name.file_uri)
                logging.error(artifact.uri)
                if cfht_name.file_uri == artifact.uri:
                    count += _do_energy(artifact, science_file, working_dir)
    logging.info('Completed ESPaDOnS energy augmentation for {}.'.format(
        observation.observation_id))
    return {'chunks': count}


def _do_energy(artifact, science_file, working_dir):
    # PD slack 08-01-20
    # espadons is a special case because using bounds allows one to
    # define "tiles" and then the SODA cutout service can extract the
    # subset of tiles that overlap the desired region. That's the best
    # that can be done because it is not possible to create a
    # CoordFunction1D to say what the wavelength of each pixel is
    #
    # If the coverage had significant gaps (eg SCUBA or SCUBA2 from
    # JCMT)  then the extra detail in bounds would enable better
    # discovery (as the gaps would be captured in the plane metadata).
    # In the case of espadons I don't think the gaps  are significant
    # (iirc, espadons is an eschelle spectrograph but I don't recall
    # whether the discontinuity between eschelle was a small gap or an
    # overlap)
    #
    # So: bounds provides more detail and it can in principle improve
    # data discovery (if gaps) and enable extraction of subsections of
    # the spectrum via the SODA service. Espadons was one of the use
    # cases that justified having bounds there

    # SF slack 08-01-20
    # We need the information that is contained in bounds. Gaps need
    # to be captured. So keep bounds. If you decide to remove range,
    # then advanced users would have to dig in the info to understand
    # range is first and last bounds.

    # read in the complete fits file, including the data
    fqn = f'{working_dir}/{science_file}'
    logging.info(f'Reading data from {fqn}.')
    hdus = ac.read_fits_data(fqn)
    wave = hdus[0].data[0, :]
    axis = Axis('WAVE', 'nm')
    coord_bounds = ac.build_chunk_energy_bounds(wave, axis)
    coord_axis = CoordAxis1D(axis=axis, bounds=coord_bounds)
    params = {'header': hdus[0].header,
              'uri': artifact.uri}
    resolving_power = main_app.get_espadons_energy_resolving_power(params)
    artifact.parts['0'].chunks[0].energy = SpectralWCS(
        coord_axis,
        specsys='TOPOCENT',
        ssysobs='TOPOCENT',
        ssyssrc='TOPOCENT',
        resolving_power=resolving_power)
    hdus.close()
    return 1
