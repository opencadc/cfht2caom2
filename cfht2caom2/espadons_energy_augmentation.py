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

import copy
import logging

from astropy.io import fits
from caom2 import Axis, Chunk, CoordAxis1D, ObservableAxis, Observation, Slice, SpectralWCS
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from cfht2caom2 import metadata as md


def visit(observation, **kwargs):
    mc.check_param(observation, Observation)

    working_dir = kwargs.get('working_directory', './')
    storage_name = kwargs.get('storage_name')
    if storage_name is None:
        raise mc.CadcException('Visitor needs a storage_name parameter.')
    clients = kwargs.get('clients')
    if clients is None:
        logging.warning(f'No clients for ESPaDoNS energy augmentation.')
    science_fqn = storage_name.get_file_fqn(working_dir)
    count = 0
    if (
        storage_name.instrument is md.Inst.ESPADONS
        and storage_name.suffix
        in [
            'i',
            'p',
        ]
    ):
        for plane in observation.planes.values():
            for artifact in plane.artifacts.values():
                if storage_name.file_uri == artifact.uri:
                    count += _do_energy(artifact, science_fqn, storage_name)
        logging.info(
            f'Completed ESPaDOnS energy augmentation for {count} chunks in {observation.observation_id}.'
        )
    return observation


def _do_energy(artifact, science_fqn, storage_name):
    # CW l826
    # Set up observable axes, row 1 is wavelength, row 2 is normalized flux, row 3 ('p' only) is Stokes spectrum

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
    logging.info(f'Reading ESPaDOnS energy data from {science_fqn}.')
    count = 0
    hdus = fits.open(science_fqn, memmap=True, lazy_load_hdus=False)
    if hdus[0].data is not None:
        wave = hdus[0].data[0, :].copy()
        del hdus[0].data
    else:
        wave = hdus[1].data[0, :].copy()
        del hdus[1].data
    hdr = hdus[0].header.copy()
    hdus.close()
    del hdus
    axis = Axis('WAVE', 'nm')
    coord_bounds = ac.build_chunk_energy_bounds(wave)
    coord_axis = CoordAxis1D(axis=axis, bounds=coord_bounds)
    resolving_power = get_energy_resolving_power(hdr)
    chunk = artifact.parts['0'].chunks[0]
    chunk.energy = SpectralWCS(
        coord_axis,
        specsys='TOPOCENT',
        ssyssrc='TOPOCENT',
        resolving_power=resolving_power,
    )
    chunk.energy_axis = 1
    chunk.naxis = hdr.get('NAXIS')
    logging.info(f'chunk.naxis is {chunk.naxis}')
    if chunk.naxis is not None and chunk.naxis == 2:
        independent_axis = Axis('WAVE', 'nm')
        independent = Slice(independent_axis, 1)
        dependent_axis = Axis('flux', 'counts')
        dependent = Slice(dependent_axis, 2)
        chunk.observable = ObservableAxis(dependent, independent)
        chunk.observable_axis = 2
        count += 1
        logging.info('Chunk 0 updated')

    if storage_name.suffix == 'p':
        if len(artifact.parts['0'].chunks) == 2:
            # replace the existing value
            artifact.parts['0'].chunks.pop(1)

        # caom2IngestEspadons.py, l863
        dependent_axis = Axis('polarized flux', 'percent')
        dependent = Slice(dependent_axis, 3)
        independent_axis = Axis('WAVE', 'nm')
        independent = Slice(independent_axis, 1)
        new_chunk = copy.deepcopy(chunk)
        new_chunk.observable = ObservableAxis(dependent, independent)
        new_chunk.energy = chunk.energy
        new_chunk._id = Chunk._gen_id()
        artifact.parts['0'].chunks.append(new_chunk)
        count += 1
        logging.info('Chunk 1 added')

    logging.info(f'Done reading energy, changed {count} chunks.')
    return count


def get_energy_resolving_power(header):
    result = None
    obstype = header.get('OBSTYPE')
    if obstype and obstype not in ['BIAS', 'DARK']:
        instmode = header.get('INSTMODE')
        if instmode is None or 'R=' not in instmode:
            # CW - Default if resolving power value not in header caom2IngestEspadons.py, l377
            result = 65000.0
        else:
            # CW - This string is already in instrument keywords but also need to extract resolving power from it:
            # 'Spectroscopy, star only, R=80,000'
            temp = instmode.split('R=')
            values = temp[1].split(',')
            if len(values) == 1:
                result = values[0]
            else:
                result = f'{values[0]}{values[1]}'
            result = mc.to_float(result)
    return result
