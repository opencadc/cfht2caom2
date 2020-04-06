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

from datetime import datetime

from caom2 import Observation, ProductType, ReleaseType, ObservationIntentType
from caom2pipe import astro_composable as ac
from caom2pipe import caom_composable as cc
from caom2pipe import manage_composable as mc
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import metadata as md

__all__ = ['visit']


MIME_TYPE = 'image/jpeg'


def visit(observation, **kwargs):
    mc.check_param(observation, Observation)

    working_dir = kwargs.get('working_directory', './')
    science_file = kwargs.get('science_file')
    if science_file is None:
        raise mc.CadcException('Visitor needs a science_file parameter.')
    cadc_client = kwargs.get('cadc_client')
    if cadc_client is None:
        logging.warning('Visitor needs a cadc_client parameter to store '
                        'previews and thumbnails.')
    stream = kwargs.get('stream')
    if stream is None:
        raise mc.CadcException('Visitor needs a stream parameter.')
    observable = kwargs.get('observable')
    if observable is None:
        raise mc.CadcException('Visitor needs a observable parameter.')

    count = 0
    for plane in observation.planes.values():
        if (plane.data_release is None or
                plane.data_release > datetime.utcnow()):
            logging.info('Plane {} is proprietary. No preview access or '
                         'thumbnail creation.'.format(plane.product_id))
            continue
        for artifact in plane.artifacts.values():
            cfht_name = cn.CFHTName(file_name=science_file,
                                    instrument=observation.instrument.name)
            if cfht_name.suffix == 'g':
                # 19-03-20 - seb -
                # for all 'g' files, use the previews from the 'p' files
                continue

            if cfht_name.file_uri == artifact.uri:
                count += _do_prev(working_dir, plane, cfht_name, cadc_client,
                                  stream, observable.metrics,
                                  observation.intent,
                                  observation.instrument.name,
                                  observation.observation_id)
                if cfht_name.suffix == 'p':
                    count += _update_g_artifact(observation, plane)
                break

    logging.info('Completed preview augmentation for {}.'.format(
        observation.observation_id))
    return {'artifacts': count}


def _do_prev(working_dir, plane, cfht_name, cadc_client,
             stream, metrics, intent, instrument, obs_id):
    science_fqn = os.path.join(working_dir, cfht_name.file_name)
    preview = cfht_name.prev
    preview_fqn = os.path.join(working_dir, preview)
    thumb = cfht_name.thumb
    thumb_fqn = os.path.join(working_dir, thumb)
    zoom = cfht_name.zoom
    zoom_fqn = os.path.join(working_dir, zoom)

    # from genWirprevperplane.py
    # if it's a datacube, just take the first slice
    # e.g. fitscopy '928690p.fits[*][*,*,1:1]' s1928690p.fits

    headers = ac.read_fits_headers(science_fqn)
    num_extensions = headers[0].get('NEXTEND')

    zoom_science_fqn = science_fqn

    # set up the correct input file - may need to use fitscopy
    rotate_param = ''
    if md.Inst(instrument) is md.Inst.WIRCAM:
        if science_fqn.endswith('.fz'):
            naxis_3 = headers[0].get('ZNAXIS3', 1)
        else:
            naxis_3 = headers[0].get('NAXIS3', 1)

        if naxis_3 != 1:
            logging.info(f'Observation {obs_id}: using first slice of '
                         f'{science_fqn}.')
            temp_science_f_name = cfht_name.file_name.replace('.fz', '_slice.fz')
            slice_cmd = f'fitscopy {cfht_name.file_name}[*][*,*,1:1,1:1] ' \
                        f'{temp_science_f_name}'
            _exec_cmd_chdir(working_dir, temp_science_f_name, slice_cmd)
            science_fqn = f'{working_dir}/{temp_science_f_name}'

        if num_extensions >= 4:
            logging.info(f'Observation {obs_id}: using slice for zoom preview '
                         f'of {science_fqn}.')
            zoom_science_f_name = cfht_name.file_name.replace(
                '.fits', '_zoom.fits')
            slice_cmd = f'fitscopy {cfht_name.file_name}[4][*,*,1:1] ' \
                        f'{zoom_science_f_name}'
            _exec_cmd_chdir(working_dir, zoom_science_f_name, slice_cmd)
    elif md.Inst(instrument) in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
        rotate_param = '-rotate 180'

    # set up the correct parameters to the ds9 command
    scope_param = 'local'
    if intent is ObservationIntentType.SCIENCE:
        scope_param = 'global'

    geometry = '256x521'
    _gen_image(science_fqn, geometry, thumb_fqn, scope_param, rotate_param)

    geometry = '1024x1024'
    if md.Inst(instrument) in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
        _gen_image(science_fqn, geometry, preview_fqn, scope_param,
                   rotate_param, mode_param='')
    else:
        _gen_image(science_fqn, geometry, preview_fqn, scope_param,
                   rotate_param)

    mosaic_param = '-fits'
    zoom_param = '1'
    scope_param = 'global'
    # set zoom parameters
    if md.Inst(instrument) is md.Inst.WIRCAM:
        pan_param = '-pan 484 -484 image'
    elif md.Inst(instrument) in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
        pan_param = '-pan -9 1780'
        rotate_param = '-rotate 180'
        if num_extensions >= 23:
            rotate_param = ''
            mosaic_param = f'-fits {zoom_science_fqn}[23]'
            zoom_science_fqn = ''
        elif num_extensions >= 14:
            mosaic_param = f'-fits {zoom_science_fqn}[14]'
            zoom_science_fqn = ''
        else:
            mosaic_param = f'-fits {zoom_science_fqn}[1]'
            zoom_science_fqn = ''
    _gen_image(zoom_science_fqn, geometry, zoom_fqn, scope_param, rotate_param,
               zoom_param, pan_param, mosaic_param)

    prev_uri = cfht_name.prev_uri
    thumb_uri = cfht_name.thumb_uri
    zoom_uri = cfht_name.zoom_uri
    _augment(plane, prev_uri, preview_fqn, ProductType.PREVIEW)
    _augment(plane, thumb_uri, thumb_fqn, ProductType.THUMBNAIL)
    _augment(plane, zoom_uri, zoom_fqn, ProductType.PREVIEW)
    if cadc_client is not None:
        _store_smalls(cadc_client, working_dir, stream, preview, thumb, zoom,
                      metrics)
    return 3


def _exec_cmd_chdir(working_dir, temp_file, cmd):
    orig_dir = os.getcwd()
    try:
        os.chdir(working_dir)
        if os.path.exists(temp_file):
            os.unlink(temp_file)
        mc.exec_cmd(cmd)
    finally:
        os.chdir(orig_dir)


def _gen_image(science_fqn, geometry, save_fqn, scope_param, rotate_param,
               zoom_param='to fit',
               pan_param='', mosaic_param='-mosaicimage iraf',
               mode_param='-mode none'):
    # 20-03-20 - seb - always use iraf - do not trust wcs coming from the data
    # acquisition. A proper one needs processing which is often not done on
    # observations.
    cmd = f'xvfb-run -a ds9 {mosaic_param} {science_fqn} ' \
          f'{pan_param} ' \
          f'-geometry {geometry} ' \
          f'{rotate_param} ' \
          f'-scale squared ' \
          f'-scale mode zscale ' \
          f'-scale scope {scope_param} ' \
          f'-scale datasec yes ' \
          f'-invert ' \
          f'{mode_param} ' \
          f'-view colorbar no ' \
          f'-zoom {zoom_param} ' \
          f'-saveimage jpeg {save_fqn} ' \
          f'-quit'
    mc.exec_cmd(cmd, timeout=900)  # wait 15 minutes till killing


def _store_smalls(cadc_client, working_directory, stream, preview_fname,
                  thumb_fname, zoom_fname, metrics):
    mc.data_put(cadc_client, working_directory, preview_fname, cn.COLLECTION,
                stream, mime_type='image/jpeg', metrics=metrics)
    mc.data_put(cadc_client, working_directory, thumb_fname, cn.COLLECTION,
                stream, mime_type='image/jpeg', metrics=metrics)
    mc.data_put(cadc_client, working_directory, zoom_fname, cn.COLLECTION,
                stream, mime_type='image/jpeg', metrics=metrics)


def _update_g_artifact(observation, plane):
    count = 0
    g_product_id = plane.product_id.replace('p', 'g')
    if g_product_id in observation.planes.keys():
        g_plane = observation.planes[g_product_id]
        for artifact in plane.artifacts.values():
            if artifact.product_type in [ProductType.THUMBNAIL,
                                         ProductType.PREVIEW]:
                features = mc.Features()
                features.supports_latest_caom = True
                artifact_copy = cc.copy_artifact(artifact, features)
                g_plane.artifacts.add(artifact_copy)
                count += 1
    return count


def _augment(plane, uri, fqn, product_type):
    temp = None
    if uri in plane.artifacts:
        temp = plane.artifacts[uri]
    plane.artifacts[uri] = mc.get_artifact_metadata(
        fqn, product_type, ReleaseType.DATA, uri, temp)
