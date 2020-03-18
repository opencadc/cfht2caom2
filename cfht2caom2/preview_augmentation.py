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

from caom2 import Observation, ProductType, ReleaseType
from caom2pipe import manage_composable as mc
from cfht2caom2 import cfht_name as cn

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
        raise mc.CadcException('Visitor needs a cadc_client parameter.')
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
            if cfht_name.file_uri == artifact.uri:
                count += _do_prev(working_dir, plane, cfht_name, cadc_client,
                                  stream, observable.metrics)
                break
    logging.info('Completed preview augmentation for {}.'.format(
        observation.observation_id))
    return {'artifacts': count}


# def _check_for_delete(file_name, uri, observable, plane):
#     """If the preview file doesn't exist, but the artifact that represents it
#     does, remove that artifact from the Observation instance."""
#     result = 0
#     if (observable.rejected.is_no_preview(
#             file_name) and uri in plane.artifacts.keys()):
#         logging.warning(
#             'Removing artifact for non-existent preview {}'.format(uri))
#         plane.artifacts.pop(uri)
#         result = 1
#     return result


# def _do_prev(file_id, science_fqn, working_dir, plane, cadc_client, metrics):
def _do_prev(working_dir, plane, cfht_name, cadc_client,
             stream, metrics):
    science_fqn = os.path.join(working_dir, cfht_name.file_name)
    preview = cfht_name.prev
    preview_fqn = os.path.join(working_dir, preview)
    thumb = cfht_name.thumb
    thumb_fqn = os.path.join(working_dir, thumb)
    zoom = cfht_name.zoom
    zoom_fqn = os.path.join(working_dir, zoom)

    if os.access(preview_fqn, 0):
        os.remove(preview_fqn)
    # prev_cmd = 'fitscut --all --autoscale=99.5 --asinh-scale --jpg --invert ' \
    #            '--compass {}'.format(science_fqn)
    # mc.exec_cmd_redirect(prev_cmd, preview_fqn)

    # cmd3="ds9 -mosaicimage wcs  %s
    #           -geometry 256x521
    #           -scale squared
    #           -scale mode zscale
    #           -scale scope global
    #           -scale datasec yes
    #           -invert
    #           -mode no
    #           -view colorbar no
    #           -zoom to fit
    #           -saveimage jpeg  %s -quit" % (fitsfile, prev1file)

    # cmd4="ds9 -mosaicimage wcs  %s
    #           -geometry 1024x1289
    #           -scale squared
    #           -scale mode zscale
    #           -scale scope global
    #           -scale datasec yes
    #           -invert
    #           -view colorbar no
    #           -zoom to fit
    #           -saveimage jpeg  %s -quit" % (fitsfile, prev2file)
    mosaic_param = 'iraf'
    geometry = '1024x1289'
    prev_cmd = f'xvfb-run ds9 -mosaicimage {mosaic_param} {science_fqn} ' \
               f'-geometry {geometry} ' \
               f'-scale squared ' \
               f'-scale mode zscale ' \
               f'-scale scope global ' \
               f'-scale datasec yes ' \
               f'-invert ' \
               f'-mode no ' \
               f'-view colorbar no ' \
               f'-zoom to fit ' \
               f'-saveimage jpeg {preview_fqn} ' \
               f'-quit'
    mc.exec_cmd(prev_cmd)

    if os.access(thumb_fqn, 0):
        os.remove(thumb_fqn)
    geometry = '256x521'
    thumb_cmd = f'xvfb-run ds9 -mosaicimage {mosaic_param} {science_fqn} ' \
                f'-geometry {geometry} ' \
                f'-scale squared ' \
                f'-scale mode zscale ' \
                f'-scale scope global ' \
                f'-scale datasec yes ' \
                f'-invert ' \
                f'-mode none ' \
                f'-view colorbar no ' \
                f'-zoom to fit ' \
                f'-saveimage jpeg {thumb_fqn} ' \
                f'-quit'
    mc.exec_cmd(thumb_cmd)

    # if os.access(zoom_fqn, 0):
    #     os.remove(zoom_fqn)
    # geometry = '1024x1289'
    # zoom_cmd = f'xvfb-run ds9 -mosaicimage {mosaic_param} {science_fqn} ' \
    #            f'-geometry {geometry} ' \
    #            f'-scale squared ' \
    #            f'-scale mode zscale ' \
    #            f'-scale scope global ' \
    #            f'-scale datasec yes ' \
    #            f'-invert ' \
    #            f'-mode none ' \
    #            f'-view colorbar no ' \
    #            f'-zoom to fit ' \
    #            f'-saveimage jpeg {zoom_fqn} ' \
    #            f'-quit'
    # mc.exec_cmd(zoom_cmd)

    prev_uri = cfht_name.prev_uri
    thumb_uri = cfht_name.thumb_uri
    zoom_uri = cfht_name.zoom_uri
    _augment(plane, prev_uri, preview_fqn, ProductType.PREVIEW)
    _augment(plane, thumb_uri, thumb_fqn, ProductType.THUMBNAIL)
    # _augment(plane, zoom_uri, zoom_fqn, ProductType.PREVIEW)
    if cadc_client is not None:
        _store_smalls(cadc_client, working_dir, stream, preview, thumb, zoom,
                      metrics)
    return 2


def _store_smalls(cadc_client, working_directory, stream, preview_fname,
                  thumb_fname, zoom_fname, metrics):
    mc.data_put(cadc_client, working_directory, preview_fname, cn.COLLECTION,
                stream, mime_type='image/jpeg', metrics=metrics)
    mc.data_put(cadc_client, working_directory, thumb_fname, cn.COLLECTION,
                stream, mime_type='image/jpeg', metrics=metrics)
    # mc.data_put(cadc_client, working_directory, zoom_fname, cn.COLLECTION,
    #             stream, mime_type='image/jpeg', metrics=metrics)


def _augment(plane, uri, fqn, product_type):
    temp = None
    if uri in plane.artifacts:
        temp = plane.artifacts[uri]
    plane.artifacts[uri] = mc.get_artifact_metadata(
        fqn, product_type, ReleaseType.DATA, uri, temp)
