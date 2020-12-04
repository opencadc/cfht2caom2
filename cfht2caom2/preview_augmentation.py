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

import aplpy
import logging
import os

import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy.table import Table
from PIL import Image

from caom2 import ProductType, ReleaseType, ObservationIntentType
from caom2pipe import astro_composable as ac
from caom2pipe import manage_composable as mc
from cfht2caom2 import cfht_name as cn
from cfht2caom2 import metadata as md

__all__ = ['visit']


class CFHTPreview(mc.PreviewVisitor):

    def __init__(self, instrument, intent, obs_type, target, **kwargs):
        super(CFHTPreview, self).__init__(cn.ARCHIVE, **kwargs)
        self._instrument = md.Inst(instrument)
        self._intent = intent
        self._obs_type = obs_type
        self._storage_name = cn.CFHTName(instrument=self._instrument,
                                         file_name=self._science_file)
        self._science_fqn = os.path.join(self._working_dir,
                                         self._storage_name.file_name)
        self._preview_fqn = os.path.join(
            self._working_dir,  self._storage_name.prev)
        self._thumb_fqn = os.path.join(
            self._working_dir, self._storage_name.thumb)
        self._zoom_fqn = os.path.join(
            self._working_dir, self._storage_name.zoom)
        self._target_name = (target.name if target is not None else
                             self._storage_name.file_id)
        self._logger = logging.getLogger(__name__)

    def generate_plots(self, obs_id):
        if (self._instrument is md.Inst.SITELLE and
                self._storage_name.suffix == 'p'):
            count = self._sitelle_calibrated_cube()
        elif (self._instrument is md.Inst.ESPADONS and
              self._storage_name.suffix in ['i', 'p']):
            count = self._do_espadons_science()
        elif (self._instrument is md.Inst.SPIROU and
              self._storage_name.suffix in ['e', 'p', 's', 't', 'v']):
            if self._storage_name.suffix == 'v':
                count = self._do_spirou_bintable()
            else:
                count = self._do_spirou_intensity_spectrum()
        else:
            count = self._do_ds9_prev(obs_id)
        return count

    def _do_espadons_science(self):
        self._logger.debug(f'Do espadons science preview augmentation with '
                           f'{self._science_fqn}')
        # from genEspaprevperplane.py

        #Polarization scale factor
        pScale=5.0

        espadons = fits.open(self._science_fqn)
        ext = 0
        try:
            ignore = espadons[0].header.get('OBJECT')
        except LookupError:
            ext = 1
            ignore = espadons[1].header.get('OBJECT')

        bzero = espadons[ext].header.get('BZERO')
        bscale = espadons[ext].header.get('BSCALE')

        if bzero is not None and bzero > 0.0:
            # wavelength array (nm)
            sw = bscale * (espadons[ext].data[0]) - bzero
            # intensity array (normalized)
            si = bscale * (espadons[ext].data[1]) - bzero
            if self._storage_name.suffix == 'p':
                # Stokes array
                sp = bscale * (espadons[ext].data[2]) - bzero
        else:
            sw = espadons[ext].data[0]  # wavelength array (nm)
            si = espadons[ext].data[1]  # intensity array (normalized)
            if self._storage_name.suffix == 'p':
                sp = espadons[ext].data[2]  # Stokes array

        espadons.close(self._science_fqn)
        self._logger.debug(f'{espadons[ext].shape}, {sw}, {si}')

        npix = sw.shape[0]

        swa = 10.0 * sw
        sia = np.arange(0., npix, 1.0)
        spa = None
        if self._storage_name.suffix == 'p':
            spa = np.arange(0., npix, 1.0)

        # determine upper/lower y limits for two planes from intensity values
        for i in range(sia.shape[0]):
            sia[i] = float(si[i])
            if self._storage_name.suffix == 'p':
                spa[i] = float(sp[i]) * pScale  # increase scale of polarization

        fig = plt.figure(figsize=(10.24, 10.24), dpi=100)
        self._subplot(fig, swa, sia, spa, 4300.0, 4600.0, 1, 4408.0, 4412.0,
                      'Stokes spectrum (x5)')
        self._subplot(fig, swa, sia, spa, 6500.0, 6750.0, 2, 6589.0, 6593.0,
                      'Stokes spectrum (x5)')
        plt.savefig(self._preview_fqn, format='jpg')
        self.add_preview(self._storage_name.prev_uri, self._storage_name.prev,
                         ProductType.PREVIEW, ReleaseType.DATA)
        self.add_to_delete(self._preview_fqn)
        count = 1
        count += self._gen_thumbnail()
        if count == 2:
            self.add_preview(self._storage_name.thumb_uri,
                             self._storage_name.thumb, ProductType.THUMBNAIL,
                             ReleaseType.META)
            self.add_to_delete(self._thumb_fqn)
        return count

    def _subplot(self, fig, swa, sia, spa, wl_low, wl_high, subplot_index,
                 text_1, text_2, text_3):
        label = f'{self._storage_name.product_id}: {self._target_name}'
        wl = swa[(swa > wl_low) & (swa < wl_high)]
        flux = sia[(swa > wl_low) & (swa < wl_high)]
        wl_sort = wl[wl.argsort()]
        flux_sort = flux[wl.argsort()]
        if self._storage_name.suffix == 'p':
            pflux = spa[(swa > wl_low) & (swa < wl_high)]
            pflux_sort = pflux[wl.argsort()]
            flux = np.append(flux, pflux)
        ymax = 1.1 * np.max(flux)
        if np.isnan(ymax):
            ymax = 1.0
        ymin = np.min([0.0, np.min(flux) - (ymax - np.max(flux))])
        if np.isnan(ymin):
            ymin = 0.0

        # pylab.subplot(2, 1, subplot_index)
        axis = fig.add_subplot(2, 1, subplot_index)
        axis.grid(True)
        axis.plot(wl_sort, flux_sort, color='k')
        if self._storage_name.suffix == 'p':
            axis.plot(wl_sort, pflux_sort, color='b')
            axis.text(text_1, (ymin + 0.02 * (ymax - ymin)), text_3, size=16,
                      color='b')
        axis.set(title=label,
                 xlabel=r'Wavelength ($\AA$)',
                 ylabel='Relative Intensity')
        axis.title.set_weight('bold')
        axis.title.set_color('m')
        axis.text(text_2, (ymin + 0.935 * (ymax - ymin)),
                  'Intensity spectrum', size=16)
        axis.set_ylim(ymin, ymax)

    def _do_ds9_prev(self, obs_id):
        """
                    256               1024                zoom
        ESPaDOnS:
        mosaic      ''                ''                  -fits
        pan         ''                ''                  ''
        rotate      ''                ''                  ''
        scale       zscale            zscale              zscale
        scope       global            global              global
        mode        -mode none        -mode none          -mode none
        zoom        to fit            to fit              1

        MegaPrime, not 'p' and 'o':
        mosaic     -mosaicimage iraf -mosaicimage iraf    -fits
        pan        ''                ''                   -pan -9 1780
        scale      zscale            zscale               zscale
        scope      local             local                global
        mode       -mode none        ''                   -mode none
        zoom       to fit            to fit               1

        MegaPrime, 'p' and 'o':
        mosaic     -mosaicimage wcs  -mosaicimage wcs     -fits
        pan        ''                ''                   -pan -9 1780
        scale      zscale            zscale               zscale
        scope      global            global               global
        mode       -mode none        ''                   -mode none
        zoom       to fit            to fit               1

        MegaPrime extensions:
        rotate[23] -rotate 180       -rotate 180          ''
        rotate[14] -rotate 180       -rotate 180          -rotate 180
        rotate[1]  -rotate 180       -rotate 180          -rotate 180

        WIRCam 'o', 'p', 'and 's':
        mosaic     -mosaicimage wcs  -mosaicimage wcs     -fits
        rotate     ''                ''                   ''
        scale      zscale            zscale               zscale
        scope      global            global               global
        mode       -mode none        ''                   -mode none
        zoom       to fit            to fit               1

        WIRCam not 'o', 'p', 'and 's':
        mosaic     -mosaicimage iraf -mosaicimage iraf    -fits
        rotate     ''                ''                   ''
        scale      zscale            zscale               zscale
        scope      local             local                global
        mode       -mode none        ''                   -mode none
        zoom       to fit            to fit               1

        WIRCam extensions:
        pan[4]     ''               ''                    -pan 484 -484
        pan[1]     ''               ''                    -pan -484 -484

        SITELLE 2D images:
        mosaic     ''               ''                    -fits
        pan        ''               ''                    -pan -512 1544
        rotate     ''               ''                    ''
        scale      zscale           zscale                zscale
        scope      global           global                global
        mode       -mode none       -mode none            -mode none
        zoom       to fit           to fit                1

        SPIRou Raw 2D:
        mosaic     ''               ''                    -fits
        pan        ''               ''                    ''
        rotate     ''               ''                    ''
        scale      zscale           zscale                zscale
        scope      global           global                global
        mode       -mode none       -mode none            -mode none
        zoom       to fit           to fit                1

        """
        self._logger.debug(f'Do ds9 preview augmentation with '
                           f'{self._science_fqn}')
        count = 0
        delete_list = []
        headers = ac.read_fits_headers(self._science_fqn)
        num_extensions = headers[0].get('NEXTEND')

        zoom_science_fqn = self._science_fqn

        # from genWirprevperplane.py
        # if it's a datacube, just take the first slice
        # e.g. fitscopy '928690p.fits[*][*,*,1:1]' s1928690p.fits

        # set up the correct input file - may need to use fitscopy
        rotate_param = ''
        scale_param = 'zscale'
        if self._instrument is md.Inst.WIRCAM:
            if self._science_fqn.endswith('.fz'):
                naxis_3 = headers[0].get('ZNAXIS3', 1)
            else:
                naxis_3 = headers[0].get('NAXIS3', 1)

            # SF - 08-04-20 - for 'g' use fitscopy, then regular ds9 for zoom
            # calibration. This is a change from guidance of 19-03-20, which was
            # to use the previews from 'p' files for the 'g' files
            if naxis_3 != 1 or self._storage_name.suffix == 'g':
                self._logger.info(f'Observation {obs_id}: using first slice '
                                  f'of {self._science_fqn}.')
                # TODO - fix this
                if self._storage_name.file_name.endswith('.fz'):
                    temp_science_f_name = self._storage_name.file_name.replace(
                        '.fz', '_slice.fz')
                elif self._storage_name.file_name.endswith('.gz'):
                    temp_science_f_name = self._storage_name.file_name.replace(
                        '.gz', '_slice.gz')
                else:
                    temp_science_f_name = self._storage_name.file_name.replace(
                        '.fits', '_slice.fits')

                slice_cmd = f'fitscopy ' \
                            f'{self._storage_name.file_name}[*][*,*,1:1,1:1] ' \
                            f'{temp_science_f_name}'
                self._exec_cmd_chdir(temp_science_f_name, slice_cmd)
                science_fqn = f'{self._working_dir}/{temp_science_f_name}'
                delete_list.append(science_fqn)

            if num_extensions >= 4:
                self._logger.info(f'Observation {obs_id}: using slice for '
                                  f'zoom preview of {self._science_fqn}.')
                zoom_science_f_name = self._storage_name.file_name.replace(
                    '.fits', '_zoom.fits')
                slice_cmd = f'fitscopy ' \
                            f'{self._storage_name.file_name}[4][*,*,1:1] ' \
                            f'{zoom_science_f_name}'
                self._exec_cmd_chdir(zoom_science_f_name, slice_cmd)
                zoom_science_fqn = f'{self._working_dir}/{zoom_science_f_name}'
                delete_list.append(zoom_science_fqn)

        elif self._instrument in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
            rotate_param = '-rotate 180'
            # SF - 09-04-20 - mosaic MEFs i.e. number of HDU > 1
            mode_param = ''
            if num_extensions > 1:
                mode_param = '-mode none'

        # SF - 08-04-20 - change to minmax for 'm' files instead of zscale
        # 'm' is equivalent to 'MASK'
        if self._storage_name.suffix == 'm' or self._obs_type == 'MASK':
            scale_param = 'minmax'

        # set up the correct parameters to the ds9 command
        scope_param = 'local'
        if (self._instrument in [md.Inst.ESPADONS, md.Inst.SITELLE,
                                 md.Inst.SPIROU] or
                self._intent is ObservationIntentType.SCIENCE):
            scope_param = 'global'

        # 20-03-20 - seb - always use iraf - do not trust wcs coming from
        # the data acquisition. A proper one needs processing which is often
        # not done on observations.
        mosaic_param = '-mosaicimage iraf'
        if self._instrument in [md.Inst.SITELLE, md.Inst.ESPADONS,
                                md.Inst.SPIROU]:
            mosaic_param = ''

        geometry = '256x521'

        count += CFHTPreview._gen_image(self._science_fqn, geometry,
                                        self._thumb_fqn, scope_param,
                                        rotate_param,
                                        mosaic_param=mosaic_param,
                                        scale_param=scale_param)
        if count == 1:
            self.add_preview(self._storage_name.thumb_uri,
                             self._storage_name.thumb, ProductType.THUMBNAIL,
                             ReleaseType.META)
            self.add_to_delete(self._thumb_fqn)

        geometry = '1024x1024'
        if self._instrument in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
            count += CFHTPreview._gen_image(self._science_fqn, geometry,
                                            self._preview_fqn,
                                            scope_param, rotate_param,
                                            mosaic_param=mosaic_param,
                                            mode_param=mode_param,
                                            scale_param=scale_param)
        else:
            count += CFHTPreview._gen_image(self._science_fqn, geometry,
                                            self._preview_fqn,
                                            scope_param, rotate_param,
                                            mosaic_param=mosaic_param,
                                            scale_param=scale_param)
        if count == 2:
            self.add_preview(self._storage_name.prev_uri,
                             self._storage_name.prev, ProductType.PREVIEW,
                             ReleaseType.DATA)
            self.add_to_delete(self._preview_fqn)

        mosaic_param = '-fits'
        zoom_param = '1'
        scope_param = 'global'
        # set zoom parameters
        if self._instrument in [md.Inst.ESPADONS, md.Inst.SPIROU]:
            pan_param = ''
        elif self._instrument is md.Inst.WIRCAM:
            pan_param = '-pan 484 -484 image'
            if self._storage_name.suffix == 'g':
                pan_param = ''
                zoom_param = 'to fit'
        elif self._instrument in [md.Inst.MEGACAM, md.Inst.MEGAPRIME]:
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
        elif self._instrument is md.Inst.SITELLE:
            pan_param = '-pan -512 1544'
        count += CFHTPreview._gen_image(zoom_science_fqn, geometry,
                                        self._zoom_fqn, scope_param,
                                        rotate_param, zoom_param,
                                        pan_param, mosaic_param=mosaic_param,
                                        scale_param=scale_param)
        CFHTPreview._gen_square(self._zoom_fqn)
        if count == 3:
            self.add_preview(self._storage_name.zoom_uri,
                             self._storage_name.zoom, ProductType.PREVIEW,
                             ReleaseType.DATA)
            self.add_to_delete(self._zoom_fqn)
        return count

    def _do_spirou_bintable(self):
        label = f'{self._storage_name.product_id}: {self._target_name}'

        df = Table.read(self._science_fqn)
        plt.figure(figsize=(10.24, 10.24), dpi=100)
        plt.plot(df['Velocity'], df['Combined'])
        plt.title(label, weight='bold', color='m')
        plt.xlabel('Radial Velocity (km/s)')
        plt.ylabel('Weighted mean echelle order')
        plt.savefig(self._preview_fqn, format='jpg')
        return self._save_figure()

    def _do_spirou_intensity_spectrum(self):

        spirou = fits.open(self._science_fqn)
        # Polarization scale factor

        if self._storage_name.suffix in ['e', 't']:
            sw2d = spirou['WaveAB'].data  # wavelength array (nm)
            si2d = spirou['FluxAB'].data  # intensity array (normalized)
            sw = np.ravel(sw2d)
            si = np.ravel(si2d)

        if self._storage_name.suffix == 'p':
            sw2d = spirou['WaveAB'].data  # wavelength array (nm)
            si2d = spirou['StokesI'].data  # intensity array (normalized)
            sp2d = spirou['Pol'].data  # Pol Stokes array
            sw = np.ravel(sw2d)
            si = np.ravel(si2d)
            sp = np.ravel(sp2d)
            pScale = 5.0 * max(si)

        if self._storage_name.suffix == 's':
            # using uniform wavelength bins
            sw = spirou[1].data.field(0)
            si = spirou[1].data.field(1)

        spirou.close(self._science_fqn)
        npix = sw.shape[0]

        swa = 10.0 * sw
        sia = np.arange(0., npix, 1.0)
        spa = None
        if self._storage_name.suffix == 'p':
            spa = np.arange(0., npix, 1.0)
        # determine upper/lower y limits for two planes from intensity values
        for i in range(sia.shape[0]):
            sia[i] = float(si[i])
            if self._storage_name.suffix == 'p':
                spa[i] = float(sp[i]) * pScale  # increase polarization scale
        fig = plt.figure(figsize=(10.24, 10.24), dpi=100)
        self._subplot(fig, swa, sia, spa, 15000.0, 15110.0, 1, 15030.0,
                      15030.0, 'Stokes spectrum')
        self._subplot(fig, swa, sia, spa, 22940.0, 23130.0, 2, 22990.0,
                      22990.0, 'Stokes spectrum')
        plt.tight_layout()
        plt.savefig(self._preview_fqn, format='jpg')
        return self._save_figure()

    def _exec_cmd_chdir(self, temp_file, cmd):
        orig_dir = os.getcwd()
        try:
            os.chdir(self._working_dir)
            if os.path.exists(temp_file):
                os.unlink(temp_file)
            mc.exec_cmd(cmd)
        finally:
            os.chdir(orig_dir)

    @staticmethod
    def _gen_image(in_science_fqn, geometry, save_fqn, scope_param,
                   rotate_param,
                   zoom_param='to fit',
                   pan_param='', mosaic_param='',
                   mode_param='-mode none', scale_param=''):
        # 20-03-20 - seb - always use iraf - do not trust wcs coming from the
        # data acquisition. A proper one needs processing which is often not
        # done on observations.
        cmd = f'xvfb-run -a ds9 {mosaic_param} {in_science_fqn} ' \
              f'{pan_param} ' \
              f'-geometry {geometry} ' \
              f'{rotate_param} ' \
              f'-scale squared ' \
              f'-scale mode {scale_param} ' \
              f'-scale scope {scope_param} ' \
              f'-scale datasec yes ' \
              f'-invert ' \
              f'{mode_param} ' \
              f'-view colorbar no ' \
              f'-zoom {zoom_param} ' \
              f'-saveimage jpeg {save_fqn} ' \
              f'-quit'
        mc.exec_cmd(cmd, timeout=900)  # wait 15 minutes till killing
        count = 0
        if os.path.exists(save_fqn):
            count = 1
        return count

    @staticmethod
    def _gen_square(f_name):
        im = Image.open(f_name)
        min_size = 1024
        extent = max(min_size, im.size[0], im.size[1])
        # fill color is white
        new_im = Image.new('RGB', (extent, extent), (255, 255, 255, 0))
        new_im.paste(im, (int((extent - im.size[0]) / 2),
                          int((extent - im.size[1]) / 2)))
        new_im.save(f_name)

    def _gen_thumbnail(self):
        self._logger.debug(f'Generating thumbnail for file '
                           f'{self._science_fqn}.')
        count = 0
        if os.path.exists(self._preview_fqn):
            thumb = image.thumbnail(self._preview_fqn, self._thumb_fqn,
                                    scale=0.25)
            if thumb is not None:
                count = 1
        else:
            self._logger.warning(f'Could not find {self._preview_fqn} for '
                                 f'thumbnail generation.')
        return count

    def _sitelle_calibrated_cube(self):
        self._logger.debug(f'Do sitelle calibrated cube preview augmentation '
                           f'with {self._science_fqn}')
        # from genSiteprevperplane.py
        sitelle = fits.open(self._science_fqn)

        # Make a RGB colour image if it's a calibrated 3D cube
        # scan through cube to look for strongest lines
        data = sitelle[0].data
        self._logger.debug(f'{data.shape}, {data.size}')

        # trim off ends to make 2048x2048
        data = data[:, 8:2056]

        # trim off 30% of spectral edges - might be noisy
        nspecaxis = data.shape[0]
        numedgechannels = int(0.15 * nspecaxis)
        self._logger.debug(f'{numedgechannels}')

        data[:numedgechannels, :, :] = 0.0
        data[(-1 * numedgechannels):, :, :] = 0.0
        nspecaxis = data.shape[0]
        nspataxis = data.shape[1] * data.shape[2]
        self._logger.debug(f'{nspecaxis}, {nspataxis}, {data.size}, '
                           f'{data.shape}')
        data2d = np.reshape(data, (nspecaxis, -1))
        self._logger.debug(f'{data2d.shape}')

        for k in range(nspecaxis):
            medianvswavenumber = np.median(data2d[k, :])
            data2d[k, :] = data2d[k, :] - medianvswavenumber
        meanbgsubvswavenumber = np.mean(data2d, axis=1)

        self._logger.debug(f'{meanbgsubvswavenumber}, '
                           f'{meanbgsubvswavenumber.shape}')
        indexmax1 = np.nanargmax(meanbgsubvswavenumber)
        self._logger.debug(f'{indexmax1}, {meanbgsubvswavenumber[indexmax1]}')

        # remove 7 channels around strongest line
        indexmax1lo = indexmax1 - 3
        indexmax1hi = indexmax1 + 3
        meanbgsubvswavenumber[indexmax1lo:indexmax1hi] = 0.0
        indexmax2 = np.nanargmax(meanbgsubvswavenumber)
        self._logger.debug(f'{indexmax2}, {meanbgsubvswavenumber[indexmax2]}')

        # remove 7 channels around second strongest line
        indexmax2lo = indexmax2 - 3
        indexmax2hi = indexmax2 + 3
        meanbgsubvswavenumber[indexmax2lo:indexmax2hi] = 0.0
        indexmax1loline = indexmax1 - 1
        indexmax1hiline = indexmax1 + 1
        indexmax2loline = indexmax2 - 1
        indexmax2hiline = indexmax2 + 1
        self._logger.debug(f'{indexmax1loline}, {indexmax1hiline}, '
                           f'{indexmax2loline}, {indexmax2hiline}')
        self._logger.debug(f'{meanbgsubvswavenumber}')

        w = np.where(meanbgsubvswavenumber > 0.0)
        self._logger.debug(f'{w[0]}')

        head = sitelle[0].header

        head['NAXIS1'] = 1024
        head['NAXIS2'] = 1024

        head256 = head
        head256['NAXIS1'] = 256
        head256['NAXIS2'] = 256

        # Make two line images in 3 different sizes
        dataline1 = data[indexmax1loline:indexmax1hiline, :, :]
        data2dline1 = np.mean(dataline1, axis=0)
        self._logger.debug(f'{data2dline1.shape}')

        dataline2 = data[indexmax2loline:indexmax2hiline, :, :]
        data2dline2 = np.mean(dataline2, axis=0)
        self._logger.debug(f'{data2dline2.shape}')

        # Make "continuum" image with two strongest lines removed in 3
        # different sizes and add this to line image so whole image not green
        datanolines = data[w[0], :, :]
        data2dcont = np.mean(datanolines, axis=0)

        data2dline1pluscont = data2dline1 + data2dcont
        data2dline2pluscont = data2dline2 + data2dcont
        self._logger.debug(f'{np.mean(data2dline1)}, '
                           f'{np.mean(data2dline1pluscont)}, '
                           f'{np.mean(data2dline2pluscont)}, '
                           f'{np.mean(data2dcont)}')

        self._create_rgb_inputs(data2dline1pluscont, head, head256,
                                'imageline1size1024.fits',
                                'imageline1size256.fits',
                                'imageline1zoom1024.fits')
        self._create_rgb_inputs(data2dline2pluscont, head, head256,
                                'imageline2size1024.fits',
                                'imageline2size256.fits',
                                'imageline2zoom1024.fits')
        self._create_rgb_inputs(data2dcont, head, head256,
                                'imagecontsize1024.fits',
                                'imagecontsize256.fits',
                                'imagecontzoom1024.fits')

        os.system("pwd")
        del data
        del datanolines
        sitelle.close(self._science_fqn)

        self._create_rgb('imageline1size1024.fits', 'imageline2size1024.fits',
                         'imagecontsize1024.fits', self._preview_fqn)
        self._create_rgb('imageline1size256.fits', 'imageline2size256.fits',
                         'imagecontsize256.fits', self._thumb_fqn)
        self._create_rgb('imageline1zoom1024.fits', 'imageline2zoom1024.fits',
                         'imagecontzoom1024.fits', self._zoom_fqn)
        self.add_to_delete('./imageline1size1024.fits')
        self.add_to_delete('./imageline1size256.fits')
        self.add_to_delete('./imageline1zoom1024.fits')
        self.add_to_delete('./imageline2size1024.fits')
        self.add_to_delete('./imageline2size256.fits')
        self.add_to_delete('./imageline2zoom1024.fits')
        self.add_to_delete('./imagecontsize1024.fits')
        self.add_to_delete('./imagecontsize256.fits')
        self.add_to_delete('./imagecontzoom1024.fits')
        self.add_preview(self._storage_name.thumb_uri,
                         self._storage_name.thumb, ProductType.THUMBNAIL,
                         ReleaseType.META)
        self.add_preview(self._storage_name.prev_uri, self._storage_name.prev,
                         ProductType.PREVIEW, ReleaseType.DATA)
        self.add_to_delete(self._thumb_fqn)
        self.add_to_delete(self._preview_fqn)
        self.add_preview(self._storage_name.zoom_uri,
                         self._storage_name.zoom, ProductType.PREVIEW,
                         ReleaseType.DATA)
        self.add_to_delete(self._zoom_fqn)
        return 3

    def _create_rgb_inputs(self, input_data, head, head256, preview_f_name,
                           thumb_f_name, zoom_f_name):
        size1024 = self._rebin_factor(input_data, (1024, 1024))
        size256 = self._rebin_factor(input_data, (256, 256))
        zoom1024 = input_data[512:1536, 512:1536]
        self._logger.debug(f'{size1024.shape}, {size256.shape}, '
                           f'{zoom1024.shape}')
        fits.writeto(preview_f_name, size1024, head, overwrite=True)
        fits.writeto(thumb_f_name, size256, head256, overwrite=True)
        fits.writeto(zoom_f_name, zoom1024, head, overwrite=True)

    @staticmethod
    def _create_rgb(line1_f_name, line2_f_name, cont_f_name, fqn):
        aplpy.make_rgb_image([line1_f_name, line2_f_name, cont_f_name], fqn,
                             stretch_r='linear', stretch_g='linear',
                             stretch_b='linear', pmax_r=99.5, pmax_g=99.5,
                             pmax_b=99.5, pmin_r=50.0, pmin_g=95.0,
                             pmin_b=50.0, embed_avm_tags=False)

    def _rebin_factor(self, a, new_shape):
        """
        Re-bin an array to a new shape.

        :param new_shape must be a factor of a.shape.
        """
        assert len(a.shape) == len(new_shape)
        assert not np.sometrue(np.mod(a.shape, new_shape))

        slices = [slice(None, None, mc.to_int(old / new)) for old, new in
                  zip(a.shape, new_shape)]
        self._logger.debug(slices)
        return a[slices]

    def _save_figure(self):
        self.add_to_delete(self._preview_fqn)
        count = 1
        self.add_preview(self._storage_name.prev_uri, self._storage_name.prev,
                         ProductType.PREVIEW, ReleaseType.DATA)
        count += self._gen_thumbnail()
        if count == 2:
            self.add_preview(self._storage_name.thumb_uri,
                             self._storage_name.thumb, ProductType.THUMBNAIL,
                             ReleaseType.META)
            self.add_to_delete(self._thumb_fqn)
        return count


def visit(observation, **kwargs):
    previewer = CFHTPreview(observation.instrument.name,
                            observation.intent, observation.type,
                            observation.target, **kwargs)
    cfht_name = cn.CFHTName(instrument=md.Inst(observation.instrument.name),
                            file_name=previewer.science_file)
    previewer.storage_name = cfht_name
    return previewer.visit(observation, cfht_name)
