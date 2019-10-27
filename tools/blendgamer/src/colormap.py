# ***************************************************************************
# This file is part of the GAMer software.
# Copyright (C) 2016-2018
# by Christopher Lee, Tom Bartol, John Moody, Rommie Amaro, J. Andrew McCammon,
#    and Michael Holst

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# ***************************************************************************

import bpy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from blendgamer.util import *

# Recommended backend if a non-interactive backend is the default
# $ pip install pyqt5
# mpl.use('QT5Agg')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# Mapping colormap_enums.colormap_enums to values
colormapDict = {
  'VIRIDIS': plt.cm.viridis,
  'PRGN' : plt.cm.PRGn,
}

class DivergingNorm(mpl.colors.Normalize):
    def __init__(self, vcenter, vmin=None, vmax=None):
        """
        Normalize data with a set center.

        Useful when mapping data with an unequal rates of change around a
        conceptual center, e.g., data that range from -2 to 4, with 0 as
        the midpoint.

        Parameters
        ----------
        vcenter : float
            The data value that defines ``0.5`` in the normalization.
        vmin : float, optional
            The data value that defines ``0.0`` in the normalization.
            Defaults to the min value of the dataset.
        vmax : float, optional
            The data value that defines ``1.0`` in the normalization.
            Defaults to the the max value of the dataset.

        Examples
        --------
        This maps data value -4000 to 0., 0 to 0.5, and +10000 to 1.0; data
        between is linearly interpolated::

            >>> import matplotlib.colors as mcolors
            >>> offset = mcolors.DivergingNorm(vmin=-4000.,
                                               vcenter=0., vmax=10000)
            >>> data = [-4000., -2000., 0., 2500., 5000., 7500., 10000.]
            >>> offset(data)
            array([0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.0])
        """

        self.vcenter = vcenter
        self.vmin = vmin
        self.vmax = vmax
        if not (vcenter and vmin and vmax) and (vcenter >= vmax or vcenter <= vmin):
            raise ValueError('vmin(%f), vcenter(%f), and vmax(%f) must be in '
                             'ascending order'%(vmin, vcenter, vmax))

    def autoscale_None(self, A):
        """
        Get vmin and vmax, and then clip at vcenter
        """
        super().autoscale_None(A)
        if self.vmin > self.vcenter:
            self.vmin = self.vcenter
        if self.vmax < self.vcenter:
            self.vmax = self.vcenter


    def __call__(self, value, clip=False):
        """
        Map value to the interval [0, 1].
        """
        result, is_scalar = self.process_value(value)
        self.autoscale_None(result)  # sets self.vmin, self.vmax if None

        if not self.vmin <= self.vcenter <= self.vmax:
            raise ValueError("vmin, vcenter, vmax must increase monotonically")
        result = np.ma.masked_array(
            np.interp(result, [self.vmin, self.vcenter, self.vmax],
                      [0, 0.5, 1.]), mask=np.ma.getmask(result))
        if is_scalar:
            result = np.atleast_1d(result)[0]
        return result


def curveToData(crv, context, report):
    """
    Helper function to take a curvature object and return smoothed data.

    :param      crv:       Curvature object
    :type       crv:       RNA object with curvature meta information
    :param      context:   Context of the scene
    :type       context:   bpy context

    :returns:   Array of converted
    :rtype:     numpy.ndarray
    """
    obj = getActiveMeshObject(report)
    bm = bmesh_from_object(obj)

    with ObjectMode():
        layer = getCurvatureLayer(obj, crv.algorithm, crv.curvatureType)

    data = np.zeros(len(obj.data.vertices), dtype=np.float)
    # Copy curvatures over
    layer.foreach_get('value', data)

    if crv.curveIter > 0:
        for i in range(0, crv.curveIter):
            tmp = np.zeros(len(obj.data.vertices), dtype=np.float)
            for v in bm.verts:
                count = 1
                tmp[v.index] = data[v.index]

                for e in v.link_edges:
                    v_other = e.other_vert(v)

                    tmp[v.index] += data[v_other.index]
                    count += 1
                tmp[v.index] /= count
            data = np.array(tmp, copy=True)
    return data


def dataToVertexColor(crv, context, report, showplot=False, saveplot=False):
    """
    Convert curvature object to colormap

    :param      crv:       Curvature object
    :type       crv:       RNA object with curvature meta information
    :param      context:   Context of the scene
    :type       context:   bpy context
    :param      report:    Reporter to return info to the UI
    :type       report:    bpy reporter
    :param      showplot:  Show the generated plots
    :type       showplot:  bool
    :param      saveplot:  Save plots to file
    :type       saveplot:  bool
    """

    data = curveToData(crv, context, report)
    cmap = colormapDict[crv.colormap]
    file_prefix = "%s_%s_m%dM%dI%dmx%0.2f%s"%(
            bpy.path.basename(bpy.context.blend_data.filepath).split('.')[0],
            context.object.name,
            crv.minCurve,
            crv.maxCurve,
            crv.curveIter,
            crv.mixpoint,
            crv.curvatureType)

    fig = plt.figure(figsize=(8,5))
    ax = fig.add_axes([0.1,0.05,0.6,0.9])

    ## This code helps make the plots easier to read...
    # tmin = np.percentile(data,1)
    # tmax = np.percentile(data,99)
    # data[data < tmin] =tmin
    # data[data > tmax] =tmax
    #
    tmin = crv.minCurve
    tmax = crv.maxCurve
    amin = np.amin(data)
    amax = np.amax(data)
    amean = np.mean(data)
    amedian = np.median(data)

    if showplot:
        ax.hist(data, bins='auto')
        ax.set_title("%s Distribution"%(file_prefix))
        ax.axvline(amin, color='r', linestyle='dashed', linewidth=1)
        ax.axvline(amax, color='r', linestyle='dashed', linewidth=1)

    extend = 'neither'
    # the tmin/tmax values are percentiles instead
    if crv.limitsArePercentiles:
        # print("Using min/max values as percentiles!")
        if tmin < 0 or tmin >= 100:
            print("Minimum percentile must be 0<=x<100. Setting to 0")
            lowerPercentile = 0
        else:
            lowerPercentile = int(tmin)

        if tmax <= 0 or tmax > 100:
            print("Maximum percentile must be 0<x<=100. Setting to 100")
            upperPercentile = 100
        else:
            upperPercentile = int(tmax)

        tmin = np.percentile(data,lowerPercentile)
        tmax = np.percentile(data,upperPercentile)
        print("Data truncated at %f and %f percentiles\n"%(lowerPercentile,upperPercentile) )

    if showplot:
        ax.axvline(tmin, color='g', linestyle='dashed', linewidth=2)
        ax.axvline(tmax, color='g', linestyle='dashed', linewidth=2)

        if tmin > amin and tmax < amax:
            extend = 'both'
        elif tmin > amin:
            extend = 'min'
        elif tmax < amax:
            extend = 'max'

    data[data < tmin] =tmin
    data[data > tmax] =tmax

    amin = np.amin(data)
    amax = np.amax(data)

    # Construct the norm and colorbar
    if amin < 0 and amax > 0:
        norm = DivergingNorm(vmin=amin, vcenter=0, vmax=amax)
        # Python 3.5 matplotlib may not support?
        # norm = mpl.colors.DivergingNorm(vmin=amin, vcenter=0, vmax=amax)
        colors_neg = cmap(np.linspace(0, crv.mixpoint, 256))
        colors_pos = cmap(np.linspace(crv.mixpoint, 1, 256))

        all_colors = np.vstack((colors_neg, colors_pos))
        curvature_map = mpl.colors.LinearSegmentedColormap.from_list('curvature_map', all_colors)
    else:
        norm = mpl.colors.Normalize(vmin=amin, vmax=amax)
        curvature_map = cmap

    # Map values to colors and add vertex color layer
    colors = curvature_map(norm(data))

    # Create view without alpha channel if Blender < 2.80
    if bpy.app.version < (2,80,0):
        colors = colors[:,:3]

    mesh = bpy.context.object.data
    vlayer = "%s%s"%(crv.algorithm, crv.curvatureType)

    if vlayer not in mesh.vertex_colors:
        if len(mesh.vertex_colors) == 8:
            report({'ERROR'}, "Maximum of 8 vertex Layers reached cannot create a new layer. Please delete a layer to continue.")
            return False
        mesh.vertex_colors.new(name=vlayer)

    color_layer = mesh.vertex_colors[vlayer]
    mesh.vertex_colors[vlayer].active = True

    mloops = np.zeros((len(mesh.loops)), dtype=np.int)
    mesh.loops.foreach_get("vertex_index", mloops)
    color_layer.data.foreach_set("color", colors[mloops].flatten())

    # Add axis for colorbar and plot it
    ax = fig.add_axes([0.75,0.05,0.05,0.9])
    cb = mpl.colorbar.ColorbarBase(ax, cmap=curvature_map, norm=norm,
                orientation='vertical')

    ticks = cb.get_ticks()
    ticks.sort()

    if amin != ticks[0]:
        ticks = np.insert(ticks, 0, amin)
    if amax != ticks[-1]:
        ticks = np.append(ticks, amax)
    cb.set_ticks(ticks)

    ticklabels = [r"{:0.1f}".format(tick) for tick in ticks]

    if extend == 'neither':
       pass
    elif extend == 'both':
        ticklabels[0] = "< " + ticklabels[0]
        ticklabels[-1] = "> " + ticklabels[-1]
    elif extend == 'max':
        ticklabels[-1] = "> " + ticklabels[-1]
    elif extend == 'min':
        ticklabels[0] = "< " + ticklabels[0]
    cb.set_ticklabels(ticklabels)
    cb.ax.tick_params(labelsize=14)
    cb.set_label("%s [$\mu m^{-1}$]"%(vlayer), size=16)

    if saveplot:
        plt.savefig(file_prefix+'.pdf', format='pdf')
    if showplot:
        plt.show()
    plt.close()
    return True


def differencePlotter(context, report, difftype='K1'):
    obj = getActiveMeshObject(report)

    with ObjectMode():
        mdbsk1 = getCurvatureLayer(obj, 'MDSB', difftype)
        jetsk1 = getCurvatureLayer(obj, 'JETS', difftype)

        mdsbk1_data = np.zeros(len(obj.data.vertices), dtype=np.float)
        mdbsk1.foreach_get('value', mdsbk1_data)

        jetsk1_data = np.zeros(len(obj.data.vertices), dtype=np.float)
        jetsk1.foreach_get('value', jetsk1_data)

    data = mdsbk1_data - jetsk1_data

    cmap = colormapDict['PRGN']
    file_prefix = "%s_difference"%(difftype)

    fig = plt.figure(figsize=(8,5))
    ax = fig.add_axes([0.1,0.05,0.6,0.9])

    amin = np.amin(data)
    amax = np.amax(data)
    amean = np.mean(data)
    amedian = np.median(data)

    plt.hist(data, bins='auto')
    plt.title("%s Distribution"%(file_prefix))
    plt.axvline(amin, color='r', linestyle='dashed', linewidth=1)
    plt.axvline(amax, color='r', linestyle='dashed', linewidth=1)


    extend = 'neither'
    tmin = amin
    tmax = amax
    # tmin = np.percentile(data,2)
    # tmax = np.percentile(data,98)

    ax.axvline(tmin, color='g', linestyle='dashed', linewidth=2)
    ax.axvline(tmax, color='g', linestyle='dashed', linewidth=2)

    if tmin > amin and tmax < amax:
        extend = 'both'
    elif tmin > amin:
        extend = 'min'
    elif tmax < amax:
        extend = 'max'

    data[data < tmin] =tmin
    data[data > tmax] =tmax

    amin = np.amin(data)
    amax = np.amax(data)

    # Construct the norm and colorbar
    if amin < 0 and amax > 0:
        norm = DivergingNorm(vmin=amin, vcenter=0, vmax=amax)
        # Python 3.5 matplotlib may not support?
        # norm = mpl.colors.DivergingNorm(vmin=amin, vcenter=0, vmax=amax)
        colors_neg = cmap(np.linspace(0, .5, 256))
        colors_pos = cmap(np.linspace(.5, 1, 256))

        all_colors = np.vstack((colors_neg, colors_pos))
        curvature_map = mpl.colors.LinearSegmentedColormap.from_list('curvature_map', all_colors)
    else:
        norm = mpl.colors.Normalize(vmin=amin, vmax=amax)
        curvature_map = cmap

    # Map values to colors and add vertex color layer
    colors = curvature_map(norm(data))

    # Create view without alpha channel if Blender < 2.80
    if bpy.app.version < (2,80,0):
        colors = colors[:,:3]

    mesh = bpy.context.object.data
    vlayer = "%s_diff"%(difftype)

    if vlayer not in mesh.vertex_colors:
        mesh.vertex_colors.new(name=vlayer)

    color_layer = mesh.vertex_colors[vlayer]
    mesh.vertex_colors[vlayer].active = True

    mloops = np.zeros((len(mesh.loops)), dtype=np.int)
    mesh.loops.foreach_get("vertex_index", mloops)
    color_layer.data.foreach_set("color", colors[mloops].flatten())

    # Add axis for colorbar and plot it
    ax = fig.add_axes([0.75,0.05,0.05,0.9])
    cb = mpl.colorbar.ColorbarBase(ax, cmap=curvature_map, norm=norm,
                orientation='vertical')

    ticks = cb.get_ticks()
    ticks.sort()

    if amin != ticks[0]:
        ticks = np.insert(ticks, 0, amin)
    if amax != ticks[-1]:
        ticks = np.append(ticks, amax)
    cb.set_ticks(ticks)

    ticklabels = [r"{:0.1f}".format(tick) for tick in ticks]

    extend = 'neither'

    if extend == 'neither':
       pass
    elif extend == 'both':
        ticklabels[0] = "< " + ticklabels[0]
        ticklabels[-1] = "> " + ticklabels[-1]
    elif extend == 'max':
        ticklabels[-1] = "> " + ticklabels[-1]
    elif extend == 'min':
        ticklabels[0] = "< " + ticklabels[0]
    cb.set_ticklabels(ticklabels)
    cb.ax.tick_params(labelsize=14)
    cb.set_label("%s [$\mu m^{-1}$]"%(vlayer), size=16)

    plt.show()
    plt.close()



def eng_notation(x,pos):
    num, power = '{:.1e}'.format(x).split('e')
    power=int(power)
    return r'${} \times 10^{{{}}}$'.format(num, power)