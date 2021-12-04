# -*- coding: utf-8 -*-
"""
NetCDFInterpolate is a python package for easy accessing VTU/PVD files as
outputed by Finite Element software like OpenGeoSys. It uses the VTK python
wrapper and linear interpolation between time steps and grid points access
any points in and and time within the simulation domain.

Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
            Distributed under a Modified BSD License.
              See accompanying file LICENSE or
              http://www.opengeosys.org/project/license

"""

# pylint: disable=C0103, R0902, R0914, R0913
import numpy as np
import pandas as pd
import netCDF4 as nc4
from scipy.interpolate import griddata
from scipy.interpolate import interp1d

class NetCDFInterpolate:
    """
    Class for interpolating HDF5 griddata

    Parameters
    ----------
    filename : `str`
    nneighbors : `int`, optional
                 default: 20
    dim : `int`, optional
          default: 3
    one_d_axis : `int`
                 between 0 and 2, default: 0
    group : `str`
                 group where to find data
    subgrou : `str`
                 subgroup of group where data is stored
    """
    def __init__(self, filename, nneighbors=20, dim=3, one_d_axis=0, two_d_planenormal=2,
            group=None, subgroup=None):
        self.fileobject = nc4.Dataset(filename)
        try:
            if subgroup is None:
                self.data = self.fileobject.groups[group]
            else:
                self.data = self.fileobject.groups[group].groups[subgroup]
        except KeyError:
            for grp in list(self.fileobject.groups):
                for subgrp in list(self.fileobject.groups[grp].groups):
                    self.data = self.fileobject.groups[grp].groups[subgrp]
        self.points = self.data.variables["geometry"][:][0]
        self.times = self.fileobject.variables["times"][:]
        self.dim = dim
        self.nneighbors = nneighbors
        self.one_d_axis=one_d_axis
        self.two_d_planenormal = two_d_planenormal
        if self.dim == 1:
            self.one_d_axis = one_d_axis
            self.points = self.points[:,one_d_axis]
        if self.dim == 2:
            self.plane = [0, 1, 2]
            self.plane.pop(two_d_planenormal)
            self.points = np.delete(self.points, two_d_planenormal, 1)


    @property
    def header(self):
        header_dict = {
                # "N Cells": [str(len(self.cell_center_points))],
                "N Points": [str(len(self.points))],
                # "N Cell Arrays": [len(self.get_cell_field_names())],
                "N Point Arrays": [len(self.get_point_field_names())],
                "X Bounds": [str(np.min(self.points[:,0])) + ", " + str(np.max(self.points[:,0]))]}
        if self.dim > 1:
            header_dict["Y Bounds"] = [str(np.min(self.points[:,1]))+" "+str(np.max(self.points[:,1]))]
        if self.dim > 2:
            header_dict["Z Bounds"] = [str(np.min(self.points[:,2]))+" "+str(np.max(self.points[:,2]))]
        df = pd.DataFrame(header_dict).T
        return df.rename({0:"Information"}, axis='columns')

    @property
    def data_arrays(self):
        pf_names = self.get_point_field_names()
        # cf_names = self.get_cell_field_names()
        data_array_dict = {}
        for name in pf_names:
            field = self.get_point_field(name)[0]
            if field.ndim == 1:
                components = 1
            else:
                components = field.shape[1]
            data_array_dict[name] = ["point", components, np.min(field), np.max(field)]
        #for name in cf_names:
        #    field = self.get_cell_field(name)
        #    data_array_dict[name] = ["cell", components, np.min(field), np.max(field)]
        df = pd.DataFrame(data_array_dict).T
        return df.rename({0:"type", 1: "components", 2: "Min", 3: "Max"}, axis='columns')

    def get_neighbors(self, points_interpol):
        """
        Method for obtaining neighbor points for interpolation.
        """
        points = self.points
        df = pd.DataFrame(points)
        neighbors = {}
        if self.dim == 1:
            return neighbors
        for i, (_, val) in enumerate(points_interpol.items()):
            if self.dim == 2:
                x, y = self.plane
                df["r_"+str(i)] = (df[x]-val[x]) * (df[x]-val[x]) + (df[y]-val[y]) * (df[y]-val[y])
            elif self.dim == 3:
                df["r_"+str(i)] = ((df[0]-val[0]) * (df[0]-val[0]) + (df[1]-val[1]) * (df[1]-val[1])
                        + (df[2]-val[2]) * (df[2]-val[2]))
            neighbors[i] = df.sort_values(by=["r_" + str(i)]).head(self.nneighbors).index
        return neighbors

    def get_nearest_points(self, points_interpol):
        """
        Return a dictionary with closest mesh points

        Parameters
        ----------
        points_interpol : `dict`
        """
        nb = self.get_neighbors(points_interpol)
        nearest = {}
        for i, (key, _) in enumerate(points_interpol.items()):
            nearest[key] = self.points[nb[i][0]]
        return nearest

    def get_nearest_indices(self, points_interpol):
        """
        Return a dictionary with closest mesh point indices

        Parameters
        ----------
        points_interpol : `dict`
        """
        nb = self.get_neighbors(points_interpol)
        nearest = {}
        for i, (key, _) in enumerate(points_interpol.items()):
            nearest[key] = nb[i][0]
        return nearest

    def get_data_scipy(self, neighbors, points_interpol, fieldname, timestep, interpolation_method="linear"):
        """
        Get interpolated data for points_interpol using neighbor points.
        """
        field = self.get_point_field(fieldname, timestep=timestep)
        points = self.points
        resp = {}
        for i, (key, val) in enumerate(points_interpol.items()):
            if self.dim == 1:
                data = pd.DataFrame(points, columns = ['x'])
                data["y"] = field
                data.sort_values(by = ['x'], inplace=True)
                data.drop_duplicates(subset=['x'], inplace=True)
                f = interp1d(data['x'], data['y'])
                resp[key] = f(val[self.one_d_axis])
            elif self.dim == 2:
                x, y = self.plane
                grid_x, grid_y = np.array([[[val[x]]],[[val[y]]]])
                resp[key] = griddata(points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y), method=interpolation_method)[0][0]
            else:
                grid_x, grid_y, grid_z = np.array([[[[val[0]]]], [[[val[1]]]], [[[val[2]]]]])
                resp[key] = griddata(points[neighbors[i]], field[neighbors[i]],
                        (grid_x, grid_y, grid_z), method=interpolation_method)[0][0][0]
        return resp

    def get_point_field(self, fieldname, timestep=None):
        """
        Return vtu point field as numpy array.
        Parameters
        ----------
        fieldname : `str`
        """
        if timestep is None:
            field = self.data.variables[fieldname][:]
        else:
            field = self.data.variables[fieldname][timestep, :]
        return field

    def get_point_field_names(self):
        """
        Get names of all point fields in the vtu file.
        """
        fieldnames = [fieldname for fieldname in self.data.variables]
        return fieldnames

    def read_data(self, fieldname, timestep, pts = None, interpolation_method="linear"):
        """
        Get data of field "fieldname" at all points specified in "pts".

        Parameters
        ----------
        fieldname : `str`
        timestep : `int`
        pts : `dict`, optional
              default: {'pt0': (0.0,0.0,0.0)}
        interpolation_method : `str`, optional
                               default: 'linear'
        """
        if pts is None:
            pts = {'pt0': (0.0,0.0,0.0)}
        resp = {}
        for pt in pts:
            if isinstance(fieldname, str):
                resp[pt] = []
            elif isinstance(fieldname, list):
                resp[pt] = {}
                for field in fieldname:
                    resp[pt][field] = []
        nb = self.get_neighbors(pts)
        if isinstance(fieldname, str):
            data = self.get_data_scipy(nb, pts, fieldname, timestep,
                        interpolation_method=interpolation_method)
            for pt in pts:
                resp[pt] = data[pt]
        elif isinstance(fieldname, list):
            data = {}
            for field in fieldname:
                data[field] = self.get_data_scipy(nb, pts, field, timestep,
                            interpolation_method=interpolation_method)
            for pt in pts:
                for field in fieldname:
                    resp[pt][field] = data[field][pt]
        return resp

    def get_set_data(self, fieldname, timestep, pointsetarray=None, interpolation_method="linear"):
        """
        Get data specified in fieldname at all points specified in "pointsetarray".
        Parameters
        ----------
        fieldname : `str`
        timestep : `int`
        pointsetarray : `list`, optional
                        default: [(0,0,0)]
        interpolation_method : `str`, optional
                               default: 'linear'
        """
        if pointsetarray is None:
            pointsetarray = [(0,0,0)]
        pts = {}
        # convert into point dictionary
        for i, entry in enumerate(pointsetarray):
            pts['pt'+str(i)] = entry
        resp = self.read_data(fieldname, timestep, pts=pts, interpolation_method=interpolation_method)
        resp_list = []
        # convert point dictionary into list
        for i, entry in enumerate(pointsetarray):
            resp_list.append(resp['pt' + str(i)])
        resp_array = np.array(resp_list)
        return resp_array

    def read_set_data(self, fieldname, time, pointsetarray = None, interpolation_method="linear"):
        """
        Get data of field "fieldname" at time "timestep" alon a given "pointsetarray".

        Parameters
        ----------
        timestep : `int`
        fieldname : `str`
        pointsetarray : `array`, optional
                        default: [(0,0,0)]
        interpolation_method : `str`
                               default: 'linear'
        """
        if pointsetarray is None:
            pointsetarray = [(0,0,0)]
        for timestep, ts in enumerate(self.times):
            if time == ts:
                field = self.get_set_data(fieldname, timestep, pointsetarray, interpolation_method=interpolation_method)
                return field
        time1 = 0.0
        time2 = 0.0
        timestep = 0
        for i, ts in enumerate(self.times):
            try:
                if ts < time < self.times[i+1]:
                    time1 = ts
                    time2 = self.times[i+1]
                    timestep = i
            except IndexError:
                print("time step is out of range")
        field1 = self.get_set_data(fieldname, timestep, pointsetarray, interpolation_method=interpolation_method)
        field2 = self.get_set_data(fieldname, timestep+1, pointsetarray, interpolation_method=interpolation_method)
        fieldslope = (field2-field1)/(time2-time1)
        field = field1 + fieldslope * (time-time1)
        return field

    def read_time_series(self, fieldname, pts=None, interpolation_method="linear"):
        """
        Return time series data of field "fieldname" at points pts.
        Also a list of fieldnames can be provided as "fieldname"

        Parameters
        ----------
        fieldname : `str`
        pts : `dict`, optional
              default: {'pt0': (0.0,0.0,0.0)}
        interpolation_method : `str`, optional
                               default: 'linear
        """
        if pts is None:
            pts = {'pt0': (0.0,0.0,0.0)}
        resp_t = {}
        for pt in pts:
            if isinstance(fieldname, str):
                resp_t[pt] = []
            elif isinstance(fieldname, list):
                resp_t[pt] = {}
                for field in fieldname:
                    resp_t[pt][field] = []
        for timestep, _ in enumerate(self.times):
            if timestep == 0:
                nb = self.get_neighbors(pts)
            if isinstance(fieldname, str):
                data = self.get_data_scipy(nb, pts, fieldname, timestep,
                            interpolation_method=interpolation_method)
                for pt in pts:
                    resp_t[pt].append(data[pt])
            elif isinstance(fieldname, list):
                data = {}
                for field in fieldname:
                    data[field] = self.get_data_scipy(nb, pts, field, timestep,
                                interpolation_method=interpolation_method)
                for pt in pts:
                    for field in fieldname:
                        resp_t[pt][field].append(data[field][pt])
        resp_t_array = {}
        for pt, field in resp_t.items():
            if isinstance(fieldname, str):
                resp_t_array[pt] = np.array(field)
            elif isinstance(fieldname, list):
                resp_t_array[pt] = {}
                for field_, fieldarray in resp_t[pt].items():
                    resp_t_array[pt][field_] = np.array(fieldarray)
        return resp_t_array

class XDMFreader:
    """
    Interface for XDMF data

    """
    def __init__(self):
        pass
