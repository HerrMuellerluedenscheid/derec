# An example from scipy cookbook demonstrating the use of numpy arrys in vtk 

import vtk
import numpy as num
from collections import defaultdict
import matplotlib.pyplot as plt
from gmtpy import GMT

def make_3D_pseudo_data():
    n_samples = 75

    D2_data = {}
    D3_data = {}

    for i in range(n_samples):
        D2_data[i] = num.array(num.random.random(n_samples))

    for i in range(n_samples):
        D2_data[i] *= i
        D3_data[i] = D2_data

    return D3_data 

def scale_2_int_perc(a):
    '''
    Scale entire matrix to range between 0 and 100 uint8.
    '''
    return num.uint8((a / num.amax(a))*1e2)


def dict_2_3d_array(ddict):
    '''
    convert dictionary to 3d numpy array. 
    '''
    cube = []
    for k_1, val_1 in ddict.iteritems():
        sheet = []
        for k_2, val_2 in val_1.iteritems():
            sheet.append(val_2) 
        cube.append(sheet)
    return num.array(cube)


def gmt_map(event_lats=None, event_lons=None, station_lats=None,
        station_lons=None, **kwargs):

    with_stations = False
    if kwargs.get('stations', False):
        with_stations = True

    with_events= False
    if kwargs.get('events', False):
        with_events = True

    gmt = GMT( config={'BASEMAP_TYPE':'fancy'} )

    lat_min = min(station_lats+ event_lats)
    lat_max = max(station_lats+ event_lats)
    lon_min = min(station_lons+ event_lons)
    lon_max = max(station_lons+ event_lons)

    lat_offset = (lat_max-lat_min)*0.3
    lon_offset = (lon_max-lon_min)*0.3

    gmt.pscoast( R='%i/%i/%i/%i'%(lon_min-lon_offset, 
                                   lon_max+lon_offset, 
                                   lat_min-lat_offset, 
                                   lat_max+lat_offset),
                J='M10c',
                B='4g4',
                D='f',
                S=(114,159,207),
                G=(233,185,110),
                W='thinnest')

    if station_lats and station_lons and with_stations:
        gmt.psxy('-St0.5c', R=True, J=True, G=(0,255,123),
                    in_columns=[station_lons, station_lats])

    if event_lats and event_lons and with_events:
        gmt.psxy('-Sa0.2c', R=True, J=True, G=(255,0,0), 
                    in_columns=[event_lons, event_lats])

    gmt.save('mapplot.pdf')


class OpticBase():
    def __init__(self, test_case):

        self.test_case=test_case 
        self.test_case.numpy_it()
        self.xdim = len(self.test_case.test_parameters[test_case.xkey])
        self.ydim = len(self.test_case.test_parameters[test_case.ykey])
        self.zdim = len(self.test_case.test_parameters[test_case.zkey])

        data=test_case.num_array[3].reshape(self.xdim,
                                            self.ydim,
                                            self.zdim)
        data = scale_2_int_perc(data)
        #self.vtkCube(data)

    def value_to_index(k, val):
        return self.dim_mapper[k].index(val)

    def add_cases(self, cases):
        self.test_tin.extend(cases)
        self.update_dimensions(cases)

 
    def plot_1d(self, fix_parameters={}):
        '''
        fix parameters is a dict, e.g. {lat:10, lon:10}
        '''
        x, y = self.test_tin.numpyrize_1d(fix_parameters)
        plt.plot(x,y, 'o')
        plt.show()


    def plot_2d(self, fix_parameters={}):
        '''
        fix parameters is a dict, e.g. {lat:10, lon:10}
        '''
        x, y = self.test_tin.numpyrize_2d(fix_parameters)
        plt.plot(x,y, 'o')
        plt.show()


    def gmt_map(self, **kwargs):
        sources_lon_lat = set()
        for s in self.test_case.sources:
            sources_lon_lat.add((s.lon, s.lat))
        sources_lon_lat = zip(*list(sources_lon_lat))

        targets_lon_lat = set()
        for t in self.test_case.targets:
            targets_lon_lat.add((t.lon, t.lat))
        targets_lon_lat = zip(*list(targets_lon_lat))
        gmt_map(sources_lon_lat[1], sources_lon_lat[0],
                targets_lon_lat[1], targets_lon_lat[0],
                    **kwargs)


    def vtkCube(self, data_matrix=None):

        # We begin by creating the data we want to render.
        # For this tutorial, we create a 3D-image containing three overlaping cubes.
        # This data can of course easily be replaced by data from a medical CT-scan or anything else three dimensional.
        # The only limit is that the data must be reduced to unsigned 8 bit or 16 bit integers.
        #data_matrix = zeros([75, 75, 75], dtype=uint8)
        #data_matrix[0:35, 0:35, 0:35] = 50
        #data_matrix[25:55, 25:55, 25:55] = 100
        #data_matrix[45:74, 45:74, 45:74] = 150

        # For VTK to be able to use the data, it must be stored as a VTK-image. This can be done by the vtkImageImport-class which
        # imports raw data and stores it.
        dataImporter = vtk.vtkImageImport()
        # The preaviusly created array is converted to a string of chars and imported.
        data_string = data_matrix.tostring()
        dataImporter.CopyImportVoidPointer(data_string, len(data_string))
        # The type of the newly imported data is set to unsigned char (uint8)
        dataImporter.SetDataScalarTypeToUnsignedChar()
        # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
        # must be told this is the case.
        dataImporter.SetNumberOfScalarComponents(1)
        # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
        # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
        # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
        # VTK complains if not both are used.
        dataImporter.SetDataExtent(0, 9, 0, 9, 0, 9)
        dataImporter.SetWholeExtent(0, 9, 0, 9, 0, 9)
        #dataImporter.SetDataExtent(0, 74, 0, 74, 0, 74)
        #dataImporter.SetWholeExtent(0, 74, 0, 74, 0, 74)

        # The following class is used to store transparencyv-values for later retrival. In our case, we want the value 0 to be
        # completly opaque whereas the three different cubes are given different transperancy-values to show how it works.
        alphaChannelFunc = vtk.vtkPiecewiseFunction()
        alphaChannelFunc.AddPoint(0, 0.6)
        alphaChannelFunc.AddPoint(33, 0.2)
        alphaChannelFunc.AddPoint(66, 0.1)
        alphaChannelFunc.AddPoint(100, 0.01)

        # Gradient opacity
        # other way: misfit 0 is anti opacity
        volumeGradientOpacity = vtk.vtkPiecewiseFunction()
        volumeGradientOpacity.AddPoint(70,   1.0)
        volumeGradientOpacity.AddPoint(50,  0.5)
        volumeGradientOpacity.AddPoint(20, 0.0)

        # This class stores color data and can create color tables from a few color points. For this demo, we want the three cubes
        # to be of the colors red green and blue.
        colorFunc = vtk.vtkColorTransferFunction()
        colorFunc.AddRGBPoint(00, 1.0, 0.0, 0.0)
        colorFunc.AddRGBPoint(30, 0.0, 1.0, 0.0)
        colorFunc.AddRGBPoint(60, 0.0, 0.0, 1.0)

        # The preavius two classes stored properties. Because we want to apply these properties to the volume we want to render,
        # we have to store them in a class that stores volume prpoperties.
        volumeProperty = vtk.vtkVolumeProperty()
        volumeProperty.SetColor(colorFunc)
        volumeProperty.SetScalarOpacity(alphaChannelFunc)
        volumeProperty.SetGradientOpacity(volumeGradientOpacity)
        volumeProperty.SetInterpolationTypeToLinear()
        volumeProperty.ShadeOff()
        volumeProperty.SetAmbient(0.1)
        volumeProperty.SetDiffuse(0.6)
        volumeProperty.SetSpecular(0.2)

        # This class describes how the volume is rendered (through ray tracing).
        compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
        # We can finally create our volume. We also have to specify the data for it, as well as how the data will be rendered.
        volumeMapper = vtk.vtkVolumeRayCastMapper()
        volumeMapper.SetVolumeRayCastFunction(compositeFunction)
        volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

        # The class vtkVolume is used to pair the preaviusly declared volume as well as the properties to be used when rendering that volume.
        volume = vtk.vtkVolume()
        volume.SetMapper(volumeMapper)
        volume.SetProperty(volumeProperty)

        # Text am Nullpunkt
        atext = vtk.vtkVectorText()
        atext.SetText("(0,0,0)")
        textMapper = vtk.vtkPolyDataMapper()
        textMapper.SetInputConnection(atext.GetOutputPort())
        textActor = vtk.vtkFollower()
        textActor.SetMapper(textMapper)
        textActor.SetScale(10, 10, 10)
        textActor.AddPosition(0, -0.1, 78)

        # Cube to give some orientation 
        # (from http://www.vtk.org/Wiki/VTK/Examples/Python/Widgets/OrientationMarkerWidget)

        axesActor = vtk.vtkAnnotatedCubeActor();
        axesActor.SetXPlusFaceText('N')
        axesActor.SetXMinusFaceText('S')
        axesActor.SetYMinusFaceText('W')
        axesActor.SetYPlusFaceText('E')
        axesActor.SetZMinusFaceText('D')
        axesActor.SetZPlusFaceText('U')
        axesActor.GetTextEdgesProperty().SetColor(1,1,0)
        axesActor.GetTextEdgesProperty().SetLineWidth(2)
        axesActor.GetCubeProperty().SetColor(0,0,1)

        # With almost everything else ready, its time to initialize the renderer and window, as well as creating a method for exiting the application
        renderer = vtk.vtkRenderer()
        renderWin = vtk.vtkRenderWindow()
        renderWin.AddRenderer(renderer)
        renderInteractor = vtk.vtkRenderWindowInteractor()
        renderInteractor.SetRenderWindow(renderWin)

        axes = vtk.vtkOrientationMarkerWidget()
        axes.SetOrientationMarker(axesActor)
        axes.SetInteractor(renderInteractor)
        axes.EnabledOn()
        axes.InteractiveOn()
        renderer.ResetCamera()

        # We add the volume to the renderer ...
        renderer.AddVolume(volume)
        # ... set background color to white ...
        renderer.SetBackground(0.7,0.7,0.7)
        # ... and set window size.
        renderWin.SetSize(400, 400)

        # Fuege Text am Nullpunkt hinzu:
        renderer.AddActor(textActor)
        
        # A simple function to be called when the user decides to quit the application.
        def exitCheck(obj, event):
            if obj.GetEventPending() != 0:
                obj.SetAbortRender(1)

        # Tell the application to use the function as an exit check.
        renderWin.AddObserver("AbortCheckEvent", exitCheck)

        renderInteractor.Initialize()
        # Because nothing will be rendered without any input, we order the first render manually before control is handed over to the main-loop.
        renderWin.Render()
        renderInteractor.Start()


if __name__=='__main__':
    data_matrix = dict_2_3d_array(make_3D_pseudo_data())

    data_matrix = scale_2_int_perc(data_matrix)
    vtkCube(data_matrix)
