# An example from scipy cookbook demonstrating the use of numpy arrys in vtk 
from derec.core import TestCaseData, TestCaseSetup, TestCase
import derec.derec_utils as du
import matplotlib
import matplotlib.transforms as transforms
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import vtk
import numpy as num
from collections import defaultdict
from gmtpy import GMT
from copy import copy


matplotlib.rcParams['font.size'] = 7

def gca_label(label_string, **kwargs):
    """
    Add target codes to top right corner in axes object.
    """
    ax = plt.gca()
    plt.text(1, 1, label_string,
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes,
                    **kwargs)


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

    def __init__(self, test_case=None, test_case_data=None):
        if test_case:
            self.test_case = test_case
            data_input = self.test_case

        elif test_case_data:
            data_input = self.test_case_data

        self.candidates = data_input.candidates
        self.references = data_input.references
        self.processed_candidates = data_input.processed_candidates
        self.processed_references = data_input.processed_references
        self.channel_map = data_input.test_case_setup.channel_map

        self.processed_candidates_lines = None
        self.processed_references_lines = None
        self.candidates_lines = None
        self.references_lines = None
        
        self.test_case_setup = data_input.test_case_setup
        self.misfit_setup = self.test_case_setup.misfit_setup
        self.targets = self.test_case_setup.targets
        self.sources = self.test_case_setup.sources
        self.reference_source = self.test_case_setup.reference_source

    def gmt_map(self, **kwargs):
        sources_lon_lat = set()
        for s in self.sources:
            sources_lon_lat.add((s.lon, s.lat))
        sources_lon_lat = zip(*list(sources_lon_lat))

        targets_lon_lat = set()
        for t in self.targets:
            targets_lon_lat.add((t.lon, t.lat))
        targets_lon_lat = zip(*list(targets_lon_lat))
        gmt_map(sources_lon_lat[1], sources_lon_lat[0],
                targets_lon_lat[1], targets_lon_lat[0],
                    **kwargs)

    @staticmethod
    def figure_dict(keys):
        figures = [plt.figure(i, facecolor='grey') for i in range(len(keys))]
        return dict(zip(keys, figures))

    def waveforms_plot(self):
        """
        plot waveforms. 
        One figure per stations. 
        Three subplots, one per channel.
        several graphs in each subplot. One graph represents one source depth
        some ideas taken from
        http://wiki.scipy.org/Cookbook/Matplotlib/MultilinePlots
        """
        
        stats = TestCase.targets_nsl_of(self.targets)
        num_stations = len(stats)

        # TODO check if static method can be
        fig_dict = self.figure_dict(stats)

        axes_dict = {}
        # overlapping axes objects <- google
        for target in self.targets:
            fig = fig_dict[target.codes[:3]]
            axes_dict[target] = fig.add_subplot(1,3,self.channel_map[target.codes[3]])

        if not self.processed_candidates_lines:
            self.processed_candidates_lines = TestCase.lines_dict(self.processed_candidates)

        for source in self.sources:
            for target in self.targets:
                axes_dict[target].add_line(self.processed_candidates_lines[source][target])

        px=140.
        y_shifts = dict(zip(self.sources, num.linspace(-px/2./72., px/2./72., 
            len(self.sources))))
        
        # vertically distribute graphs
        for s in self.sources:
            for t in self.targets:
                line = self.processed_candidates_lines[s][t]
                fig = fig_dict[t.codes[:3]]
                trans = transforms.ScaledTranslation(0, y_shifts[s], 
                                                     fig.dpi_scale_trans)

                line.set_transform(line.get_transform()+trans)

        for ax in axes_dict.values():
            ax.axes.get_yaxis().set_visible(False)
            ax.autoscale()

    def plot_marker_vs_distance(self):
        """
        Plot:
            x-axis: Distance
            y-axis: Marker tmin, and marker max
        """
        fig = plt.figure()

    def distance_sort_targets(self, ref_source, targets=[]):
        """
        return a list of targets sorted by distance to given event
        """
        targets = self.targets if not targets else targets
        return sorted(targets, key=lambda tr: tr.distance_to(ref_source))
    
    def plot_z_components(self, traces_dict, markers_dict, sources=[], targets=[]):
        """
        Plot vertical components of each station. One figure per source, one
        subplot per station.
        """
        sources = self.sources if not sources else sources
        targets = self.targets if not targets else targets

        for source in sources:
            i = 0
            fig, axs = plt.subplots(len(targets)/3, sharex=True)

            for target in [t for t in self.distance_sort_targets(source,targets=targets)\
                        if t.codes[3]=='Z']:
                    m = markers_dict[source][target]
                    c = traces_dict[source][target]
                    pt = c.pyrocko_trace()

                    axs[i].plot(pt.get_xdata(), pt.get_ydata())
                    axs[i].axvline(m.tmin, label='P')
                    axs[i].axvline(m.tmax, label='P')

                    gca_label('.'.join(target.codes[:3]))
                    axs[i].yaxis.set_major_locator(MaxNLocator(prune='lower'))

                    i+=1

            fig.subplots_adjust(hspace=0)
            plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
            fig.suptitle('Vertical (z) Components at z=%s'%source.depth)

    def stack_plot(self, sources=None, depths=None ):
        '''
        param_dict is something like {latitude:10, longitude:10}, defining the 
        area of source, that you would like to look at.

        TODO: marker_dict unused?
        '''
        #if param_dict and sources:
        #    sources = TestCase.get_sources_where(param_dict, sources)
        #
        #else:
        #    sources = self.sources
        cmap = plt.get_cmap('Paired')

        gs = gridspec.GridSpec(len(self.targets)/3,3)
        gs_dict= dict(zip(self.targets, gs))

        axes_dict = defaultdict()

        for source,t, pr_cand_line in\
            TestCase.iter_dict(self.get_processed_candidates_lines(
                         reduce=self.reference_source.time)):
            if depths and not source.depth in depths:
                continue

            ax = plt.subplot(gs_dict[t])
            pr_ref = self.processed_references[source][t]

            x_ref = pr_ref.pyrocko_trace().get_xdata()-\
                    self.reference_source.time
            y_ref = pr_ref.pyrocko_trace().get_ydata()
                
            gca_label(label_string='.'.join(t.codes), fontsize=8)

            pr_cand_line.set_label("%s m"%float(source.depth))
            pr_cand_line.set_color(cmap(self.scalez255(source.depth)))
            ax.add_line(pr_cand_line)

            p = ax.fill_between(x_ref, 0, y_ref, facecolor='grey', alpha=0.5)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            plt.tick_params(axis='both', which='major', labelsize=6)
            
            axes_dict[t] = ax

        plt.subplots_adjust(left=None,
                            bottom=None,
                            right=None,
                            top=None,
                            wspace=None,
                            hspace=0.5)

        plt.legend(loc=2, prop={'size':6})
        plt.gcf().suptitle('%s, %s'%(
            self.test_case_setup.test_parameter,
            self.test_case_setup.test_parameter_value))

        return axes_dict

    def get_candidate_line(self, source, target):
        if not self.candidates_lines:
            self.get_candidates_lines()
        return self.candidates_lines[source][target]
    
    def get_reference_line(self, source, target):
        if not self.references_lines:
            self.get_references_lines()
        return self.references_lines[source][target]

    def get_processed_candidate_line(self, source, target, **kwargs):
        if not self.processed_candidates_lines:
            self.get_processed_candidates_lines(**kwargs)
        return self.processed_candidates_lines[source][target]

    def get_processed_reference_line(self, source, target):
        if not self.processed_references_lines:
            self.get_processed_references_lines()
        return self.processed_references_lines[source][target]

    def get_processed_references_lines(self, sources=[]):
        sources=sources if sources else self.sources
        if not self.processed_references_lines:
            self.processed_references_lines=\
                            TestCase.lines_dict(self.processed_references)
        out_dict = dict((s, l) for s,l in\
                self.processed_references_lines.iteritems() if s in sources)
        return out_dict

    def get_processed_candidates_lines(self, **kwargs):
        if not self.processed_candidates_lines:
            self.processed_candidates_lines =\
                            TestCase.lines_dict(self.processed_candidates,
                                    **kwargs)
        return self.processed_candidates_lines
        
    def get_candidates_lines(self):
        if not self.candidates_lines:
            self.candidates_lines = TestCase.lines_dict(self.candidates)
        return self.candidates_lines

    def get_references_lines(self):
        if not self.references_lines:
            self.references_lines= TestCase.lines_dict(self.references)
        return self.references_lines
    
    def blinded_key(self, _key, ignore):
        """
        make a new key entry for dict, neglecting one key/value pair, defined by
        :param ignore: string
        """
        if ignore:
            key_copy = copy(_key)
            key_copy.pop(ignore)
            return key_copy
        else: 
            return key

    def blinded_compare(self, item, test_list, ignore):
        """
        Check if item is in list, by means of their __dict__ representation,
        neglecting the parameter defined by *ignore*. 
        :param ignore: string of key, to be ignored in comparison
        """
        item_dict = copy(item.__dict__)
        item_dict.pop(ignore)
        for item_i in test_list:
            item_i_dict = copy(item_i.__dict__)
            item_i_dict.pop(ignore)
            if item_i_dict==item:
                return True

        return False 

    def make_blinded_dict(self, sources, targets, the_dict, ignore=None):
        """
        make a dict with items if key sources and targets are in *sources* and
        *tagets*.
        """
        match_dict = defaultdict(dict) 
        for s,t,l in TestCase.iter_dict(the_dict):
            if self.blinded_compare(s, sources, ignore) and t in targets:
                if match_dict[self.blinded_key(s, ignore)] is not None:
                    print 'warning: blinded key already exists'
                    raise Exception
                else:
                    match_dict[self.blinded_key(s, ignore)][t] = l

    def process_compare(self, sources=None, targets=None):
        """
        4 lines plot (raw/processed reference/candidat) for given sources and
        targets.
        One plot per target.
        If omitted, use all sources and targets possible
        :param ignoe: when comparing reference sources and candidates' sources,
        ignore this parameter.
        """
        sources = self.sources if not sources else sources
        targets = self.targets if not targets else targets

        if not self.processed_references_lines:
            self.processed_references_lines=\
                            TestCase.lines_dict(self.processed_references)
        if not self.processed_candidates_lines:
            self.processed_candidates_lines =\
                            TestCase.lines_dict(self.processed_candidates)
        if not self.candidates_lines:
            self.candidates_lines = TestCase.lines_dict(self.candidates)
        if not self.references_lines:
            self.references_lines= TestCase.lines_dict(self.references)

        for line_dict in [self.processed_candidates_lines, 
                          self.processed_references_lines, 
                          self.candidates_lines, 
                          self.references_lines]:

            ignored_line_dict = self.make_blinded_dict(sources, targets,
                    line_dict, ignore=self.test_case_setup.test_parameter)

    def get_sources_where(self, param_dict):
        return TestCase.get_sources_where(param_dict, self.sources)

    def scalez255(self, z): 
        return (z-min(self.test_case_setup.depths))/max(\
                self.test_case_setup.depths)*255 


class Cube():
    def __init__(self, test_case):

        self.test_case=test_case 
        #self.test_case.numpy_it()
        #self.xdim = len(self.test_case.test_parameters[test_case.xkey])
        #self.ydim = len(self.test_case.test_parameters[test_case.ykey])
        #self.zdim = len(self.test_case.test_parameters[test_case.zkey])

        #data=test_case.num_array[3].reshape(self.xdim,
        #                                    self.ydim,
        #                                    self.zdim)
        #data = scale_2_int_perc(data)
        #self.vtkCube(data)

    def value_to_index(k, val):
        return self.dim_mapper[k].index(val)

    def add_cases(self, cases):
        self.test_tin.extend(cases)
        self.update_dimensions(cases)

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

