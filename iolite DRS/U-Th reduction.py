# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Type: DRS
#/ Name: U-Th reduction
#/ Authors: Zoe Moser
#/ Description: U-Th zircon reduction
#/ References: based on Guillong et al. 2016
#/ Version: 1.0
#/ Contact: moserz@eaps.ethz.ch

"""
Before functions are called a few additional objects are
added to the module:

data	an interface to iolite's C++ data. E.g. you can get
        existing time series data or selection groups with it
        as well as make new ones.

IoLog	an interface to iolite's logging facility. You can add
        messages with, e.g., IoLog.debug('My message')

drs	an interface to the PythonDRS C++ class in iolite from
        which some built-in features can be accessed, e.g.,
        baselineSubtract(group, channels, mask)

Qt imports can be done through 'iolite', e.g.
from iolite.QtGui import QLabel
"""

from iolite import QtGui
from iolite import QtCore
from iolite.Qt import Qt, QColor
from iolite.ui import CommonUIPyInterface as CUI
from iolite.ui import IolitePlotPyInterface as Plot
from iolite.ui import IolitePlotSettingsDialog as PlotSettings
from iolite.QtGui import QAction, QPen, QBrush, QVBoxLayout, QWidget, QHBoxLayout, QLabel
from iolite.types import Result

from time import sleep
import numpy as np
import pandas as pd
from functools import partial
from scipy.optimize import curve_fit
import time


'''
A menu to handle which Samples and Reference Materials (RMs) to use in the Abundance correction
'''
class CheckableMenu(QtGui.QMenu):
    itemsChanged = QtCore.Signal(list)

    def __init__(self, parent, combined_names):
        super().__init__(parent)
        self.selectedItems = []  # To track checked items

        # Create checkboxes for each name in the combined list
        for name in combined_names:
            a = QtGui.QWidgetAction(self)
            cb = QtGui.QCheckBox(name, self)
            cb.setStyleSheet('QCheckBox { padding-left: 5px; margin: 3px; }')
            a.setDefaultWidget(cb)
            self.addAction(a)
            cb.clicked.connect(partial(self.updateChannels, name))

        self.aboutToShow.connect(self.updateMenu)

    def updateMenu(self):
         for a in self.actions():
            cb = a.defaultWidget()
            cb.setChecked(False)
            try:
                if cb.text in self.selectedItems:
                    cb.setChecked(True)
            except Exception as e:
                print(e)
        
    def updateChannels(self, combined_names, b):
        if b:
            self.selectedItems = list(set(self.selectedItems + [combined_names]))
        else:
            self.selectedItems = list(filter(lambda rm: rm != combined_names, self.selectedItems))

        self.itemsChanged.emit(self.selectedItems)

'''
A plot to handle plotting of CaPO correction
If this is the first time the DRS script has been run (e.g.
when you first click on the DRS in iolite) the plot will be
initiated.
It also adds a Settings menu item to the context menu
'''
def showSettings():
    d = PlotSettings(PLOT)
    d.exec_()

try:
    PLOT
except:
    PLOT = Plot()
    PLOT.setAttribute(Qt.WA_DeleteOnClose)
    PLOT.setFixedSize(600,300)
    
    # Add annotation to show fit parameters
    ann = PLOT.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
    ann.visible = False

    def showSettings():
        d = PlotSettings(PLOT)
        d.exec_()

    settingsAction = QAction(PLOT.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT.contextMenu().addAction(settingsAction)

try:
    PLOT_Zr2O3
except:
    PLOT_Zr2O3 = Plot()
    PLOT_Zr2O3.setAttribute(Qt.WA_DeleteOnClose)
    PLOT_Zr2O3.setFixedSize(600,300)
    
    # Add annotation to show fit parameters
    ann2 = PLOT_Zr2O3.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
    ann2.visible = False

    def showSettings():
        d = PlotSettings(PLOT_Zr2O3)
        d.exec_()

    settingsAction = QAction(PLOT_Zr2O3.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT_Zr2O3.contextMenu().addAction(settingsAction)

try:
    PLOT_Mon228
except:
    PLOT_Mon228 = Plot()
    PLOT_Mon228.setAttribute(Qt.WA_DeleteOnClose)
    PLOT_Mon228.setFixedSize(600,300)
    
    # Add annotation to show fit parameters
    ann3 = PLOT_Mon228.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
    ann3.visible = False

    def showSettings():
        d = PlotSettings(PLOT_Mon228)
        d.exec_()

    settingsAction = QAction(PLOT_Mon228.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT_Mon228.contextMenu().addAction(settingsAction)

try:
    PLOT_Mon230
except:
    PLOT_Mon230 = Plot()
    PLOT_Mon230.setAttribute(Qt.WA_DeleteOnClose)
    PLOT_Mon230.setFixedSize(600,300)
    
    # Add annotation to show fit parameters
    ann4 = PLOT_Mon230.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
    ann4.visible = False

    def showSettings():
        d = PlotSettings(PLOT_Mon230)
        d.exec_()

    settingsAction = QAction(PLOT_Mon230.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT_Mon230.contextMenu().addAction(settingsAction)

try:
    PLOT_RSF
except:
    PLOT_RSF = Plot()
    PLOT_RSF.setAttribute(Qt.WA_DeleteOnClose)
    PLOT_RSF.setFixedSize(600,300)
    
    # Add annotation to show fit parameters
    ann5 = PLOT_RSF.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
    ann5.visible = False

    def showSettings():
        d = PlotSettings(PLOT_RSF)
        d.exec_()

    settingsAction = QAction(PLOT_RSF.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT_RSF.contextMenu().addAction(settingsAction)

try:
    PLOT_Evolution
except:
    PLOT_Evolution = Plot()
    PLOT_Evolution.setAttribute(Qt.WA_DeleteOnClose)
    PLOT_Evolution.setFixedSize(600,300)
    
    # Add annotation to show fit parameters
    ann6 = PLOT_Evolution.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
    ann6.visible = False

    def showSettings():
        d = PlotSettings(PLOT_Evolution)
        d.exec_()

    settingsAction = QAction(PLOT_Evolution.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT_Evolution.contextMenu().addAction(settingsAction)

try:
        PLOT_ppm
except:
        PLOT_ppm = Plot()
        PLOT_ppm.setAttribute(Qt.WA_DeleteOnClose)
        PLOT_ppm.setFixedSize(600,300)
        PLOT_ppm.hide()
        
        # Add annotation to show fit parameters
        ann7 = PLOT_ppm.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
        ann7.visible = False

        def showSettings():
                d = PlotSettings(PLOT_ppm)
                d.exec_()

        settingsAction = QAction(PLOT_ppm.contextMenu())
        settingsAction.setText('Settings')
        settingsAction.triggered.connect(showSettings)
        PLOT_ppm.contextMenu().addAction(settingsAction)

try:
        PLOT_Mon234
except:
        PLOT_Mon234 = Plot()
        PLOT_Mon234.setAttribute(Qt.WA_DeleteOnClose)
        PLOT_Mon234.setFixedSize(600,300)
        PLOT_Mon234.hide()
        
        # Add annotation to show fit parameters
        ann8 = PLOT_Mon234.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
        ann8.visible = False

        def showSettings():
                d = PlotSettings(PLOT_Mon234)
                d.exec_()

        settingsAction = QAction(PLOT_Mon234.contextMenu())
        settingsAction.setText('Settings')
        settingsAction.triggered.connect(showSettings)
        PLOT_Mon234.contextMenu().addAction(settingsAction)

try:
        PLOT_234_238
except:
        PLOT_234_238 = Plot()
        PLOT_234_238.setAttribute(Qt.WA_DeleteOnClose)
        PLOT_234_238.setFixedSize(600,300)
        PLOT_234_238.hide()
        
        # Add annotation to show fit parameters
        ann9 = PLOT_234_238.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
        ann9.visible = False

        def showSettings():
                d = PlotSettings(PLOT_234_238)
                d.exec_()

        settingsAction = QAction(PLOT_234_238.contextMenu())
        settingsAction.setText('Settings')
        settingsAction.triggered.connect(showSettings)
        PLOT_234_238.contextMenu().addAction(settingsAction)

PLOT_COLORS = [
    QColor(239, 71, 111),   # Carmine Red
    QColor(255, 209, 102),  # Orange Yellow
    QColor(6, 214, 160),    # Green
    QColor(17, 138, 178),   # Blue
    QColor(7, 59, 76),      # Midnight Green
    QColor(204, 245, 172),  # Light Green
    QColor(128, 138, 159),  # Roman Silver
    QColor(61, 43, 86),     # Russian Purple
    QColor(247, 169, 168),  # Pastel Pink
    QColor(207, 212, 197),  # Bone
]


def runDRS():
        drs.message("Starting baseline subtract DRS...")
        drs.progress(0)

##%     Get settings
        settings = drs.settings()
        print(settings)

        indexChannel = data.timeSeries(settings["IndexChannel"])
        maskChannel = data.timeSeries(settings["MaskChannel"])
        cutoff = settings["MaskCutoff"]
        trim = settings["MaskTrim"]
        monazite_group  = data.selectionGroup(settings["Monazite"])
        #RSF_group = data.selectionGroup(settings["RSF_on"])
        degree_of_MB_fit = settings["MB_fit"]
        MB_groups = settings["MB_groups"]
        AS228_groups = settings["Abundance228Items"]
        Result_groups = settings["Result_groups"]
        degree_of_Ab230_fit = settings["Ab230_fit"]
        degree_of_Ab228_fit = settings["Ab228_fit"]
        Zr2O3corr_logical = settings['Zr2O3corr']
        RSF_groups_plot = settings["RSF_plot"]
        RSF_outlier_logical = settings["RSFOutlier"]
        RSF_how = settings["RSF_how"]
        PlotLc = settings['plotLc1']
        RSF_fit = settings["RSF_fit"]
        InfCorr228_fit = settings['Corr228_fit']
        InfCorr_outlier_logical = settings["InfCorrOutlier"]
        how235or238 = settings["U235orU238"]
        how238_232 = settings["How238_232"]
        how230_232 = settings["How230_232"]
        fillLabels = settings["fill_Labels"]
        calc_ppm = settings["calculate_ppm"]
        ppm_fit = settings["ppm_fit"]
        ppm_on = settings["ppm_calculated_on"]
        U234secDis = settings["calculate_U234_U238"]
        U234secDis_group = settings["U234-U238_groups"]
        mask234 = settings["mask_U234_U238"]
        how234_238 = settings["How234_238"]
        corr228_int_ext = settings["Corr228_intext"]
        external_correction = settings["Corr228_external"]
        external_correction_SE = settings["Corr228_external_SE"]

##%     Constants used
        lambda_Th230 = 0.0000091577
        lambda_Th232 = 0.000000000049475
        lambda_U235 = 0.00000000098571
        lambda_U234 = 2.8226052879421155e-06
        lambda_U238 = 0.000000000155125
        scaling_factor = (17.521+1.25+0.158+0.0138)/26.28
        if RSF_how == 'based on RM U & Th concentration':
                U_RM_ppm_RSF = data.referenceMaterialData(settings["RSF_on"])["U"].value()
                Th_RM_ppm_RSF = data.referenceMaterialData(settings["RSF_on"])["Th"].value()
        U238_U235 = 137.818
        dwell_Th230_s = data.timeSeries('Th230').property('Dwell Time (ms)')/1000
        print(f'Dwell time Th230: {dwell_Th230_s*1000} ms?')


##%     Interp onto index time and baseline subtract
        drs.message("Interpolating onto index time and baseline subtracting...")

        allInputChannels = data.timeSeriesList(data.Input)
        blGrp = None

        if len(data.selectionGroupList(data.Baseline)) != 1:
                IoLog.error("There should be exactly one baseline group.")
                drs.message("DRS did not finish. Please check Messages")
                drs.progress(100)
                drs.finished()
                return
        else:
                blGrp = data.selectionGroupList(data.Baseline)[0]
        
        baselinemask = drs.createMaskFromCutoff(maskChannel,cutoff,trim)
        drs.baselineSubtract(blGrp, allInputChannels, baselinemask, 15, 35)

##%     Calculate ratios
        Th230_U238_measured = data.timeSeries('Th230_CPS').data()/data.timeSeries('U238_CPS').data()
        Th230_U238_timeSeries = data.createTimeSeries('Th230/U238', data.Intermediate, indexChannel.time(), Th230_U238_measured)

        Th232_U238_measured = data.timeSeries('Th232_CPS').data()/data.timeSeries('U238_CPS').data()
        Th232_U238_timeSeries = data.createTimeSeries('Th232/U238', data.Intermediate, indexChannel.time(), Th232_U238_measured)

        Th230_Th232_measured = data.timeSeries('Th230_CPS').data()/data.timeSeries('Th232_CPS').data()
        Th230_Th232_timeSeries = data.createTimeSeries('Th230/Th232', data.Intermediate, indexChannel.time(), Th230_Th232_measured)

        U238_U235_measured = data.timeSeries('U238_CPS').data()/data.timeSeries('U235_CPS').data()
        U238_U235_timeSeries = data.createTimeSeries('U238/U235', data.Intermediate, indexChannel.time(), U238_U235_measured)

##%     Abundance sensitivity correction of 232 on 230
        drs.message("Abundance sensitivity correction of 232 on 230...")
        drs.progress(10)

        PLOT_Mon230.clearGraphs()

        def Abundance232_on_230(sel):
                result = Result()
                result = ((data.result(sel, data.timeSeries('Th230_CPS')).value() - ((lambda_U238/lambda_Th230* (data.result(sel, data.timeSeries('U238_CPS')).value()))))/(data.result(sel, data.timeSeries('Th232_CPS')).value()))*1e+6
                return result
 
        def Abundance232_on_230_2SE(sel):
                result = Result()
                Th230_CPS = data.result(sel, data.timeSeries('Th230_CPS')).value()
                Th232_CPS = data.result(sel, data.timeSeries('Th232_CPS')).value()
                U238_CPS = data.result(sel, data.timeSeries('U238_CPS')).value()
                Th230_CPSSE = data.result(sel, data.timeSeries('Th230_CPS')).uncertaintyAs2SE()
                Th232_CPSSE = data.result(sel, data.timeSeries('Th232_CPS')).uncertaintyAs2SE()
                U238_CPSSE = data.result(sel, data.timeSeries('U238_CPS')).uncertaintyAs2SE()
                result = np.sqrt((1e+6/Th232_CPS*Th230_CPSSE)**2+(-1e+6/Th232_CPS*lambda_U238/lambda_Th230*U238_CPSSE)**2+((1e+6*Th232_CPSSE*(U238_CPS*lambda_U238/lambda_Th230-Th230_CPS))/(Th232_CPS**2))**2)
                return result       

        data.registerAssociatedResult('mass230_ppm_Th232_seq',Abundance232_on_230)  
        data.registerAssociatedResult('mass230_ppm_Th232_seq_2SE',Abundance232_on_230_2SE) 

        Th230_Th232_ppm_Monazite = []
        Th230_Th232_ppm_Monazite_2SE = []  
        start_time = []      

        sg = monazite_group
        for sel in sg.selections():
                Th230_Th232_ppm_timeseries_Monazite = data.associatedResult(sel,'mass230_ppm_Th232_seq').value()
                Th230_Th232_ppm_timeseries_Monazite_2SE = data.associatedResult(sel,'mass230_ppm_Th232_seq_2SE').value()
                start_time_sel = indexChannel.timeForSelection(sel)
                Th230_Th232_ppm_Monazite.append(Th230_Th232_ppm_timeseries_Monazite)
                Th230_Th232_ppm_Monazite_2SE.append(Th230_Th232_ppm_timeseries_Monazite_2SE)
                start_time.append(start_time_sel[0])
        
        start_times_array = np.array(start_time-indexChannel.time()[0])
        meas_array = np.array(Th230_Th232_ppm_Monazite)
        
        fit_x = np.linspace(
            start_times_array.min(),
            start_times_array.max(),
            50
        )

        if degree_of_Ab230_fit == 'Constant':
                def fitConst(x, a, b):
                        return a + b*0

                params, cov = curve_fit(fitConst, start_times_array, meas_array,  ftol=1e-5)
                params = np.mean(meas_array)
                cov = np.mean(meas_array)/np.sqrt(len(meas_array))
                fit_y = params*np.ones(len(fit_x))
                
                mass230_ppm_Th232_seq_fit = params+0*(indexChannel.time()-min(indexChannel.time()))
                mass230_ppm_Th232_seq_fit_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit)
                
                mass230_ppm_Th232_seq_fit_uncer = cov * np.ones(len(indexChannel.time()))
                mass230_ppm_Th232_seq_fit_uncer_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit_uncer', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit_uncer
                )        
        
        if degree_of_Ab230_fit == 'Linear':
                def fitLin(x, a, b):
                        return a + b*x

                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                fit_y = params[1]*fit_x + params[0]
                
                mass230_ppm_Th232_seq_fit = params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))
                mass230_ppm_Th232_seq_fit_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit)
                fit_uncertainties = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                            2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                mass230_ppm_Th232_seq_fit_uncer_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
        

        if degree_of_Ab230_fit == 'Logistic':
                def logistic_func(x, a, b, c, d):
                        return a + (b - a) / (1 + np.exp(c * (x - d)))

                # Smart initial guess without scaling x
                a0 = meas_array.min()                  # lower asymptote
                b0 = meas_array.max()                  # upper asymptote
                c0 = 10 / (start_times_array.max() - start_times_array.min())  # slope guess ~1 / range of x
                d0 = np.quantile(start_times_array, 0.75)

                p0 = [a0, b0, c0, d0]
                bounds = (
                        [0, 0, -10, start_times_array.min()],       # lower bounds
                        [10, 10, 10, start_times_array.max()]       # upper bounds
                )

                try:
                        params, cov = curve_fit(logistic_func, start_times_array, meas_array, p0=p0, bounds=bounds)
                        
                        fit_y = logistic_func(fit_x, *params)
                        
                        mass230_ppm_Th232_seq_fit = logistic_func(indexChannel.time()-min(indexChannel.time()), *params)
                        mass230_ppm_Th232_seq_fit_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit)
     
                        residuals = meas_array - logistic_func(start_times_array, *params)
                        residual_std = np.std(residuals)

                        mass230_ppm_Th232_seq_fit_uncer = np.full_like(indexChannel.time(), residual_std)
                        mass230_ppm_Th232_seq_fit_uncer_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit_uncer', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit_uncer)

                except RuntimeError as e:
                        if 'Optimal parameters not found: Number of calls' in str(e):
                                IoLog.warning('Could not find fit to data with Logistic model. Switching to Linear model.')
                                settings["MB_fit"] = 'Linear'

                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                mass230_ppm_Th232_seq_fit = params[0] + params[1] * indexChannel.time()
                                mass230_ppm_Th232_seq_fit_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit)

                                fit_uncertainties = np.sqrt(
                                        cov[0, 0] + 
                                        (indexChannel.time())**2 * cov[1, 1] +
                                        2 * indexChannel.time() * cov[0, 1]
                                )
                                mass230_ppm_Th232_seq_fit_uncer_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
                        else:
                                raise

        if degree_of_Ab230_fit == 'Polynomial':
                def fitPol1(x, a, b, c):
                        return a + b*x + c*x*x

                params, cov = curve_fit(fitPol1, start_times_array, meas_array, ftol=1e-5)

                fit_y = params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                mass230_ppm_Th232_seq_fit = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2)
                mass230_ppm_Th232_seq_fit_timeseries = data.createTimeSeries('mass230_ppm_Th232_seq_fit', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit)


                x_vals = indexChannel.time() - min(indexChannel.time())
                mass230_ppm_Th232_seq_fit_uncer = np.sqrt(
                        cov[0, 0]
                        + x_vals ** 2 * cov[1, 1]
                        + x_vals ** 4 * cov[2, 2]
                        + 2 * x_vals * cov[0, 1]
                        + 2 * x_vals ** 2 * cov[0, 2]
                        + 2 * x_vals ** 3 * cov[1, 2]
                )
                mass230_ppm_Th232_seq_fit_uncer_timeseries = data.createTimeSeries(
                        f'mass230_ppm_Th232_seq_fit_uncer', data.Intermediate, indexChannel.time(), mass230_ppm_Th232_seq_fit_uncer
                )

        g = PLOT_Mon230.addGraph()
        g.setName('Monazite')
        g.setLineStyle('lsNone')
        g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[2 % len(PLOT_COLORS)])
        g.setData(np.array(start_times_array),np.array(Th230_Th232_ppm_Monazite))    

        for x, y, y_err in zip(start_times_array, Th230_Th232_ppm_Monazite,Th230_Th232_ppm_Monazite_2SE):
                v_error = PLOT_Mon230.addGraph()
                vERR_pen = QPen(QColor('dark grey'))
                vERR_pen.setStyle(Qt.SolidLine)
                v_error.setPen(vERR_pen)
                v_error.setData([x, x], [y - y_err, y + y_err])
                v_error.removeFromLegend()  # Remove vertical error bar from legend

        g = PLOT_Mon230.addGraph()
        g.setName("Fit")
        fit_pen = QPen(QColor('black'))
        fit_pen.setStyle(Qt.DashLine)
        fit_pen.setWidth(2)
        g.pen = fit_pen
        g.setData(np.array(fit_x), np.array(fit_y))

        if degree_of_Ab230_fit == 'Linear':
                ann4.visible = True
                ann4.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                ({params[0]:.2f}±{cov[0,0]:.1e})+({params[1]:.2e}±{cov[1,1]:.1e})*t</p>'''
                
        if degree_of_Ab230_fit == 'Constant':
                ann4.visible = True
                ann4.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params:.2f}±{cov:.2f}</p>'''        
        
        if degree_of_Ab230_fit == 'Logistic':
                ann4.visible = True
                if len(params) == 4:
                        param_uncertainties = np.sqrt(np.diag(cov))
                        ann4.text = f'''
                        <p style="color:black;">
                        <b>Fit Equation:</b><br>
                        {params[0]:.2f} + ({params[1]:.2f} - {params[0]:.2f}) / 
                        (1 + exp({params[2]:.4f} × (t - {params[3]:.2f})))</p>'''
                else:
                        # fallback for linear case
                        ann4.text = f'''
                        <p style="color:black;">
                        <b>Fit Parameters (linear):</b><br>
                        ({params[0]:.2f} ± {cov[0,0]:.1e}) + ({params[1]:.2e} ± {cov[1,1]:.1e}) × t</p>'''

        if degree_of_Ab230_fit == 'Polynomial':
                ann4.visible = True
                param_uncertainties = np.sqrt(np.diag(cov))
                ann4.text = f'''
                <p style="color:black;">
                <b>Fit Equation:</b><br>
                {params[0]:.2f} + 
                {params[1]:.2e} × t + 
                {params[2]:.2e} × t²</p>'''         

        PLOT_Mon230.left().label = 'mass230 CPS ppm Th232 CPS'
        PLOT_Mon230.bottom().label = 'Time since start of session (s)'
        PLOT_Mon230.setLegendVisible(True)
        PLOT_Mon230.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)           
        PLOT_Mon230.setToolsVisible(False)
        PLOT_Mon230.rescaleAxes()
        PLOT_Mon230.replot()

##%     Abundance sensitivity correction of 228 on 230
        drs.message("Abundance sensitivity correction of 228 on 230...")
        drs.progress(30)
        PLOT_Mon228.clearGraphs()
        PLOT_Zr2O3.clearGraphs()
        
        if Zr2O3corr_logical == True:
                if corr228_int_ext == 'Add number from external file':

                        InfCorr228_fitted = external_correction+0*(indexChannel.time()-min(indexChannel.time()))
                        InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)
                        
                        InfCorr228_fitted_uncer = external_correction_SE * np.ones(len(indexChannel.time()))
                        InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), InfCorr228_fitted_uncer)        

                        median_AbundCorr228 = external_correction
                        median_AbundCorr228_SE = external_correction_SE

                        PLOT_Mon228.hide()
                        PLOT_Zr2O3.hide() 
#**                
                if corr228_int_ext == 'Th230 directly measured in Zirconblank':
                        PLOT_Mon228.hide()
                        meas_main = []
                        start_times_all = []

                        for i, gr in enumerate(AS228_groups):
                                meas = []
                                meas_2SE = []
                                start_times = []
                                sg = data.selectionGroup(gr)
                                for sel in sg.selections():
                                        ratio = data.result(sel, data.timeSeries('Th230_CPS')).value()
                                        ratio_2SE = data.result(sel, data.timeSeries('Th230_CPS')).uncertaintyAs2SE()
                                        start_time = indexChannel.timeForSelection(sel)
                                        if not np.isnan(ratio):
                                                meas.append(ratio)
                                                meas_2SE.append(ratio_2SE)
                                                meas_main.append(ratio)
                                                start_times.append(start_time[0])
                                                start_times_all.append(start_time[0]) 
                                meas = np.array(meas)
                                meas_2SE = np.array(meas_2SE)
                                start_times = np.array(start_times)
                                mean_meas = np.average(meas_main)
                                std_meas = np.std(meas_main)
                                if InfCorr_outlier_logical == True:
                                        mask = (meas >= mean_meas - 2 * std_meas) & (meas <= mean_meas + 2 * std_meas)
                                else: 
                                        mask = np.ones_like(meas, dtype=bool)
                                filtered_meas = meas[mask]
                                filtered_meas_2SE = meas_2SE[mask]
                                filtered_start_times = start_times[mask]

                                filtered_start_times = filtered_start_times-indexChannel.time()[0]
                                g = PLOT_Zr2O3.addGraph()
                                g.setName(gr)
                                g.setLineStyle('lsNone')
                                g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[i % len(PLOT_COLORS)])
                                g.setData(np.array(filtered_start_times), np.array(filtered_meas))

                                for x, y, y_err in zip(filtered_start_times, filtered_meas, filtered_meas_2SE):
                                        v_error = PLOT_Zr2O3.addGraph()
                                        vERR_pen = QPen(QColor('light grey'))
                                        vERR_pen.setStyle(Qt.SolidLine)
                                        v_error.setPen(vERR_pen)
                                        v_error.setData([x, x], [y - y_err, y + y_err])
                                        v_error.removeFromLegend()  # Remove vertical error bar from legend

                        start_times_all = start_times_all-indexChannel.time()[0]
                        meas_array = np.array(meas_main)
                        start_times_array = np.array(start_times_all)

                        if InfCorr_outlier_logical == True:
                                mask = (meas_array >= np.mean(meas_array) - 2 * np.std(meas_array)) & (meas_array <= np.mean(meas_array) + 2 * np.std(meas_array))
                        else: 
                                mask = np.ones_like(meas_array, dtype=bool)
                        meas_array = meas_array[mask]
                        start_times_array = start_times_array[mask]

                        fit_x = np.linspace(start_times_array.min(),start_times_array.max(),50)

                        if InfCorr228_fit == 'Constant':
                                def fitConst(x, a, b):
                                        return a + b*0

                                params, cov = curve_fit(fitConst, start_times_array, meas_array,  ftol=1e-5)
                                params = np.mean(meas_array)
                                cov = np.std(meas_array)/np.sqrt(len(meas_array))
                                fit_y = params*np.ones(len(fit_x))
                                
                                InfCorr228_fitted = params+0*(indexChannel.time()-min(indexChannel.time()))
                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)
                                
                                InfCorr228_fitted_uncer = cov * np.ones(len(indexChannel.time()))
                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), InfCorr228_fitted_uncer)        
                        
                        if InfCorr228_fit == 'Linear':
                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                
                                InfCorr228_fitted = params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))
                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)
                                fit_uncertainties = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                                        2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
                        
                        if InfCorr228_fit == 'Logistic':
                                def logistic_func(x, a, b, c, d):
                                        return a + (b - a) / (1 + np.exp(c * (x - d)))

                                # Smart initial guess without scaling x
                                a0 = meas_array.min()                  # lower asymptote
                                b0 = meas_array.max()                  # upper asymptote
                                c0 = 10 / (start_times_array.max() - start_times_array.min())  # slope guess ~1 / range of x
                                d0 = np.quantile(start_times_array, 0.5)

                                p0 = [a0, b0, c0, d0]
                                bounds = (
                                        [-10, -10, -10, start_times_array.min()],       # lower bounds
                                        [30, 30, 10, start_times_array.max()]       # upper bounds
                                )

                                try:
                                        params, cov = curve_fit(logistic_func, start_times_array, meas_array, p0=p0, bounds=bounds)
                                        
                                        fit_y = logistic_func(fit_x, *params)
                                        
                                        InfCorr228_fitted = logistic_func(indexChannel.time()-min(indexChannel.time()), *params)
                                        InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)

                                        residuals = meas_array - logistic_func(start_times_array, *params)
                                        residual_std = np.std(residuals)

                                        InfCorr228_fitted_uncer = np.full_like(indexChannel.time(), residual_std)
                                        InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), InfCorr228_fitted_uncer)

                                except RuntimeError as e:
                                        if 'Optimal parameters not found: Number of calls' in str(e):
                                                IoLog.warning('Could not find fit to data with Logistic model. Switching to Linear model.')
                                                settings["MB_fit"] = 'Linear'

                                                def fitLin(x, a, b):
                                                        return a + b*x

                                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                                fit_y = params[1]*fit_x + params[0]
                                                InfCorr228_fitted = params[0] + params[1] * indexChannel.time()
                                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)

                                                fit_uncertainties = np.sqrt(
                                                        cov[0, 0] + 
                                                        (indexChannel.time())**2 * cov[1, 1] +
                                                        2 * indexChannel.time() * cov[0, 1]
                                                )
                                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
                                        else:
                                                raise

                        if InfCorr228_fit == 'Polynomial':
                                def fitPol1(x, a, b, c):
                                        return a + b*x + c*x*x

                                params, cov = curve_fit(fitPol1, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                                InfCorr228_fitted = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2)
                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)


                                x_vals = indexChannel.time() - min(indexChannel.time())
                                InfCorr228_fitted_uncer = np.sqrt(
                                        cov[0, 0]
                                        + x_vals ** 2 * cov[1, 1]
                                        + x_vals ** 4 * cov[2, 2]
                                        + 2 * x_vals * cov[0, 1]
                                        + 2 * x_vals ** 2 * cov[0, 2]
                                        + 2 * x_vals ** 3 * cov[1, 2]
                                )
                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries(
                                        f'InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), InfCorr228_fitted_uncer
                                )
                        
                        median_AbundCorr228 = np.median(meas_array)
                        median_AbundCorr228_SE = (np.percentile(meas_array, 75) - np.percentile(meas_array, 25)) / np.sqrt(len(meas_array))
                        median_AbundCorr228_array = np.linspace(np.median(meas_array),np.median(meas_array),len(start_times_array))
                        
                        PLOT_Zr2O3.show()
                        g = PLOT_Zr2O3.addGraph()
                        g.setName("Fit")
                        fit_pen = QPen(QColor('black'))
                        fit_pen.setStyle(Qt.DashLine)
                        fit_pen.setWidth(2)
                        g.pen = fit_pen
                        g.setData(np.array(fit_x), np.array(fit_y))

                        if InfCorr228_fit == 'Linear':
                                ann2.visible = True
                                ann2.text = f'''
                                <p style="color:black;">
                                <b>Fit Parameters:</b><br>
                                ({params[0]:.2f}±{cov[0,0]:.1e})+({params[1]:.2e}±{cov[1,1]:.1e})*t</p>'''
                                
                        if InfCorr228_fit == 'Constant':
                                ann2.visible = True
                                ann2.text = f'''
                                <p style="color:black;">
                                <b>Fit Parameters:</b><br>
                                {params:.2f}±{cov:.2f}</p>'''        

                        if InfCorr228_fit == 'Logistic':
                                ann2.visible = True
                                if len(params) == 4:
                                        param_uncertainties = np.sqrt(np.diag(cov))
                                        ann2.text = f'''
                                        <p style="color:black;">
                                        <b>Fit Equation:</b><br>
                                        {params[0]:.2f} + ({params[1]:.2f} - {params[0]:.2f}) / 
                                        (1 + exp({params[2]:.4f} × (t - {params[3]:.2f})))</p>'''
                                else:
                                        # fallback for linear case
                                        ann2.text = f'''
                                        <p style="color:black;">
                                        <b>Fit Parameters (linear):</b><br>
                                        ({params[0]:.2f} ± {cov[0,0]:.1e}) + ({params[1]:.2e} ± {cov[1,1]:.1e}) × t</p>'''

                        if InfCorr228_fit == 'Polynomial':
                                ann2.visible = True
                                param_uncertainties = np.sqrt(np.diag(cov))
                                ann2.text = f'''
                                <p style="color:black;">
                                <b>Fit Equation:</b><br>
                                {params[0]:.2f} + 
                                {params[1]:.2e} × t + 
                                {params[2]:.2e} × t²</p>'''    

                        PLOT_Zr2O3.left().label = '230Th CPS in Zirconblank'
                        PLOT_Zr2O3.bottom().label = 'Time since start of session (s)'
                        PLOT_Zr2O3.setLegendVisible(True)
                        PLOT_Zr2O3.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
                        PLOT_Zr2O3.setToolsVisible(True)
                        PLOT_Zr2O3.rescaleAxes()
                        PLOT_Zr2O3.replot()                        
#**
                if corr228_int_ext == 'Zr90 measurements in this file':
                        def Abundance228_on_230(sel):
                                result = Result()
                                result = (data.result(sel, data.timeSeries('Zr90_CPS')).value())/(data.result(sel, data.timeSeries('Th232_CPS')).value())*1e+6
                                return result
                
                        def Abundance228_on_230_2SE(sel):
                                result = Result()
                                Zr90_CPS = data.result(sel, data.timeSeries('Zr90_CPS')).value()
                                Th232_CPS = data.result(sel, data.timeSeries('Th232_CPS')).value()
                                Zr90_CPSSE = data.result(sel, data.timeSeries('Zr90_CPS')).uncertaintyAs2SE()
                                Th232_CPSSE = data.result(sel, data.timeSeries('Th232_CPS')).uncertaintyAs2SE()
                                result = np.sqrt((1e+6 / Th232_CPS * Zr90_CPSSE)**2 + (-1e+6 * Zr90_CPS / (Th232_CPS**2) * Th232_CPSSE)**2)
                                return result       

                        data.registerAssociatedResult('mass228_ppm_Th232',Abundance228_on_230)  
                        data.registerAssociatedResult('mass228_ppm_Th232_2SE',Abundance228_on_230_2SE)   

                        Zr228_Th232_ppm_Monazite = []
                        Zr228_Th232_ppm_Monazite_2SE = [] 
                        start_time = []    

                        sg = monazite_group
                        for sel in sg.selections():
                                Zr228_Th232__ppm_timeseries_monazite = data.associatedResult(sel,'mass228_ppm_Th232').value()
                                Zr228_Th232__ppm_timeseries_monazite_2SE = data.associatedResult(sel,'mass228_ppm_Th232_2SE').value()
                                start_time_sel = indexChannel.timeForSelection(sel)
                                Zr228_Th232_ppm_Monazite.append(Zr228_Th232__ppm_timeseries_monazite)
                                Zr228_Th232_ppm_Monazite_2SE.append(Zr228_Th232__ppm_timeseries_monazite_2SE)                
                                start_time.append(start_time_sel[0])          

                        start_times_array = np.array(start_time-indexChannel.time()[0])
                        meas_array = np.array(Zr228_Th232_ppm_Monazite)
                        
                        fit_x = np.linspace(start_times_array.min(),start_times_array.max(),50)

                        if degree_of_Ab228_fit == 'Constant':
                                def fitConst(x, a, b):
                                        return a + b*0

                                params, cov = curve_fit(fitConst, start_times_array, meas_array,  ftol=1e-5)
                                params = np.mean(meas_array)
                                cov = np.std(meas_array)/np.sqrt(len(meas_array))
                                fit_y = params*np.ones(len(fit_x))
                                
                                mass228_ppm_Th232_fit = params+0*(indexChannel.time()-min(indexChannel.time()))
                                mass228_ppm_Th232_fit_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit)
                                
                                mass228_ppm_Th232_fit_uncer = cov * np.ones(len(indexChannel.time()))
                                mass228_ppm_Th232_fit_uncer_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit_uncer', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit_uncer)        
                        
                        if degree_of_Ab228_fit == 'Linear':
                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                
                                mass228_ppm_Th232_fit = params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))
                                mass228_ppm_Th232_fit_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit)
                                fit_uncertainties = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                                        2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                                mass228_ppm_Th232_fit_uncer_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)

                        if degree_of_Ab228_fit == 'Logistic':
                                def logistic_func(x, a, b, c, d):
                                        return a + (b - a) / (1 + np.exp(c * (x - d)))

                                # Smart initial guess without scaling x
                                a0 = meas_array.min()                  # lower asymptote
                                b0 = meas_array.max()                  # upper asymptote
                                c0 = 10 / (start_times_array.max() - start_times_array.min())  # slope guess ~1 / range of x
                                d0 = np.quantile(start_times_array, 0.5)

                                p0 = [a0, b0, c0, d0]
                                bounds = (
                                        [0, 0, -10, start_times_array.min()],       # lower bounds
                                        [10, 10, 10, start_times_array.max()]       # upper bounds
                                )

                                try:
                                        params, cov = curve_fit(logistic_func, start_times_array, meas_array, p0=p0, bounds=bounds)
                                        
                                        fit_y = logistic_func(fit_x, *params)
                                        
                                        mass228_ppm_Th232_fit = logistic_func(indexChannel.time()-min(indexChannel.time()), *params)
                                        mass228_ppm_Th232_fit_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit)

                                        residuals = meas_array - logistic_func(start_times_array, *params)
                                        residual_std = np.std(residuals)

                                        mass228_ppm_Th232_fit_uncer = np.full_like(indexChannel.time(), residual_std)
                                        mass228_ppm_Th232_fit_uncer_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit_uncer', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit_uncer)

                                except RuntimeError as e:
                                        if 'Optimal parameters not found: Number of calls' in str(e):
                                                IoLog.warning('Could not find fit to data with Logistic model. Switching to Linear model.')
                                                settings["MB_fit"] = 'Linear'

                                                def fitLin(x, a, b):
                                                        return a + b*x

                                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                                fit_y = params[1]*fit_x + params[0]
                                                mass228_ppm_Th232_fit = params[0] + params[1] * indexChannel.time()
                                                mass228_ppm_Th232_fit_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit)

                                                fit_uncertainties = np.sqrt(
                                                        cov[0, 0] + 
                                                        (indexChannel.time())**2 * cov[1, 1] +
                                                        2 * indexChannel.time() * cov[0, 1]
                                                )
                                                mass228_ppm_Th232_fit_uncer_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
                                        else:
                                                raise

                        if degree_of_Ab228_fit == 'Polynomial':
                                def fitPol1(x, a, b, c):
                                        return a + b*x + c*x*x

                                params, cov = curve_fit(fitPol1, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                                mass228_ppm_Th232_fit = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2)
                                mass228_ppm_Th232_fit_timeseries = data.createTimeSeries('mass228_ppm_Th232_fit', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit)


                                x_vals = indexChannel.time() - min(indexChannel.time())
                                mass228_ppm_Th232_fit_uncer = np.sqrt(
                                        cov[0, 0]
                                        + x_vals ** 2 * cov[1, 1]
                                        + x_vals ** 4 * cov[2, 2]
                                        + 2 * x_vals * cov[0, 1]
                                        + 2 * x_vals ** 2 * cov[0, 2]
                                        + 2 * x_vals ** 3 * cov[1, 2]
                                )
                                mass228_ppm_Th232_fit_uncer_timeseries = data.createTimeSeries(
                                        f'mass228_ppm_Th232_fit_uncer', data.Intermediate, indexChannel.time(), mass228_ppm_Th232_fit_uncer
                                )

                        PLOT_Mon228.show()
                        g = PLOT_Mon228.addGraph()
                        g.setName('Monazite')
                        g.setLineStyle('lsNone')
                        g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[2 % len(PLOT_COLORS)])
                        g.setData(np.array(start_times_array), np.array(Zr228_Th232_ppm_Monazite))    

                        for x, y, y_err in zip(start_times_array, Zr228_Th232_ppm_Monazite,Zr228_Th232_ppm_Monazite_2SE):
                                v_error = PLOT_Mon228.addGraph()
                                vERR_pen = QPen(QColor('dark grey'))
                                vERR_pen.setStyle(Qt.SolidLine)
                                v_error.setPen(vERR_pen)
                                v_error.setData([x, x], [y - y_err, y + y_err])
                                v_error.removeFromLegend() 

                        g = PLOT_Mon228.addGraph()
                        g.setName("Fit")
                        fit_pen = QPen(QColor('black'))
                        fit_pen.setStyle(Qt.DashLine)
                        fit_pen.setWidth(2)
                        g.pen = fit_pen
                        g.setData(np.array(fit_x), np.array(fit_y))

                        if degree_of_Ab228_fit == 'Linear':
                                ann3.visible = True
                                ann3.text = f'''
                                <p style="color:black;">
                                <b>Fit Parameters:</b><br>
                                ({params[0]:.2f}±{cov[0,0]:.1e})+({params[1]:.2e}±{cov[1,1]:.1e})*t</p>'''
                                
                        if degree_of_Ab228_fit == 'Constant':
                                ann3.visible = True
                                ann3.text = f'''
                                <p style="color:black;">
                                <b>Fit Parameters:</b><br>
                                {params:.2f}±{cov:.2f}</p>'''        

                        if degree_of_Ab228_fit == 'Logistic':
                                ann3.visible = True
                                if len(params) == 4:
                                        param_uncertainties = np.sqrt(np.diag(cov))
                                        ann3.text = f'''
                                        <p style="color:black;">
                                        <b>Fit Equation:</b><br>
                                        {params[0]:.2f} + ({params[1]:.2f} - {params[0]:.2f}) / 
                                        (1 + exp({params[2]:.4f} × (t - {params[3]:.2f})))</p>'''
                                else:
                                        # fallback for linear case
                                        ann3.text = f'''
                                        <p style="color:black;">
                                        <b>Fit Parameters (linear):</b><br>
                                        ({params[0]:.2f} ± {cov[0,0]:.1e}) + ({params[1]:.2e} ± {cov[1,1]:.1e}) × t</p>'''

                        if degree_of_Ab228_fit == 'Polynomial':
                                ann3.visible = True
                                param_uncertainties = np.sqrt(np.diag(cov))
                                ann3.text = f'''
                                <p style="color:black;">
                                <b>Fit Equation:</b><br>
                                {params[0]:.2f} + 
                                {params[1]:.2e} × t + 
                                {params[2]:.2e} × t²</p>'''   

                        PLOT_Mon228.left().label = 'mass228 CPS ppm Th232 CPS'
                        PLOT_Mon228.bottom().label = 'Time since start of session (s)'
                        PLOT_Mon228.setLegendVisible(True)
                        PLOT_Mon228.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
                        PLOT_Mon228.setToolsVisible(False)
                        PLOT_Mon228.rescaleAxes()
                        PLOT_Mon228.replot()
                        
                        def AbundCorr228(sel):
                                result = Result()
                                Zr90_CPS = data.result(sel, data.timeSeries('Zr90_CPS')).value()
                                Th232_CPS = data.result(sel, data.timeSeries('Th232_CPS')).value()
                                AbundSens228 = data.result(sel,mass228_ppm_Th232_fit_timeseries).value() 
                                result = Zr90_CPS - (AbundSens228/1e+6*Th232_CPS) 
                                return result
                        def AbundCorr228_uncer(sel):
                                result = Result()
                                Th232_CPS = data.result(sel, data.timeSeries('Th232_CPS')).value()
                                AbundSens228 = data.result(sel,mass228_ppm_Th232_fit_timeseries).value() 
                                Zr90_2SE =  data.result(sel,data.timeSeries('Zr90_CPS')).uncertaintyAs2SE() 
                                Th232_2SE = data.result(sel,data.timeSeries('Th232_CPS')).uncertaintyAs2SE() 
                                AbundSens228_2SE = 2*data.result(sel,mass228_ppm_Th232_fit_uncer_timeseries).value() 
                                result = np.sqrt((Zr90_2SE)**2+(Th232_CPS/1e+6*AbundSens228_2SE)**2+(AbundSens228/1e+6*Th232_2SE)**2)
                                return result

                        data.registerAssociatedResult('InterferCorr228',AbundCorr228)  
                        data.registerAssociatedResult('InterferCorr228_2SE',AbundCorr228_uncer)

                        meas_main = []
                        start_times_all = []

                        for i, gr in enumerate(AS228_groups):
                                meas = []
                                meas_2SE = []
                                start_times = []
                                sg = data.selectionGroup(gr)
                                for sel in sg.selections():
                                        ratio = data.associatedResult(sel,'InterferCorr228').value()
                                        ratio_2SE = data.associatedResult(sel,'InterferCorr228_2SE').value()
                                        start_time = indexChannel.timeForSelection(sel)
                                        if not np.isnan(ratio):
                                                meas.append(ratio)
                                                meas_2SE.append(ratio_2SE)
                                                meas_main.append(ratio)
                                                start_times.append(start_time[0])
                                                start_times_all.append(start_time[0]) 
                                meas = np.array(meas)
                                meas_2SE = np.array(meas_2SE)
                                start_times = np.array(start_times)
                                median_meas = np.average(meas_main)
                                std_meas = np.std(meas_main)
                                if InfCorr_outlier_logical == True:
                                        mask = (meas >= median_meas - 2 * std_meas) & (meas <= median_meas + 2 * std_meas)
                                else: 
                                        mask = np.ones_like(meas, dtype=bool)
                                filtered_meas = meas[mask]
                                filtered_meas_2SE = meas_2SE[mask]
                                filtered_start_times = start_times[mask]

                                filtered_start_times = filtered_start_times-indexChannel.time()[0]
                                g = PLOT_Zr2O3.addGraph()
                                g.setName(gr)
                                g.setLineStyle('lsNone')
                                g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[i % len(PLOT_COLORS)])
                                g.setData(np.array(filtered_start_times), np.array(filtered_meas))

                                for x, y, y_err in zip(filtered_start_times, filtered_meas, filtered_meas_2SE):
                                        v_error = PLOT_Zr2O3.addGraph()
                                        vERR_pen = QPen(QColor('light grey'))
                                        vERR_pen.setStyle(Qt.SolidLine)
                                        v_error.setPen(vERR_pen)
                                        v_error.setData([x, x], [y - y_err, y + y_err])
                                        v_error.removeFromLegend()  # Remove vertical error bar from legend

                        start_times_all = start_times_all-indexChannel.time()[0]
                        meas_array = np.array(meas_main)
                        start_times_array = np.array(start_times_all)

                        if InfCorr_outlier_logical == True:
                                mask = (meas_array >= np.mean(meas_array) - 2 * np.std(meas_array)) & (meas_array <= np.mean(meas_array) + 2 * np.std(meas_array))
                        else: 
                                mask = np.ones_like(meas_array, dtype=bool)
                        meas_array = meas_array[mask]
                        start_times_array = start_times_array[mask]

                        fit_x = np.linspace(start_times_array.min(),start_times_array.max(),50)

                        if InfCorr228_fit == 'Constant':
                                def fitConst(x, a, b):
                                        return a + b*0

                                params, cov = curve_fit(fitConst, start_times_array, meas_array,  ftol=1e-5)
                                params = np.mean(meas_array)
                                cov = np.std(meas_array)/np.sqrt(len(meas_array))
                                fit_y = params*np.ones(len(fit_x))
                                
                                InfCorr228_fitted = params+0*(indexChannel.time()-min(indexChannel.time()))
                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)
                                
                                InfCorr228_fitted_uncer = cov * np.ones(len(indexChannel.time()))
                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), InfCorr228_fitted_uncer)        
                        
                        if InfCorr228_fit == 'Linear':
                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                
                                InfCorr228_fitted = params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))
                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)
                                fit_uncertainties = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                                        2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)

                        if InfCorr228_fit == 'Logistic':
                                def logistic_func(x, a, b, c, d):
                                        return a + (b - a) / (1 + np.exp(c * (x - d)))

                                # Smart initial guess without scaling x
                                a0 = meas_array.min()                  # lower asymptote
                                b0 = meas_array.max()                  # upper asymptote
                                c0 = 10 / (start_times_array.max() - start_times_array.min())  # slope guess ~1 / range of x
                                d0 = np.quantile(start_times_array, 0.75)

                                p0 = [a0, b0, c0, d0]
                                bounds = (
                                        [-10, -10, -10, start_times_array.min()],       # lower bounds
                                        [30, 30, 10, start_times_array.max()]       # upper bounds
                                )

                                try:
                                        params, cov = curve_fit(logistic_func, start_times_array, meas_array, p0=p0, bounds=bounds)
                                        
                                        fit_y = logistic_func(fit_x, *params)
                                        
                                        InfCorr228_fitted = logistic_func(indexChannel.time()-min(indexChannel.time()), *params)
                                        InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)

                                        residuals = meas_array - logistic_func(start_times_array, *params)
                                        residual_std = np.std(residuals)

                                        InfCorr228_fitted_uncer = np.full_like(indexChannel.time(), residual_std)
                                        InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), InfCorr228_fitted_uncer)

                                except RuntimeError as e:
                                        if 'Optimal parameters not found: Number of calls' in str(e):
                                                IoLog.warning('Could not find fit to data with Logistic model. Switching to Linear model.')
                                                settings["MB_fit"] = 'Linear'

                                                def fitLin(x, a, b):
                                                        return a + b*x

                                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                                fit_y = params[1]*fit_x + params[0]
                                                InfCorr228_fitted = params[0] + params[1] * indexChannel.time()
                                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)

                                                fit_uncertainties = np.sqrt(
                                                        cov[0, 0] + 
                                                        (indexChannel.time())**2 * cov[1, 1] +
                                                        2 * indexChannel.time() * cov[0, 1]
                                                )
                                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries('InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
                                        else:
                                                raise

                        if InfCorr228_fit == 'Polynomial':
                                def fitPol1(x, a, b, c):
                                        return a + b*x + c*x*x

                                params, cov = curve_fit(fitPol1, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                                InfCorr228_fitted = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2)
                                InfCorr228_fitted_timeseries = data.createTimeSeries('InfCorr228_fit', data.Intermediate, indexChannel.time(), InfCorr228_fitted)


                                x_vals = indexChannel.time() - min(indexChannel.time())
                                InfCorr228_fitted_uncer = np.sqrt(
                                        cov[0, 0]
                                        + x_vals ** 2 * cov[1, 1]
                                        + x_vals ** 4 * cov[2, 2]
                                        + 2 * x_vals * cov[0, 1]
                                        + 2 * x_vals ** 2 * cov[0, 2]
                                        + 2 * x_vals ** 3 * cov[1, 2]
                                )
                                InfCorr228_fitted_uncer_timeseries = data.createTimeSeries(
                                        f'InfCorr228_fit_uncer', data.Intermediate, indexChannel.time(), InfCorr228_fitted_uncer
                                )

                        median_AbundCorr228 = np.median(meas_array)
                        median_AbundCorr228_SE = (np.percentile(meas_array, 75) - np.percentile(meas_array, 25)) / np.sqrt(len(meas_array))
                        median_AbundCorr228_array = np.linspace(np.median(meas_array),np.median(meas_array),len(start_times_array))
                        
                        PLOT_Zr2O3.show()
                        g = PLOT_Zr2O3.addGraph()
                        g.setName("Fit")
                        fit_pen = QPen(QColor('black'))
                        fit_pen.setStyle(Qt.DashLine)
                        fit_pen.setWidth(2)
                        g.pen = fit_pen
                        g.setData(np.array(fit_x), np.array(fit_y))

                        if InfCorr228_fit == 'Linear':
                                ann2.visible = True
                                ann2.text = f'''
                                <p style="color:black;">
                                <b>Fit Parameters:</b><br>
                                ({params[0]:.2f}±{cov[0,0]:.1e})+({params[1]:.2e}±{cov[1,1]:.1e})*t</p>'''
                                
                        if InfCorr228_fit == 'Constant':
                                ann2.visible = True
                                ann2.text = f'''
                                <p style="color:black;">
                                <b>Fit Parameters:</b><br>
                                {params:.2f}±{cov:.2f}</p>'''        
                        if InfCorr228_fit == 'Logistic':
                                ann2.visible = True
                                if len(params) == 4:
                                        param_uncertainties = np.sqrt(np.diag(cov))
                                        ann2.text = f'''
                                        <p style="color:black;">
                                        <b>Fit Equation:</b><br>
                                        {params[0]:.2f} + ({params[1]:.2f} - {params[0]:.2f}) / 
                                        (1 + exp({params[2]:.4f} × (t - {params[3]:.2f})))</p>'''
                                else:
                                        # fallback for linear case
                                        ann2.text = f'''
                                        <p style="color:black;">
                                        <b>Fit Parameters (linear):</b><br>
                                        ({params[0]:.2f} ± {cov[0,0]:.1e}) + ({params[1]:.2e} ± {cov[1,1]:.1e}) × t</p>'''

                        if InfCorr228_fit == 'Polynomial':
                                ann2.visible = True
                                param_uncertainties = np.sqrt(np.diag(cov))
                                ann2.text = f'''
                                <p style="color:black;">
                                <b>Fit Equation:</b><br>
                                {params[0]:.2f} + 
                                {params[1]:.2e} × t + 
                                {params[2]:.2e} × t²</p>'''    
                        PLOT_Zr2O3.left().label = 'Zr2O3 CPS'
                        PLOT_Zr2O3.bottom().label = 'Time since start of session (s)'
                        PLOT_Zr2O3.setLegendVisible(True)
                        PLOT_Zr2O3.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
                        PLOT_Zr2O3.setToolsVisible(True)
                        PLOT_Zr2O3.rescaleAxes()
                        PLOT_Zr2O3.replot()
        else:
                PLOT_Mon228.hide()
                PLOT_Zr2O3.hide()          


##%     Secular Equilibrium U234 and U238

        if U234secDis: 
                PLOT_Mon234.clearGraphs()
                PLOT_Mon234.show()


                def Abundance232_on_234(sel):
                        result = Result()
                        result = ((data.result(sel, data.timeSeries('U234_CPS')).value() - ((lambda_U238/lambda_U234* (data.result(sel, data.timeSeries('U238_CPS')).value()))))/(data.result(sel, data.timeSeries('Th232_CPS')).value()))*1e+6
                        return result
        
                def Abundance232_on_234_2SE(sel):
                        result = Result()
                        U234_CPS = data.result(sel, data.timeSeries('U234_CPS')).value()
                        Th232_CPS = data.result(sel, data.timeSeries('Th232_CPS')).value()
                        U238_CPS = data.result(sel, data.timeSeries('U238_CPS')).value()
                        U234_CPSSE = data.result(sel, data.timeSeries('U234_CPS')).uncertaintyAs2SE()
                        Th232_CPSSE = data.result(sel, data.timeSeries('Th232_CPS')).uncertaintyAs2SE()
                        U238_CPSSE = data.result(sel, data.timeSeries('U238_CPS')).uncertaintyAs2SE()
                        result = np.sqrt((1e+6/Th232_CPS*U234_CPSSE)**2+(-1e+6/Th232_CPS*lambda_U238/lambda_Th230*U238_CPSSE)**2+((1e+6*Th232_CPSSE*(U238_CPS*lambda_U238/lambda_Th230-U234_CPS))/(Th232_CPS**2))**2)
                        return result       

                data.registerAssociatedResult('mass234_ppm_Th232_seq',Abundance232_on_234)  
                data.registerAssociatedResult('mass234_ppm_Th232_seq_2SE',Abundance232_on_234_2SE) 

                U234_Th232_ppm_Monazite = []
                U234_Th232_ppm_Monazite_2SE = []  
                start_time = []      

                sg = monazite_group
                for sel in sg.selections():
                        U234_Th232_ppm_timeseries_Monazite = data.associatedResult(sel,'mass234_ppm_Th232_seq').value()
                        U234_Th232_ppm_timeseries_Monazite_2SE = data.associatedResult(sel,'mass234_ppm_Th232_seq_2SE').value()
                        start_time_sel = indexChannel.timeForSelection(sel)
                        U234_Th232_ppm_Monazite.append(U234_Th232_ppm_timeseries_Monazite)
                        U234_Th232_ppm_Monazite_2SE.append(U234_Th232_ppm_timeseries_Monazite_2SE)
                        start_time.append(start_time_sel[0])
                
                start_times_array = np.array(start_time-indexChannel.time()[0])
                meas_array = np.array(U234_Th232_ppm_Monazite)
                
                fit_x = np.linspace(
                start_times_array.min(),
                start_times_array.max(),
                50
                )

                if degree_of_Ab230_fit == 'Constant':
                        def fitConst(x, a, b):
                                return a + b*0

                        params, cov = curve_fit(fitConst, start_times_array, meas_array,  ftol=1e-5)
                        params = np.mean(meas_array)
                        cov = np.mean(meas_array)/np.sqrt(len(meas_array))
                        fit_y = params*np.ones(len(fit_x))
                        
                        mass234_ppm_Th232_seq_fit = params+0*(indexChannel.time()-min(indexChannel.time()))
                        mass234_ppm_Th232_seq_fit_timeseries = data.createTimeSeries('mass234_ppm_Th232_seq_fit', data.Intermediate, indexChannel.time(), mass234_ppm_Th232_seq_fit)
                        
                        mass234_ppm_Th232_seq_fit_uncer = cov * np.ones(len(indexChannel.time()))
                        mass234_ppm_Th232_seq_fit_uncer_timeseries = data.createTimeSeries('mass234_ppm_Th232_seq_fit_uncer', data.Intermediate, indexChannel.time(), mass234_ppm_Th232_seq_fit_uncer
                        )        
                
                if degree_of_Ab230_fit == 'Linear':
                        def fitLin(x, a, b):
                                return a + b*x

                        params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                        fit_y = params[1]*fit_x + params[0]
                        
                        mass234_ppm_Th232_seq_fit = params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))
                        mass234_ppm_Th232_seq_fit_timeseries = data.createTimeSeries('mass234_ppm_Th232_seq_fit', data.Intermediate, indexChannel.time(), mass234_ppm_Th232_seq_fit)
                        fit_uncertainties = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                                2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                        mass234_ppm_Th232_seq_fit_uncer_timeseries = data.createTimeSeries('mass234_ppm_Th232_seq_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
                
                g = PLOT_Mon234.addGraph()
                g.setName('Monazite')
                g.setLineStyle('lsNone')
                g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[2 % len(PLOT_COLORS)])
                g.setData(np.array(start_times_array),np.array(U234_Th232_ppm_Monazite))    

                for x, y, y_err in zip(start_times_array, U234_Th232_ppm_Monazite,U234_Th232_ppm_Monazite_2SE):
                        v_error = PLOT_Mon234.addGraph()
                        vERR_pen = QPen(QColor('dark grey'))
                        vERR_pen.setStyle(Qt.SolidLine)
                        v_error.setPen(vERR_pen)
                        v_error.setData([x, x], [y - y_err, y + y_err])
                        v_error.removeFromLegend()  # Remove vertical error bar from legend

                g = PLOT_Mon234.addGraph()
                g.setName("Fit")
                fit_pen = QPen(QColor('black'))
                fit_pen.setStyle(Qt.DashLine)
                fit_pen.setWidth(2)
                g.pen = fit_pen
                g.setData(np.array(fit_x), np.array(fit_y))

                if degree_of_Ab230_fit == 'Linear':
                        ann4.visible = True
                        ann4.text = f'''
                        <p style="color:black;">
                        <b>Fit Parameters:</b><br>
                        ({params[0]:.2f}±{cov[0,0]:.1e})+({params[1]:.2e}±{cov[1,1]:.1e})*t</p>'''
                        
                if degree_of_Ab230_fit == 'Constant':
                        ann4.visible = True
                        ann4.text = f'''
                        <p style="color:black;">
                        <b>Fit Parameters:</b><br>
                        {params:.2f}±{cov:.2f}</p>'''        
                
                PLOT_Mon234.left().label = 'mass234 CPS ppm Th232 CPS'
                PLOT_Mon234.bottom().label = 'Time since start of session (s)'
                PLOT_Mon234.setLegendVisible(True)
                PLOT_Mon234.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)           
                PLOT_Mon234.setToolsVisible(False)
                PLOT_Mon234.rescaleAxes()
                PLOT_Mon234.replot()

        else:
                PLOT_Mon234.hide()  

##%     Mass Bias Calculation
        drs.message("MB correction...")
        drs.progress(40)
        PLOT.clearGraphs()
        meas_main = []
        start_times_all = []

        for i, gr in enumerate(MB_groups):
                meas = []
                meas_2SE = []
                start_times = []
                sg = data.selectionGroup(gr)
                for sel in sg.selections():
                        ratio = data.result(sel, U238_U235_timeSeries).value()
                        ratio_2SE = data.result(sel,U238_U235_timeSeries).uncertaintyAs2SE()
                        start_time = U238_U235_timeSeries.timeForSelection(sel)
                        if not np.isnan(ratio):  # Check if the ratio and time are valid
                                meas.append(ratio)
                                meas_2SE.append(ratio_2SE)
                                meas_main.append(ratio)
                                start_times.append(start_time[0])
                                start_times_all.append(start_time[0])  # Collect the corresponding time value
                # reject outliers for plot
                meas = np.array(meas)
                meas_2SE = np.array(meas_2SE)
                start_times = np.array(start_times)
                mask = (meas >= np.percentile(meas, 25) - 1.5 * (np.percentile(meas, 75) - np.percentile(meas, 25))) & (meas <= np.percentile(meas, 75) + 1.5 * (np.percentile(meas, 75) - np.percentile(meas, 25)))
                filtered_meas = meas[mask]
                filtered_meas_2SE = meas_2SE[mask]
                filtered_start_times = start_times[mask]

                filtered_start_times = filtered_start_times-U238_U235_timeSeries.time()[0]
                g = PLOT.addGraph()
                g.setName(gr)
                g.setLineStyle('lsNone')
                g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[i % len(PLOT_COLORS)])
                g.setData(np.array(filtered_start_times), np.array(filtered_meas))                

                for x, y, y_err in zip(filtered_start_times, filtered_meas, filtered_meas_2SE):
                        v_error = PLOT.addGraph()
                        vERR_pen = QPen(QColor('light grey'))
                        vERR_pen.setStyle(Qt.SolidLine)
                        v_error.setPen(vERR_pen)
                        v_error.setData([x, x], [y - y_err, y + y_err])
                        v_error.removeFromLegend()  # Remove vertical error bar from legend

        PLOT.left().label = 'Measured Value'
        PLOT.bottom().label = 'Time'
        PLOT.setLegendVisible(True)
        PLOT.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
        PLOT.setToolsVisible(True)
        PLOT.rescaleAxes()
        PLOT.replot()


        # reject outliers for fit
        mask_MB = (meas_main >= np.percentile(meas_main, 25) - 1.5 * (np.percentile(meas_main, 75) - np.percentile(meas_main, 25))) & (meas_main <= np.percentile(meas_main, 75) + 1.5 * (np.percentile(meas_main, 75) - np.percentile(meas_main, 25)))
        meas_main = np.array(meas_main)
        start_times_all = np.array(start_times_all)
        meas_main_filtered = meas_main[mask_MB]
        start_times_all_filtered = start_times_all[mask_MB]
        start_times_all_filtered = start_times_all_filtered-U238_U235_timeSeries.time()[0]

        meas_array = np.array(meas_main_filtered)
        start_times_array = np.array(start_times_all_filtered)
   
        fit_x = np.linspace(
            start_times_array.min(),
            start_times_array.max(),
            50
        )

        if degree_of_MB_fit == 'Constant':
                def fitConst(x, a, b):
                        return a + b*0

                params, cov = curve_fit(fitConst, start_times_array, meas_array,  ftol=1e-5)

                fit_y = params[0]*np.ones(len(fit_x))
                Mass_Bias_per_mass = ((params[0]+0*(indexChannel.time()-min(indexChannel.time())))-U238_U235)/(U238_U235*(238-235))
                Mass_Bias_per_mass_timeSeries = data.createTimeSeries('MassBias', data.Intermediate, indexChannel.time(), Mass_Bias_per_mass)
        if degree_of_MB_fit == 'Linear':
                def fitLin(x, a, b):
                        return a + b*x

                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                fit_y = params[1]*fit_x + params[0]
                Mass_Bias_per_mass = ((params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))-U238_U235)/(U238_U235*(238-235))
                Mass_Bias_per_mass_timeSeries = data.createTimeSeries('MassBias', data.Intermediate, indexChannel.time(), Mass_Bias_per_mass)

        if degree_of_MB_fit == 'Exponential':
                def fitExp(x, a, b, c):
                        return a + b * np.exp(-c * x)
                initial_guess = [np.mean(meas_array), 0, 0.01]
                try:
                        params, cov = curve_fit(fitExp, start_times_array, meas_array, p0=initial_guess)
                        fit_y = params[0] + params[1] * np.exp(-params[2] * fit_x)
                        Mass_Bias_per_mass = ((params[0]+params[1]*np.exp(-params[2]*(indexChannel.time()-min(indexChannel.time()))))-U238_U235)/(U238_U235*(238-235))
                        Mass_Bias_per_mass_timeSeries = data.createTimeSeries('MassBias', data.Intermediate, indexChannel.time(), Mass_Bias_per_mass)

                except RuntimeError as e:
                        if 'Optimal parameters not found: Number of calls' in str(e):

                                IoLog.warning('Could not find fit to data with Exp model. Switching to Linear model.')

                                settings["MB_fit"] = 'Linear'

                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)
                                
                                fit_y = params[1]*fit_x + params[0]
                                Mass_Bias_per_mass = ((params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))-U238_U235)/(U238_U235*(238-235))
                                Mass_Bias_per_mass_timeSeries = data.createTimeSeries('MassBias', data.Intermediate, indexChannel.time(), Mass_Bias_per_mass)

                        else:
                                raise

        if degree_of_MB_fit == 'Logarithmic':
                def fitLog(x, a, b, c):
                        return a + b * np.log(x + c)  # Ensure x + c is > 0
                
                initial_guess = [np.mean(meas_array), 1, 2000]
                try:    
                        params, cov = curve_fit(fitLog, start_times_array, meas_array, p0=initial_guess)
                        fit_y = params[0] + params[1] * np.log(fit_x + params[2])
                        Mass_Bias_per_mass = ((params[0]+params[1]*np.log(params[2]+(indexChannel.time()-min(indexChannel.time()))))-U238_U235)/(U238_U235*(238-235))
                        Mass_Bias_per_mass_timeSeries = data.createTimeSeries('MassBias', data.Intermediate, indexChannel.time(), Mass_Bias_per_mass)

                except RuntimeError as e:
                        if 'Optimal parameters not found: Number of calls' in str(e):

                                IoLog.warning('Could not find fit to data with Exp model. Switching to Linear model.')

                                settings["MB_fit"] = 'Linear'

                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                Mass_Bias_per_mass = ((params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))-U238_U235)/(U238_U235*(238-235))
                                Mass_Bias_per_mass_timeSeries = data.createTimeSeries('MassBias', data.Intermediate, indexChannel.time(), Mass_Bias_per_mass)                        
                        else:
                                raise
        if degree_of_MB_fit == 'Polynomial':
                def fitPol1(x, a, b, c):
                        return a + b*x + c*x*x

                params, cov = curve_fit(fitPol1, start_times_array, meas_array, ftol=1e-5)

                fit_y = params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                Mass_Bias_per_mass = ((params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2)-U238_U235)/(U238_U235*(238-235))
                Mass_Bias_per_mass_timeSeries = data.createTimeSeries('MassBias', data.Intermediate, indexChannel.time(), Mass_Bias_per_mass)

        g = PLOT.addGraph()
        g.setName("Fit")
        fit_pen = QPen(QColor('black'))
        fit_pen.setStyle(Qt.DashLine)
        fit_pen.setWidth(2)
        g.pen = fit_pen
        g.setData(fit_x, fit_y)

        ann.visible = True
        if degree_of_MB_fit == 'Constant':
                ann.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f}</p>'''
        elif len(params) == 2:
                ann.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f} + {params[1]:.3e}*t</p>'''
        elif degree_of_MB_fit == 'Exponential': 
                ann.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f} + {params[1]:.2f} * exp(-{params[2]:.3e}*t)</p>'''
        elif degree_of_MB_fit == 'Logarithmic': 
                ann.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f} + {params[1]:.2f} * ln(t + {params[2]:.2f})</p>'''
        elif degree_of_MB_fit == 'Polynomial':
                ann.visible = True
                ann.text = f'''
                <p style="color:black;">
                <b>Fit Equation:</b><br>
                {params[0]:.2f} + 
                {params[1]:.2e} × t + 
                {params[2]:.2e} × t²</p>'''    
        PLOT.left().label = '238U/235U'
        PLOT.bottom().label = 'Time since start of session (s)'
        PLOT.setToolsVisible(True)
        PLOT.rescaleAxes()
        PLOT.replot()


##%     Th230 (and U234) corrected        
        if Zr2O3corr_logical == True:
                Th230_CPS_corr = data.timeSeries('Th230_CPS').data() - InfCorr228_fitted*scaling_factor - (mass230_ppm_Th232_seq_fit_timeseries.data()*data.timeSeries('Th232_CPS').data()/1e+6)
                #Th230_CPS_corr = data.timeSeries('Th230_CPS').data() - median_AbundCorr228*scaling_factor - (mass230_ppm_Th232_seq_fit_timeseries.data()*data.timeSeries('Th232_CPS').data()/1e+6)

        else:
                Th230_CPS_corr = data.timeSeries('Th230_CPS').data() - (mass230_ppm_Th232_seq_fit_timeseries.data()*data.timeSeries('Th232_CPS').data()/1e+6)
        Th230_CPS_corr_timeseries = data.createTimeSeries('Th230_CPS_corr', data.Intermediate, indexChannel.time(), Th230_CPS_corr)

        if how235or238 == 'based on 235U*137.818':
                def Th230_U238_activity_corr_AR(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('Th230_CPS_corr')).value()/(data.result(sel,data.timeSeries('U235_CPS')).value()*U238_U235)*(lambda_Th230/lambda_U238)
                        return result
                def Th230_U238_activity_corr_uncer_AR(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(Th230)/(U238)_ROI').value()*np.sqrt((data.result(sel,data.timeSeries('Th230_CPS_corr')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('Th230_CPS_corr')).value())**2+(data.result(sel,data.timeSeries('U235_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('U235_CPS')).value())**2)
                        return result
        
        if how235or238 == 'based on 238U':
                def Th230_U238_activity_corr_AR(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('Th230_CPS_corr')).value()/(data.result(sel,data.timeSeries('U238_CPS')).value())*(lambda_Th230/lambda_U238)
                        return result
                def Th230_U238_activity_corr_uncer_AR(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(Th230)/(U238)_ROI').value()*np.sqrt((data.result(sel,data.timeSeries('Th230_CPS_corr')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('Th230_CPS_corr')).value())**2+(data.result(sel,data.timeSeries('U238_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('U238_CPS')).value())**2)
                        return result

        data.registerAssociatedResult('(Th230)/(U238)_ROI',Th230_U238_activity_corr_AR) 
        data.registerAssociatedResult('(Th230)/(U238)_2SE_ROI',Th230_U238_activity_corr_uncer_AR) 

        if U234secDis: 
                U234_CPS_corr = data.timeSeries('U234_CPS').data() - (mass234_ppm_Th232_seq_fit_timeseries.data()*data.timeSeries('Th232_CPS').data()/1e+6)
                U234_CPS_corr_timeseries = data.createTimeSeries('U234_CPS_corr', data.Intermediate, indexChannel.time(), U234_CPS_corr)
       
##%     RSF 
        drs.message("RSF...")
        drs.progress(50)

        PLOT_RSF.clearGraphs()
        RSF_230_238 = []
        RSF_final = []
        RSF_final_2SE = []
        start_time = []
        for i, gr in enumerate(RSF_groups_plot):
                sg = data.selectionGroup(gr)

                RSF_RM = []
                RSF_RM_2SE = []
                start_times_all = []
                if RSF_how == 'based on RM U & Th concentration':
                        try:
                                U_RM_ppm_RSF = data.referenceMaterialData(gr)["U"].value()
                                Th_RM_ppm_RSF = data.referenceMaterialData(gr)["Th"].value()

                        except:
                                print(f"Skipping group {gr}: U or Th value not found.")
                                continue

                        RSF = (data.timeSeries('U238_CPS').data()/U_RM_ppm_RSF)/(data.timeSeries('Th232_CPS').data()/Th_RM_ppm_RSF)
                        RSF_timeseries = data.createTimeSeries(f'RSF_{settings["RSF_on"]}', data.Intermediate, indexChannel.time(), RSF)
                for sel in sg.selections():
                        if RSF_how == 'based on secular equilibrium':
                                RSF_sel_230_238 = data.associatedResult(sel,'(Th230)/(U238)_ROI').value()
                                RSF_timeseries_RM = 1/data.associatedResult(sel,'(Th230)/(U238)_ROI').value()
                                RSF_timeseries_RM_2SE = (1/data.associatedResult(sel,'(Th230)/(U238)_ROI').value())*data.associatedResult(sel,'(Th230)/(U238)_2SE_ROI').value()/(data.associatedResult(sel,'(Th230)/(U238)_ROI').value())
                                start_times_RM = indexChannel.timeForSelection(sel)
                                RSF_230_238.append(RSF_sel_230_238)

                        else:
                                RSF_timeseries_RM = data.result(sel, RSF_timeseries).value()
                                RSF_timeseries_RM_2SE = data.result(sel, RSF_timeseries).uncertaintyAs2SE()                               
                                start_times_RM = RSF_timeseries.timeForSelection(sel)
                        
                        RSF_RM.append(RSF_timeseries_RM)
                        RSF_RM_2SE.append(RSF_timeseries_RM_2SE)
                        start_times_all.append(start_times_RM[0])

                RSF_RM = np.array(RSF_RM)
                RSF_RM_2SE = np.array(RSF_RM_2SE)
                start_times_all = np.array(start_times_all)
                
                if RSF_outlier_logical == True:
                        mask_RSF = (RSF_RM >= np.percentile(RSF_RM, 25) - 1.5 * (np.percentile(RSF_RM, 75) - np.percentile(RSF_RM, 25))) & (RSF_RM <= np.percentile(RSF_RM, 75) + 1.5 * (np.percentile(RSF_RM, 75) - np.percentile(RSF_RM, 25)))
                else:
                        mask_RSF = np.ones_like(RSF_RM, dtype=bool)
                RSF_RM_filtered = RSF_RM[mask_RSF]
                RSF_RM_2SE_filtered = RSF_RM_2SE[mask_RSF]
                start_times_all_filtered = start_times_all[mask_RSF]
                start_times_all_filtered = start_times_all_filtered-indexChannel.time()[0]
                
                if RSF_how == 'based on RM U & Th concentration':
                        if gr == settings["RSF_on"]:
                                RSF_final = np.array(RSF_RM_filtered)
                                RSF_final_2SE = np.array(RSF_RM_2SE_filtered)
                                start_time = np.array(start_times_all_filtered)
                if RSF_how == 'based on secular equilibrium':
                        RSF_final = np.concatenate((RSF_final, RSF_RM_filtered), axis=0) if len(RSF_final) > 0 else np.array(RSF_RM_filtered)
                        RSF_final_2SE = np.concatenate((RSF_final_2SE, RSF_RM_2SE_filtered), axis=0) if len(RSF_final_2SE) > 0 else np.array(RSF_RM_2SE_filtered)
                        start_time = np.concatenate((start_time, start_times_all_filtered), axis=0) if len(start_time) > 0 else np.array(start_times_all_filtered)                      

                g = PLOT_RSF.addGraph()
                g.setName(gr)
                g.setLineStyle('lsNone')
                g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[i % len(PLOT_COLORS)])
                g.setData(np.array(start_times_all_filtered), np.array(RSF_RM_filtered))   
                
                for x, y, y_err in zip(start_times_all_filtered, RSF_RM_filtered, RSF_RM_2SE_filtered):
                        v_error = PLOT_RSF.addGraph()
                        vERR_pen = QPen(QColor('dark grey'))
                        vERR_pen.setStyle(Qt.SolidLine)
                        v_error.setPen(vERR_pen)
                        v_error.setData([x, x], [y - y_err, y + y_err])
                        v_error.removeFromLegend() 
                
                if RSF_how == 'based on RM U & Th concentration':
                        data.removeTimeSeries(f'RSF_{settings["RSF_on"]}')

        start_time = np.array(start_time)
        RSF_final = np.array(RSF_final)
        RSF_final_2SE = np.array(RSF_final_2SE)
        fit_x = np.linspace(start_time.min(),start_time.max(),50)

        # Scaling of values to ensure secular equilibrium
        if RSF_how == 'based on secular equilibrium':
               RSF_final = 1/(np.average(RSF_final)*np.average(1/RSF_final))*RSF_final

        if RSF_fit == 'Constant':
                params = np.mean(RSF_final)
                cov = np.std(RSF_final)/np.sqrt(len(RSF_final))
                fit_y = params*np.ones(len(fit_x))
                
                RSF_fitted = params+0*(indexChannel.time()-min(indexChannel.time()))
                RSF_fit_timeseries = data.createTimeSeries('RSF_fit', data.Intermediate, indexChannel.time(), RSF_fitted)
                
                RSF_fit_uncer = cov * np.ones(len(indexChannel.time()))
                RSF_fit_uncer_timeseries = data.createTimeSeries('RSF_fit_uncer', data.Intermediate, indexChannel.time(), RSF_fit_uncer)        
        
        if RSF_fit == 'Linear':
                def fitLin(x, a, b):
                        return a + b*x

                params, cov = curve_fit(fitLin, start_time, RSF_final, ftol=1e-5)

                fit_y = params[1]*fit_x + params[0]
                RSF_fitted = params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))
                RSF_fit_timeseries = data.createTimeSeries('RSF_fit', data.Intermediate, indexChannel.time(), RSF_fitted)
                fit_uncertainties = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                        2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                RSF_fit_uncer_timeseries = data.createTimeSeries('RSF_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)

        if RSF_fit == 'Logistic':
                def logistic_func(x, a, b, c, d):
                        return a + (b - a) / (1 + np.exp(c * (x - d)))

                # Smart initial guess without scaling x
                a0 = RSF_final.min()                  # lower asymptote
                b0 = RSF_final.max()                  # upper asymptote
                c0 = 10 / (start_time.max() - start_time.min())  # slope guess ~1 / range of x
                d0 = np.quantile(start_time, 0.5)
                p0 = [a0, b0, c0, d0]
                bounds = (
                        [0, 0, -10, start_time.min()],       # lower bounds
                        [20, 10, 10, start_time.max()]       # upper bounds
                )

                try:
                        params, cov = curve_fit(logistic_func, start_time, RSF_final, p0=p0, bounds=bounds)
                        
                        fit_y = logistic_func(fit_x, *params)
                        
                        RSF_fitted = logistic_func(indexChannel.time()-min(indexChannel.time()), *params)
                        RSF_fit_timeseries = data.createTimeSeries('RSF_fit', data.Intermediate, indexChannel.time(), RSF_fitted)

                        residuals = RSF_final - logistic_func(start_time, *params)
                        residual_std = np.std(residuals)

                        RSF_fit_uncer = np.full_like(indexChannel.time(), residual_std)
                        RSF_fit_uncer_timeseries = data.createTimeSeries('RSF_fit_uncer', data.Intermediate, indexChannel.time(), RSF_fit_uncer)

                except RuntimeError as e:
                        if 'Optimal parameters not found: Number of calls' in str(e):
                                IoLog.warning('Could not find fit to data with Logistic model. Switching to Linear model.')
                                settings["MB_fit"] = 'Linear'

                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_time, RSF_final, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                RSF_fitted = params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))
                                RSF_fit_timeseries = data.createTimeSeries('RSF_fit', data.Intermediate, indexChannel.time(), RSF_fitted)
                                fit_uncertainties = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                                        2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                                RSF_fit_uncer_timeseries = data.createTimeSeries('RSF_fit_uncer', data.Intermediate, indexChannel.time(), fit_uncertainties)
                        else:
                                raise

        if RSF_fit == 'Polynomial':
                def fitPol1(x, a, b, c):
                        return a + b*x + c*x*x

                params, cov = curve_fit(fitPol1, start_time, RSF_final, ftol=1e-5)

                fit_y = params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                RSF_fitted = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2)
                RSF_fit_timeseries = data.createTimeSeries('RSF_fit', data.Intermediate, indexChannel.time(), RSF_fitted)


                x_vals = indexChannel.time() - min(indexChannel.time())
                RSF_fit_uncer = np.sqrt(
                        cov[0, 0]
                        + x_vals ** 2 * cov[1, 1]
                        + x_vals ** 4 * cov[2, 2]
                        + 2 * x_vals * cov[0, 1]
                        + 2 * x_vals ** 2 * cov[0, 2]
                        + 2 * x_vals ** 3 * cov[1, 2]
                )
                RSF_fit_uncer_timeseries = data.createTimeSeries(
                        f'RSF_fit_uncer', data.Intermediate, indexChannel.time(), RSF_fit_uncer
                )

        if RSF_fit == 'Polynomial3':
                def fitLin(x, a, b, c, d):
                        return a + b*x + c*x*x + d*x*x*x

                params, cov = curve_fit(fitLin, start_time, RSF_final, ftol=1e-5)

                fit_y = params[3]*fit_x*fit_x*fit_x + params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                RSF_fitted = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2+params[3]*(indexChannel.time()-min(indexChannel.time()))**3)
                RSF_fit_timeseries = data.createTimeSeries(f'RSF_fit', data.Intermediate, indexChannel.time(), RSF_fitted)


                x_vals = indexChannel.time() - min(indexChannel.time())
                x1 = x_vals
                x2 = x_vals**2
                x3 = x_vals**3

                RSF_fit_uncer = np.sqrt(
                cov[0, 0] +
                cov[1, 1] * x1**2 +
                cov[2, 2] * x2**2 +
                cov[3, 3] * x3**2 +
                2 * cov[0, 1] * x1 +
                2 * cov[0, 2] * x2 +
                2 * cov[0, 3] * x3 +
                2 * cov[1, 2] * x1 * x2 +
                2 * cov[1, 3] * x1 * x3 +
                2 * cov[2, 3] * x2 * x3
                )
                RSF_fit_uncer_timeseries = data.createTimeSeries(
                        f'RSF_fit_uncer', data.Intermediate, indexChannel.time(), RSF_fit_uncer
                )

        g = PLOT_RSF.addGraph()
        if RSF_how == 'based on RM U & Th concentration':
                RSF_average_for = settings["RSF_on"]
                g.setName(f"Fit for {RSF_average_for}")
        if RSF_how == 'based on secular equilibrium':
                g.setName(f"Gobal fit")
        fit_pen = QPen(QColor('black'))
        fit_pen.setStyle(Qt.DashLine)
        fit_pen.setWidth(2)
        g.pen = fit_pen
        g.setData(np.array(fit_x), np.array(fit_y))

        if RSF_fit == 'Linear':
                ann5.visible = True
                ann5.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                ({params[0]:.2f}±{cov[0,0]:.1e})+({params[1]:.2e}±{cov[1,1]:.1e})*t</p>'''
                
        if RSF_fit == 'Constant':
                ann5.visible = True
                ann5.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params:.2f}±{cov:.2f}</p>'''        
        if RSF_fit == 'Logistic':
                ann5.visible = True
                if len(params) == 4:
                        param_uncertainties = np.sqrt(np.diag(cov))
                        ann5.text = f'''
                        <p style="color:black;">
                        <b>Fit Equation:</b><br>
                        {params[0]:.2f} + ({params[1]:.2f} - {params[0]:.2f}) / 
                        (1 + exp({params[2]:.4f} × (t - {params[3]:.2f})))</p>'''
                else:
                        # fallback for linear case
                        ann5.text = f'''
                        <p style="color:black;">
                        <b>Fit Parameters (linear):</b><br>
                        ({params[0]:.2f} ± {cov[0,0]:.1e}) + ({params[1]:.2e} ± {cov[1,1]:.1e}) × t</p>'''

        if RSF_fit == 'Polynomial':
                ann5.visible = True
                param_uncertainties = np.sqrt(np.diag(cov))
                ann5.text = f'''
                <p style="color:black;">
                <b>Fit Equation:</b><br>
                {params[0]:.2f} + 
                {params[1]:.2e} × t + 
                {params[2]:.2e} × t²</p>'''          

        if RSF_fit == 'Polynomial3':
                ann5.visible = True
                param_uncertainties = np.sqrt(np.diag(cov))
                ann5.text = f'''
                <p style="color:black;">
                <b>Fit Equation:</b><br>
                {params[0]:.2f} + 
                {params[1]:.2e} × t + 
                {params[2]:.2e} × t² +
                {params[3]:.2e} × t³</p>'''                    


        PLOT_RSF.left().label = 'RSF'
        PLOT_RSF.bottom().label = 'Time since start of session (s)'
        PLOT_RSF.setLegendVisible(True)
        PLOT_RSF.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)    
        PLOT_RSF.setToolsVisible(True)
        PLOT_RSF.rescaleAxes()
        PLOT_RSF.replot()

##%     Lc check
        drs.message("Analyze Detection Limit...")
        drs.progress(60)
##%     Length of Actual signal
        def Signal_s(sel):
                result = Result()
                start_time_sig_sel = Th230_CPS_corr_timeseries.timeForSelection(sel)
                result = len(start_time_sig_sel)*dwell_Th230_s
                return result
 
##%     Length of Background signal
        def BG_s(sel):
                result = Result()
                start_time_sig = Th230_CPS_corr_timeseries.timeForSelection(sel)[0]
                
                BaselineGroup = data.selectionGroupNames(data.Baseline)
                duration_BG_ms = []
                Th230_BG = []
                start_time_BG = []
                for i,gr in enumerate(BaselineGroup):
                        bg = data.selectionGroup(gr)
                        for sel_bg in bg.selections():
                                Th230_BG_sel = data.result(sel_bg,data.timeSeries('Th230')).value() 
                                start_time_BG_sel = data.timeSeries('Th230').timeForSelection(sel_bg)
                                duration_bg_ms_sel = len(start_time_BG_sel)*dwell_Th230_s

                                duration_BG_ms.append(duration_bg_ms_sel)
                                Th230_BG.append(Th230_BG_sel)
                                start_time_BG.append(start_time_BG_sel[0])

                backgrounds = sorted(zip(start_time_BG, duration_BG_ms), key=lambda x: x[0])
                matched_bg = None
                for bg_start, bg_duration in backgrounds:
                        if bg_start < start_time_sig:
                                matched_bg = bg_duration  # Store current background as the candidate
                        else:
                                break  # Stop when background start time is later than signal start time

                if matched_bg:
                        result = matched_bg
                else:
                        result = None

                return result
        
##%     Th230 Background signal        
        def BG_Th230(sel):
                result = Result()
                start_time_sig = Th230_CPS_corr_timeseries.timeForSelection(sel)[0]
                
                BaselineGroup = data.selectionGroupNames(data.Baseline)
                duration_BG_ms = []
                Th230_BG = []
                start_time_BG = []

                for i,gr in enumerate(BaselineGroup):
                        bg = data.selectionGroup(gr)
                        for sel_bg in bg.selections():
                                Th230_BG_sel = data.result(sel_bg,data.timeSeries('Th230')).value() 
                                start_time_BG_sel = data.timeSeries('Th230').timeForSelection(sel_bg)
                                duration_bg_ms_sel = len(start_time_BG_sel)*dwell_Th230_s

                                duration_BG_ms.append(duration_bg_ms_sel)
                                Th230_BG.append(Th230_BG_sel)
                                start_time_BG.append(start_time_BG_sel[0])

                backgrounds = sorted(zip(start_time_BG, Th230_BG), key=lambda x: x[0])
                matched_bg = None                
                for bg_start, bg_Th230 in backgrounds:
                        if bg_start < start_time_sig:
                                matched_bg = bg_Th230  # Store current background as the candidate
                        else:
                                break  # Stop when background start time is later than signal start time

                if matched_bg:
                        result = matched_bg
                else:
                        result = None

                return result

##%     Th230 Background signal uncertainty                
        def BG_1SE_Th230(sel):
                result = Result()
                start_time_sig = Th230_CPS_corr_timeseries.timeForSelection(sel)[0]
                
                BaselineGroup = data.selectionGroupNames(data.Baseline)
                duration_BG_ms = []
                Th230_BG_1SE = []
                start_time_BG = []
                matched_bg = None
                for i,gr in enumerate(BaselineGroup):
                        bg = data.selectionGroup(gr)
                        for sel_bg in bg.selections():
                                Th230_BG_sel = 0.5*data.result(sel_bg,data.timeSeries('Th230')).uncertaintyAs2SE() 
                                start_time_BG_sel = data.timeSeries('Th230').timeForSelection(sel_bg)
                                duration_bg_ms_sel = len(start_time_BG_sel)*dwell_Th230_s

                                duration_BG_ms.append(duration_bg_ms_sel)
                                Th230_BG_1SE.append(Th230_BG_sel)
                                start_time_BG.append(start_time_BG_sel[0])

                backgrounds = sorted(zip(start_time_BG, Th230_BG_1SE), key=lambda x: x[0])
                for bg_start, bg_Th230 in backgrounds:
                        if bg_start < start_time_sig:
                                matched_bg = bg_Th230  # Store current background as the candidate
                        else:
                                break  # Stop when background start time is later than signal start time

                if matched_bg:
                        result = matched_bg
                else:
                        result = None

                return result
        
##%     Uncertainty on Th230 corr (gem. Guillong)
        def Th230_corr_uncer(sel):
                result = Result()
                Th230corr = data.result(sel,data.timeSeries('Th230_CPS_corr')).value() 
                Abund230 = data.result(sel,mass230_ppm_Th232_seq_fit_timeseries).value()
                Abund230_unc = data.result(sel,mass230_ppm_Th232_seq_fit_uncer_timeseries).value()
                duration_sig = data.associatedResult(sel,'Signal_Th230_s').value()
                duration_bg = data.associatedResult(sel,'BG_Th230_s').value()
                result_bg = data.associatedResult(sel,'BG_Th230').value()
                if Zr2O3corr_logical == True:
                        if corr228_int_ext == 'Th230 directly measured in Zirconblank':
                                result = (np.sqrt((((np.sqrt(Th230corr* duration_sig) / duration_sig) +
                                (np.sqrt(result_bg * duration_bg) / duration_bg)) 
                                / Th230corr)**2 + 
                                (Abund230_unc / Abund230)**2)) * Th230corr   

                        else:
                                abund228 = data.result(sel,InfCorr228_fitted_timeseries).value()
                                abund228_uncer = data.result(sel,InfCorr228_fitted_uncer_timeseries).value()
                                result = (np.sqrt((((np.sqrt(Th230corr* duration_sig) / duration_sig) + 
                                                (np.sqrt(result_bg * duration_bg) / duration_bg)) 
                                                / Th230corr)**2 + 
                                                (abund228_uncer / abund228)**2 + 
                                                (Abund230_unc /Abund230)**2)) * Th230corr           
                                 
                else:
                        result = (np.sqrt((((np.sqrt(Th230corr* duration_sig) / duration_sig) +
                                        (np.sqrt(result_bg * duration_bg) / duration_bg)) 
                                        / Th230corr)**2 + 
                                        (Abund230_unc / Abund230)**2)) * Th230corr    
                         
                return result

##%     Limit of detection for Th230
        def Lc_Th230(sel):
                result = Result()
                duration_sig = data.associatedResult(sel,'Signal_Th230_s').value()
                duration_bg = data.associatedResult(sel,'BG_Th230_s').value()
                result_bg = data.associatedResult(sel,'BG_Th230').value()  

                result = (3.29*np.sqrt(result_bg*duration_sig*(1+duration_sig/duration_bg))+2.71)/duration_sig

                return result
        
##%     Lc check
        def Lc_check(sel):
                result = Result()
                Th230_corr_uncer = data.associatedResult(sel,'Th230_corr_uncer').value()
                Lc_Th230 = data.associatedResult(sel,'Lc_Th230').value()                
                Th230corr = data.result(sel,data.timeSeries('Th230_CPS_corr')).value() 

                if (Th230corr ) > (Lc_Th230):

                        result = 1
                else:
                        result = 0
                return result

##%     Lc check
        def Lc_check_safe(sel):
                result = Result()
                Th230_corr_uncer = data.associatedResult(sel,'Th230_corr_uncer').value()
                Lc_Th230 = data.associatedResult(sel,'Lc_Th230').value()                
                Th230corr = data.result(sel,data.timeSeries('Th230_CPS_corr')).value() 

                if (Th230corr - 1*Th230_corr_uncer) > (2*Lc_Th230):

                        result = 1
                else:
                        result = 0
                return result

        data.registerAssociatedResult('Signal_Th230_s',Signal_s)  
        data.registerAssociatedResult('BG_Th230_s',BG_s) 
        data.registerAssociatedResult('BG_Th230',BG_Th230)
        data.registerAssociatedResult('BG_1SE_Th230',BG_1SE_Th230)    
        data.registerAssociatedResult('Th230_corr_uncer',Th230_corr_uncer)     
        data.registerAssociatedResult('Lc_Th230',Lc_Th230) 
        data.registerAssociatedResult('Lc_check',Lc_check)
        data.registerAssociatedResult('Lc_check_conservative',Lc_check_safe) 
 

##%     Calculating resulting ratios
        drs.message("Calculate results...")
        drs.progress(80)
        Th232_U238_activity_RSF_corr = (data.timeSeries('Th232_CPS').data()/data.timeSeries('U238_CPS').data())*RSF_fit_timeseries.data()*(lambda_Th232/lambda_U238)
        Th232_U238_activity_RSF_corr_timeseries = data.createTimeSeries('(Th232)/(U238)_RSF_corr', data.Intermediate, indexChannel.time(), Th232_U238_activity_RSF_corr)
        
        U238_Th232_activity_RSF_corr = 1/data.timeSeries('(Th232)/(U238)_RSF_corr').data()
        U238_Th232_activity_RSF_corr_timeseries = data.createTimeSeries('(U238)/(Th232)_RSF_corr', data.Output, indexChannel.time(), U238_Th232_activity_RSF_corr) 
        
        Th230_Th232_activity_MB_corr = data.timeSeries('Th230_CPS_corr').data()/(data.timeSeries('Th232_CPS').data()-data.timeSeries('Th232_CPS').data()*2*data.timeSeries('MassBias').data())*(lambda_Th230/lambda_Th232)
        Th230_Th232_activity_MB_corr_timeseries = data.createTimeSeries('(Th230)/(Th232)_MB_corr', data.Output, indexChannel.time(), Th230_Th232_activity_MB_corr)

        if how235or238 == 'based on 235U*137.818':
                Th230_U238_activity_RSF_corr = data.timeSeries('Th230_CPS_corr').data()/(data.timeSeries('U235_CPS').data()*U238_U235)*RSF_fit_timeseries.data()*(lambda_Th230/lambda_U238)
        if how235or238 == 'based on 238U':
                Th230_U238_activity_RSF_corr = data.timeSeries('Th230_CPS_corr').data()/(data.timeSeries('U238_CPS').data())*RSF_fit_timeseries.data()*(lambda_Th230/lambda_U238)
        
        Th230_U238_activity_RSF_corr_timeseries = data.createTimeSeries('(Th230)/(U238)_RSF_corr', data.Output, indexChannel.time(), Th230_U238_activity_RSF_corr)
        
        if U234secDis:                 
                U234_U238_activity_MB_corr = data.timeSeries('U234_CPS_corr').data()/(data.timeSeries('U238_CPS').data()-data.timeSeries('U238_CPS').data()*4*data.timeSeries('MassBias').data())*(lambda_U234/lambda_U238)
                U234_U238_activity_MB_corr_timeseries = data.createTimeSeries('(U234)/(U238)_MB_corr', data.Output, indexChannel.time(), U234_U238_activity_MB_corr)

        #Calculating resulting ratios as ROI's
        def U238_Th232_activity_RSF_corr_AR(sel):
                result= Result()
                result = data.result(sel,data.timeSeries('U238_CPS')).value()/(data.result(sel,data.timeSeries('Th232_CPS')).value())/data.result(sel,data.timeSeries('RSF_fit')).value()*(lambda_U238/lambda_Th232)
                return result
        def U238_Th232_activity_RSF_corr_uncer_AR(sel):
                result= Result()
                result = data.associatedResult(sel,'(U238)/(Th232)_RSF_corr_ROI').value()*np.sqrt((data.result(sel,data.timeSeries('Th232_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('Th232_CPS')).value())**2+(data.result(sel,data.timeSeries('U238_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('U238_CPS')).value())**2+(2*data.result(sel,data.timeSeries('RSF_fit_uncer')).value()/data.result(sel,data.timeSeries('RSF_fit')).value())**2)
                return result

        data.registerAssociatedResult('(U238)/(Th232)_RSF_corr_ROI',U238_Th232_activity_RSF_corr_AR) 
        data.registerAssociatedResult('(U238)/(Th232)_RSF_corr_2SE_ROI',U238_Th232_activity_RSF_corr_uncer_AR)

        def Th230_Th232_activity_MB_corr_AR(sel):
                result= Result()
                result = data.result(sel,data.timeSeries('Th230_CPS_corr')).value()/(data.result(sel,data.timeSeries('Th232_CPS')).value()-data.result(sel,data.timeSeries('Th232_CPS')).value()*2*data.result(sel,data.timeSeries('MassBias')).value())*(lambda_Th230/lambda_Th232)
                return result
        def Th230_Th232_activity_MB_corr_uncer_AR(sel):
                result= Result()
                result = data.associatedResult(sel,'(Th230)/(Th232)_MB_corr_ROI').value()*np.sqrt((data.result(sel,data.timeSeries('Th230_CPS_corr')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('Th230_CPS_corr')).value())**2+(data.result(sel,data.timeSeries('Th232_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('Th232_CPS')).value())**2)
                return result

        data.registerAssociatedResult('(Th230)/(Th232)_MB_corr_ROI',Th230_Th232_activity_MB_corr_AR) 
        data.registerAssociatedResult('(Th230)/(Th232)_MB_corr_2SE_ROI',Th230_Th232_activity_MB_corr_uncer_AR) 

        if how235or238 == 'based on 238U':
                def Th230_U238_activity_RSF_corr_AR(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('Th230_CPS_corr')).value()/(data.result(sel,data.timeSeries('U238_CPS')).value())*data.result(sel,data.timeSeries('RSF_fit')).value()*(lambda_Th230/lambda_U238)
                        return result
                def Th230_U238_activity_RSF_corr_uncer_AR(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(Th230)/(U238)_RSF_corr_ROI').value()*np.sqrt((data.result(sel,data.timeSeries('Th230_CPS_corr')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('Th230_CPS_corr')).value())**2+(data.result(sel,data.timeSeries('U238_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('U238_CPS')).value())**2+(2*data.result(sel,data.timeSeries('RSF_fit_uncer')).value()/data.result(sel,data.timeSeries('RSF_fit')).value())**2)
                        return result
        if how235or238 == 'based on 235U*137.818':
                def Th230_U238_activity_RSF_corr_AR(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('Th230_CPS_corr')).value()/(data.result(sel,data.timeSeries('U235_CPS')).value()*U238_U235)*data.result(sel,data.timeSeries('RSF_fit')).value()*(lambda_Th230/lambda_U238)
                        return result
                def Th230_U238_activity_RSF_corr_uncer_AR(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(Th230)/(U238)_RSF_corr_ROI').value()*np.sqrt((data.result(sel,data.timeSeries('Th230_CPS_corr')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('Th230_CPS_corr')).value())**2+(data.result(sel,data.timeSeries('U235_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('U235_CPS')).value())**2+(2*data.result(sel,data.timeSeries('RSF_fit_uncer')).value()/data.result(sel,data.timeSeries('RSF_fit')).value())**2)
                        return result
               
        data.registerAssociatedResult('(Th230)/(U238)_RSF_corr_ROI',Th230_U238_activity_RSF_corr_AR) 
        data.registerAssociatedResult('(Th230)/(U238)_RSF_corr_2SE_ROI',Th230_U238_activity_RSF_corr_uncer_AR) 

        if U234secDis:
                def U234_U238_activity_MB_corr_AR(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('U234_CPS_corr')).value()/(data.result(sel,data.timeSeries('U238_CPS')).value()-data.result(sel,data.timeSeries('U238_CPS')).value()*4*data.result(sel,data.timeSeries('MassBias')).value())*(lambda_U234/lambda_U238)
                        return result
                def U234_U238_activity_MB_corr_uncer_AR(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(U234)/(U238)_MB_corr_ROI').value()*np.sqrt((data.result(sel,data.timeSeries('U234_CPS_corr')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('U234_CPS_corr')).value())**2+(data.result(sel,data.timeSeries('U238_CPS')).uncertaintyAs2SE()/data.result(sel,data.timeSeries('U238_CPS')).value())**2)
                        return result
               
                data.registerAssociatedResult('(U234)/(U238)_MB_corr_ROI',U234_U238_activity_MB_corr_AR) 
                data.registerAssociatedResult('(U234)/(U238)_MB_corr_2SE_ROI',U234_U238_activity_MB_corr_uncer_AR) 

                def sec_check_234(sel):
                        result = Result()                        
                        if how234_238 == 'Ratio of Intensities':
                                U234_U238 = data.associatedResult(sel,'(U234)/(U238)_MB_corr_ROI').value()
                                U234_U238_2SE = data.associatedResult(sel,'(U234)/(U238)_MB_corr_2SE_ROI').value()                
                        if how234_238 == 'Mean of Ratios':
                                U234_U238 = data.result(sel,data.timeSeries('(U234)/(U238)_MB_corr')).value()
                                U234_U238_2SE = data.result(sel,data.timeSeries('(U234)/(U238)_MB_corr')).uncertaintyAs2SE()                                  
                        if (U234_U238 - U234_U238_2SE) < 1 and (U234_U238 + U234_U238_2SE) > 1:
                                result = 1
                        else:
                                result = 0
                        
                        return result
                
                data.registerAssociatedResult('Check234_238',sec_check_234) 
        
        corr_Th230_percent = (data.timeSeries('Th230_CPS').data()-data.timeSeries('Th230_CPS_corr').data())/data.timeSeries('Th230_CPS').data()*100
        corr_Th230_percent_timeseries = data.createTimeSeries('Th230_corr_percent', data.Intermediate, indexChannel.time(), corr_Th230_percent)       

        PLOT_Evolution.clearGraphs()
        U238_Th232_showall = []
        Th230_Th232_showall = []
        Th230_U238_showall = []

        for i, gr in enumerate(Result_groups):
                U238_Th232_sel_gr = []
                U238_Th232_sel_2SE_gr = []
                Th230_Th232_sel_gr = []
                Th230_Th232_sel_2SE_gr = []                
                sg = data.selectionGroup(gr)
                for sel in sg.selections():
                        if PlotLc == True:
                                Lc_check_val = data.associatedResult(sel, 'Lc_check').value()
                        else:
                                Lc_check_val = 1
                        if U234secDis:
                                if mask234:
                                        Sec_check_val = data.associatedResult(sel, 'Check234_238').value()
                                else:
                                        Sec_check_val = 1
                                mask_val = Sec_check_val*Lc_check_val
                        else:
                                mask_val = Lc_check_val
                        mask = mask_val == 1

                        if mask:  # Only proceed if Lc_check is 1
                                U238_Th232_sel = data.result(sel, U238_Th232_activity_RSF_corr_timeseries).value()
                                U238_Th232_sel_2SE = data.result(sel, U238_Th232_activity_RSF_corr_timeseries).uncertaintyAs2SE()
                                Th230_Th232_sel = data.result(sel, Th230_Th232_activity_MB_corr_timeseries).value()
                                Th230_Th232_sel_2SE = data.result(sel, Th230_Th232_activity_MB_corr_timeseries).uncertaintyAs2SE()
                                Th230_U238_sel = data.result(sel, Th230_U238_activity_RSF_corr_timeseries).value()

                                U238_Th232_selAR = data.associatedResult(sel, '(U238)/(Th232)_RSF_corr_ROI').value()
                                U238_Th232_sel_2SEAR = data.associatedResult(sel, '(U238)/(Th232)_RSF_corr_2SE_ROI').value()
                                Th230_Th232_selAR = data.associatedResult(sel, '(Th230)/(Th232)_MB_corr_ROI').value()
                                Th230_Th232_sel_2SEAR = data.associatedResult(sel, '(Th230)/(Th232)_MB_corr_2SE_ROI').value()
                                Th230_U238_selAR = data.associatedResult(sel, '(Th230)/(U238)_RSF_corr_ROI').value()

                                Th230_U238_showall.append(Th230_U238_selAR)

                                if not np.isnan(U238_Th232_sel):

                                        if how238_232 == 'Ratio of Intensities':
                                                U238_Th232_sel_gr.append(U238_Th232_selAR)
                                                U238_Th232_sel_2SE_gr.append(U238_Th232_sel_2SEAR)
                                                U238_Th232_showall.append(U238_Th232_selAR)
                                        if how238_232 == 'Mean of Ratios':
                                                U238_Th232_sel_gr.append(U238_Th232_sel)
                                                U238_Th232_sel_2SE_gr.append(U238_Th232_sel_2SE)
                                                U238_Th232_showall.append(U238_Th232_sel)
                                        if how230_232 == 'Ratio of Intensities':
                                                Th230_Th232_sel_gr.append(Th230_Th232_selAR)
                                                Th230_Th232_sel_2SE_gr.append(Th230_Th232_sel_2SEAR)
                                                Th230_Th232_showall.append(Th230_Th232_selAR)
                                        if how230_232 == 'Mean of Ratios':
                                                Th230_Th232_sel_gr.append(Th230_Th232_sel)
                                                Th230_Th232_sel_2SE_gr.append(Th230_Th232_sel_2SE)
                                                Th230_Th232_showall.append(Th230_Th232_sel)

                U238_Th232_sel_gr = np.array(U238_Th232_sel_gr)
                U238_Th232_sel_2SE_gr = np.array(U238_Th232_sel_2SE_gr)
                Th230_Th232_sel_gr = np.array(Th230_Th232_sel_gr)
                Th230_Th232_sel_2SE_gr = np.array(Th230_Th232_sel_2SE_gr)                
                g = PLOT_Evolution.addGraph()
                g.setName(gr)
                g.setLineStyle('lsNone')
                g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[i % len(PLOT_COLORS)])
                g.setData(U238_Th232_sel_gr, Th230_Th232_sel_gr)

                for x, y, x_err, y_err in zip(U238_Th232_sel_gr, Th230_Th232_sel_gr, U238_Th232_sel_2SE_gr, Th230_Th232_sel_2SE_gr):
                        h_error = PLOT_Evolution.addGraph()
                        hERR_pen = QPen(QColor('dark grey'))
                        hERR_pen.setStyle(Qt.SolidLine)
                        h_error.setPen(hERR_pen)
                        h_error.setData([x - x_err, x + x_err], [y, y])
                        h_error.removeFromLegend() 

                        v_error = PLOT_Evolution.addGraph()
                        vERR_pen = QPen(QColor('dark grey'))
                        vERR_pen.setStyle(Qt.SolidLine)
                        v_error.setPen(vERR_pen)
                        v_error.setData([x, x], [y - y_err, y + y_err])
                        v_error.removeFromLegend() 
        U238_Th232_showall = np.array(U238_Th232_showall)
        Th230_Th232_showall = np.array(Th230_Th232_showall)
        Th230_U238_showall = np.array(Th230_U238_showall)
        Th230_U238_mean = np.mean(Th230_U238_showall)

        max_sec_equil = max([max(U238_Th232_showall),max(Th230_Th232_showall)])
        g = PLOT_Evolution.addGraph()
        g.setName("Secular Equilibrium")
        fit_pen = QPen(QColor('black'))
        fit_pen.setStyle(Qt.DashLine)
        g.pen = fit_pen
        g.setData(np.array([0,max_sec_equil]), np.array([0,max_sec_equil]))

        def fitLin(x, a, b):
                return a + b*x

        params, cov = curve_fit(fitLin, U238_Th232_showall, Th230_Th232_showall, ftol=1e-5)

        n_samples = len(U238_Th232_showall)

        if params[1]<1:
                fit_x = np.linspace((-params[0]/(params[1]-1)),max_sec_equil,50)
                fit_y = params[0] + params[1]*fit_x
                g = PLOT_Evolution.addGraph()
                g.setName(f"Global isochron")
                fit_pen_iso = QPen(QColor('grey'))
                fit_pen_iso.setStyle(Qt.DashDotLine)
                fit_pen_iso.setWidth(2)
                g.pen = fit_pen_iso
                g.setData(np.array(fit_x), np.array(fit_y))   

                age = -np.log(1-params[1])/lambda_Th230/1000
                ann6.visible = True
                ann6.text = f'''
                <p style="color:black;">
                <b>Average (230Th)/(238U):</b> {Th230_U238_mean:.4f}<br>
                <b>Global Isochron age (n={n_samples}):</b> {age:.2f} ka</p>'''
        else:
                ann6.visible = True
                ann6.text = f'''
                <p style="color:black;">
                <b>Average (230Th)/(238U):</b> {Th230_U238_mean:.4f}<br>
                <b>Global Isochron age (n={n_samples}):</b> secular equilibrium</p>'''

        PLOT_Evolution.left().label = '(230Th)/(232Th)'
        PLOT_Evolution.bottom().label = '(238U)/(232Th)'
        PLOT_Evolution.setLegendVisible(True)
        PLOT_Evolution.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
        PLOT_Evolution.setToolsVisible(True)
        PLOT_Evolution.rescaleAxes()
        PLOT_Evolution.replot()

##%     Calculating Secular Equilibrium U234 and U238

        if U234secDis:                 
                PLOT_234_238.clearGraphs()
                PLOT_234_238.show()

                meas_main = []
                start_times_all = []

                for i, gr in enumerate(U234secDis_group):
                        meas = []
                        meas_2SE = []
                        start_times = []
                        sg = data.selectionGroup(gr)
                        for sel in sg.selections():
                                if PlotLc == True:
                                        Lc_check_val = data.associatedResult(sel, 'Lc_check').value()
                                else:
                                        Lc_check_val = 1
                                mask = Lc_check_val == 1
                                
                                if mask:
                                        if how234_238 == 'Ratio of Intensities':
                                                ratio = data.associatedResult(sel,'(U234)/(U238)_MB_corr_ROI').value()
                                                ratio_2SE = data.associatedResult(sel,'(U234)/(U238)_MB_corr_2SE_ROI').value()
                                        if how234_238 == 'Mean of Ratios':
                                                ratio = data.result(sel, U234_U238_activity_MB_corr_timeseries).value()
                                                ratio_2SE = data.result(sel,U234_U238_activity_MB_corr_timeseries).uncertaintyAs2SE()                                                
                                        start_time = indexChannel.timeForSelection(sel)
                                        if not np.isnan(ratio):  # Check if the ratio and time are valid
                                                meas.append(ratio)
                                                meas_2SE.append(ratio_2SE)
                                                meas_main.append(ratio)
                                                start_times.append(start_time[0])
                                                start_times_all.append(start_time[0])  # Collect the corresponding time value
                        # reject outliers for plot
                        meas = np.array(meas)
                        meas_2SE = np.array(meas_2SE)
                        start_times = np.array(start_times)
                        mask = (meas >= 0)
                        filtered_meas = meas[mask]
                        filtered_meas_2SE = meas_2SE[mask]
                        filtered_start_times = start_times[mask]

                        filtered_start_times = filtered_start_times-indexChannel.time()[0]
                        g = PLOT_234_238.addGraph()
                        g.setName(gr)
                        g.setLineStyle('lsNone')
                        g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[i % len(PLOT_COLORS)])
                        g.setData(np.array(filtered_start_times), np.array(filtered_meas))                

                        for x, y, y_err in zip(filtered_start_times, filtered_meas, filtered_meas_2SE):
                                v_error = PLOT_234_238.addGraph()
                                vERR_pen = QPen(QColor('light grey'))
                                vERR_pen.setStyle(Qt.SolidLine)
                                v_error.setPen(vERR_pen)
                                v_error.setData([x, x], [y - y_err, y + y_err])
                                v_error.removeFromLegend()  # Remove vertical error bar from legend
                
                mean_234_238 = np.mean(meas)
                g = PLOT_234_238.addGraph()
                g.setName("Secular Equilibrium")
                fit_pen = QPen(QColor('black'))
                fit_pen.setStyle(Qt.DashLine)
                g.pen = fit_pen
                g.setData(np.array([(min(start_times_all)-indexChannel.time()[0]),(max(start_times_all)-indexChannel.time()[0])]), np.array([1,1]))

                f = PLOT_234_238.addGraph()
                f.setName("Average")
                fit_pen = QPen(QColor('black'))
                fit_pen.setStyle(Qt.DashLine)
                fit_pen.setWidth(2)
                f.pen = fit_pen
                f.setData(np.array([(min(start_times_all)-indexChannel.time()[0]),(max(start_times_all)-indexChannel.time()[0])]), np.array([mean_234_238,mean_234_238]))

                ann9.visible = True
                ann9.text = f'''
                <p style="color:black;">
                <b>Average (234U)/(238U):</b> {mean_234_238:.4f}</p>'''

                PLOT_234_238.left().label = '(U234)/(U238)'
                PLOT_234_238.bottom().label = 'Time'
                PLOT_234_238.setLegendVisible(True)
                PLOT_234_238.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
                PLOT_234_238.setToolsVisible(True)
                PLOT_234_238.rescaleAxes()
                PLOT_234_238.replot()
        else:
                PLOT_234_238.hide()

##%     ppm U238 and Th232

        if calc_ppm: 
                PLOT_ppm.clearGraphs()
                PLOT_ppm.show()
                U_RM_ppm_NIST = data.referenceMaterialData(ppm_on)["U"].value()
                Th_RM_ppm_NIST = data.referenceMaterialData(ppm_on)["Th"].value()

                Th232_cps_NIST = []
                Th232_cps_2SE_NIST = []
                U238_cps_NIST = []
                U238_cps_2SE_NIST = []
                start_time_NIST = []

                nist = data.selectionGroup(ppm_on)
                for sel_nist in nist.selections():
                        Th232_cps = data.result(sel_nist,data.timeSeries('Th232_CPS')).value() 
                        Th232_cps_2SE = data.result(sel_nist,data.timeSeries('Th232_CPS')).uncertaintyAs2SE() 
                        U238_cps = data.result(sel_nist,data.timeSeries('U238_CPS')).value() 
                        U238_cps_2SE = data.result(sel_nist,data.timeSeries('U238_CPS')).uncertaintyAs2SE() 
                        start_time_sel = data.timeSeries('Th232_CPS').timeForSelection(sel_nist)

                        Th232_cps_NIST.append(Th232_cps)
                        Th232_cps_2SE_NIST.append(Th232_cps_2SE)
                        U238_cps_NIST.append(U238_cps)
                        U238_cps_2SE_NIST.append(U238_cps_2SE)
                        start_time_NIST.append(start_time_sel[0])

                Th232_cps_NIST = np.array(Th232_cps_NIST)
                Th232_cps_2SE_NIST = np.array(Th232_cps_2SE_NIST)
                U238_cps_NIST = np.array(U238_cps_NIST)
                U238_cps_2SE_NIST = np.array(U238_cps_2SE_NIST)                
                start_time_NIST = start_time_NIST-indexChannel.time()[0]
                start_time_NIST = np.array(start_time_NIST)    

                fit_x = np.linspace(
                        start_time_NIST.min(),
                        start_time_NIST.max(),
                        50
                        )

                for i,element,NIST_conc in zip([Th232_cps_NIST,U238_cps_NIST],['Th232','U238'],[Th_RM_ppm_NIST,U_RM_ppm_NIST]):

                        if ppm_fit == 'Constant':
                                def fitConst(x, a, b):
                                        return a + b*0

                                params, cov = curve_fit(fitConst, start_time_NIST, i,  ftol=1e-5)
                                params = np.mean(i)
                                cov = np.std(i)/np.sqrt(len(i))

                                fit_y = params*np.ones(len(fit_x))
                                ppm_fit_NIST = (params+0*(indexChannel.time()-min(indexChannel.time())))
                                ppm_fit_NIST_timeSeries = data.createTimeSeries(f'ppm_fit_{element}', data.Intermediate, indexChannel.time(), ppm_fit_NIST)

                                ppm_fit_NIST_uncer = cov * np.ones(len(indexChannel.time()))
                                ppm_fit_NIST_uncer_timeseries = data.createTimeSeries(f'ppm_fit_{element}_uncer', data.Intermediate, indexChannel.time(), ppm_fit_NIST_uncer)

                                ppm_element = data.timeSeries(f'{element}_CPS').data() * NIST_conc / data.timeSeries(f'ppm_fit_{element}').data()
                                ppm_element_timeseries = data.createTimeSeries(f'{element}_ppm',data.Intermediate,indexChannel.time(), ppm_element)

                        if ppm_fit == 'Linear':
                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_time_NIST, i, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                ppm_fit_NIST = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))
                                ppm_fit_NIST_timeSeries = data.createTimeSeries(f'ppm_fit_{element}', data.Intermediate, indexChannel.time(), ppm_fit_NIST)

                                ppm_fit_NIST_uncer = np.sqrt(cov[0, 0] + (indexChannel.time() - min(indexChannel.time()))**2 * cov[1, 1] + 
                                        2 * (indexChannel.time() - min(indexChannel.time())) * cov[0, 1])
                                ppm_fit_NIST_uncer_timeseries = data.createTimeSeries(f'ppm_fit_{element}_uncer', data.Intermediate, indexChannel.time(), ppm_fit_NIST_uncer)

                                ppm_element = data.timeSeries(f'{element}_CPS').data() * NIST_conc / data.timeSeries(f'ppm_fit_{element}').data()
                                ppm_element_timeseries = data.createTimeSeries(f'{element}_ppm',data.Intermediate,indexChannel.time(), ppm_element) 

                        if ppm_fit == 'Polynomial':
                                def fitLin(x, a, b, c):
                                        return a + b*x + c*x*x

                                params, cov = curve_fit(fitLin, start_time_NIST, i, ftol=1e-5)

                                fit_y = params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                                ppm_fit_NIST = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2)
                                ppm_fit_NIST_timeSeries = data.createTimeSeries(f'ppm_fit_{element}', data.Intermediate, indexChannel.time(), ppm_fit_NIST)


                                x_vals = indexChannel.time() - min(indexChannel.time())
                                ppm_fit_NIST_uncer = np.sqrt(
                                        cov[0, 0]
                                        + x_vals ** 2 * cov[1, 1]
                                        + x_vals ** 4 * cov[2, 2]
                                        + 2 * x_vals * cov[0, 1]
                                        + 2 * x_vals ** 2 * cov[0, 2]
                                        + 2 * x_vals ** 3 * cov[1, 2]
                                )
                                ppm_fit_NIST_uncer_timeseries = data.createTimeSeries(
                                        f'ppm_fit_{element}_uncer', data.Intermediate, indexChannel.time(), ppm_fit_NIST_uncer
                                )

                                ppm_element = data.timeSeries(f'{element}_CPS').data() * NIST_conc / data.timeSeries(f'ppm_fit_{element}').data()
                                ppm_element_timeseries = data.createTimeSeries(f'{element}_ppm',data.Intermediate,indexChannel.time(), ppm_element)
                        if ppm_fit == 'Exponential':
                                def fitExp(x, a, b, c):
                                        return a + b * np.exp(-c * x)
                                initial_guess = [np.mean(i), 0, 0.01]
                                try:
                                        params, cov = curve_fit(fitExp, start_time_NIST, i, p0=initial_guess)
                                        fit_y = params[0] + params[1] * np.exp(-params[2] * fit_x)

                                        ppm_fit_NIST = (params[0]+params[1]*np.exp((-indexChannel.time()-min(indexChannel.time()))*params[2]))
                                        ppm_fit_NIST_timeSeries = data.createTimeSeries(f'ppm_fit_{element}', data.Intermediate, indexChannel.time(), ppm_fit_NIST)

                                        x_vals = indexChannel.time() - min(indexChannel.time())

                                        residuals = i - fitExp(start_time_NIST, *params)
                                        residual_std = np.std(residuals)
                                        ppm_fit_NIST_uncer = np.full_like(indexChannel.time(), residual_std)

                                        ppm_fit_NIST_uncer_timeseries = data.createTimeSeries(
                                        f'ppm_fit_{element}_uncer', data.Intermediate, indexChannel.time(), ppm_fit_NIST_uncer
                                        )

                                        ppm_element = data.timeSeries(f'{element}_CPS').data() * NIST_conc / data.timeSeries(f'ppm_fit_{element}').data()
                                        ppm_element_timeseries = data.createTimeSeries(f'{element}_ppm',data.Intermediate,indexChannel.time(), ppm_element)

                                except RuntimeError as e:
                                        if 'Optimal parameters not found: Number of calls' in str(e):

                                                IoLog.warning('Could not find fit to data with Exp model. Switching to Linear model.')

                                                settings["MB_fit"] = 'Linear'

                                                def fitLin(x, a, b):
                                                        return a + b*x

                                                params, cov = curve_fit(fitLin, start_time_NIST, i, ftol=1e-5)
                                                
                                                fit_y = params[1]*fit_x + params[0]
                                                ppm_fit_NIST = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))
                                                ppm_fit_NIST_timeSeries = data.createTimeSeries(f'ppm_fit_{element}', data.Intermediate, indexChannel.time(), ppm_fit_NIST)
                                                ppm_element = data.timeSeries(f'{element}_CPS').data() * NIST_conc / data.timeSeries(f'ppm_fit_{element}').data()
                                                ppm_element_timeseries = data.createTimeSeries(f'{element}_ppm',data.Intermediate,indexChannel.time(), ppm_element)                        

                                        else:
                                                raise
                        if ppm_fit == 'Polynomial3':
                                def fitLin(x, a, b, c, d):
                                        return a + b*x + c*x*x + d*x*x*x

                                params, cov = curve_fit(fitLin, start_time_NIST, i, ftol=1e-5)

                                fit_y = params[3]*fit_x*fit_x*fit_x + params[2]*fit_x*fit_x + params[1]*fit_x + params[0]
                                ppm_fit_NIST = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time()))+params[2]*(indexChannel.time()-min(indexChannel.time()))**2+params[3]*(indexChannel.time()-min(indexChannel.time()))**3)
                                ppm_fit_NIST_timeSeries = data.createTimeSeries(f'ppm_fit_{element}', data.Intermediate, indexChannel.time(), ppm_fit_NIST)


                                x_vals = indexChannel.time() - min(indexChannel.time())
                                x1 = x_vals
                                x2 = x_vals**2
                                x3 = x_vals**3

                                ppm_fit_NIST_uncer = np.sqrt(
                                cov[0, 0] +
                                cov[1, 1] * x1**2 +
                                cov[2, 2] * x2**2 +
                                cov[3, 3] * x3**2 +
                                2 * cov[0, 1] * x1 +
                                2 * cov[0, 2] * x2 +
                                2 * cov[0, 3] * x3 +
                                2 * cov[1, 2] * x1 * x2 +
                                2 * cov[1, 3] * x1 * x3 +
                                2 * cov[2, 3] * x2 * x3
                                )
                                ppm_fit_NIST_uncer_timeseries = data.createTimeSeries(
                                        f'ppm_fit_{element}_uncer', data.Intermediate, indexChannel.time(), ppm_fit_NIST_uncer
                                )

                                ppm_element = data.timeSeries(f'{element}_CPS').data() * NIST_conc / data.timeSeries(f'ppm_fit_{element}').data()
                                ppm_element_timeseries = data.createTimeSeries(f'{element}_ppm',data.Intermediate,indexChannel.time(), ppm_element)

                        g = PLOT_ppm.addGraph()
                        g.setName(f"Fit {element}")
                        fit_pen = QPen(QColor('grey'))
                        fit_pen.setStyle(Qt.DashDotLine)
                        fit_pen.setWidth(2)
                        g.pen = fit_pen
                        g.setData(np.array(fit_x), np.array(fit_y))  

                g = PLOT_ppm.addGraph()
                g.setName(f'Th232 CPS {ppm_on}')
                g.setLineStyle('lsNone')
                g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[1 % len(PLOT_COLORS)])
                g.setData(start_time_NIST, Th232_cps_NIST)
                p = PLOT_ppm.addGraph()
                p.setName(f'U238 CPS {ppm_on}')
                p.setLineStyle('lsNone')
                p.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[2 % len(PLOT_COLORS)])
                p.setData(start_time_NIST, U238_cps_NIST)    

                for x, y, y_err in zip(start_time_NIST,Th232_cps_NIST,Th232_cps_2SE_NIST):
                        v_error = PLOT_ppm.addGraph()
                        vERR_pen = QPen(QColor('dark grey'))
                        vERR_pen.setStyle(Qt.SolidLine)
                        v_error.setPen(vERR_pen)
                        v_error.setData([x, x], [y - y_err, y + y_err])
                        v_error.removeFromLegend() 
                for x, y, y_err in zip(start_time_NIST,U238_cps_NIST,U238_cps_2SE_NIST):
                        v_error = PLOT_ppm.addGraph()
                        vERR_pen = QPen(QColor('dark grey'))
                        vERR_pen.setStyle(Qt.SolidLine)
                        v_error.setPen(vERR_pen)
                        v_error.setData([x, x], [y - y_err, y + y_err])
                        v_error.removeFromLegend() 

                PLOT_ppm.left().label = 'CPS'
                PLOT_ppm.bottom().label = 'Time since start of session (s)'
                PLOT_ppm.setLegendVisible(True)
                PLOT_ppm.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
                PLOT_ppm.setToolsVisible(True)
                PLOT_ppm.rescaleAxes()
                PLOT_ppm.replot()
        else:
                PLOT_ppm.hide()
                                   
##%     Fill empty labels      
        if fillLabels == True:
                BaselineGroup = data.selectionGroupNames(data.Baseline)
                duration_BG_ms = []
                Th230_BG = []
                start_time_BG = []

                # Step 1: Gather background data
                for gr in BaselineGroup:
                        bg = data.selectionGroup(gr)
                        for sel_bg in bg.selections():
                                Th230_BG_sel = data.result(sel_bg, data.timeSeries('Th230')).value()
                                start_time_BG_sel = data.timeSeries('Th230').timeForSelection(sel_bg)
                                duration_bg_ms_sel = len(start_time_BG_sel) * dwell_Th230_s

                                duration_BG_ms.append(duration_bg_ms_sel)
                                Th230_BG.append(Th230_BG_sel)
                                start_time_BG.append(start_time_BG_sel[0])

                # Combine background data into a sorted list by start time
                backgrounds = sorted(zip(start_time_BG, Th230_BG, duration_BG_ms), key=lambda x: x[0])

                # Step 2: Iterate through combined selection groups
                combined_names = data.selectionGroupNames(data.ReferenceMaterial) + data.selectionGroupNames(data.Sample)
                unmatched_selections = []  # Track selections without names

                # Step 3: Assign names based on matching background and subsequent selection
                for combgr in combined_names:
                        combgrsg = data.selectionGroup(combgr)
                        for sel in combgrsg.selections():
                                if not sel.property('Rep Rate'):
                                        sel.name = ''
                                if sel.name == "Linked selection":
                                        sel.name = ''
                                if sel.name:  # Skip if already named
                                        continue
                                
                                # Get the start time of the current selection
                                sel_start_time = data.timeSeries('Th230').timeForSelection(sel)[0]
                                # Match the closest background based on start time
                                matched_bg = None
                                for bg_start, bg_Th230, bg_duration in backgrounds:
                                        if bg_start < sel_start_time:
                                                matched_bg = (bg_start, bg_Th230, bg_duration)
                                        else:
                                                break  # Stop when background start time exceeds selection start time

                                if matched_bg:
                                        bg_start, bg_Th230, bg_duration = matched_bg

                                        # Find the selection with a start time slightly later than the matched background
                                        matched_name = None
                                        check_group = data.selectionGroup(combgr)
                                        for sel_check in check_group.selections():
                                                check_start_time = data.timeSeries('Th230').timeForSelection(sel_check)[0]
                                                if check_start_time > bg_start and sel_check.name:  # Look for the first valid name
                                                        matched_name = sel_check.name
                                                        break
                                                if matched_name:
                                                        break

                                        # Assign the matched name to the current selection
                                        if matched_name:
                                                sel.name = matched_name
                                        else:
                                                sel.name = "Unmatched Selection" 
                                else:
                                        sel.name = "Unmatched Selection" 
                                

##%     Final output at the end of the output file:                                        
        if how238_232 == 'Ratio of Intensities':
                def U238_Th232_final(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(U238)/(Th232)_RSF_corr_ROI').value()
                        return result
                def U238_Th232_final_uncer(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(U238)/(Th232)_RSF_corr_2SE_ROI').value()/2
                        return result
        if how238_232 == 'Mean of Ratios':
                def U238_Th232_final(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('(U238)/(Th232)_RSF_corr')).value()
                        return result
                def U238_Th232_final_uncer(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('(U238)/(Th232)_RSF_corr')).uncertaintyAs2SE()/2
                        return result
        data.registerAssociatedResult('final_(U238)/(Th232)',U238_Th232_final) 
        data.registerAssociatedResult('final_(U238)/(Th232)_1sigma',U238_Th232_final_uncer) 
        
        if how230_232 == 'Ratio of Intensities':
                def Th230_Th232_final(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(Th230)/(Th232)_MB_corr_ROI').value()
                        return result
                def Th230_Th232_final_uncer(sel):
                        result= Result()
                        result = data.associatedResult(sel,'(Th230)/(Th232)_MB_corr_2SE_ROI').value()/2
                        return result
        if how230_232 == 'Mean of Ratios':
                def Th230_Th232_final(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('(Th230)/(Th232)_MB_corr')).value()
                        return result
                def Th230_Th232_final_uncer(sel):
                        result= Result()
                        result = data.result(sel,data.timeSeries('(Th230)/(Th232)_MB_corr')).uncertaintyAs2SE()/2
                        return result          
                     
        data.registerAssociatedResult('final_(Th230)/(Th232)',Th230_Th232_final) 
        data.registerAssociatedResult('final_(Th230)/(Th232)_1sigma',Th230_Th232_final_uncer) 

##%     Setup index time
        drs.message("Finished!")
        drs.progress(100)
        drs.finished()

def settingsWidget():
        widget = QtGui.QWidget()
        formLayout = QtGui.QFormLayout()
        widget.setLayout(formLayout)

        timeSeriesNames = data.timeSeriesNames(data.Input)
        defaultChannelName = ""
        if timeSeriesNames:
                defaultChannelName = timeSeriesNames[0]

        RMNames = data.selectionGroupNames(data.ReferenceMaterial)
        SampleNames = data.selectionGroupNames(data.Sample)
        combined_names = RMNames + SampleNames

        defaultMonazitelName = next((ch for ch in combined_names if "onazit" in ch.lower()), None)
        defaultRSFName = next((ch for ch in RMNames if "91500" in ch.lower()), None)

        drs.setSetting("IndexChannel", 'TotalBeam')
        drs.setSetting("MaskChannel", 'TotalBeam')
        drs.setSetting("MaskCutoff", 10000.0)
        drs.setSetting("MaskTrim", 0.0)
        drs.setSetting("Monazite", defaultMonazitelName)
        drs.setSetting("Ab230_fit",'Linear')
        drs.setSetting("Ab228_fit",'Linear')
        drs.setSetting('Zr2O3corr',True)
        drs.setSetting("Abundance228Items", None)
        drs.setSetting("RSF_on", defaultRSFName)
        drs.setSetting("RSF_plot",None)
        drs.setSetting("RSFOutlier",True)
        drs.setSetting("MB_groups", None)
        drs.setSetting("MB_fit", 'Linear')
        drs.setSetting("Result_groups", None)
        drs.setSetting("RSF_how", 'based on RM U & Th concentration')
        drs.setSetting("plotLc1",True)
        drs.setSetting("RSF_fit",'Linear')
        drs.setSetting("Corr228_fit",'Linear')
        drs.setSetting("InfCorrOutlier",True)
        drs.setSetting("U235orU238",'based on 235U*137.818')
        drs.setSetting("How238_232",'Mean of Ratios')
        drs.setSetting("How230_232",'Ratio of Intensities')
        drs.setSetting("fill_Labels",False)
        drs.setSetting("calculate_ppm",False)
        drs.setSetting("ppm_fit",'Linear')
        drs.setSetting("ppm_calculated_on",defaultRSFName)
        drs.setSetting("calculate_U234_U238",False)
        drs.setSetting("U234-U238_groups",None)
        drs.setSetting("mask_U234_U238",False)
        drs.setSetting("How234_238",'Mean of Ratios')
        drs.setSetting("Corr228_intext", 'Zr90 measurements in this file')
        drs.setSetting("Corr228_external", None)
        drs.setSetting("Corr228_external_SE", None)

        settings = drs.settings()

##%     Set Title
        Title_Label = QtGui.QLabel("<h1>U-Th disequillibrium dating DRS<\h1>")
        formLayout.addRow(Title_Label)
        Subtitle_Label = QtGui.QLabel("Appropriate statistics: Session - Preferences - Stats - Baseline/Normal stats calculation: Mean no outlier reject")
        formLayout.addRow(Subtitle_Label)
        Subsubtitle_Label = QtGui.QLabel("Suggested Baseline Spline: StepForward")
        formLayout.addRow(Subsubtitle_Label)

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

##%     Define Index for Baselinesubtract
        Baseline_Label = QtGui.QLabel("<h2>Baseline Subtract<\h2>")
        formLayout.addRow(Baseline_Label)

        indexComboBox = QtGui.QComboBox(widget)
        indexComboBox.addItems(timeSeriesNames)
        indexComboBox.setCurrentText(settings["IndexChannel"])
        indexComboBox.textActivated.connect(lambda t: drs.setSetting("IndexChannel", t))
        formLayout.addRow("Index channel", indexComboBox)
        
        maskComboBox = QtGui.QComboBox(widget)
        maskComboBox.addItems(data.timeSeriesNames(data.Input))
        maskComboBox.setCurrentText(settings["MaskChannel"])
        maskComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MaskChannel", t))
        formLayout.addRow("Mask channel", maskComboBox)

        maskLineEdit = QtGui.QLineEdit(widget)
        maskLineEdit.setText(settings["MaskCutoff"])
        maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))
        formLayout.addRow("Mask cutoff", maskLineEdit)

        maskTrimLineEdit = QtGui.QLineEdit(widget)
        maskTrimLineEdit.setText(settings["MaskTrim"])
        maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))
        formLayout.addRow("Mask trim", maskTrimLineEdit)

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

##%     Choosing Group for Abundance Correction of 232 on 230
        Abundance230_Label = QtGui.QLabel("<h2>Abundance Correction of mass 232 on mass 230<\h2>")
        formLayout.addRow(Abundance230_Label)

        Abundance230ComboBox = QtGui.QComboBox(widget)
        SampleNames = data.selectionGroupNames(data.Sample)
        Abundance230ComboBox.addItems(SampleNames)
        Abundance230ComboBox.setCurrentText(settings["Monazite"])
        Abundance230ComboBox.textActivated.connect(lambda t: drs.setSetting("Monazite", t))
        formLayout.addRow("Select Monazite", Abundance230ComboBox)

        def updateRMCombo():
                SampleNames = data.selectionGroupNames(data.Sample)
                Abundance230ComboBox.clear()
                Abundance230ComboBox.addItems(SampleNames)

        data.selectionGroupsChanged.connect(updateRMCombo)

        Ab230_fit_ComboBox = QtGui.QComboBox(widget)
        Ab230_fit_ComboBox.addItems(['Constant','Linear', 'Logistic','Polynomial'])
        Ab230_fit_ComboBox.setCurrentText(settings["Ab230_fit"])
        Ab230_fit_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("Ab230_fit", t))
        formLayout.addRow("Fit of Abundance Sensitivity", Ab230_fit_ComboBox)

        formLayout.addRow('Abundance sensitivity of \nMonazite on mass 230', PLOT_Mon230)

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

##%     Choosing Groups for Abundance correction of 228 on 230
        Abundance228_Label = QtGui.QLabel("<h2>Interference Correction of mass 228 on mass 230 (only needed for zircons due to Zr-Oxide)<\h2>")
        formLayout.addRow(Abundance228_Label)

        ZrcorrCheckBox = QtGui.QCheckBox(widget)
        ZrcorrCheckBox.setChecked(settings["Zr2O3corr"])
        ZrcorrCheckBox.toggled.connect(lambda t: drs.setSetting("Zr2O3corr", bool(t)))
        formLayout.addRow("Apply Correction", ZrcorrCheckBox)

        Corr228_intext_ComboBox = QtGui.QComboBox(widget)
        Corr228_intext_ComboBox.addItems(['Zr90 measurements in this file','Th230 directly measured in Zirconblank','Add number from external file'])
        Corr228_intext_ComboBox.setCurrentText(settings["Corr228_intext"])
        Corr228_intext_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("Corr228_intext", t))
        formLayout.addRow("How to correct?", Corr228_intext_ComboBox)

        Corr228LineEdit = QtGui.QLineEdit(widget)
        Corr228LineEdit.setText(settings["Corr228_external"])
        Corr228LineEdit.textChanged.connect(lambda t: drs.setSetting("Corr228_external", float(t)))
        formLayout.addRow("Interference of 228 on 230 (only if external value)", Corr228LineEdit)

        Corr228SELineEdit = QtGui.QLineEdit(widget)
        Corr228SELineEdit.setText(settings["Corr228_external_SE"])
        Corr228SELineEdit.textChanged.connect(lambda t: drs.setSetting("Corr228_external_SE", float(t)))
        formLayout.addRow("SE of interference of 228 on 230 (only if external value)", Corr228SELineEdit)

        Ab228_fit_ComboBox = QtGui.QComboBox(widget)
        Ab228_fit_ComboBox.addItems(['Constant','Linear','Logistic','Polynomial'])
        Ab228_fit_ComboBox.setCurrentText(settings["Ab228_fit"])
        Ab228_fit_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("Ab228_fit", t))
        formLayout.addRow("Fit of Abundance Sensitivity", Ab228_fit_ComboBox)

        formLayout.addRow('Abundance sensitivity of \nMonazite on mass 228', PLOT_Mon228)

        Abundance228ExtButton = QtGui.QToolButton(widget)
        CombinedMenu = CheckableMenu(Abundance228ExtButton, combined_names)

        def updateSelected228Groups(selected):
                if selected:
                        selectedGroups228Label.setText(f"Selected Groups: {', '.join(selected)}")
                else:
                        selectedGroups228Label.setText("Selected Groups: None")

        CombinedMenu.itemsChanged.connect(lambda selected: [drs.setSetting("Abundance228Items", selected), updateSelected228Groups(selected)])
        Abundance228ExtButton.setMenu(CombinedMenu)
        Abundance228ExtButton.setPopupMode(QtGui.QToolButton.InstantPopup)
        Abundance228ExtButton.setIcon(CUI().icon('trophy'))
        Abundance228ExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        formLayout.addRow("Select Zircons", Abundance228ExtButton)
        selectedGroups228Label = QtGui.QLabel("Selected Groups: None")
        formLayout.addRow(selectedGroups228Label)

        Corr228_fit_ComboBox = QtGui.QComboBox(widget)
        Corr228_fit_ComboBox.addItems(['Constant','Linear','Logistic','Polynomial'])
        Corr228_fit_ComboBox.setCurrentText(settings["Corr228_fit"])
        Corr228_fit_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("Corr228_fit", t))
        formLayout.addRow("Fit of Inferrence Correction", Corr228_fit_ComboBox)

        InfCorrOutlierCheckBox = QtGui.QCheckBox(widget)
        InfCorrOutlierCheckBox.setChecked(settings["InfCorrOutlier"])
        InfCorrOutlierCheckBox.toggled.connect(lambda t: drs.setSetting("InfCorrOutlier", bool(t)))
        formLayout.addRow("Outlier Rejection", InfCorrOutlierCheckBox)

        formLayout.addRow('Interference correction on \nmass 230 due to Zr2O3 (mass228)', PLOT_Zr2O3)

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

        def toggleZrcorrWidgets(checked):
                Corr228_intext_ComboBox.setVisible(checked)
                Corr228LineEdit.setVisible(checked)
                Corr228SELineEdit.setVisible(checked)
                Ab228_fit_ComboBox.setVisible(checked)
                PLOT_Mon228.setVisible(checked)
                Abundance228ExtButton.setVisible(checked)
                Corr228_fit_ComboBox.setVisible(checked)
                InfCorrOutlierCheckBox.setVisible(checked)
                PLOT_Zr2O3.setVisible(checked)
        ZrcorrCheckBox.toggled.connect(toggleZrcorrWidgets)
        toggleZrcorrWidgets(ZrcorrCheckBox.isChecked())

##%     Choosing Group for the Relative Sensitivity Factor
        RSF_Label = QtGui.QLabel("<h2>Relative Sensitivity Factor<\h2>")
        formLayout.addRow(RSF_Label)

        HowRSFComboBox = QtGui.QComboBox(widget)
        HowRSFComboBox.addItems(['based on RM U & Th concentration','based on secular equilibrium'])
        HowRSFComboBox.setCurrentText(settings["RSF_how"])
        HowRSFComboBox.textActivated.connect(lambda t: drs.setSetting("RSF_how", t))
        formLayout.addRow("How to calculate RSF", HowRSFComboBox)

        def updatehowRMCombo():
                HowRSFComboBox.clear()
                HowRSFComboBox.addItems(['based on RM U & Th concentration','based on secular equilibrium'])

        data.selectionGroupsChanged.connect(updatehowRMCombo)

        If_plot_Label = QtGui.QLabel("If based on RM U & Th concentration, choose RM to chompare (U, Th need to be available in RM file). If based on secular equilibrium, choose samples in secular equilibrium.")
        formLayout.addRow(If_plot_Label)

        RSFPlotExtButton = QtGui.QToolButton(widget)
        Menu = CheckableMenu(RSFPlotExtButton, combined_names)
        Menu.itemsChanged.connect(lambda selected: [drs.setSetting("RSF_plot", selected), updateSelectedRSFPlotGroups(selected)])
        RSFPlotExtButton.setMenu(Menu)
        RSFPlotExtButton.setPopupMode(QtGui.QToolButton.InstantPopup)
        RSFPlotExtButton.setIcon(CUI().icon('trophy'))
        RSFPlotExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        formLayout.addRow("Which materials to plot", RSFPlotExtButton)

        selectedGroupsRSFPlotLabel = QtGui.QLabel("Selected Groups: None")
        formLayout.addRow(selectedGroupsRSFPlotLabel)

        If_Label = QtGui.QLabel("If based on Matrix Match, choose RM based on which you want to calculate RSF")
        formLayout.addRow(If_Label)

        def updateSelectedRSFPlotGroups(selected):
                if selected:
                        selectedGroupsRSFPlotLabel.setText(f"Selected Groups: {', '.join(selected)}")
                else:
                        selectedGroupsRSFPlotLabel.setText("Selected Groups: None")

        RSFComboBox = QtGui.QComboBox(widget)
        RSFComboBox.addItems(RMNames)
        RSFComboBox.setCurrentText(settings["RSF_on"])
        RSFComboBox.textActivated.connect(lambda t: drs.setSetting("RSF_on", t))
        formLayout.addRow("Which RM material to use", RSFComboBox)

        def updateRMCombo():
                RSFComboBox.clear()
                RSFComboBox.addItems(RMNames)

        data.selectionGroupsChanged.connect(updateRMCombo)

        RSF_fit_ComboBox = QtGui.QComboBox(widget)
        RSF_fit_ComboBox.addItems(['Constant','Linear','Logistic','Polynomial', 'Polynomial3'])
        RSF_fit_ComboBox.setCurrentText(settings["RSF_fit"])
        RSF_fit_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("RSF_fit", t))
        formLayout.addRow("Fit of RSF", RSF_fit_ComboBox)

        RSFOutlierCheckBox = QtGui.QCheckBox(widget)
        RSFOutlierCheckBox.setChecked(settings["RSFOutlier"])
        RSFOutlierCheckBox.toggled.connect(lambda t: drs.setSetting("RSFOutlier", bool(t)))
        formLayout.addRow("Outlier Rejection", RSFOutlierCheckBox)

        formLayout.addRow('U/Th RSF', PLOT_RSF)

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

##%     Choosing groups and degree of fit for Mass Bias correction
        MB_Label = QtGui.QLabel("<h2>Mass Bias correction<\h2>")
        formLayout.addRow(MB_Label)
        MB_explain_Label = QtGui.QLabel("Assuming linear relationship between masses, calculated based on 238U/235U = 137.818, and applied on 230Th/232Th.")
        formLayout.addRow(MB_explain_Label)

        def updateSelectedMBGroups(selected):
                if selected:
                        selectedGroupsMBLabel.setText(f"Selected Groups: {', '.join(selected)}")
                else:
                        selectedGroupsMBLabel.setText("Selected Groups: None")

        MBExtButton = QtGui.QToolButton(widget)
        CombinedMenu = CheckableMenu(MBExtButton, combined_names)
        CombinedMenu.itemsChanged.connect(lambda selected: [drs.setSetting("MB_groups", selected), updateSelectedMBGroups(selected)])
        MBExtButton.setMenu(CombinedMenu)
        MBExtButton.setPopupMode(QtGui.QToolButton.InstantPopup)
        MBExtButton.setIcon(CUI().icon('trophy'))
        MBExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        formLayout.addRow("Select only natural Samples", MBExtButton)

        selectedGroupsMBLabel = QtGui.QLabel("Selected Groups: None")
        formLayout.addRow(selectedGroupsMBLabel)

        MB_fit_ComboBox = QtGui.QComboBox(widget)
        MB_fit_ComboBox.addItems(['Constant','Linear', 'Exponential', 'Logarithmic','Polynomial'])
        MB_fit_ComboBox.setCurrentText(settings["MB_fit"])
        MB_fit_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MB_fit", t))
        formLayout.addRow("Fit of Mass Bias", MB_fit_ComboBox)

        formLayout.addRow('Plot of Mass Bias', PLOT)
        
        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

##%     Displaying results
        Results_Label = QtGui.QLabel("<h2>Results<\h2>")
        formLayout.addRow(Results_Label)

        FillLabelsCheckBox = QtGui.QCheckBox(widget)
        FillLabelsCheckBox.setChecked(settings["fill_Labels"])
        FillLabelsCheckBox.toggled.connect(lambda t: drs.setSetting("fill_Labels", bool(t)))
        formLayout.addRow("Fill in empty labels", FillLabelsCheckBox)

        def updateSelectedResultsGroups(selected):
                if selected:
                        selectedGroupsResultsLabel.setText(f"Selected Groups: {', '.join(selected)}")
                else:
                        selectedGroupsResultsLabel.setText("Selected Groups: None")

        ResultsExtButton = QtGui.QToolButton(widget)
        CombinedMenu = CheckableMenu(ResultsExtButton, combined_names)
        CombinedMenu.itemsChanged.connect(lambda selected: [drs.setSetting("Result_groups", selected), updateSelectedResultsGroups(selected)])
        ResultsExtButton.setMenu(CombinedMenu)
        ResultsExtButton.setPopupMode(QtGui.QToolButton.InstantPopup)
        ResultsExtButton.setIcon(CUI().icon('trophy'))
        ResultsExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        formLayout.addRow("Select groups to show result", ResultsExtButton)

        selectedGroupsResultsLabel = QtGui.QLabel("Selected Groups: None")
        formLayout.addRow(selectedGroupsResultsLabel)

        How230_232_ComboBox = QtGui.QComboBox(widget)
        How230_232_ComboBox.addItems(['Ratio of Intensities','Mean of Ratios'])
        How230_232_ComboBox.setCurrentText(settings["How230_232"])
        How230_232_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("How230_232", t))
        formLayout.addRow("How to calculate (230Th)/(232Th)", How230_232_ComboBox)

        How238_232_ComboBox = QtGui.QComboBox(widget)
        How238_232_ComboBox.addItems(['Ratio of Intensities','Mean of Ratios'])
        How238_232_ComboBox.setCurrentText(settings["How238_232"])
        How238_232_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("How238_232", t))
        formLayout.addRow("How to calculate (238U)/(232Th)", How238_232_ComboBox)        

        U235orU238_ComboBox = QtGui.QComboBox(widget)
        U235orU238_ComboBox.addItems(['based on 238U','based on 235U*137.818'])
        U235orU238_ComboBox.setCurrentText(settings["U235orU238"])
        U235orU238_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("U235orU238", t))
        formLayout.addRow("How to calculate (230Th)/(238U)", U235orU238_ComboBox)

        ResultLcCheckBox = QtGui.QCheckBox(widget)
        ResultLcCheckBox.setChecked(settings["plotLc1"])
        ResultLcCheckBox.toggled.connect(lambda t: drs.setSetting("plotLc1", bool(t)))
        formLayout.addRow("Plot only samples above detection limit", ResultLcCheckBox)

        formLayout.addRow('Evolution Plot', PLOT_Evolution)
        
        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

##%     Calculating ppm
        ppm_Label = QtGui.QLabel("<h2>Th ppm & U ppm based on Referenz Material<\h2>")
        formLayout.addRow(ppm_Label)
        
        ppmCheckBox = QtGui.QCheckBox(widget)
        ppmCheckBox.setChecked(settings["calculate_ppm"])
        ppmCheckBox.toggled.connect(lambda t: drs.setSetting("calculate_ppm", bool(t)))
        formLayout.addRow("Calculate Th & U ppm", ppmCheckBox)

        ppm_fit_ComboBox = QtGui.QComboBox(widget)
        ppm_fit_ComboBox.addItems(['Constant','Linear','Polynomial','Polynomial3','Exponential'])
        ppm_fit_ComboBox.setCurrentText(settings["MB_fit"])
        ppm_fit_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("ppm_fit", t))
        formLayout.addRow("Fit of CPS", ppm_fit_ComboBox)

        ppmComboBox = QtGui.QComboBox(widget)
        ppmComboBox.addItems(RMNames)
        ppmComboBox.setCurrentText(settings["ppm_calculated_on"])
        ppmComboBox.textActivated.connect(lambda t: drs.setSetting("ppm_calculated_on", t))
        formLayout.addRow("Which RM material to use", ppmComboBox)

        def updateRMCombo():
                ppmComboBox.clear()
                ppmComboBox.addItems(RMNames)

        data.selectionGroupsChanged.connect(updateRMCombo)

        formLayout.addRow('Simplified ppm calculation', PLOT_ppm)

        def toggleppmWidgets(checked):
                ppm_fit_ComboBox.setVisible(checked)
                ppmComboBox.setVisible(checked)
                PLOT_ppm.setVisible(checked)
        ppmCheckBox.toggled.connect(toggleppmWidgets)
        toggleppmWidgets(ppmCheckBox.isChecked())

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)


##%     Calculating SecEquilibrium 234/238
        SecEquil234_Label = QtGui.QLabel("<h2>Check for U234-U238 secular equilibrium<\h2>")
        formLayout.addRow(SecEquil234_Label)
        SecEquil234_explain_Label = QtGui.QLabel("If U234 has been measured, the secular equiblibrium between U234 and U238 can be tested.")
        formLayout.addRow(SecEquil234_explain_Label)
        
        SecEquil234CheckBox = QtGui.QCheckBox(widget)
        SecEquil234CheckBox.setChecked(settings["calculate_U234_U238"])
        SecEquil234CheckBox.toggled.connect(lambda t: drs.setSetting("calculate_U234_U238", bool(t)))
        formLayout.addRow("Check (U234)/(U238)", SecEquil234CheckBox)
        
        def updateSelectedSecGroups(selected):
                if selected:
                        selectedGroupsSecLabel.setText(f"Selected Groups: {', '.join(selected)}")
                else:
                        selectedGroupsSecLabel.setText("Selected Groups: None")

        SecExtButton = QtGui.QToolButton(widget)
        CombinedMenu = CheckableMenu(SecExtButton, combined_names)
        CombinedMenu.itemsChanged.connect(lambda selected: [drs.setSetting("U234-U238_groups", selected), updateSelectedSecGroups(selected)])
        SecExtButton.setMenu(CombinedMenu)
        SecExtButton.setPopupMode(QtGui.QToolButton.InstantPopup)
        SecExtButton.setIcon(CUI().icon('trophy'))
        SecExtButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        formLayout.addRow("Select groups to show result", SecExtButton)

        selectedGroupsSecLabel = QtGui.QLabel("Selected Groups: None")
        formLayout.addRow(selectedGroupsSecLabel)

        How234_238_ComboBox = QtGui.QComboBox(widget)
        How234_238_ComboBox.addItems(['Ratio of Intensities','Mean of Ratios'])
        How234_238_ComboBox.setCurrentText(settings["How234_238"])
        How234_238_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("How234_238", t))
        formLayout.addRow("How to calculate (234U)/(238U)", How234_238_ComboBox)

        Mask234CheckBox = QtGui.QCheckBox(widget)
        Mask234CheckBox.setChecked(settings["mask_U234_U238"])
        Mask234CheckBox.toggled.connect(lambda t: drs.setSetting("mask_U234_U238", bool(t)))
        formLayout.addRow("Plot only samples in secular equilibrium on Evolution Plot", Mask234CheckBox)

        formLayout.addRow('Abundance sensitivity of Monazite on mass 234', PLOT_Mon234)
        formLayout.addRow('Check for secular equilibrium between U234 and U238', PLOT_234_238)

        def toggleSecEquilWidgets(checked):
                SecExtButton.setVisible(checked)
                selectedGroupsSecLabel.setVisible(checked)
                How234_238_ComboBox.setVisible(checked)
                Mask234CheckBox.setVisible(checked)
                PLOT_Mon234.setVisible(checked)
                PLOT_234_238.setVisible(checked)
        SecEquil234CheckBox.toggled.connect(toggleSecEquilWidgets)
        toggleSecEquilWidgets(SecEquil234CheckBox.isChecked())

        drs.setSettingsWidget(widget)
