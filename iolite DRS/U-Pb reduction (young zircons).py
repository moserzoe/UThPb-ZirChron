# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Type: DRS
#/ Name: U-Pb reduction (young zircons)
#/ Authors: Zoe Moser
#/ Description: 
#/ References:
#/ Version:
#/ Contact:

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
from iolite.QtCore import QDateTime

import numpy as np
from scipy.optimize import curve_fit

try:
    PLOT_RS
except:
    PLOT_RS = Plot()
    PLOT_RS.setAttribute(Qt.WA_DeleteOnClose)
    PLOT_RS.setFixedSize(600,300)
    
    # Add annotation to show fit parameters
    ann2 = PLOT_RS.annotate('', 0.99, 0.95, 'ptAxisRectRatio', Qt.AlignRight | Qt.AlignBottom)
    ann2.visible = False

    def showSettings():
        d = PlotSettings(PLOT_RS)
        d.exec_()

    settingsAction = QAction(PLOT_RS.contextMenu())
    settingsAction.setText('Settings')
    settingsAction.triggered.connect(showSettings)
    PLOT_RS.contextMenu().addAction(settingsAction)

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
        """
        This method will be called by iolite when the user clicks
        Crunch Data in the DRS window or as part of a processing
        template. It should transform 'input' data into 'output'
        data using the provided settings.

        DRS progress can be updated via the 'message' and 'progress'
        signals. These will be displayed in the iolite interface.

        When finished, the 'finished' signal should be emitted.

        As an example, we will do baseline subtraction of all
        input channels using a DRS helper function.
        """

        drs.message("Starting baseline subtract DRS...")
        drs.progress(0)

        # Get settings
        settings = drs.settings()
        print(settings)

        indexChannel = data.timeSeries(settings["IndexChannel"])
        maskChannel = data.timeSeries(settings["MaskChannel"])
        cutoff = settings["MaskCutoff"]
        trim = settings["MaskTrim"]
        zircon_RM = data.selectionGroup(settings["RM_zircon"])
        RS_fit = settings["RS_fit"]

        # Create debug messages for the settings being used
        IoLog.debug("indexChannelName = %s" % indexChannel.name)
        IoLog.debug("maskChannelName = %s" % maskChannel.name)
        IoLog.debug("maskCutoff = %f" % cutoff)
        IoLog.debug("maskTrim = %f" % trim)

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

##% ROI and MOR calculations
        MOR_206_238 = data.timeSeries('Pb206_CPS').data()/data.timeSeries('U238_CPS').data()
        MOR_206_238_timeSeries = data.createTimeSeries('Pb206/U238_MOR', data.Intermediate, indexChannel.time(), MOR_206_238)

        MOR_238_206 = data.timeSeries('U238_CPS').data()/data.timeSeries('Pb206_CPS').data()
        MOR_238_206_timeSeries = data.createTimeSeries('U238/Pb206_MOR', data.Intermediate, indexChannel.time(), MOR_238_206)

        MOR_207_206 = data.timeSeries('Pb207_CPS').data()/data.timeSeries('Pb206_CPS').data()
        MOR_207_206_timeSeries = data.createTimeSeries('Pb207/Pb206_MOR', data.Intermediate, indexChannel.time(), MOR_207_206)

        def ROI_206_238(sel):
                U238 = data.result(sel,data.timeSeries('U238_CPS')).value() 
                Pb206 = data.result(sel,data.timeSeries('Pb206_CPS')).value()     
                result  = Pb206/U238           
                return result
        
        '''
        def ROI_206_238_2SE(sel):
                U238 = data.result(sel,data.timeSeries('U238_CPS')).value() 
                Pb206 = data.result(sel,data.timeSeries('Pb206_CPS')).value()     
                U238_2SE = data.result(sel,data.timeSeries('U238_CPS')).uncertaintyAs2SE() 
                Pb206_2SE = data.result(sel,data.timeSeries('Pb206_CPS')).uncertaintyAs2SE() 
                Pb206_U238 = data.associatedResult(sel,'Pb206/U238_ROI').value() 
                result  =  Pb206_U238*np.sqrt((U238_2SE/U238)**2+(Pb206_2SE/Pb206)**2)        
                return result
        '''
        def ROI_206_238_2SE(sel):
                result = data.result(sel,data.timeSeries('Pb206/U238_MOR')).uncertaintyAs2SE()
                return result
                
        def ROI_238_206(sel):
                U238 = data.result(sel,data.timeSeries('U238_CPS')).value() 
                Pb206 = data.result(sel,data.timeSeries('Pb206_CPS')).value()     
                result  = U238/Pb206           
                return result
        
        '''
        def ROI_238_206_2SE(sel):
                U238 = data.result(sel,data.timeSeries('U238_CPS')).value() 
                Pb206 = data.result(sel,data.timeSeries('Pb206_CPS')).value()     
                U238_2SE = data.result(sel,data.timeSeries('U238_CPS')).uncertaintyAs2SE() 
                Pb206_2SE = data.result(sel,data.timeSeries('Pb206_CPS')).uncertaintyAs2SE() 
                U238_Pb206 = data.associatedResult(sel,'U238/Pb206_ROI').value() 
                result  =  U238_Pb206*np.sqrt((U238_2SE/U238)**2+(Pb206_2SE/Pb206)**2)        
                return result
        '''
        def ROI_238_206_2SE(sel):
                result = data.result(sel,data.timeSeries('U238/Pb206_MOR')).uncertaintyAs2SE()
                return result
        
        def ROI_207_206(sel):
                Pb207 = data.result(sel,data.timeSeries('Pb207_CPS')).value() 
                Pb206 = data.result(sel,data.timeSeries('Pb206_CPS')).value()     
                result  = Pb207/Pb206           
                return result
        
        '''
        def ROI_207_206_2SE(sel):
                Pb207 = data.result(sel,data.timeSeries('Pb207_CPS')).value() 
                Pb206 = data.result(sel,data.timeSeries('Pb206_CPS')).value()     
                Pb207_2SE = data.result(sel,data.timeSeries('Pb207_CPS')).uncertaintyAs2SE() 
                Pb206_2SE = data.result(sel,data.timeSeries('Pb206_CPS')).uncertaintyAs2SE() 
                Pb207_Pb206 = data.associatedResult(sel,'Pb207/Pb206_ROI').value() 
                result  =  Pb207_Pb206*np.sqrt((Pb207_2SE/Pb207)**2+(Pb206_2SE/Pb206)**2)        
                return result        
        '''
        def ROI_207_206_2SE(sel):
                result = data.result(sel,data.timeSeries('Pb207/Pb206_MOR')).uncertaintyAs2SE()
                return result
        
        data.registerAssociatedResult('Pb206/U238_ROI',ROI_206_238)
        data.registerAssociatedResult('Pb206/U238_ROI_2SE',ROI_206_238_2SE) 
        data.registerAssociatedResult('U238/Pb206_ROI',ROI_238_206)
        data.registerAssociatedResult('U238/Pb206_ROI_2SE',ROI_238_206_2SE)    
        data.registerAssociatedResult('Pb207/Pb206_ROI',ROI_207_206)
        data.registerAssociatedResult('Pb207/Pb206_ROI_2SE',ROI_207_206_2SE)   

        def MOR_ROI_207_206(sel):
                Pb207_Pb206_MOR = data.result(sel,data.timeSeries('Pb207/Pb206_MOR')).value()
                Pb207_Pb206_ROI = data.associatedResult(sel,'Pb207/Pb206_ROI').value() 
                result  =   Pb207_Pb206_MOR/Pb207_Pb206_ROI      
                return result          
        
        data.registerAssociatedResult('Pb207/Pb206_MOR/ROI',MOR_ROI_207_206)

##% calculating the corrections 
        RMNames = data.selectionGroupNames(data.ReferenceMaterial)
        #group91500 = data.selectionGroup(next((ch for ch in RMNames if "91500" in ch.lower()), None))
        #BaselineGroup = data.selectionGroupNames(data.Baseline)[0]

        meas = []
        meas_2SE = []
        start_times = []
        for sel in zircon_RM.selections():
                sel_start_time = data.timeSeries('TotalBeam').timeForSelection(sel)
                RSF_sel_238_206 = data.associatedResult(sel,'U238/Pb206_ROI').value()
                RSF_sel_238_206_2SE = data.associatedResult(sel,'U238/Pb206_ROI_2SE').value()                
                meas.append(RSF_sel_238_206)
                meas_2SE.append(RSF_sel_238_206_2SE)
                start_times.append(sel_start_time[0])
        
        meas_array = np.array(meas)
        meas_2SE_array = np.array(meas_2SE)
        start_times_array = np.array(start_times)-indexChannel.time()[0]

        PLOT_RS.clearGraphs()
        g = PLOT_RS.addGraph()
        g.setName(f'{settings["RM_zircon"]} uncorr')
        g.setLineStyle('lsNone')
        g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[1 % len(PLOT_COLORS)])
        g.setData(start_times_array, meas_array)   
        
        for x, y, y_err in zip(start_times_array, meas_array, meas_2SE_array):
                v_error = PLOT_RS.addGraph()
                vERR_pen = QPen(QColor('dark grey'))
                vERR_pen.setStyle(Qt.SolidLine)
                v_error.setPen(vERR_pen)
                v_error.setData([x, x], [y - y_err, y + y_err])
                v_error.removeFromLegend() 

        PLOT_RS.left().label = '238/206'
        PLOT_RS.bottom().label = 'Time since start of session (s)'
        PLOT_RS.setLegendVisible(True)
        PLOT_RS.setLegendAlignment(Qt.AlignTop | Qt.AlignLeft)   
        PLOT_RS.setToolsVisible(True)
        PLOT_RS.rescaleAxes()
        PLOT_RS.replot()


        fit_x = np.linspace(
            start_times_array.min(),
            start_times_array.max(),
            50
        )

        if RS_fit == 'Constant':
                def fitConst(x, a, b):
                        return a + b*0

                params, cov = curve_fit(fitConst, start_times_array, meas_array,  ftol=1e-5)

                fit_y = params[0]*np.ones(len(fit_x))
                RS_fit_fit = (params[0]+0*(indexChannel.time()-min(indexChannel.time())))
                RS_fit_fit_timeSeries = data.createTimeSeries('RS_correction', data.Intermediate, indexChannel.time(), RS_fit_fit)
        if RS_fit == 'Linear':
                def fitLin(x, a, b):
                        return a + b*x

                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                fit_y = params[1]*fit_x + params[0]
                RS_fit_fit = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))
                RS_fit_fit_timeSeries = data.createTimeSeries('RS_correction', data.Intermediate, indexChannel.time(), RS_fit_fit)
        if RS_fit == 'Exponential':
                def fitExp(x, a, b, c):
                        return a + b * np.exp(-c * x)
                initial_guess = [np.mean(meas_array), 0, 0.01]
                try:
                        params, cov = curve_fit(fitExp, start_times_array, meas_array, p0=initial_guess)
                        fit_y = params[0] + params[1] * np.exp(-params[2] * fit_x)
                        RS_fit_fit = (params[0]+params[1]*np.exp(-params[2]*(indexChannel.time()-min(indexChannel.time()))))
                        RS_fit_fit_timeSeries = data.createTimeSeries('RS_correction', data.Intermediate, indexChannel.time(), RS_fit_fit)

                except RuntimeError as e:
                        if 'Optimal parameters not found: Number of calls' in str(e):

                                IoLog.warning('Could not find fit to data with Exp model. Switching to Linear model.')

                                settings["MB_fit"] = 'Linear'

                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)
                                
                                fit_y = params[1]*fit_x + params[0]
                                RS_fit_fit = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))
                                RS_fit_fit_timeSeries = data.createTimeSeries('RS_correction', data.Intermediate, indexChannel.time(), RS_fit_fit)

                        else:
                                raise

        if RS_fit == 'Logarithmic':
                def fitLog(x, a, b, c):
                        return a + b * np.log(x + c)  # Ensure x + c is > 0
                
                initial_guess = [np.mean(meas_array), 1, 2000]
                try:    
                        params, cov = curve_fit(fitLog, start_times_array, meas_array, p0=initial_guess)
                        fit_y = params[0] + params[1] * np.log(fit_x + params[2])
                        RS_fit_fit = (params[0]+params[1]*np.log(params[2]+(indexChannel.time()-min(indexChannel.time()))))
                        RS_fit_fit_timeSeries = data.createTimeSeries('RS_correction', data.Intermediate, indexChannel.time(), RS_fit_fit)

                except RuntimeError as e:
                        if 'Optimal parameters not found: Number of calls' in str(e):

                                IoLog.warning('Could not find fit to data with Exp model. Switching to Linear model.')

                                settings["MB_fit"] = 'Linear'

                                def fitLin(x, a, b):
                                        return a + b*x

                                params, cov = curve_fit(fitLin, start_times_array, meas_array, ftol=1e-5)

                                fit_y = params[1]*fit_x + params[0]
                                RS_fit_fit = (params[0]+params[1]*(indexChannel.time()-min(indexChannel.time())))
                                RS_fit_fit_timeSeries = data.createTimeSeries('RS_correction', data.Intermediate, indexChannel.time(), RS_fit_fit)                        
                        else:
                                raise
        
        g = PLOT_RS.addGraph()
        g.setName("Fit")
        fit_pen = QPen(QColor('black'))
        fit_pen.setStyle(Qt.DashLine)
        fit_pen.setWidth(2)
        g.pen = fit_pen
        g.setData(fit_x, fit_y)

        ann2.visible = True
        if RS_fit == 'Constant':
                ann2.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f}</p>'''
        elif RS_fit == 'Linear':
                ann2.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f} + {params[1]:.3e}*t</p>'''
        elif RS_fit == 'Exponential': 
                ann2.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f} + {params[1]:.2f} * exp(-{params[2]:.3e}*t)</p>'''
        elif RS_fit == 'Logarithmic': 
                ann2.text = f'''
                <p style="color:black;">
                <b>Fit Parameters:</b><br>
                {params[0]:.2f} + {params[1]:.2f} * ln(t + {params[2]:.2f})</p>'''

        PLOT_RS.setToolsVisible(True)
        PLOT_RS.rescaleAxes()
        PLOT_RS.replot()

        soll_0638_91500 = data.referenceMaterialData(settings["RM_zircon"])["206Pb/238U"].value() #0.17917
        soll_3806_91500 = 1/soll_0638_91500
        soll_0706_91500 = data.referenceMaterialData(settings["RM_zircon"])["207Pb/206Pb"].value() #0.07488
        ist_0638_91500 = []
        ist_3806_91500 = []
        ist_0706_91500 = []  
        corr_ist_0638_91500 = []
        corr_ist_3806_91500 = []
        start_time = []
        for sel in zircon_RM.selections():
                ist_0638_91500_sel = data.associatedResult(sel,'Pb206/U238_ROI').value()
                ist_3806_91500_sel = data.associatedResult(sel,'U238/Pb206_ROI').value()
                ist_0706_91500_sel = data.associatedResult(sel,'Pb207/Pb206_ROI').value()
                start_time_sel = indexChannel.timeForSelection(sel)
                ist_0638_91500.append(ist_0638_91500_sel)
                ist_3806_91500.append(ist_3806_91500_sel)
                ist_0706_91500.append(ist_0706_91500_sel)
                start_time.append(start_time_sel[0])
                corr_ist_0638_91500_sel = ist_0638_91500_sel/(1/data.result(sel,data.timeSeries('RS_correction')).value())
                corr_ist_0638_91500.append(corr_ist_0638_91500_sel)
                corr_ist_3806_91500_sel = ist_3806_91500_sel/(data.result(sel,data.timeSeries('RS_correction')).value())
                corr_ist_3806_91500.append(corr_ist_3806_91500_sel)

        corr_ist_0638_91500_norm_mean = np.array(corr_ist_0638_91500)*np.mean(ist_0638_91500)
        corr_ist_3806_91500_norm_mean = np.array(corr_ist_3806_91500)*np.mean(ist_3806_91500)

        corr_MB = 1-(soll_0706_91500/np.mean(ist_0706_91500))
        corr_MB_2SE  = 2*np.std(ist_0706_91500)/np.sqrt(len(ist_0706_91500))
        corr_DF_0638 = 1-(soll_0638_91500/np.mean(corr_ist_0638_91500_norm_mean))
        corr_DF_0638_2SE = 2*np.std(corr_ist_0638_91500_norm_mean)/np.sqrt(len(corr_ist_0638_91500_norm_mean))
        corr_DF_3806 = 1-(soll_3806_91500/np.mean(corr_ist_3806_91500_norm_mean))
        corr_DF_3806_2SE = 2*np.std(corr_ist_3806_91500_norm_mean)/np.sqrt(len(corr_ist_3806_91500_norm_mean))
        print(f'corr_MB: {corr_MB} pm {corr_MB_2SE}, corr_DF_0638: {corr_DF_0638} pm {corr_DF_0638_2SE},corr_DF_3806: {corr_DF_3806} pm {corr_DF_3806_2SE}')


        def RS_corr_3806(sel):
                result= Result()
                result = data.associatedResult(sel,'U238/Pb206_ROI').value()/(data.result(sel,data.timeSeries('RS_correction')).value())*np.mean(ist_3806_91500)
                return result
        
        def RS_corr_0638(sel):
                result= Result()
                result = data.associatedResult(sel,'Pb206/U238_ROI').value()/(1/data.result(sel,data.timeSeries('RS_correction')).value())*np.mean(ist_0638_91500)
                return result
        
        data.registerAssociatedResult('Pb206/U238_ROI_RS_corr', RS_corr_0638)
        data.registerAssociatedResult('U238/Pb206_ROI_RS_corr', RS_corr_3806)
        
##% matching backgrounds
        BaselineGroup = data.selectionGroupNames(data.Baseline)
        start_time_BG = []

        for gr in BaselineGroup:
                bg = data.selectionGroup(gr)
                for sel_bg in bg.selections():
                        start_time_BG_sel = data.timeSeries('TotalBeam').timeForSelection(sel_bg)
                        if start_time_BG_sel is not None and len(start_time_BG_sel) > 0:
                                start_time_BG.append(start_time_BG_sel[0])

                backgrounds = sorted(start_time_BG)

        def matched_background_start_time(sel):
                sel_start_time = data.timeSeries('TotalBeam').timeForSelection(sel)
                if sel_start_time is None or len(sel_start_time) == 0:
                        return None  # If selection has no time, return None
                
                sel_start_time = sel_start_time[0]  # Get the first time value
                
                matched_bg = None
                for bg_start in backgrounds:
                        if bg_start < sel_start_time:
                                matched_bg = bg_start        
                        else:
                                break  # Stop when we exceed the selection start time
                
                return matched_bg if matched_bg is not None else np.nan  # Return NaN if no match

        data.registerAssociatedResult('BG_StartTime', matched_background_start_time)

        def SG_StartTime(sel):
                signal_starttime = indexChannel.timeForSelection(sel)[0]
                result  =  signal_starttime   
                return result

        data.registerAssociatedResult('SG_StartTime',SG_StartTime)

        def SG_Duration(sel):
                signal_duration= sel.property("Duration (s)")
                result  =  signal_duration   
                return result

        data.registerAssociatedResult('SG_Duration',SG_Duration)

        def SG_after_BG(sel):
                signal_start = data.associatedResult(sel,'SG_StartTime').value() 
                BG_start = data.associatedResult(sel,'BG_StartTime').value() 
                result  =  signal_start-BG_start
                return result

        data.registerAssociatedResult('SG_after_BG',SG_after_BG)
        
##% Applying corrections
        BaselineGroup = data.selectionGroupNames(data.Baseline)[0]
        def ROI_206_238_DF_corr_complex(sel):
                if sel in data.selectionGroup(BaselineGroup).selections():
                        return NaN  
                
                Pb206_U238 = data.associatedResult(sel,'Pb206/U238_ROI_RS_corr').value() 

                if f'{sel.name}_91500' in data.selectionGroupNames(data.Sample):
                        data.removeSelectionGroup(f'{sel.name}_91500')

                data.createSelectionGroup(f'{sel.name}_91500', data.Sample)
                sel_group = data.selectionGroup(f'{sel.name}_91500')   
                SG_after_BG = data.associatedResult(sel, 'SG_after_BG').value()
                SG_Duration = data.associatedResult(sel, 'SG_Duration').value()

                soll_0638_91500 = data.referenceMaterialData(settings["RM_zircon"])["206Pb/238U"].value() #0.17917
                ist_0638_91500 = []
                ist_3806_91500 = []
                ist_0706_91500 = []  
                corr_ist_0638_91500 = []
                for sel_91500 in zircon_RM.selections():
                        x_BG_StartTimes_sel = data.associatedResult(sel_91500,'BG_StartTime').value()
                        x_SG_StartTimes_sel = x_BG_StartTimes_sel+SG_after_BG
                        x_SG_EndTimes_sel = x_SG_StartTimes_sel+SG_Duration
                        
                        start_qdatetime = QDateTime.fromMSecsSinceEpoch(x_SG_StartTimes_sel*1000)
                        end_qdatetime = QDateTime.fromMSecsSinceEpoch(x_SG_EndTimes_sel*1000)
                        data.createSelection(sel_group,start_qdatetime, end_qdatetime,'')

                for sel_new in  sel_group.selections():       
                        ist_0638_91500_sel = data.associatedResult(sel_new,'Pb206/U238_ROI_RS_corr').value()
                        ist_3806_91500_sel = data.associatedResult(sel_new,'U238/Pb206_ROI_RS_corr').value()
                        ist_0706_91500_sel = data.associatedResult(sel_new,'Pb207/Pb206_ROI').value()
                        ist_0638_91500.append(ist_0638_91500_sel)
                        ist_3806_91500.append(ist_3806_91500_sel)
                        ist_0706_91500.append(ist_0706_91500_sel)

                corr_DF_0638 = 1-(soll_0638_91500/np.mean(ist_0638_91500))             
                
                data.removeSelectionGroup(f'{sel.name}_91500')

                result  =  Pb206_U238 + Pb206_U238*(-corr_DF_0638)  
                return result        

        def ROI_238_206_DF_corr_complex(sel):
                ROI_206_238_DF = data.associatedResult(sel,'Pb206/U238_ROI_DFcorr').value() 
                result = 1/ ROI_206_238_DF
                return result             
        
        def ROI_206_238_2SE_DF_corr(sel):
                Pb206_U238 = data.associatedResult(sel,'Pb206/U238_ROI').value() 
                Pb206_U238_2SE = data.associatedResult(sel,'Pb206/U238_ROI_2SE').value() 
                #result = Pb206_U238_2SE
                result = np.sqrt(((1-corr_DF_0638)*Pb206_U238_2SE)**2+(Pb206_U238*corr_DF_0638_2SE)**2)
                return result   
        
        def ROI_238_206_2SE_DF_corr(sel):
                U238_Pb206 = data.associatedResult(sel,'U238/Pb206_ROI').value() 
                U238_Pb206_2SE = data.associatedResult(sel,'U238/Pb206_ROI_2SE').value() 
                #result = U238_Pb206_2SE
                result = np.sqrt(((1-corr_DF_3806)*U238_Pb206_2SE)**2+(-U238_Pb206*corr_DF_3806_2SE)**2)
                return result

        def ROI_207_206_MB_corr(sel):
                Pb207_Pb206 = data.associatedResult(sel,'Pb207/Pb206_ROI').value() 
   
                result  =  Pb207_Pb206 + Pb207_Pb206*(-corr_MB)          
                return result
        
        def ROI_207_206_2SE_MB_corr(sel):
                Pb207_Pb206 = data.associatedResult(sel,'Pb207/Pb206_ROI').value()
                Pb207_Pb206_2SE = data.associatedResult(sel,'Pb207/Pb206_ROI_2SE').value()
                result = Pb207_Pb206_2SE
                result  =  np.sqrt(((1-corr_MB)*Pb207_Pb206_2SE)**2+(-Pb207_Pb206*corr_MB_2SE)**2)
                return result      

        data.registerAssociatedResult('Pb206/U238_ROI_DFcorr',ROI_206_238_DF_corr_complex)
        data.registerAssociatedResult('Pb206/U238_ROI_DFcorr_2SE',ROI_206_238_2SE_DF_corr)   
        data.registerAssociatedResult('U238/Pb206_ROI_DFcorr',ROI_238_206_DF_corr_complex)  
        data.registerAssociatedResult('U238/Pb206_ROI_DFcorr_2SE',ROI_238_206_2SE_DF_corr)          
        data.registerAssociatedResult('Pb207/Pb206_ROI_MBcorr',ROI_207_206_MB_corr)
        data.registerAssociatedResult('Pb207/Pb206_ROI_MBcorr_2SE',ROI_207_206_2SE_MB_corr)  
        
        ###
        meas = []
        meas_2SE = []
        start_times = []
        for sel in zircon_RM.selections():
                sel_start_time = data.timeSeries('TotalBeam').timeForSelection(sel)
                RSF_sel_238_206 = data.associatedResult(sel,'U238/Pb206_ROI_RS_corr').value()
                RSF_sel_238_206_2SE = data.associatedResult(sel,'U238/Pb206_ROI_DFcorr_2SE').value()                
                meas.append(RSF_sel_238_206)
                meas_2SE.append(RSF_sel_238_206_2SE)
                start_times.append(sel_start_time[0])
        
        meas_array = np.array(meas)
        meas_2SE_array = np.array(meas_2SE)
        start_times_array = np.array(start_times)-indexChannel.time()[0]

        g = PLOT_RS.addGraph()
        g.setName(f'{settings["RM_zircon"]} RS corr')
        g.setLineStyle('lsNone')
        g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[2 % len(PLOT_COLORS)])
        g.setData(start_times_array, meas_array)   
        
        for x, y, y_err in zip(start_times_array, meas_array, meas_2SE_array):
                v_error = PLOT_RS.addGraph()
                vERR_pen = QPen(QColor('dark grey'))
                vERR_pen.setStyle(Qt.SolidLine)
                v_error.setPen(vERR_pen)
                v_error.setData([x, x], [y - y_err, y + y_err])
                v_error.removeFromLegend()   
        PLOT_RS.setToolsVisible(True)
        PLOT_RS.rescaleAxes()
        PLOT_RS.replot()  
        ###
        '''
        meas = []
        meas_2SE = []
        start_times = []
        for sel in zircon_RM.selections():
                sel_start_time = data.timeSeries('TotalBeam').timeForSelection(sel)
                RSF_sel_238_206 = data.associatedResult(sel,'U238/Pb206_ROI_DFcorr').value()
                RSF_sel_238_206_2SE = data.associatedResult(sel,'U238/Pb206_ROI_DFcorr_2SE').value()                
                meas.append(RSF_sel_238_206)
                meas_2SE.append(RSF_sel_238_206_2SE)
                start_times.append(sel_start_time[0])
        
        meas_array = np.array(meas)
        meas_2SE_array = np.array(meas_2SE)
        start_times_array = np.array(start_times)-indexChannel.time()[0]

        g = PLOT_RS.addGraph()
        g.setName(f'{settings["RM_zircon"]} DF corr')
        g.setLineStyle('lsNone')
        g.setScatterStyle('ssDisc', 6.0, PLOT_COLORS[3 % len(PLOT_COLORS)])
        g.setData(start_times_array, meas_array)   
        
        for x, y, y_err in zip(start_times_array, meas_array, meas_2SE_array):
                v_error = PLOT_RS.addGraph()
                vERR_pen = QPen(QColor('dark grey'))
                vERR_pen.setStyle(Qt.SolidLine)
                v_error.setPen(vERR_pen)
                v_error.setData([x, x], [y - y_err, y + y_err])
                v_error.removeFromLegend()   
        
        PLOT_RS.setToolsVisible(True)
        PLOT_RS.rescaleAxes()
        PLOT_RS.replot()     
        '''
##% End
        drs.message("Finished!")
        drs.progress(100)
        drs.finished()


def settingsWidget():
        """
        This function puts together a user interface to configure the DRS.

        It is important to have the last line of this function call:
        drs.setSettingsWidget(widget)
        """

        widget = QtGui.QWidget()
        formLayout = QtGui.QFormLayout()
        widget.setLayout(formLayout)
        
        RMNames = data.selectionGroupNames(data.ReferenceMaterial)
        defaultRSFName = next((ch for ch in RMNames if "91500" in ch.lower()), None)

        drs.setSetting("IndexChannel", 'TotalBeam')
        drs.setSetting("MaskChannel", 'TotalBeam')
        drs.setSetting("MaskCutoff", 100000.0)
        drs.setSetting("MaskTrim", 0.0)
        drs.setSetting("RM_zircon", defaultRSFName)
        drs.setSetting("RS_fit", 'Constant')

        settings = drs.settings()
        
        Title_Label = QtGui.QLabel("<h1>U-Pb dating of young zircons DRS<\h1>")
        formLayout.addRow(Title_Label)
        Subtitle_Label = QtGui.QLabel("Appropriate statistics: Session - Preferences - Stats - Baseline/Normal stats calculation: Mean no outlier reject")
        formLayout.addRow(Subtitle_Label)
        Subsubtitle_Label = QtGui.QLabel("Suggested Baseline Spline: StepForward")
        formLayout.addRow(Subsubtitle_Label)

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

        Baseline_Label = QtGui.QLabel("<h2>Baseline Subtract<\h2>")
        formLayout.addRow(Baseline_Label)

        indexComboBox = QtGui.QComboBox(widget)
        indexComboBox.addItems(data.timeSeriesNames(data.Input))
        indexComboBox.setCurrentText(settings["IndexChannel"])
        indexComboBox.currentTextChanged.connect(lambda t: drs.setSetting("IndexChannel", t))

        maskComboBox = QtGui.QComboBox(widget)
        maskComboBox.addItems(data.timeSeriesNames(data.Input))
        maskComboBox.setCurrentText(settings["MaskChannel"])
        maskComboBox.currentTextChanged.connect(lambda t: drs.setSetting("MaskChannel", t))

        maskLineEdit = QtGui.QLineEdit(widget)
        maskLineEdit.setText(settings["MaskCutoff"])
        maskLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskCutoff", float(t)))

        maskTrimLineEdit = QtGui.QLineEdit(widget)
        maskTrimLineEdit.setText(settings["MaskTrim"])
        maskTrimLineEdit.textChanged.connect(lambda t: drs.setSetting("MaskTrim", float(t)))

        formLayout.addRow("Index channel", indexComboBox)
        formLayout.addRow("Mask channel", maskComboBox)
        formLayout.addRow("Mask cutoff", maskLineEdit)
        formLayout.addRow("Mask trim", maskTrimLineEdit)

        verticalSpacer2 = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum)
        formLayout.addItem(verticalSpacer2)

        Correction_Label = QtGui.QLabel("<h2>Corrections<\h2>")
        formLayout.addRow(Correction_Label)
        explain1_Label = QtGui.QLabel("1. Ratios of intensities (ROI's) are calculated.")
        formLayout.addRow(explain1_Label)
        explain2_Label = QtGui.QLabel("2. Correction factors for Downhole Fractionation and Mass Bias based on known ratios from reference zircon are calculated.")
        formLayout.addRow(explain2_Label)

        ComboBox = QtGui.QComboBox(widget)
        ComboBox.addItems(RMNames)
        ComboBox.setCurrentText(settings["RM_zircon"])
        ComboBox.textActivated.connect(lambda t: drs.setSetting("RM_zircon", t))
        formLayout.addRow("Which reference zircon to use", ComboBox)

        def updateRMCombo():
                ComboBox.clear()
                ComboBox.addItems(RMNames)

        data.selectionGroupsChanged.connect(updateRMCombo)

        RS_fit_ComboBox = QtGui.QComboBox(widget)
        RS_fit_ComboBox.addItems(['Constant','Linear','Logarithmic','Exponential'])
        RS_fit_ComboBox.setCurrentText(settings["RS_fit"])
        RS_fit_ComboBox.currentTextChanged.connect(lambda t: drs.setSetting("RS_fit", t))
        formLayout.addRow("Temporal fit for\nDownhole fractiontion", RS_fit_ComboBox)

        formLayout.addRow('Relative correction\nbetween U and Pb', PLOT_RS)

        explain3_Label = QtGui.QLabel("3. Downhole Fractionation correction is applied to U238/Pb206 and Pb206/U238.")
        formLayout.addRow(explain3_Label)
        explain4_Label = QtGui.QLabel("4. Mass Bias correction is applied to Pb207/Pb206.")
        formLayout.addRow(explain4_Label)

        drs.setSettingsWidget(widget)
