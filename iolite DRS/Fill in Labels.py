# A python-based data reduction scheme for iolite 4 starts with some metadata
#/ Type: DRS
#/ Name: Fill in Labels
#/ Authors: Zoe Moser
#/ Description: Add missing labels for additional selections.
#/ References:
#/ Version: 1.0
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
from time import sleep
import numpy as np


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


##%     Fill empty labels      
        BaselineGroup = data.selectionGroupNames(data.Baseline)
        start_time_BG = []

        # Step 1: Gather background data
        for gr in BaselineGroup:
                bg = data.selectionGroup(gr)
                for sel_bg in bg.selections():
                        start_time_BG_sel = data.timeSeries('TotalBeam').timeForSelection(sel_bg)

                        start_time_BG.append(start_time_BG_sel[0])

        # Combine background data into a sorted list by start time
        backgrounds = sorted(zip(start_time_BG), key=lambda x: x[0])

        # Step 2: Iterate through combined selection groups
        combined_names = data.selectionGroupNames(data.ReferenceMaterial) + data.selectionGroupNames(data.Sample)
        unmatched_selections = []  # Track selections without names

        # Step 3: Assign names based on matching background and subsequent selection
        for combgr in combined_names:
                combgrsg = data.selectionGroup(combgr)
                for sel in combgrsg.selections():
                        if not sel.property('Rep Rate'):
                                sel.name = ''
                        if sel.name.strip().lower() == "linked selection":
                                sel.name = ''
                        if sel.name:  # Skip if already named
                                continue
                        
                        # Get the start time of the current selection
                        sel_start_time = data.timeSeries('TotalBeam').timeForSelection(sel)[0]
                        # Match the closest background based on start time
                        matched_bg = None
                        for bg_start in backgrounds:
                                if bg_start < sel_start_time:
                                        matched_bg = (bg_start)
                                else:
                                        break  # Stop when background start time exceeds selection start time

                        if matched_bg:
                                bg_start= matched_bg

                                # Find the selection with a start time slightly later than the matched background
                                matched_name = None
                                check_group = data.selectionGroup(combgr)
                                for sel_check in check_group.selections():
                                        check_start_time = data.timeSeries('TotalBeam').timeForSelection(sel_check)[0]
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

        drs.setSettingsWidget(widget)
