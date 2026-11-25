
#by Henry Schumacher
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import time
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import os
import sys
import json
import uuid
import h5py
import math
import xraydb
import plotly
import base64
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
import numpy as np
import pandas as pd
# import pyxray as xy
import odrpack as odr
import seaborn as sb
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column
from bokeh.models import Slider, HoverTool, BoxZoomTool, ResetTool, PanTool, SaveTool, FileInput, Div, ColumnDataSource, Paragraph
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.io import show

import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
from getmac import get_mac_address as gma
from itertools import chain
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea, VPacker
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
from colors import load_colors
from PIXE_functions import *
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#

from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times"],
    "text.usetex": True,
    "font.size": 8,
    "pgf.rcfonts": False
})


plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.preamble": "\n".join([
          r'\usepackage{amsmath}',
     ]),
})

#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
color_schemes = load_colors()
#-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-#
#http://localhost:5006/PIXE_interactive_plot
#bokeh serve --show .\PIXE_interactive_plot.py 


'''
sys.argv[1] sollte der Dateiname für das Spektrum sein.
'''


source = ColumnDataSource(data=dict(x=[], y=[]))
# -------------------------------------------------
# Plot
# -------------------------------------------------
plot = figure(
    title="VSPC Data",
    height=400,
    width=800,
    tools="pan,wheel_zoom,box_zoom,undo,redo,reset",
    sizing_mode="stretch_width"
)

plot.line(
    x="x",
    y="y",
    source=source,
    line_width=2,
    color=color_schemes['c_five2'][1]
)

output = Div(text='Uplode a .vspc file.',)
file_input = FileInput(accept='.vspc')


def load_file(attr, old, new):

    if not new:
        plot.title.text = f"Nothing loaded yet."
        return

    try:
        # Decode uploaded file
        file_bytes = base64.b64decode(new)

        # Custom unpacking
        json_data = json.loads(file_bytes.decode('utf-8'))
        # json_data = read_json_formatted_file(file_bytes, encoding="utf-8")

        # Update plot
        source.data = {
            "x": np.arange(0,len(json_data['RawData'][:-1]),1),
            "y": json_data['RawData'][:-1]
        }

        plot.title.text = f"Loaded: {json_data['MeasurementInfo']['Time']['Start']}"

    except Exception as e:
        plot.title.text = f"Error loading file: {e}"




file_input.on_change("value", load_file)

plot.add_tools(HoverTool(
            tooltips=[("x", "@x"), ("y", "@y")]
        ))

layout = column(
    output,
    file_input,
    plot,
    sizing_mode="stretch_width"
)

curdoc().add_root(layout)
curdoc().title = "Raw Data Plot"
