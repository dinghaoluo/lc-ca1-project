# -*- coding: utf-8 -*-
'''
Created on Mon 14 Apr 14:15:45 2025
Updated on Fri 18 Apr 14:50:12 2025:
    - improved reference-image rendering and ROI-overlay switching
    - added segmentation controls and loading of saved ROI dictionaries
    - adjusted ROI colours and added persistent fixed ROIs
Updated on Thu 27 Nov 2025:
    - added remove-fixed and clear-board controls
    - loaded ROI dictionaries are fixed automatically
    - space toggles the ROI overlay
Updated on Mon May 11 2026:
    - reorganised the interface around grouped controls, status feedback, and shortcuts
    - added light/dark themes and selection utilities

FibreSegger v5.0

same input/output paths, ROI dict format, and segmentation workflow as v4.

@author: Dinghao Luo
'''

#%% imports
import colorsys
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Qt5Agg')

import cv2
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.ndimage import median_filter
from skimage.measure import regionprops
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor, QFont, QIcon, QKeySequence, QPalette, QTextCursor
from PyQt5.QtWidgets import (
    QAbstractSpinBox,
    QApplication,
    QButtonGroup,
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QPlainTextEdit,
    QScrollArea,
    QShortcut,
    QSpinBox,
    QSizePolicy,
    QStatusBar,
    QStyle,
    QToolButton,
    QVBoxLayout,
    QWidget,
)


PARAMETER_SPECS = [
    {
        'name': 'solidity min',
        'kind': 'float',
        'default': 0.1,
        'minimum': 0.0,
        'maximum': 1.0,
        'step': 0.01,
        'decimals': 2,
    },
    {
        'name': 'eccentricity min',
        'kind': 'float',
        'default': 0.75,
        'minimum': 0.0,
        'maximum': 1.0,
        'step': 0.01,
        'decimals': 2,
    },
    {
        'name': 'thinness max',
        'kind': 'float',
        'default': 0.8,
        'minimum': 0.0,
        'maximum': 1.0,
        'step': 0.01,
        'decimals': 2,
    },
    {
        'name': 'tophat kernel',
        'kind': 'int',
        'default': 11,
        'minimum': 1,
        'maximum': 99,
        'step': 2,
    },
    {
        'name': 'clahe clip',
        'kind': 'float',
        'default': 2.0,
        'minimum': 0.1,
        'maximum': 10.0,
        'step': 0.1,
        'decimals': 2,
    },
    {
        'name': 'MSER max variation',
        'kind': 'float',
        'default': 1.2,
        'minimum': 0.0,
        'maximum': 10.0,
        'step': 0.05,
        'decimals': 2,
    },
    {
        'name': 'MSER delta',
        'kind': 'int',
        'default': 5,
        'minimum': 1,
        'maximum': 50,
        'step': 1,
    },
    {
        'name': 'MSER min area',
        'kind': 'int',
        'default': 30,
        'minimum': 1,
        'maximum': 5000,
        'step': 1,
    },
    {
        'name': 'MSER max area',
        'kind': 'int',
        'default': 15000,
        'minimum': 10,
        'maximum': 100000,
        'step': 100,
    },
    {
        'name': 'aspect ratio min',
        'kind': 'float',
        'default': 1.2,
        'minimum': 1.0,
        'maximum': 10.0,
        'step': 0.05,
        'decimals': 2,
    },
    {
        'name': 'clip-percentile',
        'kind': 'float',
        'default': 99.0,
        'minimum': 0.0,
        'maximum': 100.0,
        'step': 0.1,
        'decimals': 2,
    },
    {
        'name': 'area min',
        'kind': 'int',
        'default': 30,
        'minimum': 1,
        'maximum': 5000,
        'step': 1,
    },
    {
        'name': 'MSER threshold',
        'kind': 'float',
        'default': 85.0,
        'minimum': 0.0,
        'maximum': 100.0,
        'step': 0.1,
        'decimals': 2,
    },
]


PARAMETER_TOOLTIPS = {
    'solidity min': 'minimum filled-area ratio accepted for candidate ROIs.',
    'eccentricity min': 'minimum elongation; higher values favour fibre-like shapes.',
    'thinness max': 'upper compactness filter; lower values favour thin objects.',
    'tophat kernel': 'background-removal kernel size before MSER detection.',
    'clahe clip': 'local contrast equalisation strength.',
    'MSER max variation': 'MSER stability filter; lower values are stricter.',
    'MSER delta': 'MSER intensity step size.',
    'MSER min area': 'smallest MSER component area kept, in pixels.',
    'MSER max area': 'largest MSER component area kept, in pixels.',
    'aspect ratio min': 'minimum long-axis to short-axis ratio.',
    'clip-percentile': 'upper intensity percentile used before normalising.',
    'area min': 'smallest final ROI area kept, in pixels.',
    'MSER threshold': 'intensity percentile used to keep bright MSER candidates.',
}


def generate_distinct_colours(n):
    '''
    Generate visually distinct colours.
    '''
    colours = []
    if n <= 0:
        return colours

    hue = 0.0
    golden_ratio = 0.61803398875
    for _ in range(n):
        hue = (hue + golden_ratio) % 1.0
        rgb = colorsys.hsv_to_rgb(hue, 0.68, 0.92)
        colours.append(rgb)
    return colours


class OutputStream:
    def __init__(self, write_func):
        self.write_func = write_func

    def write(self, text):
        self.write_func(text)

    def flush(self):
        pass


class ZoomableCanvas(FigureCanvas):
    def __init__(self, figure, ax):
        super().__init__(figure)
        self.ax = ax
        self.zoom = 1.0
        self._base_xlim = self.ax.get_xlim()
        self._base_ylim = self.ax.get_ylim()
        self._drag_start_pos = None
        self._drag_start_xlim = None
        self._drag_start_ylim = None

    def mousePressEvent(self, event):
        if event.button() == Qt.RightButton:
            self._drag_start_pos = event.pos()
            self._drag_start_xlim = self.ax.get_xlim()
            self._drag_start_ylim = self.ax.get_ylim()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.RightButton and self._drag_start_pos is not None:
            dx = event.pos().x() - self._drag_start_pos.x()
            dy = event.pos().y() - self._drag_start_pos.y()

            trans = self.ax.transData.inverted()
            dx_data = trans.transform((0, 0))[0] - trans.transform((dx, 0))[0]
            dy_data = trans.transform((0, dy))[1] - trans.transform((0, 0))[1]

            new_xlim = (
                self._drag_start_xlim[0] + dx_data,
                self._drag_start_xlim[1] + dx_data,
            )
            new_ylim = (
                self._drag_start_ylim[0] + dy_data,
                self._drag_start_ylim[1] + dy_data,
            )

            self.ax.set_xlim(new_xlim)
            self.ax.set_ylim(new_ylim)
            self.draw()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.RightButton:
            self._drag_start_pos = None
        else:
            super().mouseReleaseEvent(event)

    def reset_view(self):
        if self.ax.images:
            img = self.ax.images[0]
            self.ax.set_xlim(0, img.get_array().shape[1])
            self.ax.set_ylim(img.get_array().shape[0], 0)
            self.draw()

    def wheelEvent(self, event):
        if self.ax.images:
            xmouse = event.position().x()
            ymouse = event.position().y()
            xdata, ydata = self.ax.transData.inverted().transform((xmouse, ymouse))

            cur_xlim = self.ax.get_xlim()
            cur_ylim = self.ax.get_ylim()
            zoom_factor = 1 / 1.2 if event.angleDelta().y() > 0 else 1.2

            xleft = xdata - (xdata - cur_xlim[0]) * zoom_factor
            xright = xdata + (cur_xlim[1] - xdata) * zoom_factor
            ytop = ydata - (ydata - cur_ylim[0]) * zoom_factor
            ybottom = ydata + (cur_ylim[1] - ydata) * zoom_factor

            self.ax.set_xlim([xleft, xright])
            self.ax.set_ylim([ytop, ybottom])
            self.draw()


def enhance_contrast_u8(img, tophat_kernel=11, clahe_clip=2.0):
    '''
    Apply a white top-hat and CLAHE enhancement.
    '''
    img = img.astype(np.float32)
    lo, hi = np.percentile(img, 1), np.percentile(img, 99)
    if hi <= lo:
        hi = lo + 1.0
    img01 = np.clip((img - lo) / (hi - lo), 0, 1)

    k = int(tophat_kernel)
    if k % 2 == 0:
        k += 1
    se = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    tophat = cv2.morphologyEx((img01 * 255).astype(np.uint8), cv2.MORPH_TOPHAT, se)

    clahe = cv2.createCLAHE(clipLimit=float(clahe_clip), tileGridSize=(8, 8))
    return clahe.apply(tophat)


class ROIEditor(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('FibreSegger v5.0')
        self.showMaximized()
        icon_path = Path(__file__).resolve().with_name('fibre-segmenter.ico')
        if icon_path.exists():
            self.setWindowIcon(QIcon(str(icon_path)))

        self.ref_image = None
        self.ref_image_path = None
        self.recname = ''
        self.roi_dict = {}
        self.selected = set()
        self.labelled = None
        self.fixed_ids = set()
        self.undo_stack = []
        self.show_overlay = True
        self.last_dir = os.getcwd()
        self.param_widgets = {}
        self._shortcuts = []
        self.dark_mode = False

        self._build_ui()
        self._register_shortcuts()
        sys.stdout = OutputStream(self.append_output)
        sys.stderr = OutputStream(self.append_output)
        self.refresh_status()

    def _theme(self):
        if self.dark_mode:
            return {
                'window': '#181a1b',
                'surface': '#242729',
                'surface_alt': '#2c3033',
                'surface_hover': '#34393d',
                'surface_strong': '#41484d',
                'border': '#4b5359',
                'border_strong': '#7b858c',
                'text': '#f1f3f4',
                'muted': '#bdc3c7',
                'primary': '#69a8d7',
                'primary_hover': '#82b8df',
                'primary_text': '#10202d',
                'danger_bg': '#3b2628',
                'danger_hover': '#503235',
                'danger_border': '#83565b',
                'danger_text': '#efc9cd',
                'selection': '#294f6e',
                'canvas': '#1c2022',
                }
        return {
            'window': '#dcdcdc',
            'surface': '#ffffff',
            'surface_alt': '#f0f0f0',
            'surface_hover': '#e8e8e8',
            'surface_strong': '#d4d4d4',
            'border': '#c6c6c6',
            'border_strong': '#8c8c8c',
            'text': '#000000',
            'muted': '#4f4f4f',
            'primary': '#0063b7',
            'primary_hover': '#004f92',
            'primary_text': '#ffffff',
            'danger_bg': '#f8e8e8',
            'danger_hover': '#f1d6d6',
            'danger_border': '#c78f8f',
            'danger_text': '#812f35',
            'selection': '#c6def1',
            'canvas': '#ffffff',
            }

    def _apply_palette_and_style(self):
        theme = self._theme()
        palette = self.palette()
        palette.setColor(QPalette.Window, QColor(theme['window']))
        palette.setColor(QPalette.Base, QColor(theme['surface_alt']))
        palette.setColor(QPalette.Button, QColor(theme['surface']))
        palette.setColor(QPalette.ButtonText, QColor(theme['text']))
        palette.setColor(QPalette.Text, QColor(theme['text']))
        palette.setColor(QPalette.WindowText, QColor(theme['text']))
        self.setPalette(palette)

        self.setStyleSheet(
            f'''
            QMainWindow {{
                background: {theme['window']};
            }}
            QGroupBox {{
                border: 1px solid {theme['border']};
                border-radius: 10px;
                margin-top: 16px;
                padding: 12px 10px 10px 10px;
                background: {theme['surface']};
                font-weight: 600;
                color: {theme['text']};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 12px;
                top: 1px;
                padding: 0 5px;
                background: {theme['window']};
            }}
            QPushButton, QToolButton {{
                border: 1px solid {theme['border']};
                border-radius: 7px;
                padding: 5px 10px;
                background: {theme['surface']};
                color: {theme['text']};
                min-height: 26px;
                font-weight: 600;
            }}
            QPushButton:hover, QToolButton:hover {{
                background: {theme['surface_hover']};
                border-color: {theme['border_strong']};
            }}
            QPushButton:pressed, QToolButton:pressed {{
                background: {theme['surface_strong']};
            }}
            QToolButton:checked {{
                background: {theme['primary']};
                border-color: {theme['primary']};
                color: {theme['primary_text']};
                font-weight: 600;
            }}
            QToolButton:checked:hover {{
                background: {theme['primary_hover']};
                border-color: {theme['primary_hover']};
            }}
            QPushButton[role="primary"] {{
                background: {theme['primary']};
                border-color: {theme['primary']};
                color: {theme['primary_text']};
                font-weight: 600;
            }}
            QPushButton[role="primary"]:hover {{
                background: {theme['primary_hover']};
                border-color: {theme['primary_hover']};
            }}
            QPushButton[role="danger"] {{
                background: {theme['danger_bg']};
                border-color: {theme['danger_border']};
                color: {theme['danger_text']};
            }}
            QPushButton[role="danger"]:hover {{
                background: {theme['danger_hover']};
                border-color: {theme['danger_text']};
            }}
            QDoubleSpinBox, QSpinBox {{
                border: 1px solid {theme['border']};
                border-radius: 7px;
                padding: 4px 8px;
                background: {theme['surface_alt']};
                color: {theme['text']};
                selection-background-color: {theme['selection']};
                min-height: 24px;
            }}
            QCheckBox {{
                spacing: 8px;
                color: {theme['text']};
                font-weight: 600;
            }}
            QCheckBox:hover {{
                color: {theme['primary']};
            }}
            QCheckBox::indicator {{
                width: 15px;
                height: 15px;
                border: 1px solid {theme['border_strong']};
                border-radius: 3px;
                background: {theme['surface']};
            }}
            QCheckBox::indicator:hover {{
                border-color: {theme['primary']};
                background: {theme['surface_hover']};
            }}
            QCheckBox::indicator:checked {{
                background: {theme['primary']};
                border-color: {theme['primary']};
            }}
            QTextEdit, QPlainTextEdit {{
                border: 1px solid {theme['border']};
                border-radius: 7px;
                background: {theme['surface']};
                color: {theme['text']};
                font-family: Consolas, 'Courier New', monospace;
                font-size: 9pt;
            }}
            QLabel {{
                color: {theme['text']};
            }}
            QLabel#stateSummary {{
                color: {theme['text']};
                font-weight: 600;
            }}
            QLabel#subtleLabel {{
                color: {theme['muted']};
            }}
            QLabel#panelValue {{
                color: {theme['text']};
                font-weight: 600;
                padding: 2px 0;
            }}
            QFrame#canvasFrame {{
                border: 1px solid {theme['border']};
                border-radius: 10px;
                background: {theme['canvas']};
            }}
            QScrollArea {{
                border: 0;
                background: transparent;
            }}
            QScrollArea QWidget {{
                background: transparent;
            }}
            QToolTip {{
                border: 1px solid {theme['border_strong']};
                border-radius: 6px;
                padding: 5px 7px;
                background: {theme['surface']};
                color: {theme['text']};
            }}
            QStatusBar {{
                background: {theme['surface']};
                color: {theme['muted']};
            }}
            QScrollBar:vertical {{
                border: none;
                background: {theme['surface']};
                width: 11px;
                margin: 2px 1px 2px 1px;
                border-radius: 4px;
            }}
            QScrollBar::handle:vertical {{
                background: {theme['border_strong']};
                min-height: 24px;
                border-radius: 4px;
            }}
            QScrollBar::handle:vertical:hover {{
                background: {theme['primary']};
            }}
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
                height: 0px;
            }}
            '''
        )

    def _make_button(self, text, callback, role='secondary', tooltip=None, icon=None):
        button = QPushButton(text)
        button.setProperty('role', role)
        button.clicked.connect(callback)
        button.setToolTip(tooltip or text)
        button.setMinimumHeight(28)
        button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        return button

    def _make_param_widget(self, spec):
        if spec['kind'] == 'int':
            widget = QSpinBox()
            widget.setRange(int(spec['minimum']), int(spec['maximum']))
            widget.setSingleStep(int(spec['step']))
            widget.setValue(int(spec['default']))
        else:
            widget = QDoubleSpinBox()
            widget.setDecimals(spec['decimals'])
            widget.setRange(float(spec['minimum']), float(spec['maximum']))
            widget.setSingleStep(float(spec['step']))
            widget.setValue(float(spec['default']))
        widget.setKeyboardTracking(False)
        widget.setAlignment(Qt.AlignRight)
        widget.setButtonSymbols(QAbstractSpinBox.NoButtons)
        widget.setMinimumWidth(88)
        widget.setMaximumWidth(104)
        widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        return widget

    def _get_params(self):
        params = {}
        for spec in PARAMETER_SPECS:
            widget = self.param_widgets[spec['name']]
            value = widget.value()
            params[spec['name']] = int(value) if spec['kind'] == 'int' else float(value)
        return params

    def _register_shortcuts(self):
        shortcuts = [
            ('Ctrl+O', self.load_image),
            ('Ctrl+L', self.load_roi_dict),
            ('Ctrl+S', self.save_roi_dict),
            ('Ctrl+R', self.run_segmentation),
        ]
        for sequence, slot in shortcuts:
            shortcut = QShortcut(QKeySequence(sequence), self)
            shortcut.activated.connect(slot)
            self._shortcuts.append(shortcut)

    def _build_ui(self):
        self._apply_palette_and_style()

        self.fig = Figure(dpi=100, facecolor=self._theme()['canvas'])
        self.ax = self.fig.add_subplot(111)
        self.ax.axis('off')
        self.canvas = ZoomableCanvas(self.fig, self.ax)
        self.canvas.setMinimumSize(720, 720)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.state_label = QLabel('No image loaded')
        self.state_label.setObjectName('stateSummary')
        self.state_label.setWordWrap(True)

        self.image_label = QLabel('Reference image: -')
        self.image_label.setObjectName('panelValue')
        self.image_label.setWordWrap(True)

        self.selection_label = QLabel('Selected ROIs: -')
        self.selection_label.setObjectName('panelValue')
        self.selection_label.setWordWrap(True)

        self.fixed_label = QLabel('Fixed ROIs: -')
        self.fixed_label.setObjectName('panelValue')
        self.fixed_label.setWordWrap(True)

        # left column
        session_box = QGroupBox('Session')
        session_layout = QVBoxLayout(session_box)

        session_buttons = QGridLayout()
        session_buttons.setHorizontalSpacing(8)
        session_buttons.setVerticalSpacing(8)
        session_buttons.addWidget(
            self._make_button(
                'Open Image',
                self.load_image,
                'primary',
                'Open a reference image and load the matching ROI dictionary if present.',
                self.style().standardIcon(QStyle.SP_DirOpenIcon),
            ),
            0,
            0,
        )
        session_buttons.addWidget(
            self._make_button(
                'Load ROI',
                self.load_roi_dict,
                tooltip='Load an ROI dictionary without changing the current image.',
                icon=self.style().standardIcon(QStyle.SP_FileDialogDetailedView),
            ),
            0,
            1,
        )
        session_buttons.addWidget(
            self._make_button(
                'Segment',
                self.run_segmentation,
                'primary',
                'Run the segmentation pipeline with the current parameters.',
                self.style().standardIcon(QStyle.SP_MediaPlay),
            ),
            1,
            0,
        )
        session_buttons.addWidget(
            self._make_button(
                'Save ROI',
                self.save_roi_dict,
                'primary',
                'Write the current ROI dictionary next to the reference image.',
                self.style().standardIcon(QStyle.SP_DialogSaveButton),
            ),
            1,
            1,
        )
        session_buttons.addWidget(
            self._make_button(
                'Reset View',
                self.reset_view,
                tooltip='Reset the canvas zoom and pan.',
                icon=self.style().standardIcon(QStyle.SP_BrowserReload),
            ),
            2,
            0,
        )
        session_buttons.addWidget(
            self._make_button(
                'Defaults',
                self.reset_parameters,
                tooltip='Restore the default segmentation parameters.',
                icon=self.style().standardIcon(QStyle.SP_BrowserReload),
            ),
            2,
            1,
        )

        session_layout.addLayout(session_buttons)
        session_layout.addWidget(self.image_label)

        params_box = QGroupBox('Segmentation parameters')
        params_scroll = QScrollArea()
        params_scroll.setWidgetResizable(True)
        params_scroll.setFrameShape(QFrame.NoFrame)
        params_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        params_container = QWidget()
        params_layout = QVBoxLayout(params_container)
        params_layout.setContentsMargins(6, 6, 6, 6)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignLeft)
        form.setFormAlignment(Qt.AlignTop)
        form.setVerticalSpacing(8)
        form.setHorizontalSpacing(18)
        for spec in PARAMETER_SPECS:
            label = QLabel(spec['name'])
            tooltip = PARAMETER_TOOLTIPS.get(spec['name'], '')
            if tooltip:
                label.setToolTip(tooltip)
            widget = self._make_param_widget(spec)
            if tooltip:
                widget.setToolTip(tooltip)
            self.param_widgets[spec['name']] = widget
            form.addRow(label, widget)

        params_layout.addLayout(form)
        params_layout.addStretch(1)
        params_scroll.setWidget(params_container)

        params_box_layout = QVBoxLayout(params_box)
        params_box_layout.addWidget(params_scroll)

        left_layout = QVBoxLayout()
        left_layout.setSpacing(12)
        left_layout.addWidget(session_box)
        left_layout.addWidget(params_box, 1)

        left_widget = QWidget()
        left_widget.setLayout(left_layout)
        left_widget.setMinimumWidth(380)

        # centre column
        canvas_frame = QFrame()
        canvas_frame.setObjectName('canvasFrame')
        canvas_layout = QVBoxLayout(canvas_frame)
        canvas_layout.setContentsMargins(10, 10, 10, 10)
        canvas_layout.addWidget(self.canvas)

        # right column
        state_box = QGroupBox('State')
        state_layout = QVBoxLayout(state_box)
        state_layout.addWidget(self.state_label)
        state_layout.addWidget(self.selection_label)
        state_layout.addWidget(self.fixed_label)

        selection_box = QGroupBox('Selection')
        selection_layout = QGridLayout(selection_box)
        selection_layout.setHorizontalSpacing(8)
        selection_layout.setVerticalSpacing(8)
        selection_layout.addWidget(
            self._make_button('Select All', self.select_all, tooltip='Select every ROI currently on the canvas.'),
            0,
            0,
        )
        selection_layout.addWidget(
            self._make_button('Invert', self.invert_selection, tooltip='Invert the current ROI selection.'),
            0,
            1,
        )
        selection_layout.addWidget(
            self._make_button(
                'clear selection',
                self.clear_selection,
                tooltip='Clear the current ROI selection.',
                icon=self.style().standardIcon(QStyle.SP_DialogCancelButton),
            ),
            1,
            0,
            1,
            2,
        )

        actions_box = QGroupBox('Actions')
        actions_layout = QGridLayout(actions_box)
        actions_layout.setHorizontalSpacing(8)
        actions_layout.setVerticalSpacing(8)
        actions_layout.addWidget(
            self._make_button(
                'Delete',
                self.delete_selected,
                'danger',
                'Delete the selected ROIs from the board.',
                self.style().standardIcon(QStyle.SP_TrashIcon),
            ),
            0,
            0,
        )
        actions_layout.addWidget(
            self._make_button('Merge', self.merge_selected, tooltip='Merge all selected ROIs into one.'),
            0,
            1,
        )
        actions_layout.addWidget(
            self._make_button('Fix', self.fix_selected, tooltip='Keep the selected ROIs during re-segmentation.'),
            1,
            0,
        )
        actions_layout.addWidget(
            self._make_button(
                'Unfix',
                self.remove_fixed,
                'danger',
                'Remove every fixed ROI from the board.',
                self.style().standardIcon(QStyle.SP_DialogCancelButton),
            ),
            1,
            1,
        )
        actions_layout.addWidget(
            self._make_button(
                'Clear All',
                self.clear_board,
                'danger',
                'Clear the board while keeping fixed ROIs if any are present.',
                self.style().standardIcon(QStyle.SP_TrashIcon),
            ),
            2,
            0,
        )
        actions_layout.addWidget(
            self._make_button(
                'Undo',
                self.undo_action,
                tooltip='Restore the previous board state.',
                icon=self.style().standardIcon(QStyle.SP_ArrowBack),
            ),
            2,
            1,
        )

        view_box = QGroupBox('View')
        view_layout = QGridLayout(view_box)
        view_layout.setHorizontalSpacing(8)
        view_layout.setVerticalSpacing(8)
        self.overlay_toggle = QCheckBox('Overlay')
        self.overlay_toggle.setChecked(True)
        self.overlay_toggle.stateChanged.connect(self.toggle_overlay)
        self.theme_group = QButtonGroup(self)
        self.theme_group.setExclusive(True)
        self.light_theme_toggle = QToolButton()
        self.light_theme_toggle.setText('Light')
        self.light_theme_toggle.setCheckable(True)
        self.light_theme_toggle.setChecked(True)
        self.light_theme_toggle.setToolTip('use the light UI theme.')
        self.light_theme_toggle.toggled.connect(self._set_light_theme)
        self.dark_theme_toggle = QToolButton()
        self.dark_theme_toggle.setText('Dark')
        self.dark_theme_toggle.setCheckable(True)
        self.dark_theme_toggle.setToolTip('use the dark UI theme.')
        self.dark_theme_toggle.toggled.connect(self._set_dark_theme)
        self.theme_group.addButton(self.light_theme_toggle)
        self.theme_group.addButton(self.dark_theme_toggle)
        view_layout.addWidget(self.overlay_toggle, 0, 0, 1, 2)
        view_layout.addWidget(self.light_theme_toggle, 1, 0)
        view_layout.addWidget(self.dark_theme_toggle, 1, 1)

        log_box = QGroupBox('run log')
        log_layout = QVBoxLayout(log_box)
        self.output_box = QPlainTextEdit()
        self.output_box.setReadOnly(True)
        self.output_box.setMinimumHeight(240)
        self.output_box.setLineWrapMode(QPlainTextEdit.WidgetWidth)
        self.output_box.document().setMaximumBlockCount(1500)
        log_layout.addWidget(self.output_box)

        right_layout = QVBoxLayout()
        right_layout.setSpacing(12)
        right_layout.addWidget(state_box)
        right_layout.addWidget(selection_box)
        right_layout.addWidget(actions_box)
        right_layout.addWidget(view_box)
        right_layout.addWidget(log_box, 1)

        right_widget = QWidget()
        right_widget.setLayout(right_layout)
        right_widget.setMinimumWidth(390)

        # main layout
        full_layout = QHBoxLayout()
        full_layout.setContentsMargins(16, 16, 16, 16)
        full_layout.setSpacing(16)
        full_layout.addWidget(left_widget, 1)
        full_layout.addWidget(right_widget, 1)
        full_layout.addWidget(canvas_frame, 4)

        container = QWidget()
        container.setLayout(full_layout)
        self.setCentralWidget(container)

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar_info = QLabel()
        self.status_bar_info.setObjectName('subtleLabel')
        self.status_bar.addPermanentWidget(self.status_bar_info, 1)

        self.canvas.mpl_connect('button_press_event', self.on_click)

    def append_output(self, text):
        self.output_box.moveCursor(QTextCursor.End)
        self.output_box.insertPlainText(text)
        self.output_box.ensureCursorVisible()

    def refresh_status(self, message=''):
        image_name = self.recname if self.ref_image is not None and self.recname else 'No image loaded'
        roi_count = len(self.roi_dict) if self.roi_dict is not None else 0
        selected_count = len(self.selected) if self.selected is not None else 0
        fixed_count = len(self.fixed_ids) if self.fixed_ids is not None else 0

        self.state_label.setText(
            f'ROIs: {roi_count} | Selected: {selected_count} | Fixed: {fixed_count}'
        )
        self.image_label.setText(f'Reference image: {image_name}')
        self.selection_label.setText(f'Selected ROIs: {self._format_id_list(self.selected)}')
        self.fixed_label.setText(f'Fixed ROIs: {self._format_id_list(self.fixed_ids)}')
        self.status_bar_info.setText(
            f'Image: {image_name} | ROIs: {roi_count} | Selected: {selected_count} | Fixed: {fixed_count}'
        )
        if message:
            self.statusBar().showMessage(message, 4000)

    def _format_id_list(self, ids):
        if not ids:
            return '-'
        ordered = sorted(ids)
        preview = ', '.join(str(x) for x in ordered[:4])
        if len(ordered) > 4:
            preview += f' +{len(ordered) - 4}'
        return preview

    def reset_parameters(self):
        for spec in PARAMETER_SPECS:
            self.param_widgets[spec['name']].setValue(spec['default'])
        self.refresh_status('Parameters reset')

    def _current_save_path(self):
        if not self.ref_image_path or not self.recname:
            return None
        save_dir = os.path.dirname(self.ref_image_path)
        return os.path.join(save_dir, f'{self.recname}_ROI_dict.npy')

    def push_undo_state(self):
        if self.labelled is None:
            return
        roi_dict_copy = {
            key: {'xpix': value['xpix'].copy(), 'ypix': value['ypix'].copy()}
            for key, value in self.roi_dict.items()
        }
        fixed_ids_copy = set(self.fixed_ids)
        self.undo_stack.append((self.labelled.copy(), roi_dict_copy, fixed_ids_copy))

    def select_all(self):
        if self.roi_dict:
            self.selected = set(self.roi_dict.keys())
            self.plot_image(preserve_view=True)
            self.refresh_status('Selection updated')

    def clear_selection(self):
        if self.selected:
            self.selected.clear()
            self.plot_image(preserve_view=True)
            self.refresh_status('Selection cleared')

    def invert_selection(self):
        if self.roi_dict:
            all_ids = set(self.roi_dict.keys())
            self.selected = all_ids - self.selected
            self.plot_image(preserve_view=True)
            self.refresh_status('Selection inverted')

    def delete_selected(self):
        if self.labelled is None or not self.selected:
            return
        self.push_undo_state()
        for roi_id in self.selected:
            self.labelled[self.labelled == roi_id] = 0
            self.fixed_ids.discard(roi_id)
        self.selected.clear()
        self.update_roi_dict()
        self.plot_image(preserve_view=True)
        self.refresh_status('Deleted selected ROIs')

    def fix_selected(self):
        if self.labelled is None or not self.selected:
            return
        newly_fixed = []
        for roi_id in self.selected:
            if roi_id in self.roi_dict and roi_id not in self.fixed_ids:
                self.fixed_ids.add(roi_id)
                newly_fixed.append(roi_id)
        if newly_fixed:
            print(f'fixed ROI(s): {newly_fixed} (will persist across segmentation)')
            self.plot_image(preserve_view=True)
            self.refresh_status('Fixed selected ROIs')

    def remove_fixed(self):
        if self.labelled is None or not self.fixed_ids:
            return
        self.push_undo_state()
        for fid in self.fixed_ids:
            self.labelled[self.labelled == fid] = 0
        self.fixed_ids.clear()
        self.selected.clear()
        self.update_roi_dict()
        self.plot_image(preserve_view=True)
        print('removed all fixed ROIs')
        self.refresh_status('Removed fixed ROIs')

    def clear_board(self):
        if self.labelled is None:
            return
        self.push_undo_state()
        if not self.fixed_ids:
            self.labelled = np.zeros_like(self.labelled, dtype=np.int32)
            self.roi_dict = {}
            self.selected.clear()
            self.fixed_ids.clear()
            self.plot_image(preserve_view=True)
            print('cleared all ROIs (no fixed ROIs present)')
            self.refresh_status('Board cleared')
            return

        new_labelled = np.zeros_like(self.labelled, dtype=np.int32)
        for fid in self.fixed_ids:
            new_labelled[self.labelled == fid] = fid
        self.labelled = new_labelled
        self.update_roi_dict()
        self.fixed_ids = {fid for fid in self.fixed_ids if fid in self.roi_dict}
        self.selected.clear()
        self.plot_image(preserve_view=True)
        print(f'cleared non-fixed ROIs; kept {len(self.fixed_ids)} fixed ROI(s)')
        self.refresh_status('Board cleared')

    def merge_selected(self):
        if self.labelled is None or len(self.selected) < 2:
            return
        self.push_undo_state()
        target_id = min(self.selected)
        for roi_id in self.selected:
            if roi_id != target_id:
                self.labelled[self.labelled == roi_id] = target_id
        if any(roi_id in self.fixed_ids for roi_id in self.selected):
            self.fixed_ids.add(target_id)
        self.fixed_ids = {
            fid for fid in self.fixed_ids
            if fid == target_id or fid not in self.selected
        }
        self.selected = {target_id}
        self.update_roi_dict()
        self.plot_image(preserve_view=True)
        self.refresh_status('Merged selected ROIs')

    def undo_action(self):
        if self.undo_stack:
            self.labelled, self.roi_dict, self.fixed_ids = self.undo_stack.pop()
            self.selected.clear()
            self.plot_image(preserve_view=True)
            self.refresh_status('Undo')

    def toggle_overlay(self):
        self.show_overlay = self.overlay_toggle.isChecked()
        self.plot_image(preserve_view=True)
        self.refresh_status('Overlay updated')

    def _set_light_theme(self, checked):
        if not checked:
            return
        self.dark_mode = False
        self._apply_palette_and_style()
        self.plot_image(preserve_view=True)
        self.refresh_status('Theme updated')

    def _set_dark_theme(self, checked):
        if not checked:
            return
        self.dark_mode = True
        self._apply_palette_and_style()
        self.plot_image(preserve_view=True)
        self.refresh_status('Theme updated')

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Space:
            self.overlay_toggle.setChecked(not self.overlay_toggle.isChecked())
            return
        if event.key() == Qt.Key_Delete:
            self.delete_selected()
            return
        if event.key() == Qt.Key_Backspace:
            self.undo_action()
            return
        if event.key() == Qt.Key_M and (event.modifiers() & Qt.ControlModifier):
            self.merge_selected()
            return
        if event.key() == Qt.Key_F and (event.modifiers() & Qt.ControlModifier):
            if event.modifiers() & Qt.ShiftModifier:
                self.remove_fixed()
            else:
                self.fix_selected()
            return
        super().keyPressEvent(event)

    def load_image(self):
        start_dir = self.last_dir if os.path.isdir(self.last_dir) else os.getcwd()
        path, _ = QFileDialog.getOpenFileName(
            self,
            'Open NPY Image',
            start_dir,
            'NumPy image (*.npy)',
        )
        if not path:
            return

        self.last_dir = os.path.dirname(path)
        self.ref_image_path = path
        self.recname = os.path.basename(path).split('_ref_mat')[0]
        self.ref_image = np.load(path)

        self.roi_dict.clear()
        self.labelled = None
        self.selected.clear()
        self.fixed_ids.clear()
        self.undo_stack.clear()

        self.plot_image()
        self.canvas.reset_view()
        print(f'{self.recname} loaded')
        self.refresh_status('Image loaded')

        roi_dict_path = os.path.join(os.path.dirname(path), f'{self.recname}_ROI_dict.npy')
        if os.path.exists(roi_dict_path):
            try:
                roi_dict = np.load(roi_dict_path, allow_pickle=True).item()
                if isinstance(roi_dict, dict):
                    self.roi_dict = roi_dict
                    self.labelled = np.zeros_like(self.ref_image, dtype=np.int32)
                    for roi_id, coords in self.roi_dict.items():
                        self.labelled[coords['ypix'], coords['xpix']] = roi_id
                    self.selected.clear()
                    self.undo_stack.clear()
                    self.fixed_ids = set(self.roi_dict.keys())
                    self.plot_image()
                    self.canvas.reset_view()
                    print(f'ROI dict loaded automatically from {roi_dict_path}')
                    self.refresh_status('ROI dict loaded')
                    return
            except Exception as e:
                print(f'failed to load ROI dict: {e}')
                self.refresh_status('ROI dict load failed')
                return

        self.run_segmentation()

    def load_roi_dict(self):
        if self.ref_image is None:
            print('please load a reference image first')
            self.refresh_status('Load a reference image first')
            return

        start_dir = self.last_dir if os.path.isdir(self.last_dir) else os.getcwd()
        path, _ = QFileDialog.getOpenFileName(
            self,
            'Load ROI Dict',
            start_dir,
            'NumPy dict (*.npy)',
        )
        if not path:
            return

        self.last_dir = os.path.dirname(path)
        roi_dict = np.load(path, allow_pickle=True).item()
        if not isinstance(roi_dict, dict):
            print('invalid file format.')
            self.refresh_status('Invalid ROI dict')
            return

        self.roi_dict = roi_dict
        self.labelled = np.zeros_like(self.ref_image, dtype=np.int32)
        for roi_id, coords in self.roi_dict.items():
            self.labelled[coords['ypix'], coords['xpix']] = roi_id
        self.selected.clear()
        self.fixed_ids.clear()
        self.undo_stack.clear()
        self.plot_image()
        self.canvas.reset_view()
        print(f'ROI dict loaded from {path}')
        self.refresh_status('ROI dict loaded')

    def save_roi_dict(self):
        save_path = self._current_save_path()
        if save_path is None:
            print('reference image not loaded, cannot auto-save.')
            self.refresh_status('Load a reference image first')
            return
        np.save(save_path, self.roi_dict)
        print(f'saved: {save_path}')
        self.refresh_status('ROI dict saved')

    def on_click(self, event):
        if event.inaxes != self.ax or self.labelled is None:
            return
        if event.xdata is None or event.ydata is None:
            return
        x, y = int(event.xdata), int(event.ydata)
        if y < 0 or x < 0 or y >= self.labelled.shape[0] or x >= self.labelled.shape[1]:
            return
        roi_id = self.labelled[y, x]
        if roi_id > 0:
            modifiers = Qt.NoModifier
            if getattr(event, 'guiEvent', None) is not None:
                modifiers = event.guiEvent.modifiers()
            if modifiers & Qt.ShiftModifier:
                if roi_id in self.selected:
                    self.selected.remove(roi_id)
                else:
                    self.selected.add(roi_id)
            else:
                self.selected = {roi_id}
            self.plot_image(preserve_view=True)
            self.refresh_status('Selection updated')

    def plot_image(self, preserve_view=False):
        params = self._get_params()

        if preserve_view:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()

        self.ax.clear()
        self.fig.set_facecolor(self._theme()['canvas'])
        self.ax.set_facecolor(self._theme()['canvas'])
        self.ax.axis('off')

        if self.ref_image is None:
            self.canvas.draw()
            self.refresh_status()
            return

        lo = np.percentile(self.ref_image, 1)
        hi = np.percentile(self.ref_image, params['clip-percentile'])
        self.ax.imshow(self.ref_image, cmap='gray', vmin=lo, vmax=hi)

        if self.labelled is not None and self.show_overlay:
            overlay = np.zeros((*self.labelled.shape, 4))
            ids = np.unique(self.labelled)
            ids = ids[ids > 0]
            colours = generate_distinct_colours(len(ids))
            for idx, roi_id in enumerate(ids):
                mask = self.labelled == roi_id
                overlay[mask, :3] = colours[idx]
                alpha = 0.5
                if roi_id in self.fixed_ids:
                    alpha = 0.8
                overlay[mask, 3] = alpha
            self.ax.imshow(overlay)
            for roi_id in self.selected:
                coords = np.column_stack(np.where(self.labelled == roi_id))
                self.ax.plot(coords[:, 1], coords[:, 0], 'c.', markersize=1)

        if preserve_view:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)

        self.canvas.draw()
        self.refresh_status()

    def run_segmentation(self):
        if self.ref_image is None:
            self.refresh_status('Load a reference image first')
            return

        p = self._get_params()

        if self.labelled is not None:
            self.push_undo_state()

        fixed_rois = {}
        if self.roi_dict and self.fixed_ids:
            for fid in list(self.fixed_ids):
                coords = self.roi_dict.get(fid, None)
                if coords is not None:
                    fixed_rois[fid] = {
                        'xpix': coords['xpix'].copy(),
                        'ypix': coords['ypix'].copy(),
                    }
                else:
                    self.fixed_ids.discard(fid)

        ref_f = median_filter(self.ref_image, size=(3, 3))
        img_u8 = enhance_contrast_u8(
            ref_f,
            tophat_kernel=p['tophat kernel'],
            clahe_clip=p['clahe clip'],
        )

        thr = np.percentile(img_u8, p['MSER threshold'])
        soft = img_u8.copy()
        soft[soft < thr] = 0

        delta = int(max(1, round(p['MSER delta'])))
        min_area = int(max(5, round(p['MSER min area'])))
        max_area = int(max(min_area + 1, round(p['MSER max area'])))
        mser = cv2.MSER_create(delta, min_area, max_area)
        mser.setMaxVariation(p['MSER max variation'])

        regions, _ = mser.detectRegions(soft)

        lab = np.zeros_like(img_u8, dtype=np.int32)
        for i, reg in enumerate(regions):
            lab[reg[:, 1], reg[:, 0]] = i + 1

        props = regionprops(lab)
        seg_rois = []
        for r in props:
            area = r.area
            if area < p['area min']:
                continue

            ecc = r.eccentricity if np.isfinite(r.eccentricity) else 0.0
            sol = r.solidity if np.isfinite(r.solidity) else 0.0
            maj = r.major_axis_length
            minax = r.minor_axis_length if r.minor_axis_length > 1e-6 else 1e-6
            ar = maj / minax
            perim = r.perimeter if r.perimeter > 1e-6 else 1e-6
            thin = 4 * np.pi * area / (perim ** 2)

            if sol < p['solidity min']:
                continue
            if ecc < p['eccentricity min']:
                continue
            if ar < p['aspect ratio min']:
                continue
            if thin > p['thinness max']:
                continue

            ypix, xpix = r.coords[:, 0], r.coords[:, 1]
            seg_rois.append({'xpix': xpix, 'ypix': ypix})

        self.labelled = np.zeros_like(self.ref_image, dtype=np.int32)
        self.roi_dict = {}
        new_fixed_ids = set()
        roi_id = 1

        for _, coords in fixed_rois.items():
            ypix = coords['ypix']
            xpix = coords['xpix']
            self.roi_dict[roi_id] = {'xpix': xpix, 'ypix': ypix}
            self.labelled[ypix, xpix] = roi_id
            new_fixed_ids.add(roi_id)
            roi_id += 1

        for roi in seg_rois:
            ypix = roi['ypix']
            xpix = roi['xpix']
            if np.any(self.labelled[ypix, xpix] > 0):
                continue
            self.roi_dict[roi_id] = {'xpix': xpix, 'ypix': ypix}
            self.labelled[ypix, xpix] = roi_id
            roi_id += 1

        self.fixed_ids = new_fixed_ids
        self.selected.clear()
        self.plot_image()
        self.canvas.reset_view()
        print(
            f'MSER regions: {len(regions)} | kept ROIs: {len(self.roi_dict)} '
            f'(including {len(self.fixed_ids)} fixed)'
        )
        self.refresh_status('Segmentation complete')

    def update_roi_dict(self):
        if self.labelled is None:
            self.roi_dict = {}
            self.fixed_ids = set()
            return

        props = regionprops(self.labelled)
        self.roi_dict = {}
        for region in props:
            ypix, xpix = region.coords[:, 0], region.coords[:, 1]
            self.roi_dict[region.label] = {'xpix': xpix, 'ypix': ypix}
        self.fixed_ids = {fid for fid in self.fixed_ids if fid in self.roi_dict}

    def reset_view(self):
        self.canvas.reset_view()
        self.refresh_status('View reset')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    app.setFont(QFont('Segoe UI', 9))
    editor = ROIEditor()
    editor.show()
    sys.exit(app.exec_())
