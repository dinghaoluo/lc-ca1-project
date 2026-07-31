'''
Created on 12 May 2026
Modified on 13 May 2026 while wiring the MSER controls into the workbench
Modified on 21 May 2026 to separate labelling, training, prediction, and editing
Modified on 2 June 2026 during the first command-line and diagnostics pass
Modified on 23 June 2026 to bring the MSER editor and model workflow together
fibre-sight workbench for labelling, training, prediction, and ROI curation

MSER labels and model predictions return to the same editable xpix/ypix ROI
dictionary so hand curation stays in the loop

@author: Dinghao Luo
'''


#%% imports
from pathlib import Path
import sys

import matplotlib

matplotlib.use('Qt5Agg')

import numpy as np
from matplotlib import font_manager
from matplotlib.figure import Figure
from PyQt5.QtCore import QProcess, Qt
from PyQt5.QtGui import QColor, QCursor, QFont, QFontDatabase, QIcon, QKeySequence, QPalette, QTextCursor
from PyQt5.QtWidgets import (
    QAbstractSpinBox,
    QApplication,
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPushButton,
    QPlainTextEdit,
    QShortcut,
    QSizePolicy,
    QSplitter,
    QSpinBox,
    QTabWidget,
    QToolTip,
    QVBoxLayout,
    QWidget,
    )

if __package__ in {None, ''}:
    repo_root = Path(__file__).resolve().parents[3]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from preprocessing.imaging.fibre_sight._repo import default_output_root, default_source_root, get_repo_root
from preprocessing.imaging.fibre_sight.api import AxonROIPredictor, get_default_checkpoint, get_model_entry
from preprocessing.imaging.fibre_sight.config import load_config, save_config
from preprocessing.imaging.fibre_sight.gui_canvas import (
    ZoomableCanvas,
    generate_distinct_colours,
    normalise_for_display,
    squeeze_image,
    )
from preprocessing.imaging.fibre_sight.mser_segmenter import (
    PARAMETER_SPECS,
    PARAMETER_TOOLTIPS,
    run_mser_segmentation,
    )
from preprocessing.imaging.fibre_sight.postprocess import probability_to_roi_dict
from preprocessing.imaging.fibre_sight.roi_io import labels_to_roi_dict, load_roi_dict, roi_dict_to_label, save_roi_dict


REPO_ROOT = get_repo_root()
SCRIPT_DIR = Path(__file__).resolve().parent
APP_ICON_PATH = SCRIPT_DIR / 'assets' / 'fibresight_icon.ico'
MONONOKI_FONT_DIR = SCRIPT_DIR / 'assets' / 'fonts' / 'mononoki'
MONONOKI_FONT_FAMILY = 'mononoki'
CORE_SEGMENT_PARAMETERS = {
    'MSER threshold',
    'MSER min area',
    'MSER max area',
    'area min',
    'eccentricity min',
    'aspect ratio min',
    'tophat kernel',
    }


#%% helpers
def load_gui_font(point_size=9):
    for font_path in sorted(MONONOKI_FONT_DIR.glob('*.ttf')):
        QFontDatabase.addApplicationFont(str(font_path))
        font_manager.fontManager.addfont(str(font_path))

    matplotlib.rcParams['font.family'] = MONONOKI_FONT_FAMILY
    matplotlib.rcParams['font.monospace'] = [MONONOKI_FONT_FAMILY, 'Consolas', 'Courier New']
    font = QFont(MONONOKI_FONT_FAMILY, point_size)
    font.setStyleHint(QFont.Monospace)
    return font


#%% main window
class FibreSightWorkbench(QMainWindow):
    def __init__(self):
        super().__init__()
        app = QApplication.instance()
        if app is not None:
            app.setFont(load_gui_font(9))
        self.setWindowTitle('fibre-sight')
        self.setWindowIcon(QIcon(str(APP_ICON_PATH)))
        self.resize(1500, 950)
        self.showMaximized()

        self.ref_image = None
        self.image_path = None
        self.recname = None
        self.roi_dict = {}
        self.labelled = None
        self.selected = set()
        self.fixed_ids = set()
        self.undo_stack = []
        self.probability = None
        self.predictor = None
        self.process = None
        self.current_process_name = None
        self.pending_checkpoint_path = None
        self.last_saved_model_path = None
        self.dark_mode = False

        self._build_widgets()
        self._build_layout()
        self._apply_palette_and_style()
        self._connect_shortcuts()
        self.plot_image()
        self.refresh_status('ready')

    #%% setup
    def _build_widgets(self):
        self.fig = Figure(dpi=100, facecolor=self._theme()['canvas'])
        self.ax = self.fig.add_subplot(111)
        self.ax.axis('off')
        self.canvas = ZoomableCanvas(self.fig, self.ax)
        self.canvas.setMinimumSize(760, 760)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.mpl_connect('button_press_event', self.on_click)

        self.state_label = QLabel('no image loaded')
        self.state_label.setObjectName('stateSummary')
        self.roi_label = QLabel('ROIs: 0 | selected: 0')
        self.roi_label.setObjectName('panelValue')
        self.model_label = QLabel('model: not loaded')
        self.model_label.setObjectName('panelValue')
        self.fixed_label = QLabel('fixed: 0')
        self.fixed_label.setObjectName('panelValue')

        self.tabs = QTabWidget()
        self.no_model_tab = QWidget()
        self.model_tab = QWidget()
        self.tabs.addTab(self.no_model_tab, 'MSER + training')
        self.tabs.addTab(self.model_tab, 'predict + curate')
        self.tabs.setTabToolTip(0, 'make hand-corrected MSER labels and train from them')
        self.tabs.setTabToolTip(1, 'predict ROIs, correct them, and save the next labels')
        self.tabs.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.tabs.tabBar().setExpanding(False)
        self.tabs.tabBar().setUsesScrollButtons(False)
        self.tabs.currentChanged.connect(lambda _: self.refresh_status())
        self._build_training_widgets()
        self._build_prediction_widgets()
        self._build_segment_widgets()
        self._build_editing_widgets()
        self.output_box = self.make_log_box()

        self.roi_overlay_check = QCheckBox('ROI on')
        self.roi_overlay_check.setChecked(True)
        self.roi_overlay_check.setToolTip('show or hide ROI outlines on the raw image')
        self.roi_overlay_check.stateChanged.connect(self.set_roi_overlay_visible)
        self.dark_mode_check = QCheckBox('dark mode')
        self.dark_mode_check.setToolTip('switch between light and dark UI themes')
        self.dark_mode_check.stateChanged.connect(self.set_dark_mode)

    def _build_training_widgets(self):
        self.source_root_line = QLineEdit(str(default_source_root()))
        self.manifest_line = QLineEdit(str(default_output_root() / 'manifests' / 'training_manifest.csv'))
        self.config_line = QLineEdit(
            str(SCRIPT_DIR / 'configs' / 'dlight_hpc_lc_opto_unet.yaml')
            )
        self.run_name_line = QLineEdit('fibre_sight_unet')
        self.source_root_line.textChanged.connect(lambda _: self.refresh_status())
        self.manifest_line.textChanged.connect(lambda _: self.refresh_status())
        self.config_line.textChanged.connect(lambda _: self.refresh_status())
        self.run_name_line.textChanged.connect(lambda _: self.refresh_status())

        self.val_fraction_spin = QDoubleSpinBox()
        self.val_fraction_spin.setRange(0, 0.5)
        self.val_fraction_spin.setSingleStep(0.01)
        self.val_fraction_spin.setDecimals(2)
        self.val_fraction_spin.setValue(0.15)
        self.prepare_value_control(self.val_fraction_spin)

        self.test_fraction_spin = QDoubleSpinBox()
        self.test_fraction_spin.setRange(0, 0.5)
        self.test_fraction_spin.setSingleStep(0.01)
        self.test_fraction_spin.setDecimals(2)
        self.test_fraction_spin.setValue(0.15)
        self.prepare_value_control(self.test_fraction_spin)

        self.epochs_spin = QSpinBox()
        self.epochs_spin.setRange(1, 500)
        self.epochs_spin.setValue(80)
        self.prepare_value_control(self.epochs_spin)

        self.build_manifest_button = QPushButton('scan labelled sessions')
        self.train_model_button = QPushButton('train model')
        self.evaluate_model_button = QPushButton('score model')
        self.stop_process_button = QPushButton('stop process')
        self.inspect_manifest_button = QPushButton('dataset summary')
        self.preview_training_button = QPushButton('preview labels')
        self.preview_predictions_button = QPushButton('preview predictions')
        self.training_diagnostics_button = QPushButton('diagnostics')
        self.training_diagnostics_popup = None
        self.training_advanced_button = QPushButton('advanced options')
        self.training_advanced_popup = None
        self.set_button_role(self.train_model_button, 'primary')
        self.set_button_role(self.stop_process_button, 'danger')
        self.set_button_role(self.training_diagnostics_button, 'quiet')
        self.set_button_role(self.training_advanced_button, 'quiet')
        self.set_button_role(self.inspect_manifest_button, 'quiet')
        self.set_button_role(self.evaluate_model_button, 'quiet')
        self.set_button_role(self.preview_training_button, 'quiet')
        self.set_button_role(self.preview_predictions_button, 'quiet')
        self.evaluate_model_button.setToolTip('score the current trained model on held-out labelled sessions')
        self.preview_predictions_button.setToolTip('save example overlays comparing model ROIs with held-out labels')
        self.training_diagnostics_button.setToolTip('open model scoring and preview tools')
        self.training_advanced_button.setToolTip('open dataset table and training split controls')

        self.build_manifest_button.clicked.connect(self.build_manifest)
        self.train_model_button.clicked.connect(self.train_model)
        self.evaluate_model_button.clicked.connect(self.evaluate_model)
        self.inspect_manifest_button.clicked.connect(self.inspect_manifest)
        self.preview_training_button.clicked.connect(self.preview_training_labels)
        self.preview_predictions_button.clicked.connect(self.preview_model_predictions)
        self.training_diagnostics_button.clicked.connect(self.show_training_diagnostics_popup)
        self.training_advanced_button.clicked.connect(self.show_training_advanced_popup)
        self.stop_process_button.clicked.connect(self.stop_process)

    def _build_prediction_widgets(self):
        self.image_line = QLineEdit()
        self.checkpoint_line = QLineEdit(str(get_default_checkpoint()))
        self.image_line.textChanged.connect(lambda _: self.refresh_status())
        self.checkpoint_line.textChanged.connect(lambda _: self.refresh_status())

        model_entry = get_model_entry()
        # these stayed visible because they were the useful controls during tuning
        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setRange(0.01, 0.99)
        self.threshold_spin.setSingleStep(0.01)
        self.threshold_spin.setDecimals(2)
        self.threshold_spin.setValue(float(model_entry['threshold']))
        self.prepare_value_control(self.threshold_spin)

        self.min_size_spin = QSpinBox()
        self.min_size_spin.setRange(1, 100000)
        self.min_size_spin.setValue(int(model_entry['min_size']))
        self.prepare_value_control(self.min_size_spin)

        self.show_probability_check = QCheckBox('show model confidence')
        self.show_probability_check.setToolTip('show the model confidence map behind the ROI outlines')
        self.show_probability_check.stateChanged.connect(lambda _: self.plot_image(preserve_view=True))

        self.load_image_button = QPushButton('load chan2 image')
        self.load_model_button = QPushButton('load trained model')
        self.predict_button = QPushButton('predict')
        self.set_button_role(self.load_image_button, 'primary')
        self.set_button_role(self.load_model_button, 'quiet')

        for widget, tip in [
            (self.image_line, 'channel-2 reference image used for ROI prediction'),
            (self.checkpoint_line, 'saved trained model file (.pt)'),
            (self.threshold_spin, 'higher values keep only stronger model detections'),
            (self.min_size_spin, 'smallest ROI size to keep'),
            ]:
            widget.setToolTip(tip)

        self.threshold_spin.editingFinished.connect(self.apply_probability_threshold_from_controls)
        self.min_size_spin.editingFinished.connect(self.apply_probability_threshold_from_controls)
        self.load_image_button.clicked.connect(self.load_channel_image)
        self.load_model_button.clicked.connect(self.load_model)
        self.predict_button.clicked.connect(self.predict_rois)


    def _build_segment_widgets(self):
        self.segment_param_widgets = {}
        for spec in PARAMETER_SPECS:
            self.segment_param_widgets[spec['name']] = self._make_segment_param_widget(spec)

        self.segment_button = QPushButton('segment')
        self.reset_segment_button = QPushButton('reset params')
        self.advanced_segment_button = QPushButton('advanced parameters')
        self.fix_selected_button = QPushButton('fix selected')
        self.unfix_selected_button = QPushButton('unfix selected')
        self.clear_fixed_button = QPushButton('clear fixed')
        self.clear_unfixed_button = QPushButton('clear unfixed')
        self.segment_load_image_button = QPushButton('load chan2 image')
        self.segment_load_roi_button = QPushButton('load ROI dict')
        self.segment_save_roi_button = QPushButton('save ROI dict')
        self.segment_advanced_popup = None

        self.set_button_role(self.segment_button, 'primary')
        self.set_button_role(self.segment_load_image_button, 'primary')
        self.set_button_role(self.segment_load_roi_button, 'quiet')
        self.set_button_role(self.segment_save_roi_button, 'primary')
        self.set_button_role(self.reset_segment_button, 'quiet')
        self.set_button_role(self.advanced_segment_button, 'quiet')
        self.set_button_role(self.fix_selected_button, 'quiet')
        self.set_button_role(self.unfix_selected_button, 'quiet')
        self.set_button_role(self.clear_fixed_button, 'quiet')
        self.set_button_role(self.clear_unfixed_button, 'danger')

        self.segment_button.setToolTip('run MSER segmentation with the current parameter values')
        self.reset_segment_button.setToolTip('restore the default MSER parameter values')
        self.fix_selected_button.setToolTip('keep selected ROIs when segmentation is run again')
        self.unfix_selected_button.setToolTip('allow selected ROIs to change during segmentation')
        self.clear_fixed_button.setToolTip('remove all fixed marks without deleting ROIs')
        self.clear_unfixed_button.setToolTip('remove all ROIs except fixed ROIs')
        self.advanced_segment_button.setToolTip('open less commonly changed MSER parameters')

        self.segment_button.clicked.connect(self.segment_rois)
        self.reset_segment_button.clicked.connect(self.reset_segment_parameters)
        self.fix_selected_button.clicked.connect(self.fix_selected)
        self.unfix_selected_button.clicked.connect(self.unfix_selected)
        self.clear_fixed_button.clicked.connect(self.clear_fixed)
        self.clear_unfixed_button.clicked.connect(self.clear_unfixed)
        self.advanced_segment_button.clicked.connect(self.show_segment_advanced_popup)
        self.segment_load_image_button.clicked.connect(self.load_channel_image)
        self.segment_load_roi_button.clicked.connect(self.load_roi_file)
        self.segment_save_roi_button.clicked.connect(self.save_roi_file)

    def _build_editing_widgets(self):
        self.curate_buttons = {
            'select_all': [],
            'delete': [],
            'merge': [],
            'undo': [],
            'reset_view': [],
            'save_roi': [],
        }
        self.segment_curate_buttons = self._make_curate_button_set()
        self.predict_curate_buttons = self._make_curate_button_set(save_label='save ROI dict')

    def _make_curate_button_set(self, save_label=None):
        buttons = {
            'select_all': QPushButton('select all'),
            'delete': QPushButton('delete selected'),
            'merge': QPushButton('merge selected'),
            'undo': QPushButton('undo'),
            'reset_view': QPushButton('reset view'),
        }
        if save_label is not None:
            buttons['save_roi'] = QPushButton(save_label)
            self.set_button_role(buttons['save_roi'], 'primary')

        self.set_button_role(buttons['delete'], 'danger')
        self.set_button_role(buttons['select_all'], 'quiet')
        self.set_button_role(buttons['merge'], 'quiet')
        self.set_button_role(buttons['undo'], 'quiet')
        self.set_button_role(buttons['reset_view'], 'quiet')

        buttons['select_all'].clicked.connect(self.select_all)
        buttons['delete'].clicked.connect(self.delete_selected)
        buttons['merge'].clicked.connect(self.merge_selected)
        buttons['undo'].clicked.connect(self.undo)
        buttons['reset_view'].clicked.connect(self.reset_view)
        if 'save_roi' in buttons:
            buttons['save_roi'].clicked.connect(self.save_roi_file)

        for name, button in buttons.items():
            if name in self.curate_buttons:
                self.curate_buttons[name].append(button)

        return buttons

    def _build_layout(self):
        canvas_layout = QVBoxLayout()
        canvas_layout.setContentsMargins(8, 8, 8, 8)
        canvas_layout.addWidget(self.canvas)
        canvas_frame = QFrame()
        canvas_frame.setObjectName('canvasFrame')
        canvas_frame.setLayout(canvas_layout)

        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(10)
        right_layout.addWidget(canvas_frame, 1)
        right_layout.addWidget(self.output_box, 0)
        right_pane = QWidget()
        right_pane.setObjectName('mainPane')
        right_pane.setLayout(right_layout)

        status_card = QFrame()
        status_card.setObjectName('statusCard')
        status_card.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        status_card.setMaximumHeight(132)
        status_layout = QGridLayout(status_card)
        status_layout.setContentsMargins(10, 8, 10, 8)
        status_layout.setHorizontalSpacing(12)
        status_layout.setVerticalSpacing(5)
        status_layout.addWidget(self.state_label, 0, 0, 1, 2)
        status_layout.addWidget(self.roi_label, 1, 0)
        status_layout.addWidget(self.fixed_label, 1, 1)
        status_layout.addWidget(self.model_label, 2, 0, 1, 2)

        view_options = QWidget()
        view_options_layout = QHBoxLayout(view_options)
        view_options_layout.setContentsMargins(0, 0, 0, 0)
        view_options_layout.setSpacing(12)
        view_options_layout.addWidget(self.roi_overlay_check)
        view_options_layout.addWidget(self.dark_mode_check)
        view_options_layout.addStretch(1)
        status_layout.addWidget(view_options, 3, 0, 1, 2)

        self._layout_no_model_tab()
        self._layout_model_tab()

        controls_layout = QVBoxLayout()
        controls_layout.setContentsMargins(0, 0, 0, 0)
        controls_layout.setSpacing(10)
        controls_layout.addWidget(status_card)
        controls_layout.addWidget(self.tabs, 1)

        controls_widget = QWidget()
        controls_widget.setObjectName('controlsPanel')
        controls_widget.setLayout(controls_layout)
        controls_widget.setMinimumWidth(380)
        controls_widget.setMaximumWidth(440)

        splitter = QSplitter(Qt.Horizontal)
        splitter.setObjectName('mainSplitter')
        splitter.addWidget(controls_widget)
        splitter.addWidget(right_pane)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setChildrenCollapsible(False)
        splitter.setSizes([390, 1250])

        main_layout = QHBoxLayout()
        main_layout.setContentsMargins(14, 14, 14, 14)
        main_layout.setSpacing(0)
        main_layout.addWidget(splitter)

        container = QWidget()
        container.setObjectName('centralWidget')
        container.setLayout(main_layout)
        self.setCentralWidget(container)

    def _layout_no_model_tab(self):
        layout = QVBoxLayout(self.no_model_tab)
        layout.setSpacing(6)
        layout.setContentsMargins(8, 8, 8, 8)
        self._layout_mser_section(layout)
        self._layout_training_section(layout)
        layout.addStretch(1)

    def _layout_model_tab(self):
        layout = QVBoxLayout(self.model_tab)
        layout.setSpacing(6)
        layout.setContentsMargins(8, 8, 8, 8)
        self._layout_prediction_section(layout)
        self._layout_prediction_curation_section(layout)
        layout.addStretch(1)

    def _layout_mser_section(self, layout):
        layout.addWidget(self.make_section_label('MSER + curation'))

        io_buttons = QGridLayout()
        io_buttons.setHorizontalSpacing(6)
        io_buttons.setVerticalSpacing(6)
        io_buttons.addWidget(self.segment_load_image_button, 0, 0)
        io_buttons.addWidget(self.segment_load_roi_button, 0, 1)
        io_buttons.addWidget(self.segment_save_roi_button, 0, 2)
        layout.addLayout(io_buttons)

        main_specs = [spec for spec in PARAMETER_SPECS if spec['name'] in CORE_SEGMENT_PARAMETERS]
        param_grid = QGridLayout()
        param_grid.setHorizontalSpacing(8)
        param_grid.setVerticalSpacing(4)
        advanced_row = 0
        for idx, spec in enumerate(main_specs):
            tooltip = PARAMETER_TOOLTIPS.get(spec['name'], '')
            widget = self.segment_param_widgets[spec['name']]
            widget.setToolTip(tooltip)
            row = idx // 2
            col = (idx % 2) * 2
            param_grid.addWidget(self.make_form_label(spec['name'], tooltip), row, col)
            param_grid.addWidget(widget, row, col + 1)
            if spec['name'] == 'MSER threshold':
                advanced_row = row
        param_grid.addWidget(self.advanced_segment_button, advanced_row, 2, 1, 2)
        layout.addLayout(param_grid)

        buttons = QGridLayout()
        buttons.setHorizontalSpacing(6)
        buttons.setVerticalSpacing(6)
        buttons.addWidget(self.segment_button, 0, 0)
        buttons.addWidget(self.reset_segment_button, 0, 1)
        buttons.addWidget(self.fix_selected_button, 1, 0)
        buttons.addWidget(self.unfix_selected_button, 1, 1)
        buttons.addWidget(self.segment_curate_buttons['select_all'], 1, 2)
        buttons.addWidget(self.segment_curate_buttons['delete'], 2, 0)
        buttons.addWidget(self.segment_curate_buttons['merge'], 2, 1)
        buttons.addWidget(self.segment_curate_buttons['undo'], 2, 2)
        buttons.addWidget(self.clear_fixed_button, 3, 0)
        buttons.addWidget(self.clear_unfixed_button, 3, 1)
        buttons.addWidget(self.segment_curate_buttons['reset_view'], 3, 2)
        layout.addLayout(buttons)

    def _layout_training_section(self, layout):
        layout.addWidget(self.make_section_label('Train model'))

        form = QFormLayout()
        form.setVerticalSpacing(8)
        form.addRow(self.make_form_label('labelled sessions', 'folder containing processed sessions for training'), self._path_row(self.source_root_line, self.browse_source_root))
        form.addRow(self.make_form_label('model name', 'name used for the output training folder'), self.run_name_line)
        form.addRow(self.make_form_label('training length', 'full passes through the training set'), self.epochs_spin)
        layout.addLayout(form)

        buttons = QGridLayout()
        buttons.setHorizontalSpacing(6)
        buttons.setVerticalSpacing(6)
        buttons.addWidget(self.build_manifest_button, 0, 0)
        buttons.addWidget(self.train_model_button, 0, 1)
        buttons.addWidget(self.stop_process_button, 1, 0)
        buttons.addWidget(self.training_diagnostics_button, 1, 1)
        buttons.addWidget(self.training_advanced_button, 2, 0, 1, 2)
        layout.addLayout(buttons)

    def _layout_prediction_section(self, layout):
        layout.addWidget(self.make_section_label('Predict'))

        form = QFormLayout()
        form.setVerticalSpacing(6)
        form.addRow(self.make_form_label('chan2 image', 'channel-2 reference image used for ROI prediction'), self._path_row(self.image_line, self.browse_image))
        form.addRow(self.make_form_label('trained model', 'saved trained model file'), self._path_row(self.checkpoint_line, self.browse_checkpoint))
        form.addRow(self.make_form_label('strictness', 'higher values keep only stronger model detections'), self.threshold_spin)
        form.addRow(self.make_form_label('min ROI size', 'smallest ROI size to keep'), self.min_size_spin)
        layout.addLayout(form)

        buttons = QGridLayout()
        buttons.setHorizontalSpacing(6)
        buttons.setVerticalSpacing(6)
        buttons.addWidget(self.load_image_button, 0, 0)
        buttons.addWidget(self.load_model_button, 0, 1)
        buttons.addWidget(self.predict_button, 0, 2)
        layout.addLayout(buttons)

        options = QHBoxLayout()
        options.addWidget(self.show_probability_check)
        options.addStretch(1)
        layout.addLayout(options)

    def _layout_prediction_curation_section(self, layout):
        layout.addWidget(self.make_section_label('Curate ROIs'))

        buttons = self.predict_curate_buttons
        grid = QGridLayout()
        grid.setHorizontalSpacing(6)
        grid.setVerticalSpacing(6)
        grid.addWidget(buttons['select_all'], 0, 0)
        grid.addWidget(buttons['delete'], 0, 1)
        grid.addWidget(buttons['merge'], 0, 2)
        grid.addWidget(buttons['undo'], 1, 0)
        grid.addWidget(buttons['reset_view'], 1, 1)
        grid.addWidget(buttons['save_roi'], 1, 2)
        layout.addLayout(grid)


    def show_segment_advanced_popup(self):
        if self.segment_advanced_popup is None:
            self.segment_advanced_popup = QFrame(self, Qt.Popup)
            self.segment_advanced_popup.setObjectName('advancedPopup')
            layout = QFormLayout(self.segment_advanced_popup)
            layout.setContentsMargins(12, 10, 12, 10)
            layout.setHorizontalSpacing(10)
            layout.setVerticalSpacing(7)
            for spec in PARAMETER_SPECS:
                if spec['name'] in CORE_SEGMENT_PARAMETERS:
                    continue
                tooltip = PARAMETER_TOOLTIPS.get(spec['name'], '')
                widget = self.segment_param_widgets[spec['name']]
                widget.setToolTip(tooltip)
                layout.addRow(self.make_form_label(spec['name'], tooltip), widget)

        pos = self.advanced_segment_button.mapToGlobal(self.advanced_segment_button.rect().bottomLeft())
        self.segment_advanced_popup.move(pos)
        self.segment_advanced_popup.show()
        self.segment_advanced_popup.raise_()

    def show_training_diagnostics_popup(self):
        if self.training_diagnostics_popup is None:
            self.training_diagnostics_popup = QFrame(self, Qt.Popup)
            self.training_diagnostics_popup.setObjectName('advancedPopup')
            layout = QGridLayout(self.training_diagnostics_popup)
            layout.setContentsMargins(12, 10, 12, 10)
            layout.setHorizontalSpacing(6)
            layout.setVerticalSpacing(6)
            layout.addWidget(self.inspect_manifest_button, 0, 0)
            layout.addWidget(self.evaluate_model_button, 0, 1)
            layout.addWidget(self.preview_training_button, 1, 0)
            layout.addWidget(self.preview_predictions_button, 1, 1)

        pos = self.training_diagnostics_button.mapToGlobal(self.training_diagnostics_button.rect().bottomLeft())
        self.training_diagnostics_popup.move(pos)
        self.training_diagnostics_popup.show()
        self.training_diagnostics_popup.raise_()

    def show_training_advanced_popup(self):
        if self.training_advanced_popup is None:
            self.training_advanced_popup = QFrame(self, Qt.Popup)
            self.training_advanced_popup.setObjectName('advancedPopup')
            layout = QFormLayout(self.training_advanced_popup)
            layout.setContentsMargins(12, 10, 12, 10)
            layout.setHorizontalSpacing(10)
            layout.setVerticalSpacing(7)
            layout.addRow(self.make_form_label('dataset table', 'CSV index of labelled images and ROI dicts'), self._path_row(self.manifest_line, self.browse_manifest_out))
            layout.addRow(self.make_form_label('training recipe', 'YAML settings used for model training'), self._path_row(self.config_line, self.browse_config))
            layout.addRow(self.make_form_label('validation split', 'fraction used for tuning during training'), self.val_fraction_spin)
            layout.addRow(self.make_form_label('test split', 'held-out fraction used for scoring'), self.test_fraction_spin)

        pos = self.training_advanced_button.mapToGlobal(self.training_advanced_button.rect().bottomLeft())
        self.training_advanced_popup.move(pos)
        self.training_advanced_popup.show()
        self.training_advanced_popup.raise_()

    def _path_row(self, line_edit, browse_slot):
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        browse = QPushButton('browse')
        self.set_button_role(browse, 'small')
        browse.clicked.connect(browse_slot)
        layout.addWidget(line_edit, 1)
        layout.addWidget(browse)
        return row


    def make_section_label(self, text):
        label = QLabel(text)
        label.setObjectName('sectionLabel')
        label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        return label

    def make_log_box(self):
        box = QPlainTextEdit()
        box.setObjectName('logBox')
        box.setReadOnly(True)
        box.setMaximumBlockCount(3000)
        box.setMinimumHeight(112)
        box.setMaximumHeight(170)
        box.setPlaceholderText('console')
        return box

    def make_form_label(self, text, tooltip):
        label = QLabel(text)
        label.setToolTip(tooltip.rstrip('.'))
        return label

    @staticmethod
    def set_button_role(button, role):
        if button.property('role') == role:
            return
        button.setProperty('role', role)
        button.style().unpolish(button)
        button.style().polish(button)

    @staticmethod
    def prepare_value_control(control):
        control.setButtonSymbols(QAbstractSpinBox.NoButtons)
        control.setKeyboardTracking(False)

    def _make_segment_param_widget(self, spec):
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
        self.prepare_value_control(widget)
        widget.setAlignment(Qt.AlignRight)
        widget.setMinimumWidth(58)
        widget.setMaximumWidth(74)
        return widget

    def _connect_shortcuts(self):
        shortcuts = [
            ('Ctrl+O', self.load_channel_image),
            ('Ctrl+M', self.load_model),
            ('Ctrl+P', self.predict_rois),
            ('Ctrl+S', self.save_roi_file),
            ('Ctrl+R', self.segment_rois),
            ('Ctrl+F', self.fix_selected),
            ('Ctrl+Shift+F', self.clear_fixed),
            ('Delete', self.delete_selected),
            ('Ctrl+Z', self.undo),
            ('Backspace', self.undo),
            ]
        for sequence, slot in shortcuts:
            shortcut = QShortcut(QKeySequence(sequence), self)
            shortcut.activated.connect(slot)

    def _theme(self):
        if self.dark_mode:
            return {
                'window': '#151b17',
                'surface': '#1e2721',
                'surface_alt': '#242f28',
                'surface_strong': '#34433a',
                'surface_hover': '#2c3931',
                'border': '#3c4b42',
                'border_strong': '#66786d',
                'text': '#edf3ef',
                'muted': '#b7c4bb',
                'primary': '#79b28b',
                'primary_hover': '#92c5a1',
                'primary_text': '#102017',
                'danger_bg': '#3b2528',
                'danger_hover': '#503034',
                'danger_border': '#83545b',
                'danger_text': '#f0c6cb',
                'selection': '#365c42',
                'canvas': '#182019',
                }

        return {
            'window': '#f3f6f4',
            'surface': '#ffffff',
            'surface_alt': '#f7faf8',
            'surface_strong': '#dbe7df',
            'surface_hover': '#edf5f0',
            'border': '#cbd8d0',
            'border_strong': '#8aa093',
            'text': '#17211b',
            'muted': '#516157',
            'primary': '#2f6f4e',
            'primary_hover': '#25593f',
            'primary_text': '#ffffff',
            'danger_bg': '#f6e9eb',
            'danger_hover': '#efd8dc',
            'danger_border': '#c88e98',
            'danger_text': '#8f3544',
            'selection': '#c9ead2',
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

        tab_width = 138

        self.setStyleSheet(
            f'''
            QMainWindow, QWidget#centralWidget {{
                background: {theme['window']};
                font-family: 'mononoki';
            }}
            QFrame#canvasFrame {{
                border: 1px solid {theme['border']};
                border-radius: 8px;
                background: {theme['canvas']};
            }}
            QFrame#statusCard, QFrame#advancedPopup {{
                border: 1px solid {theme['border']};
                border-radius: 8px;
                background: {theme['surface']};
            }}

            QSplitter::handle {{
                background: transparent;
                width: 10px;
            }}
            QSplitter::handle:hover {{
                background: {theme['surface_hover']};
                border-radius: 5px;
            }}
            QWidget#controlsPanel, QWidget#mainPane {{
                background: transparent;
            }}
            QTabWidget::pane {{
                border: 1px solid {theme['border']};
                border-radius: 8px;
                background: {theme['surface']};
                top: -2px;
            }}
            QTabBar::tab {{
                background: {theme['surface_alt']};
                border: 1px solid {theme['border']};
                border-bottom: none;
                padding: 6px 6px;
                margin-right: 3px;
                border-top-left-radius: 7px;
                border-top-right-radius: 7px;
                color: {theme['muted']};
                font-weight: 600;
                min-width: 0px;
                width: {tab_width}px;
            }}
            QTabBar::tab:hover {{
                background: {theme['surface_hover']};
                color: {theme['text']};
                border-color: {theme['border_strong']};
            }}
            QTabBar::tab:selected {{
                background: {theme['surface']};
                color: {theme['primary']};
                border-color: {theme['border_strong']};
                border-top: 2px solid {theme['primary']};
            }}
            QLabel {{
                color: {theme['text']};
            }}
            QLabel#stateSummary {{
                color: {theme['text']};
                font-weight: 700;
                font-size: 10pt;
                padding: 0;
            }}
            QLabel#panelValue {{
                color: {theme['muted']};
                font-weight: 600;
                padding: 0;
            }}
            QLabel#sectionLabel {{
                color: {theme['text']};
                font-weight: 700;
                font-size: 10pt;
                padding: 4px 0 1px 0;
            }}

            QLineEdit, QDoubleSpinBox, QSpinBox {{
                border: 1px solid {theme['border']};
                border-radius: 6px;
                padding: 3px 7px;
                background: {theme['surface_alt']};
                color: {theme['text']};
                selection-background-color: {theme['selection']};
                min-height: 21px;
            }}
            QLineEdit:hover, QDoubleSpinBox:hover, QSpinBox:hover {{
                border-color: {theme['border_strong']};
                background: {theme['surface']};
            }}
            QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus {{
                border-color: {theme['primary']};
                background: {theme['surface']};
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
                border-radius: 4px;
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
            QPushButton {{
                border: 1px solid {theme['border']};
                border-radius: 6px;
                padding: 4px 8px;
                background: {theme['surface']};
                color: {theme['text']};
                min-height: 22px;
                font-weight: 600;
            }}
            QPushButton:hover {{
                background: {theme['surface_hover']};
                border-color: {theme['border_strong']};
            }}
            QPushButton:pressed {{
                background: {theme['surface_strong']};
            }}
            QPushButton[role='primary'] {{
                background: {theme['primary']};
                border-color: {theme['primary']};
                color: {theme['primary_text']};
            }}
            QPushButton[role='primary']:hover {{
                background: {theme['primary_hover']};
                border-color: {theme['primary_hover']};
            }}
            QPushButton[role='danger'] {{
                background: {theme['danger_bg']};
                border-color: {theme['danger_border']};
                color: {theme['danger_text']};
            }}
            QPushButton[role='danger']:hover {{
                background: {theme['danger_hover']};
                border-color: {theme['danger_text']};
            }}
            QPushButton[role='quiet'] {{
                background: transparent;
                border-color: transparent;
                color: {theme['muted']};
            }}
            QPushButton[role='quiet']:hover {{
                background: {theme['surface_hover']};
                border-color: {theme['border']};
                color: {theme['text']};
            }}
            QPushButton[role='small'] {{
                min-height: 20px;
                padding: 2px 7px;
                color: {theme['muted']};
            }}
            QPushButton:disabled {{
                background: {theme['surface_alt']};
                border-color: {theme['border']};
                color: {theme['muted']};
            }}
            QPlainTextEdit#logBox {{
                border: 1px solid {theme['border']};
                border-radius: 8px;
                background: {theme['surface']};
                color: {theme['text']};
                font-family: 'mononoki', Consolas, 'Courier New', monospace;
                font-size: 10pt;
                padding: 8px;
                selection-background-color: {theme['selection']};
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
                border-top: 1px solid {theme['border']};
            }}
            '''
            )
        self.refresh_status()

    def set_roi_overlay_visible(self, state):
        self.roi_overlay_check.setText('ROI on' if state else 'ROI off')
        self.plot_image(preserve_view=True)

    def set_dark_mode(self, state):
        self.dark_mode = bool(state)
        self._apply_palette_and_style()
        self.plot_image(preserve_view=True)


    #%% browse
    def browse_source_root(self):
        path = QFileDialog.getExistingDirectory(self, 'select training source root', self.source_root_line.text())
        if path:
            self.source_root_line.setText(path)
            self.refresh_status('source root selected')

    def browse_manifest_out(self):
        path, _ = QFileDialog.getSaveFileName(
            self,
            'select dataset table output',
            self.manifest_line.text(),
            'CSV files (*.csv)',
            )
        if path:
            self.manifest_line.setText(path)
            self.refresh_status('dataset table path selected')

    def browse_config(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            'select training recipe',
            self.config_line.text(),
            'YAML files (*.yaml *.yml)',
            )
        if path:
            self.config_line.setText(path)
            self.refresh_status('training recipe selected')

    def browse_image(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            'select chan2 reference image',
            str(REPO_ROOT),
            'NumPy images (*.npy)',
            )
        if path:
            self.image_line.setText(path)
            self.refresh_status('chan2 path selected')

    def browse_checkpoint(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            'select trained model',
            str(default_output_root() / 'runs'),
            'Trained models (*.pt)',
            )
        if path:
            self.checkpoint_line.setText(path)
            self.refresh_status('trained model selected')

    #%% training processes
    def build_manifest(self):
        source_root = Path(self.source_root_line.text().strip())
        if not source_root.exists():
            self.print_log(f'labelled sessions folder does not exist: {source_root}')
            self.refresh_status('labelled sessions missing')
            return

        args = [
            '--source-root', str(source_root),
            '--out', self.manifest_line.text().strip(),
            '--val-fraction', str(self.val_fraction_spin.value()),
            '--test-fraction', str(self.test_fraction_spin.value()),
            ]
        self.start_process('manifest', SCRIPT_DIR / 'build_manifest.py', args)

    def train_model(self):
        manifest = Path(self.manifest_line.text().strip())
        config = Path(self.config_line.text().strip())
        if not manifest.exists():
            source_root = Path(self.source_root_line.text().strip())
            if not source_root.exists():
                self.print_log(f'labelled sessions folder does not exist: {source_root}')
                self.refresh_status('labelled sessions missing')
                return
            self.print_log('dataset table not found; scanning labelled sessions first')
            self.build_manifest()
            self.print_log('scan started; click train model again after it finishes')
            return
        if not config.exists():
            self.print_log(f'training recipe does not exist: {config}')
            self.refresh_status('training recipe missing')
            return

        try:
            config_path = self.write_training_config()
        except Exception as exc:
            self.print_log(f'failed to write training recipe: {exc}')
            return

        run_name = self.run_name_line.text().strip()
        self.pending_checkpoint_path = default_output_root() / 'runs' / run_name / 'best.pt'
        args = ['--config', str(config_path)]
        self.start_process('train', SCRIPT_DIR / 'train_unet.py', args)

    def evaluate_model(self):
        model_path = self.current_model_path()
        if not model_path:
            self.print_log('please train or load a model first')
            return
        if not Path(model_path).exists():
            self.print_log(f'trained model does not exist: {model_path}')
            self.refresh_status('model file missing')
            return
        manifest = Path(self.manifest_line.text().strip())
        if not manifest.exists():
            self.print_log(f'dataset table does not exist yet: {manifest}')
            self.refresh_status('dataset table missing')
            return

        args = [
            '--manifest', str(manifest),
            '--checkpoint', model_path,
            '--split', 'test',
            '--threshold', str(self.threshold_spin.value()),
            '--min-size', str(self.min_size_spin.value()),
            '--tta',
            ]
        self.start_process('evaluate', SCRIPT_DIR / 'evaluate.py', args)

    def preview_training_labels(self):
        manifest = Path(self.manifest_line.text().strip())
        if not manifest.exists():
            self.print_log(f'dataset table does not exist yet: {manifest}')
            self.refresh_status('dataset table missing')
            return

        args = [
            '--manifest', str(manifest),
            '--split', 'train',
            '--n', '6',
            ]
        self.start_process('preview labels', SCRIPT_DIR / 'plot_training_data.py', args)

    def preview_model_predictions(self):
        manifest = Path(self.manifest_line.text().strip())
        model_path_text = self.current_model_path()
        if not manifest.exists():
            self.print_log(f'dataset table does not exist yet: {manifest}')
            self.refresh_status('dataset table missing')
            return
        if not model_path_text:
            self.print_log('please train or load a model first')
            self.refresh_status('model missing')
            return
        model_path = Path(model_path_text)
        if not model_path.exists():
            self.print_log(f'trained model does not exist: {model_path}')
            self.refresh_status('model file missing')
            return

        args = [
            '--manifest', str(manifest),
            '--split', 'test',
            '--checkpoint', str(model_path),
            '--threshold', str(self.threshold_spin.value()),
            '--min-size', str(self.min_size_spin.value()),
            '--n', '4',
            ]
        self.start_process('preview predictions', SCRIPT_DIR / 'plot_model_diagnostics.py', args)

    def inspect_manifest(self):
        path = Path(self.manifest_line.text().strip())
        if not path.exists():
            self.print_log(f'dataset table not found: {path}')
            self.refresh_status('dataset table missing')
            return

        try:
            import csv

            with open(path, 'r', newline='') as f:
                rows = list(csv.DictReader(f))
        except Exception as exc:
            self.print_log(f'failed to read dataset table: {exc}')
            self.refresh_status('dataset table read failed')
            return

        included = [
            row for row in rows
            if str(row.get('included', '')).lower() in {'true', '1', 'yes', 'y'}
            ]
        split_counts = {}
        for row in included:
            split = row.get('split', '') or 'unsplit'
            split_counts[split] = split_counts.get(split, 0) + 1

        excluded_reasons = {}
        for row in rows:
            if str(row.get('included', '')).lower() in {'true', '1', 'yes', 'y'}:
                continue
            reason = row.get('exclusion_reason', '') or 'not included'
            excluded_reasons[reason] = excluded_reasons.get(reason, 0) + 1

        roi_total = sum(int(float(row.get('roi_count') or 0)) for row in included)
        self.print_log(f'\ndataset table: {path}')
        self.print_log(f'total sessions: {len(rows)}')
        self.print_log(f'included sessions: {len(included)} | total ROIs: {roi_total}')
        self.print_log(f'splits: {split_counts if split_counts else {}}')
        if excluded_reasons:
            self.print_log(f'excluded: {excluded_reasons}')
        self.refresh_status('dataset summary ready')

    def write_training_config(self):
        config = load_config(self.config_line.text())
        run_name = self.run_name_line.text().strip()
        if not run_name:
            raise ValueError('run name cannot be empty')

        # Leave the baseline YAML alone while trying controls here; save this recipe beside the run.
        config.setdefault('data', {})
        config.setdefault('train', {})
        config.setdefault('postprocess', {})

        config['data']['manifest'] = self.path_for_config(self.manifest_line.text())
        config['train']['run_name'] = run_name
        config['train']['epochs'] = int(self.epochs_spin.value())
        config['postprocess']['threshold'] = float(self.threshold_spin.value())
        config['postprocess']['min_size'] = int(self.min_size_spin.value())

        out_path = default_output_root() / 'gui_configs' / f'{run_name}.yaml'
        save_config(config, out_path)
        self.print_log(f'training recipe written to {out_path}')
        return out_path

    def start_process(self, process_name, script_path, args):
        if self.process is not None and self.process.state() != QProcess.NotRunning:
            self.print_log('another process is already running')
            return

        self.current_process_name = process_name
        self.process = QProcess(self)
        self.process.setWorkingDirectory(str(REPO_ROOT))
        self.process.setProgram(sys.executable)
        self.process.setArguments([str(script_path)] + args)
        self.process.readyReadStandardOutput.connect(self.read_process_stdout)
        self.process.readyReadStandardError.connect(self.read_process_stderr)
        self.process.finished.connect(self.process_finished)

        command_text = ' '.join([sys.executable, str(script_path)] + args)
        self.print_log(f'\n$ {command_text}')
        self.process.start()
        # The log keeps the command and output; the status bar only confirms that it started.
        self.refresh_status('process started')

    def stop_process(self):
        if self.process is None or self.process.state() == QProcess.NotRunning:
            self.print_log('no process is running')
            return

        self.process.kill()
        self.print_log('process stopped')

    def read_process_stdout(self):
        text = bytes(self.process.readAllStandardOutput()).decode(errors='replace')
        self.print_log(text, end='')

    def read_process_stderr(self):
        text = bytes(self.process.readAllStandardError()).decode(errors='replace')
        self.print_log(text, end='')

    def process_finished(self, exit_code, exit_status):
        self.print_log(f'\nprocess finished: exit code {exit_code}')
        if (
            self.current_process_name == 'train' and
            self.pending_checkpoint_path is not None and
            self.pending_checkpoint_path.exists()
        ):
            self.checkpoint_line.setText(str(self.pending_checkpoint_path))
            self.last_saved_model_path = self.pending_checkpoint_path
            self.print_log(f'trained model ready: {self.pending_checkpoint_path}')

        self.current_process_name = None
        self.process = None
        self.refresh_status('process finished')

    #%% model prediction
    def load_channel_image(self):
        path = self.image_line.text().strip()
        if not path:
            self.browse_image()
            path = self.image_line.text().strip()
        if not path:
            return

        try:
            image = squeeze_image(np.load(path))
        except Exception as exc:
            self.print_log(f'failed to load image: {exc}')
            return

        self.ref_image = image
        self.image_path = Path(path)
        self.recname = self.get_recname(self.image_path)
        self.roi_dict = {}
        self.labelled = np.zeros_like(self.ref_image, dtype=np.int32)
        self.selected.clear()
        self.fixed_ids.clear()
        self.undo_stack.clear()
        self.probability = None

        default_roi = self.default_roi_path()
        if default_roi.exists():
            self.print_log(f'found existing ROI dict: {default_roi}')
            self.print_log('load an ROI dict to review or continue from it')
        self.plot_image()
        self.canvas.reset_view()
        self.refresh_status('chan2 loaded')

    def load_model(self):
        try:
            self.predictor = self.make_predictor()
            self.predictor.load()
            self.checkpoint_line.setText(str(self.predictor.checkpoint_path))
        except Exception as exc:
            self.print_log(f'failed to load model: {exc}')
            self.predictor = None
            self.refresh_status('model load failed')
            return

        self.model_label.setText(f'model: {self.predictor.checkpoint_path.name}')
        self.last_saved_model_path = self.predictor.checkpoint_path
        self.print_log(f'model loaded from {self.predictor.checkpoint_path}')
        self.refresh_status('model loaded')

    def predict_rois(self):
        if self.ref_image is None:
            self.print_log('please load a chan2 image first')
            return

        if self.predictor is None:
            self.load_model()
        if self.predictor is None:
            return

        self.predictor.threshold = float(self.threshold_spin.value())
        self.predictor.min_size = int(self.min_size_spin.value())
        self.predictor.tta = True

        try:
            prediction = self.predictor.predict_image(self.ref_image)
        except Exception as exc:
            self.print_log(f'prediction failed: {exc}')
            self.refresh_status('prediction failed')
            return

        self.push_undo_state()
        # prediction returns to the same editable state as an MSER or loaded ROI dict
        self.roi_dict = prediction.roi_dict
        self.labelled = prediction.labelled
        self.probability = prediction.probability
        self.selected.clear()
        self.fixed_ids.clear()

        self.plot_image()
        self.canvas.reset_view()
        self.print_log(
            f'predicted {len(self.roi_dict)} ROIs '
            f'(strictness {prediction.threshold:.2f}, min size {prediction.min_size})'
            )
        self.refresh_status('prediction complete')

    def apply_probability_threshold_from_controls(self):
        if self.probability is None:
            return
        self.apply_probability_threshold()

    def apply_probability_threshold(self):
        if self.probability is None:
            return

        # reuse the confidence map here; changing strictness should not rerun
        # the model while I am deciding which faint fibres to keep
        self.roi_dict, self.labelled = probability_to_roi_dict(
            self.probability,
            threshold=float(self.threshold_spin.value()),
            min_size=int(self.min_size_spin.value()),
            )
        self.selected.clear()
        self.fixed_ids.clear()
        self.plot_image(preserve_view=True)
        self.print_log(f'strictness applied: {len(self.roi_dict)} ROIs')
        self.refresh_status('strictness applied')

    def make_predictor(self):
        checkpoint = self.checkpoint_line.text().strip()
        checkpoint_path = Path(checkpoint) if checkpoint else None
        return AxonROIPredictor(
            checkpoint_path=checkpoint_path,
            threshold=float(self.threshold_spin.value()),
            min_size=int(self.min_size_spin.value()),
            tta=True,
            )

    def current_model_path(self):
        if self.predictor is not None:
            return str(self.predictor.checkpoint_path)
        if self.last_saved_model_path is not None:
            return str(self.last_saved_model_path)
        return self.checkpoint_line.text().strip()

    #%% ROI io
    def load_roi_file(self):
        if self.ref_image is None:
            self.print_log('please load a chan2 image first')
            return

        path, _ = QFileDialog.getOpenFileName(
            self,
            'load ROI dict',
            str(self.image_path.parent if self.image_path else REPO_ROOT),
            'NumPy dict (*.npy)',
            )
        if not path:
            return

        try:
            roi_dict = load_roi_dict(path)
            labelled, _, _ = roi_dict_to_label(roi_dict, self.ref_image.shape)
        except Exception as exc:
            self.print_log(f'failed to load ROI dict: {exc}')
            return

        self.push_undo_state()
        self.roi_dict = labels_to_roi_dict(labelled)
        self.labelled = labelled
        self.selected.clear()
        self.fixed_ids.clear()
        self.probability = None
        self.plot_image()
        self.refresh_status('ROI dict loaded')

    def save_roi_file(self):
        if self.ref_image is None:
            self.print_log('please load a chan2 image first')
            return

        self.update_roi_dict()
        out_path = self.default_roi_path()
        save_roi_dict(self.roi_dict, out_path)
        self.print_log(f'saved ROI dict to {out_path}')
        self.refresh_status('ROI dict saved')
        QToolTip.showText(QCursor.pos(), 'ROI dict saved', self, self.rect(), 1600)

    def default_roi_path(self):
        if self.image_path is None:
            return default_output_root() / 'predicted_ROI_dict.npy'
        return self.image_path.parent / f'{self.recname}_ROI_dict.npy'


    #%% MSER segmentation
    def reset_segment_parameters(self):
        for spec in PARAMETER_SPECS:
            self.segment_param_widgets[spec['name']].setValue(spec['default'])
        self.refresh_status('parameters reset')

    def get_segment_params(self):
        params = {}
        for spec in PARAMETER_SPECS:
            value = self.segment_param_widgets[spec['name']].value()
            params[spec['name']] = int(value) if spec['kind'] == 'int' else float(value)
        return params

    def get_fixed_roi_dict(self):
        if not self.roi_dict or not self.fixed_ids:
            return {}
        return {
            roi_id: {
                'xpix': self.roi_dict[roi_id]['xpix'].copy(),
                'ypix': self.roi_dict[roi_id]['ypix'].copy(),
            }
            for roi_id in sorted(self.fixed_ids)
            if roi_id in self.roi_dict
        }

    def segment_rois(self):
        if self.ref_image is None:
            self.print_log('please load a chan2 image first')
            return

        if self.labelled is not None:
            self.push_undo_state()

        # MSER still helps before a model exists and on images that need hand repair
        try:
            roi_dict, labelled, fixed_ids, stats = run_mser_segmentation(
                self.ref_image,
                self.get_segment_params(),
                fixed_rois=self.get_fixed_roi_dict(),
            )
        except Exception as exc:
            self.print_log(f'segmentation failed: {exc}')
            self.refresh_status('segmentation failed')
            return

        self.roi_dict = roi_dict
        self.labelled = labelled
        self.fixed_ids = fixed_ids
        self.selected.clear()
        self.probability = None
        self.plot_image()
        self.canvas.reset_view()
        self.print_log(
            f"MSER regions: {stats['MSER regions']} | kept ROIs: {stats['kept ROIs']} "
            f"(fixed {stats['fixed ROIs']})"
            )
        self.refresh_status('segmentation complete')

    def fix_selected(self):
        if not self.selected:
            return
        self.fixed_ids.update(roi_id for roi_id in self.selected if roi_id in self.roi_dict)
        self.plot_image(preserve_view=True)
        self.refresh_status('selected ROIs fixed')

    def unfix_selected(self):
        if not self.selected:
            return
        self.fixed_ids.difference_update(self.selected)
        self.plot_image(preserve_view=True)
        self.refresh_status('selected ROIs unfixed')

    def clear_fixed(self):
        if not self.fixed_ids:
            return
        self.fixed_ids.clear()
        self.plot_image(preserve_view=True)
        self.refresh_status('fixed marks cleared')

    def clear_unfixed(self):
        if self.labelled is None:
            return
        self.push_undo_state()
        if self.fixed_ids:
            keep_mask = np.isin(self.labelled, list(self.fixed_ids))
            self.labelled[~keep_mask] = 0
        else:
            self.labelled = np.zeros_like(self.labelled, dtype=np.int32)
        self.selected.clear()
        self.compact_labels()
        self.update_roi_dict()
        self.plot_image(preserve_view=True)
        self.refresh_status('unfixed ROIs cleared')

    #%% selection and pruning
    def on_click(self, event):
        if self.ref_image is None:
            if event.inaxes == self.ax:
                self.load_channel_image()
            return
        if event.inaxes != self.ax or self.labelled is None:
            return
        if event.xdata is None or event.ydata is None:
            return

        xpix = int(round(event.xdata))
        ypix = int(round(event.ydata))
        if not self.in_bounds(ypix, xpix):
            return

        roi_id = int(self.labelled[ypix, xpix])
        if roi_id <= 0:
            return

        modifiers = Qt.NoModifier
        if getattr(event, 'guiEvent', None) is not None:
            modifiers = event.guiEvent.modifiers()

        if modifiers & (Qt.ControlModifier | Qt.ShiftModifier):
            if roi_id in self.selected:
                self.selected.remove(roi_id)
            else:
                self.selected.add(roi_id)
        else:
            self.selected = {roi_id}

        self.plot_image(preserve_view=True)
        self.refresh_status('selection updated')

    def select_all(self):
        if self.labelled is None:
            return
        self.selected = {int(label_id) for label_id in np.unique(self.labelled) if label_id > 0}
        self.plot_image(preserve_view=True)
        self.refresh_status('all ROIs selected')

    def invert_selection(self):
        if self.labelled is None:
            return
        ids = {int(label_id) for label_id in np.unique(self.labelled) if label_id > 0}
        self.selected = ids - self.selected
        self.plot_image(preserve_view=True)
        self.refresh_status('selection inverted')

    def clear_selection(self):
        self.selected.clear()
        self.plot_image(preserve_view=True)
        self.refresh_status('selection cleared')

    def delete_selected(self):
        if self.labelled is None or not self.selected:
            return

        self.push_undo_state()
        mask = np.isin(self.labelled, list(self.selected))
        self.labelled[mask] = 0
        self.compact_labels()
        self.selected.clear()
        self.update_roi_dict()
        self.plot_image(preserve_view=True)
        self.refresh_status('selected ROIs deleted')

    def merge_selected(self):
        if self.labelled is None or len(self.selected) < 2:
            return

        self.push_undo_state()
        keep_id = min(self.selected)
        keep_fixed = bool(self.fixed_ids.intersection(self.selected))
        mask = np.isin(self.labelled, list(self.selected))
        self.labelled[mask] = keep_id
        if keep_fixed:
            self.fixed_ids.add(keep_id)
        self.compact_labels()
        self.selected.clear()
        self.update_roi_dict()
        self.plot_image(preserve_view=True)
        self.refresh_status('selected ROIs merged')

    def push_undo_state(self):
        if self.labelled is None:
            return
        self.undo_stack.append((self.labelled.copy(), set(self.selected), set(self.fixed_ids)))
        if len(self.undo_stack) > 30:
            self.undo_stack.pop(0)

    def undo(self):
        if not self.undo_stack:
            self.print_log('nothing to undo')
            return

        state = self.undo_stack.pop()
        if len(state) == 2:
            self.labelled, self.selected = state
            self.fixed_ids.clear()
        else:
            self.labelled, self.selected, self.fixed_ids = state
        self.update_roi_dict()
        self.plot_image(preserve_view=True)
        self.refresh_status('undo')

    #%% plotting
    def plot_image(self, preserve_view=False):
        xlim, ylim = self.get_current_view()
        theme = self._theme()
        self.fig.set_facecolor(theme['canvas'])
        self.ax.clear()
        self.ax.set_facecolor(theme['canvas'])
        self.ax.axis('off')

        if self.ref_image is None:
            self.canvas.setCursor(Qt.PointingHandCursor)
            self.canvas.setToolTip('click to load a chan2 image')
            self.ax.text(
                0.5,
                0.53,
                'open a chan2 image to start',
                transform=self.ax.transAxes,
                ha='center',
                va='center',
                color=theme['muted'],
                fontsize=13,
                fontweight='bold',
                fontfamily=MONONOKI_FONT_FAMILY,
                )
            self.ax.text(
                0.5,
                0.47,
                'click anywhere on this canvas to browse',
                transform=self.ax.transAxes,
                ha='center',
                va='center',
                color=theme['muted'],
                fontsize=9,
                alpha=0.8,
                fontfamily=MONONOKI_FONT_FAMILY,
                )
            self.canvas.draw_idle()
            return

        self.canvas.setCursor(Qt.CrossCursor)
        self.canvas.setToolTip('click an ROI to select; Ctrl/Shift-click to multi-select; right-drag to pan; scroll to zoom')
        base = normalise_for_display(self.ref_image)
        self.ax.imshow(base, cmap='gray', interpolation='nearest')

        if self.show_probability_check.isChecked() and self.probability is not None:
            self.plot_probability()
        if self.roi_overlay_check.isChecked():
            self.plot_roi_overlay()

        if preserve_view and xlim is not None:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)

        self.canvas.draw_idle()
        self.refresh_status()

    def plot_probability(self):
        probability = np.asarray(self.probability, dtype=np.float32)
        rgba = matplotlib.colormaps['magma'](np.clip(probability, 0, 1))
        rgba[..., 3] = np.clip(probability, 0, 1) * 0.35
        self.ax.imshow(rgba, interpolation='nearest')

    def plot_roi_overlay(self):
        if self.labelled is None:
            return

        ids = [int(label_id) for label_id in np.unique(self.labelled) if label_id > 0]
        if not ids:
            return

        overlay = np.zeros((*self.labelled.shape, 4), dtype=np.float32)
        colours = generate_distinct_colours(len(ids))
        for colour, roi_id in zip(colours, ids):
            mask = self.labelled == roi_id
            overlay[mask, :3] = colour
            if roi_id in self.selected:
                overlay[mask, 3] = 0.82
            elif roi_id in self.fixed_ids:
                overlay[mask, 3] = 0.66
            else:
                overlay[mask, 3] = 0.46

        self.ax.imshow(overlay, interpolation='nearest')
        for roi_id in self.selected:
            ypix, xpix = np.where(self.labelled == roi_id)
            if len(xpix) > 0:
                self.ax.plot(xpix, ypix, 'c.', markersize=1)

    def get_current_view(self):
        if not self.ax.images:
            return None, None
        return self.ax.get_xlim(), self.ax.get_ylim()

    def reset_view(self):
        self.canvas.reset_view()
        self.refresh_status('view reset')

    #%% state helpers
    def update_roi_dict(self):
        if self.labelled is None:
            self.roi_dict = {}
            self.fixed_ids.clear()
            return
        self.roi_dict = labels_to_roi_dict(self.labelled)
        self.fixed_ids = {roi_id for roi_id in self.fixed_ids if roi_id in self.roi_dict}

    def compact_labels(self):
        if self.labelled is None:
            return

        out = np.zeros_like(self.labelled, dtype=np.int32)
        remap = {}
        next_id = 1
        for label_id in sorted(np.unique(self.labelled)):
            if label_id == 0:
                continue
            out[self.labelled == label_id] = next_id
            remap[int(label_id)] = next_id
            next_id += 1
        self.labelled = out
        self.fixed_ids = {remap[roi_id] for roi_id in self.fixed_ids if roi_id in remap}

    def in_bounds(self, ypix, xpix):
        return (
            0 <= ypix < self.labelled.shape[0] and
            0 <= xpix < self.labelled.shape[1]
            )

    def refresh_status(self, message=None):
        image_name = self.image_path.name if self.image_path else 'no image'
        roi_count = len(self.roi_dict) if self.roi_dict is not None else 0
        selected_count = len(self.selected)
        fixed_count = len(self.fixed_ids)
        self.state_label.setText(f'image: {image_name}')
        self.roi_label.setText(f'ROIs: {roi_count} | selected: {selected_count}')
        self.fixed_label.setText(f'fixed: {fixed_count}')
        self.update_workflow_state()
        if message:
            # Keep routine confirmation off the canvas while I am selecting and fixing ROIs.
            self.statusBar().showMessage(message, 4000)

    def update_workflow_state(self):
        image_ready = self.ref_image is not None
        checkpoint_text = self.checkpoint_line.text().strip()
        checkpoint_ready = bool(checkpoint_text) and Path(checkpoint_text).exists()
        model_ready = bool(self.current_model_path()) and Path(self.current_model_path()).exists()
        model_loaded = self.predictor is not None
        manifest_text = self.manifest_line.text().strip()
        manifest_ready = bool(manifest_text) and Path(manifest_text).exists()
        config_text = self.config_line.text().strip()
        config_ready = bool(config_text) and Path(config_text).exists()
        source_text = self.source_root_line.text().strip()
        source_ready = bool(source_text) and Path(source_text).exists()
        has_rois = self.labelled is not None and bool(np.any(self.labelled > 0))
        process_running = self.process is not None and self.process.state() != QProcess.NotRunning

        self.build_manifest_button.setEnabled(source_ready and not process_running)
        self.inspect_manifest_button.setEnabled(manifest_ready and not process_running)
        self.train_model_button.setEnabled(source_ready and config_ready and bool(self.run_name_line.text().strip()) and not process_running)
        self.evaluate_model_button.setEnabled(manifest_ready and model_ready and not process_running)
        self.preview_training_button.setEnabled(manifest_ready and not process_running)
        self.preview_predictions_button.setEnabled(manifest_ready and model_ready and not process_running)
        self.stop_process_button.setEnabled(process_running)

        self.load_model_button.setEnabled(checkpoint_ready and not process_running)
        self.predict_button.setEnabled(image_ready and model_loaded and not process_running)
        self.set_button_role(self.predict_button, 'primary' if model_loaded else '')

        self.segment_button.setEnabled(image_ready)
        self.segment_load_roi_button.setEnabled(image_ready)
        self.segment_save_roi_button.setEnabled(image_ready and has_rois)
        self.fix_selected_button.setEnabled(bool(self.selected))
        self.unfix_selected_button.setEnabled(bool(self.selected))
        self.clear_fixed_button.setEnabled(bool(self.fixed_ids))
        self.clear_unfixed_button.setEnabled(has_rois)

        for button in self.curate_buttons['select_all']:
            button.setEnabled(has_rois)

        for button in self.curate_buttons['delete']:
            button.setEnabled(bool(self.selected))
        for button in self.curate_buttons['merge']:
            button.setEnabled(len(self.selected) >= 2)
        for button in self.curate_buttons['undo']:
            button.setEnabled(bool(self.undo_stack))

        for button in self.curate_buttons['save_roi']:
            button.setEnabled(image_ready and has_rois)

    def print_log(self, text, end='\n'):
        if end:
            text = f'{text}{end}'
        box = getattr(self, 'output_box', None)
        if box is None:
            return
        box.moveCursor(QTextCursor.End)
        box.insertPlainText(text)
        box.moveCursor(QTextCursor.End)

    @staticmethod
    def get_recname(path):
        name = Path(path).name
        if '_ref_mat' in name:
            return name.split('_ref_mat')[0]
        return Path(path).stem

    @staticmethod
    def path_for_config(path):
        path = Path(path)
        try:
            return str(path.resolve().relative_to(Path(REPO_ROOT).resolve()))
        except ValueError:
            return str(path)


#%% entry point
def main():
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    app.setStyle('Fusion')
    app.setFont(load_gui_font(9))
    app.setWindowIcon(QIcon(str(APP_ICON_PATH)))
    window = FibreSightWorkbench()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
