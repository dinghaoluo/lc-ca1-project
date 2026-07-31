# -*- coding: utf-8 -*-
'''
Created on Fri Nov  8 12:03:35 2024
Modified to support targeted conversion, 12 Nov 2024, Dinghao Luo
Modified on 31 March 2026: class-based GUI rewrite

convert greyscale TIFF stacks to video through a GUI

@author: Dinghao Luo
'''

#%% imports
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import queue
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import cv2
import numpy as np
from PIL import Image, ImageEnhance, ImageTk
import tifffile
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import gaussian_kernel_unity
import project_paths as pp


#%% constants
DEFAULT_OUTPUT_DIR = pp.MOVIE_FIGURES_STEM
CANVAS_SIZE        = 520

try:
    RESAMPLE_BILINEAR = Image.Resampling.BILINEAR
except AttributeError:
    RESAMPLE_BILINEAR = Image.BILINEAR

CODECS_BY_FORMAT = {
    '.avi': ['HFYU', 'MJPG', 'XVID'],
    '.mp4': ['mp4v'],
    '.mov': ['mp4v', 'MJPG'],
    '.mkv': ['FFV1', 'MJPG', 'XVID'],
}

CONTRAST_MODES = {
    'full stack range': None,
    '0.5-99.5 percentile': (0.5, 99.5),
    '1-99 percentile': (1.0, 99.0),
    '2-98 percentile': (2.0, 98.0),
}


#%% data structures
@dataclass
class DisplayTransform:
    scale: float = 1.0
    offset_x: int = 0
    offset_y: int = 0
    image_width: int = 0
    image_height: int = 0

@dataclass
class MovieOptions:
    in_path: Path
    out_dir: Path
    out_stem: str
    extension: str
    codec: str
    fps: float
    start_frame: int
    end_frame: int
    crop: tuple[int, int, int, int] | None
    smooth: bool
    sigma: float
    contrast_mode: str
    max_output_width: int


#%% movie helpers
def read_tif_stack(path):
    stack = tifffile.imread(path)
    stack = np.asarray(stack)

    if stack.ndim == 2:
        stack = stack[np.newaxis, :, :]
    if stack.ndim != 3:
        raise ValueError(
            f'Expected a 2D frame or 3D greyscale stack, got shape {stack.shape}.'
        )
    return stack

def normalise_frame(frame, vmin, vmax):
    if vmax <= vmin:
        return np.zeros(frame.shape, dtype=np.uint8)
    frame = np.clip(frame, vmin, vmax)
    return ((frame - vmin) / (vmax - vmin) * 255).astype(np.uint8)

def resize_movie_frame(frame, max_output_width):
    if max_output_width <= 0 or frame.shape[1] <= max_output_width:
        return frame

    scale = max_output_width / frame.shape[1]
    new_height = max(1, int(round(frame.shape[0] * scale)))
    return cv2.resize(
        frame,
        (max_output_width, new_height),
        interpolation=cv2.INTER_AREA,
    )

def write_movie(options, preloaded_data=None, stop_event=None, log=None, progress=None):
    log = log or print

    movie = preloaded_data if preloaded_data is not None else read_tif_stack(options.in_path)
    frame_count, frame_height, frame_width = movie.shape
    start_frame = max(0, min(options.start_frame, frame_count - 1))
    end_frame = max(start_frame, min(options.end_frame, frame_count - 1))

    movie = movie[start_frame:end_frame + 1, :, :]
    if options.crop is not None:
        x, y, width, height = options.crop
        _, frame_height, frame_width = movie.shape
        if x < 0 or y < 0 or x + width > frame_width or y + height > frame_height:
            raise ValueError(
                'Crop exceeds frame bounds: '
                f'x={x}, y={y}, width={width}, height={height}, '
                f'frame={frame_width}x{frame_height}.'
            )
        movie = movie[:, y:y + height, x:x + width]

    if options.smooth:
        log(f'applying temporal smoothing, sigma={options.sigma:g}')
        kernel = gaussian_kernel_unity(options.sigma)
        pad_width = len(kernel) // 2
        movie_padded = np.pad(
            movie,
            ((pad_width, pad_width), (0, 0), (0, 0)),
            mode='reflect',
        )
        movie = np.apply_along_axis(
            lambda x: np.convolve(x, kernel, mode='same'),
            axis=0,
            arr=movie_padded,
        )[pad_width:-pad_width, :, :]

    percentiles = CONTRAST_MODES[options.contrast_mode]
    if percentiles is None:
        vmin = float(movie.min())
        vmax = float(movie.max())
    else:
        vmin, vmax = np.percentile(movie, percentiles)
        vmin = float(vmin)
        vmax = float(vmax)

    out_dir = options.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_stem = options.out_stem.strip()
    out_stem = ''.join(
        char if char not in r'<>:"/\|?*' else '_'
        for char in out_stem
        )
    out_path = out_dir / f'{out_stem}{options.extension}'

    first_frame = resize_movie_frame(
        normalise_frame(movie[0, :, :], vmin, vmax),
        options.max_output_width,
    )
    out_height, out_width = first_frame.shape
    write_color = options.extension.lower() in {'.mp4', '.mov'}

    fourcc = cv2.VideoWriter_fourcc(*options.codec)
    writer = cv2.VideoWriter(
        str(out_path),
        fourcc,
        options.fps,
        (out_width, out_height),
        isColor=write_color,
    )

    if not writer.isOpened():
        raise RuntimeError(
            f'Could not open video writer for {out_path}. '
            f'Try another codec or output format.'
        )

    log(f'writing {out_path}')
    log(f'frames {start_frame}-{end_frame}; size {out_width}x{out_height}; fps {options.fps:g}')

    try:
        total = movie.shape[0]
        for frame_index in tqdm(range(total), file=sys.stdout, ncols=50):
            if stop_event is not None and stop_event.is_set():
                log('conversion halted before completion')
                break

            frame = normalise_frame(movie[frame_index, :, :], vmin, vmax)
            frame = resize_movie_frame(frame, options.max_output_width)
            if write_color:
                frame = cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR)
            writer.write(frame)

            if progress is not None:
                progress(frame_index + 1, total)
    finally:
        writer.release()

    if stop_event is not None and stop_event.is_set():
        return None

    log('video saved successfully')
    return out_path


#%% gui
class TifMovieConverterApp:
    def __init__(self, root):
        self.root = root
        self.root.title('TIFF movie converter')
        self.root.minsize(1120, 760)

        self.stack = None
        self.tif_min = 0.0
        self.tif_max = 1.0
        self.display_transform = DisplayTransform()
        self.preview_photo = None
        self.crop_drag_start = None
        self.crop_canvas_rect = None
        self.worker = None
        self.stop_event = threading.Event()
        self.messages = queue.Queue()

        self._build_vars()
        self._build_styles()
        self._build_layout()
        self._poll_messages()

    def _build_vars(self):
        self.in_path_var = tk.StringVar()
        self.out_dir_var = tk.StringVar(value=str(DEFAULT_OUTPUT_DIR))
        self.out_stem_var = tk.StringVar()
        self.fps_var = tk.StringVar(value='30')
        self.extension_var = tk.StringVar(value='.mp4')
        self.codec_var = tk.StringVar(value=CODECS_BY_FORMAT['.mp4'][0])
        self.contrast_var = tk.StringVar(value='full stack range')
        self.smooth_var = tk.BooleanVar(value=True)
        self.sigma_var = tk.StringVar(value='3')
        self.max_width_var = tk.StringVar(value='0')
        self.frame_var = tk.DoubleVar(value=0)
        self.start_frame_var = tk.StringVar(value='0')
        self.end_frame_var = tk.StringVar(value='0')
        self.crop_vars = [tk.StringVar(value='0') for _ in range(4)]
        self.status_var = tk.StringVar(value='Select a TIFF stack to begin.')
        self.frame_label_var = tk.StringVar(value='frame 0')
        self.summary_var = tk.StringVar(value='No TIFF loaded.')

    def _build_styles(self):
        style = ttk.Style(self.root)
        for theme in ('vista', 'xpnative', 'clam'):
            if theme in style.theme_names():
                style.theme_use(theme)
                break
        style.configure('Header.TLabel', font=('Segoe UI', 13, 'bold'))
        style.configure('Section.TLabelframe.Label', font=('Segoe UI', 10, 'bold'))
        style.configure('Accent.TButton', font=('Segoe UI', 10, 'bold'))

    def _build_layout(self):
        self.root.columnconfigure(0, weight=0)
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)

        control_panel = ttk.Frame(self.root, padding=(16, 14, 12, 14))
        control_panel.grid(row=0, column=0, sticky='nsw')
        control_panel.columnconfigure(0, weight=1)

        preview_panel = ttk.Frame(self.root, padding=(8, 14, 16, 14))
        preview_panel.grid(row=0, column=1, sticky='nsew')
        preview_panel.columnconfigure(0, weight=1)
        preview_panel.rowconfigure(1, weight=1)

        title = ttk.Label(control_panel, text='TIFF movie converter', style='Header.TLabel')
        title.grid(row=0, column=0, sticky='w')

        input_frame = ttk.Labelframe(
            control_panel,
            text='Input and output',
            style='Section.TLabelframe',
            padding=10,
        )
        input_frame.grid(row=1, column=0, sticky='ew', pady=(12, 8))
        input_frame.columnconfigure(1, weight=1)

        ttk.Label(input_frame, text='TIFF stack').grid(row=0, column=0, sticky='w', pady=3)
        ttk.Entry(input_frame, textvariable=self.in_path_var, width=54).grid(
            row=0,
            column=1,
            sticky='ew',
            padx=(8, 6),
            pady=3,
        )
        ttk.Button(input_frame, text='Browse', command=self.browse_input).grid(
            row=0,
            column=2,
            sticky='e',
            pady=3,
        )

        ttk.Label(input_frame, text='Output folder').grid(row=1, column=0, sticky='w', pady=3)
        ttk.Entry(input_frame, textvariable=self.out_dir_var, width=54).grid(
            row=1,
            column=1,
            sticky='ew',
            padx=(8, 6),
            pady=3,
        )
        ttk.Button(input_frame, text='Browse', command=self.browse_output).grid(
            row=1,
            column=2,
            sticky='e',
            pady=3,
        )

        ttk.Label(input_frame, text='Movie name').grid(row=2, column=0, sticky='w', pady=3)
        ttk.Entry(input_frame, textvariable=self.out_stem_var).grid(
            row=2,
            column=1,
            columnspan=2,
            sticky='ew',
            padx=(8, 0),
            pady=3,
        )

        export_frame = ttk.Labelframe(
            control_panel,
            text='Export',
            style='Section.TLabelframe',
            padding=10,
        )
        export_frame.grid(row=2, column=0, sticky='ew', pady=8)
        export_frame.columnconfigure(1, weight=1)
        export_frame.columnconfigure(3, weight=1)

        ttk.Label(export_frame, text='Format').grid(row=0, column=0, sticky='w', pady=3)
        self.format_box = ttk.Combobox(
            export_frame,
            textvariable=self.extension_var,
            values=list(CODECS_BY_FORMAT),
            state='readonly',
            width=10,
        )
        self.format_box.grid(row=0, column=1, sticky='w', padx=(8, 14), pady=3)
        self.format_box.bind('<<ComboboxSelected>>', self.on_format_changed)

        ttk.Label(export_frame, text='Codec').grid(row=0, column=2, sticky='w', pady=3)
        self.codec_box = ttk.Combobox(
            export_frame,
            textvariable=self.codec_var,
            values=CODECS_BY_FORMAT[self.extension_var.get()],
            state='readonly',
            width=10,
        )
        self.codec_box.grid(row=0, column=3, sticky='w', padx=(8, 0), pady=3)

        ttk.Label(export_frame, text='Frame rate').grid(row=1, column=0, sticky='w', pady=3)
        ttk.Entry(export_frame, textvariable=self.fps_var, width=10).grid(
            row=1,
            column=1,
            sticky='w',
            padx=(8, 14),
            pady=3,
        )

        ttk.Label(export_frame, text='Max width').grid(row=1, column=2, sticky='w', pady=3)
        ttk.Entry(export_frame, textvariable=self.max_width_var, width=10).grid(
            row=1,
            column=3,
            sticky='w',
            padx=(8, 0),
            pady=3,
        )

        ttk.Label(export_frame, text='Contrast').grid(row=2, column=0, sticky='w', pady=3)
        ttk.Combobox(
            export_frame,
            textvariable=self.contrast_var,
            values=list(CONTRAST_MODES),
            state='readonly',
            width=22,
        ).grid(row=2, column=1, columnspan=3, sticky='w', padx=(8, 0), pady=3)

        ttk.Checkbutton(
            export_frame,
            text='Temporal smoothing',
            variable=self.smooth_var,
        ).grid(row=3, column=0, columnspan=2, sticky='w', pady=(6, 3))

        ttk.Label(export_frame, text='Sigma').grid(row=3, column=2, sticky='w', pady=(6, 3))
        ttk.Entry(export_frame, textvariable=self.sigma_var, width=10).grid(
            row=3,
            column=3,
            sticky='w',
            padx=(8, 0),
            pady=(6, 3),
        )

        frame_range = ttk.Labelframe(
            control_panel,
            text='Frame range',
            style='Section.TLabelframe',
            padding=10,
        )
        frame_range.grid(row=3, column=0, sticky='ew', pady=8)
        frame_range.columnconfigure(1, weight=1)

        ttk.Label(frame_range, text='Start').grid(row=0, column=0, sticky='w', pady=3)
        ttk.Entry(frame_range, textvariable=self.start_frame_var, width=10).grid(
            row=0,
            column=1,
            sticky='w',
            padx=(8, 10),
            pady=3,
        )
        ttk.Button(frame_range, text='Set current', command=self.set_start_frame).grid(
            row=0,
            column=2,
            sticky='e',
            pady=3,
        )

        ttk.Label(frame_range, text='End').grid(row=1, column=0, sticky='w', pady=3)
        ttk.Entry(frame_range, textvariable=self.end_frame_var, width=10).grid(
            row=1,
            column=1,
            sticky='w',
            padx=(8, 10),
            pady=3,
        )
        ttk.Button(frame_range, text='Set current', command=self.set_end_frame).grid(
            row=1,
            column=2,
            sticky='e',
            pady=3,
        )

        crop_frame = ttk.Labelframe(
            control_panel,
            text='Crop',
            style='Section.TLabelframe',
            padding=10,
        )
        crop_frame.grid(row=4, column=0, sticky='ew', pady=8)

        for column, label in enumerate(('x', 'y', 'width', 'height')):
            ttk.Label(crop_frame, text=label).grid(row=0, column=column, sticky='w')
            entry = ttk.Entry(crop_frame, textvariable=self.crop_vars[column], width=8)
            entry.grid(row=1, column=column, sticky='w', padx=(0, 8), pady=(2, 6))
            entry.bind('<FocusOut>', lambda _event: self.redraw_crop_from_entries())
            entry.bind('<Return>', lambda _event: self.redraw_crop_from_entries())

        ttk.Button(crop_frame, text='Clear crop', command=self.clear_crop).grid(
            row=1,
            column=4,
            sticky='w',
            padx=(8, 0),
            pady=(2, 6),
        )

        run_frame = ttk.Frame(control_panel)
        run_frame.grid(row=5, column=0, sticky='ew', pady=(12, 6))
        run_frame.columnconfigure(0, weight=1)
        run_frame.columnconfigure(1, weight=1)

        self.run_button = ttk.Button(
            run_frame,
            text='Export movie',
            style='Accent.TButton',
            command=self.start_export,
        )
        self.run_button.grid(row=0, column=0, sticky='ew', padx=(0, 6))

        self.cancel_button = ttk.Button(
            run_frame,
            text='Cancel',
            command=self.cancel_export,
            state='disabled',
        )
        self.cancel_button.grid(row=0, column=1, sticky='ew', padx=(6, 0))

        self.progress = ttk.Progressbar(control_panel, mode='determinate')
        self.progress.grid(row=6, column=0, sticky='ew', pady=(4, 6))
        ttk.Label(control_panel, textvariable=self.status_var).grid(row=7, column=0, sticky='w')

        log_frame = ttk.Labelframe(
            control_panel,
            text='Log',
            style='Section.TLabelframe',
            padding=(8, 8, 8, 8),
        )
        log_frame.grid(row=8, column=0, sticky='nsew', pady=(8, 0))
        control_panel.rowconfigure(8, weight=1)
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        self.log_text = tk.Text(log_frame, height=4, width=60, wrap=tk.WORD)
        self.log_text.grid(row=0, column=0, sticky='nsew')
        log_scroll = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        log_scroll.grid(row=0, column=1, sticky='ns')
        self.log_text.configure(yscrollcommand=log_scroll.set)

        ttk.Label(preview_panel, textvariable=self.summary_var, style='Header.TLabel').grid(
            row=0,
            column=0,
            sticky='w',
            pady=(0, 8),
        )

        self.canvas = tk.Canvas(
            preview_panel,
            width=CANVAS_SIZE,
            height=CANVAS_SIZE,
            bg='#20242a',
            highlightthickness=1,
            highlightbackground='#626975',
        )
        self.canvas.grid(row=1, column=0, sticky='nsew')
        self.canvas.bind('<Button-1>', self.start_crop)
        self.canvas.bind('<B1-Motion>', self.update_crop)
        self.canvas.bind('<ButtonRelease-1>', self.finish_crop)

        preview_controls = ttk.Frame(preview_panel)
        preview_controls.grid(row=2, column=0, sticky='ew', pady=(10, 0))
        preview_controls.columnconfigure(1, weight=1)

        ttk.Label(preview_controls, textvariable=self.frame_label_var).grid(
            row=0,
            column=0,
            sticky='w',
            padx=(0, 8),
        )
        self.frame_slider = ttk.Scale(
            preview_controls,
            from_=0,
            to=0,
            variable=self.frame_var,
            orient='horizontal',
            command=self.on_frame_slider,
        )
        self.frame_slider.grid(row=0, column=1, sticky='ew')

        ttk.Button(preview_controls, text='Reload', command=self.reload_input).grid(
            row=0,
            column=2,
            sticky='e',
            padx=(8, 0),
        )

    def browse_input(self):
        path = filedialog.askopenfilename(
            filetypes=[
                ('TIFF files', '*.tif *.tiff'),
                ('All files', '*.*'),
            ]
        )
        if not path:
            return
        self.in_path_var.set(path)
        if not self.out_stem_var.get().strip():
            self.out_stem_var.set(Path(path).stem)
        self.load_tif_path(path)

    def browse_output(self):
        path = filedialog.askdirectory(initialdir=self.out_dir_var.get() or str(DEFAULT_OUTPUT_DIR))
        if path:
            self.out_dir_var.set(path)

    def reload_input(self):
        path = self.in_path_var.get().strip()
        if path:
            self.load_tif_path(path)

    def load_tif_path(self, path):
        try:
            stack = read_tif_stack(path)
        except Exception as exc:
            messagebox.showerror('Could not load TIFF', str(exc))
            self.log(f'load failed: {exc}')
            return

        self.stack = stack
        self.tif_min = float(stack.min())
        self.tif_max = float(stack.max())
        frame_count, frame_height, frame_width = stack.shape

        self.frame_slider.configure(from_=0, to=max(0, frame_count - 1))
        self.frame_var.set(0)
        self.start_frame_var.set('0')
        self.end_frame_var.set(str(frame_count - 1))
        self.clear_crop(log_message=False)
        self.summary_var.set(f'{Path(path).name}: {frame_count} frames, {frame_width}x{frame_height}')
        self.status_var.set('TIFF loaded.')
        self.log(f'loaded {path}')
        self.show_frame(0)

    def on_format_changed(self, _event=None):
        codecs = CODECS_BY_FORMAT[self.extension_var.get()]
        self.codec_box.configure(values=codecs)
        self.codec_var.set(codecs[0])

    def current_frame_index(self):
        if self.stack is None:
            return 0
        return int(round(float(self.frame_var.get())))

    def on_frame_slider(self, _value=None):
        if self.stack is None:
            return
        self.show_frame(self.current_frame_index())

    def set_start_frame(self):
        self.start_frame_var.set(str(self.current_frame_index()))

    def set_end_frame(self):
        self.end_frame_var.set(str(self.current_frame_index()))

    def clear_crop(self, log_message=True):
        for var in self.crop_vars:
            var.set('0')
        self.crop_canvas_rect = None
        self.canvas.delete('crop_rect')
        if log_message:
            self.log('crop cleared')

    def parse_int(self, value, field_name, min_value=None):
        try:
            parsed = int(float(value))
        except ValueError as exc:
            raise ValueError(f'{field_name} must be numeric.') from exc
        if min_value is not None and parsed < min_value:
            raise ValueError(f'{field_name} must be at least {min_value}.')
        return parsed

    def parse_float(self, value, field_name, min_value=None):
        try:
            parsed = float(value)
        except ValueError as exc:
            raise ValueError(f'{field_name} must be numeric.') from exc
        if min_value is not None and parsed < min_value:
            raise ValueError(f'{field_name} must be at least {min_value}.')
        return parsed

    def parse_crop(self):
        values = [
            self.parse_int(var.get(), label, min_value=0)
            for var, label in zip(self.crop_vars, ('crop x', 'crop y', 'crop width', 'crop height'))
        ]
        if not any(values) or values[2] == 0 or values[3] == 0:
            return None
        return tuple(values)

    def parse_options(self):
        if self.stack is None:
            raise ValueError('Load a TIFF stack before exporting.')

        in_path = Path(self.in_path_var.get().strip())
        if not in_path.exists():
            raise ValueError('Input TIFF path does not exist.')

        out_stem = (self.out_stem_var.get() or in_path.stem).strip()
        out_stem = ''.join(
            char if char not in r'<>:"/\|?*' else '_'
            for char in out_stem
            )
        if not out_stem:
            raise ValueError('Movie name is empty.')

        frame_count = self.stack.shape[0]
        start_frame = self.parse_int(self.start_frame_var.get(), 'start frame', min_value=0)
        end_frame = self.parse_int(self.end_frame_var.get(), 'end frame', min_value=0)
        start_frame = max(0, min(start_frame, frame_count - 1))
        end_frame = max(start_frame, min(end_frame, frame_count - 1))

        return MovieOptions(
            in_path=in_path,
            out_dir=Path(self.out_dir_var.get().strip() or DEFAULT_OUTPUT_DIR),
            out_stem=out_stem,
            extension=self.extension_var.get(),
            codec=self.codec_var.get(),
            fps=self.parse_float(self.fps_var.get(), 'frame rate', min_value=0.1),
            start_frame=start_frame,
            end_frame=end_frame,
            crop=self.parse_crop(),
            smooth=self.smooth_var.get(),
            sigma=self.parse_float(self.sigma_var.get(), 'sigma', min_value=0.1),
            contrast_mode=self.contrast_var.get(),
            max_output_width=self.parse_int(self.max_width_var.get(), 'max width', min_value=0),
        )

    def start_export(self):
        if self.worker is not None and self.worker.is_alive():
            return

        try:
            options = self.parse_options()
        except ValueError as exc:
            messagebox.showerror('Cannot export movie', str(exc))
            return

        self.stop_event.clear()
        self.progress.configure(value=0, maximum=1)
        self.status_var.set('Exporting movie...')
        self.run_button.configure(state='disabled')
        self.cancel_button.configure(state='normal')
        self.log('starting export')

        # keep conversion off the Tk main thread so the GUI does not freeze
        self.worker = threading.Thread(
            target=self._export_worker,
            args=(options,),
            daemon=True,
        )
        self.worker.start()

    def _export_worker(self, options):
        try:
            out_path = write_movie(
                options,
                preloaded_data=self.stack,
                stop_event=self.stop_event,
                log=lambda msg: self.messages.put(('log', msg)),
                progress=lambda done, total: self.messages.put(('progress', done, total)),
            )
            if out_path is None:
                self.messages.put(('status', 'Export halted.'))
            else:
                self.messages.put(('status', f'saved: {out_path}'))
        except Exception as exc:
            self.messages.put(('error', str(exc)))
        finally:
            self.messages.put(('done', None))

    def cancel_export(self):
        self.stop_event.set()
        self.status_var.set('Cancelling after current frame...')
        self.log('cancel requested')

    def _poll_messages(self):
        while True:
            try:
                message = self.messages.get_nowait()
            except queue.Empty:
                break

            kind = message[0]
            if kind == 'log':
                self.log(message[1])
            elif kind == 'progress':
                _, done, total = message
                self.progress.configure(maximum=total, value=done)
                self.status_var.set(f'Exporting frame {done}/{total}')
            elif kind == 'status':
                self.status_var.set(message[1])
                self.log(message[1])
            elif kind == 'error':
                self.status_var.set('Export failed.')
                self.log(f'error: {message[1]}')
                messagebox.showerror('Export failed', message[1])
            elif kind == 'done':
                self.run_button.configure(state='normal')
                self.cancel_button.configure(state='disabled')

        self.root.after(80, self._poll_messages)

    def log(self, message):
        self.log_text.insert(tk.END, f'{message}\n')
        self.log_text.yview(tk.END)

    def show_frame(self, frame_index):
        if self.stack is None:
            return

        frame_index = max(0, min(frame_index, self.stack.shape[0] - 1))
        self.frame_label_var.set(f'frame {frame_index}')
        frame = normalise_frame(self.stack[frame_index, :, :], self.tif_min, self.tif_max)
        image = Image.fromarray(frame)

        image_width, image_height = image.size
        scale = min(CANVAS_SIZE / image_width, CANVAS_SIZE / image_height)
        display_width = max(1, int(round(image_width * scale)))
        display_height = max(1, int(round(image_height * scale)))
        offset_x = (CANVAS_SIZE - display_width) // 2
        offset_y = (CANVAS_SIZE - display_height) // 2

        if scale != 1.0:
            image = image.resize((display_width, display_height), RESAMPLE_BILINEAR)

        image = ImageEnhance.Contrast(image).enhance(1.25)
        self.preview_photo = ImageTk.PhotoImage(image)
        self.canvas.delete('all')
        self.canvas.create_image(offset_x, offset_y, anchor=tk.NW, image=self.preview_photo)

        self.display_transform = DisplayTransform(
            scale=scale,
            offset_x=offset_x,
            offset_y=offset_y,
            image_width=image_width,
            image_height=image_height,
        )
        self.redraw_crop_from_entries()

    def canvas_to_image(self, x, y):
        transform = self.display_transform
        image_x = int(round((x - transform.offset_x) / transform.scale))
        image_y = int(round((y - transform.offset_y) / transform.scale))
        image_x = max(0, min(transform.image_width, image_x))
        image_y = max(0, min(transform.image_height, image_y))
        return image_x, image_y

    def image_to_canvas(self, x, y):
        transform = self.display_transform
        canvas_x = int(round(x * transform.scale + transform.offset_x))
        canvas_y = int(round(y * transform.scale + transform.offset_y))
        return canvas_x, canvas_y

    def start_crop(self, event):
        if self.stack is None:
            return
        self.crop_drag_start = self.canvas_to_image(event.x, event.y)
        self.canvas.delete('crop_rect')

    def update_crop(self, event):
        if self.crop_drag_start is None:
            return

        x1, y1 = self.crop_drag_start
        x2, y2 = self.canvas_to_image(event.x, event.y)
        canvas_x1, canvas_y1 = self.image_to_canvas(x1, y1)
        canvas_x2, canvas_y2 = self.image_to_canvas(x2, y2)

        self.canvas.delete('crop_rect')
        self.canvas.create_rectangle(
            canvas_x1,
            canvas_y1,
            canvas_x2,
            canvas_y2,
            outline='#ff4d4d',
            width=2,
            tags='crop_rect',
        )

    def finish_crop(self, event):
        if self.crop_drag_start is None:
            return

        x1, y1 = self.crop_drag_start
        x2, y2 = self.canvas_to_image(event.x, event.y)
        self.crop_drag_start = None

        x = min(x1, x2)
        y = min(y1, y2)
        width = abs(x2 - x1)
        height = abs(y2 - y1)

        for var, value in zip(self.crop_vars, (x, y, width, height)):
            var.set(str(value))
        self.redraw_crop_from_entries()

    def redraw_crop_from_entries(self):
        self.canvas.delete('crop_rect')
        if self.stack is None:
            return

        try:
            crop = self.parse_crop()
        except ValueError:
            return
        if crop is None:
            return

        x, y, width, height = crop
        x1, y1 = self.image_to_canvas(x, y)
        x2, y2 = self.image_to_canvas(x + width, y + height)
        self.canvas.create_rectangle(
            x1,
            y1,
            x2,
            y2,
            outline='#ff4d4d',
            width=2,
            tags='crop_rect',
        )


#%% main
def main():
    root = tk.Tk()
    TifMovieConverterApp(root)
    root.mainloop()

if __name__ == '__main__':
    main()
