'''
Created on 12 May 2026
zooming, panning, and ROI display helpers for the workbench

@author: Dinghao Luo
'''

#%% imports
import colorsys

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import Qt


#%% colours
def generate_distinct_colours(n):
    colours = []
    if n <= 0:
        return colours

    hue = 0.0
    # walking around the hue circle keeps neighbouring ROI labels distinct
    golden_ratio = 0.61803398875
    for _ in range(n):
        hue = (hue + golden_ratio) % 1.0
        colours.append(colorsys.hsv_to_rgb(hue, 0.68, 0.92))
    return colours


def normalise_for_display(image):
    image = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(image)
    if not np.any(finite):
        return np.zeros_like(image, dtype=np.float32)

    low, high = np.nanpercentile(image[finite], [1, 99.7])
    if high <= low:
        return np.zeros_like(image, dtype=np.float32)

    image = (image - low) / (high - low)
    return np.clip(image, 0, 1)


def squeeze_image(image):
    image = np.asarray(image)
    if image.ndim > 2:
        image = np.squeeze(image)
    if image.ndim != 2:
        raise ValueError(f'expected a 2D image, got shape {image.shape}')
    return image


#%% canvas
class ZoomableCanvas(FigureCanvas):
    def __init__(self, figure, ax):
        super().__init__(figure)
        self.ax = ax
        self._drag_start_pos = None
        self._drag_start_xlim = None
        self._drag_start_ylim = None

    def mousePressEvent(self, event):
        if event.button() == Qt.RightButton:
            self._drag_start_pos = event.pos()
            self._drag_start_xlim = self.ax.get_xlim()
            self._drag_start_ylim = self.ax.get_ylim()
            return

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if event.buttons() & Qt.RightButton and self._drag_start_pos is not None:
            self.pan_view(event)
            return

        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.RightButton:
            self._drag_start_pos = None
            return

        super().mouseReleaseEvent(event)

    def wheelEvent(self, event):
        if not self.ax.images:
            return

        x_mouse = event.position().x()
        y_mouse = event.position().y()
        x_data, y_data = self.ax.transData.inverted().transform((x_mouse, y_mouse))

        scale = 0.85 if event.angleDelta().y() > 0 else 1.18
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        new_width = (xlim[1] - xlim[0]) * scale
        new_height = (ylim[1] - ylim[0]) * scale
        rel_x = (x_data - xlim[0]) / (xlim[1] - xlim[0])
        rel_y = (y_data - ylim[0]) / (ylim[1] - ylim[0])

        self.ax.set_xlim(x_data - new_width * rel_x, x_data + new_width * (1 - rel_x))
        self.ax.set_ylim(y_data - new_height * rel_y, y_data + new_height * (1 - rel_y))
        self.draw_idle()

    def pan_view(self, event):
        dx = event.pos().x() - self._drag_start_pos.x()
        dy = event.pos().y() - self._drag_start_pos.y()

        trans = self.ax.transData.inverted()
        dx_data = trans.transform((0, 0))[0] - trans.transform((dx, 0))[0]
        dy_data = trans.transform((0, dy))[1] - trans.transform((0, 0))[1]

        self.ax.set_xlim(
            self._drag_start_xlim[0] + dx_data,
            self._drag_start_xlim[1] + dx_data,
            )
        self.ax.set_ylim(
            self._drag_start_ylim[0] + dy_data,
            self._drag_start_ylim[1] + dy_data,
            )
        self.draw_idle()

    def reset_view(self):
        if not self.ax.images:
            return

        image = self.ax.images[0].get_array()
        self.ax.set_xlim(0, image.shape[1])
        self.ax.set_ylim(image.shape[0], 0)
        self.draw_idle()
