# -*- coding: utf-8 -*-
'''
Created on Mon Apr 15 16:50:53 2024

dF/F, reference images and optogenetic timing for the imaging extractions
modified: added GPU acceleration using cupy, 1 Nov 2024 Dinghao
modified: added opto pulse parsing for dLight, GRABNE and nLight extraction

@author: Dinghao Luo
@contributor: Jingyu Cao
'''

#%% imports
import numpy as np
import scipy.ndimage
from pathlib import Path
from tqdm import tqdm
import gc
import matplotlib.pyplot as plt

from behaviour_functions import get_next_line

# -------------------------------
# baseline calculation functions
# -------------------------------
def convolve_gaussian(
        arr,
        sigma,
        t_axis=0,
        GPU_AVAILABLE=False
        ):
    '''
    convolve an array with a gaussian kernel along a specified axis.

    parameters:
    - arr: numpy array
        input array to convolve (1D or multi-dimensional).
    - sigma: float
        standard deviation of the gaussian kernel.
    - t_axis: int, default=0
        axis along which convolution is performed.
    - GPU_AVAILABLE: bool, default=False
        whether to perform convolution on GPU using cupy.

    returns:
    - convolved array: numpy array
        gaussian-convolved array with same shape as input.
    '''
    if GPU_AVAILABLE:
        import cupyx.scipy.ndimage as cpximg
        # arr_gpu = cp.array(arr)
        # we assume that the array in on GPU already
        return cpximg.gaussian_filter1d(arr, sigma, axis=t_axis, mode='reflect')
    else:
        return scipy.ndimage.gaussian_filter1d(arr, sigma, axis=t_axis, mode='reflect')

def rolling_min(
        arr,
        win,
        t_axis=0,
        GPU_AVAILABLE=False
        ):
    '''
    calculate rolling minimum along a specified axis of an array.

    parameters:
    - arr: numpy array
        input array (1D or multi-dimensional).
    - win: int
        size of the rolling window.
    - t_axis: int, default=0
        axis along which rolling minimum is computed.
    - GPU_AVAILABLE: bool, default=False
        whether to perform computation on GPU using cupy.

    returns:
    - minimum-filtered array: numpy array
        rolling minimum array with same shape as input.
    '''
    if GPU_AVAILABLE:
        import cupyx.scipy.ndimage as cpximg
        # arr_gpu = cp.array(arr)
        return cpximg.minimum_filter1d(arr, size=win, axis=t_axis, mode='reflect')
    else:
        return scipy.ndimage.minimum_filter1d(arr, size=win, axis=t_axis, mode='reflect')

def rolling_max(
        arr,
        win,
        t_axis=0,
        GPU_AVAILABLE=False
        ):
    '''
    calculate rolling maximum along a specified axis of an array.

    parameters:
    - arr: numpy array
        input array (1D or multi-dimensional).
    - win: int
        size of the rolling window.
    - t_axis: int, default=0
        axis along which rolling maximum is computed.
    - GPU_AVAILABLE: bool, default=False
        whether to perform computation on GPU using cupy.

    returns:
    - maximum-filtered array: numpy array
        rolling maximum array with same shape as input.
    '''
    if GPU_AVAILABLE:
        import cupyx.scipy.ndimage as cpximg
        # arr_gpu = cp.array(arr)
        return cpximg.maximum_filter1d(arr, size=win, axis=t_axis, mode='reflect')
    else:
        return scipy.ndimage.maximum_filter1d(arr, size=win, axis=t_axis, mode='reflect')

def rolling_percentile(
        arr,
        window,
        pct,
        t_axis,
        GPU_AVAILABLE,
        progress_desc=None,
        ):
    '''
    compute rolling percentile along t_axis using a centred window.

    parameters:
    - arr: np.ndarray or cp.ndarray
    - window: int, number of frames
    - pct: float, percentile (e.g. 10)
    - t_axis: int, time axis
    - GPU_AVAILABLE: bool

    returns:
    - out: same shape as arr, rolling percentile baseline
    '''
    if GPU_AVAILABLE:
        import cupy as cp
        xp = cp
    else:
        xp = np

    pad = window // 2
    T = arr.shape[t_axis]
    label = progress_desc or 'rolling percentile baseline'

    # pad reflectively along the time axis
    pad_width = [(0,0)] * arr.ndim
    pad_width[t_axis] = (pad, pad)
    if progress_desc:
        print(f'{label}: reflective padding...', flush=True)
    arr_pad = xp.pad(arr, pad_width, mode='reflect')
    if GPU_AVAILABLE:
        xp.cuda.Stream.null.synchronize()
    if progress_desc:
        print(f'{label}: allocating output...', flush=True)

    # output
    out = xp.empty_like(arr)

    # move t_axis to last for easier slicing
    arr_pad = xp.moveaxis(arr_pad, t_axis, -1)
    out_t   = xp.moveaxis(out, t_axis, -1)

    frame_iter = range(T)
    if progress_desc:
        print(f'{label}: starting {T}-frame percentile loop...', flush=True)
        frame_iter = tqdm(frame_iter, desc=progress_desc, unit='frame')

    for i in frame_iter:
        w = arr_pad[..., i:i+window]
        out_t[..., i] = xp.percentile(w, pct, axis=-1)

    # restore original axis order
    out = xp.moveaxis(out_t, -1, t_axis)
    return out
# -----------------------------------
# baseline calculation functions end
# -----------------------------------

# --------------------------
# dFF calculation functions
# --------------------------
def calculate_dFF(
        F_array,
        t_axis,
        sigma=300,
        GPU_AVAILABLE=False,
        CHUNK=False,
        chunk_size=2000
        ):
    '''
    calculate dF/F for fluorescence traces using gaussian smoothing and rolling min-max baseline.

    parameters:
    - F_array: np.ndarray
        fluorescence trace array, can be 1D or ND (time must be along t_axis).
    - t_axis: int
        axis corresponding to time.
    - sigma: int, default=300
        standard deviation for the gaussian smoothing filter.
    - GPU_AVAILABLE: bool, default=False
        whether to use CuPy to accelerate filtering operations.
    - CHUNK: bool, default=False
        whether to apply the dF/F calculation in memory-safe chunks.
    - chunk_size: int, default=10000
        number of timepoints per chunk if CHUNK is True.

    returns:
    - dFF: np.ndarray
        array of same shape as F_array, containing the computed dF/F values.
    '''
    window = sigma * 6
    T = F_array.shape[t_axis]

    if not CHUNK:
        # full-array processing
        if GPU_AVAILABLE:
            import cupy as cp
            F_array = cp.array(F_array, dtype=cp.float32)
        else:
            F_array = F_array.astype(np.float32, copy=False)

        baseline = convolve_gaussian(F_array, sigma, t_axis, GPU_AVAILABLE)
        baseline = rolling_min(baseline, window, t_axis, GPU_AVAILABLE)
        baseline = rolling_max(baseline, window, t_axis, GPU_AVAILABLE)
        dFF = (F_array - baseline) / baseline

        return dFF.get() if GPU_AVAILABLE else dFF

    # chunked processing
    pad = window // 2
    slices = []

    if GPU_AVAILABLE:
        import cupy as cp

        for start in tqdm(range(0, T, chunk_size),
                          desc='chunked dFF calculation on GPU...'):
            chunk_start = max(0, start - pad)
            chunk_end = min(T, start + chunk_size + pad)

            # extract and move to GPU
            slicer = [slice(None)] * F_array.ndim
            slicer[t_axis] = slice(chunk_start, chunk_end)
            chunk_gpu = cp.array(F_array[tuple(slicer)], dtype=cp.float32)

            baseline = convolve_gaussian(chunk_gpu, sigma, t_axis, GPU_AVAILABLE=True)
            baseline = rolling_min(baseline, window, t_axis, GPU_AVAILABLE=True)
            baseline = rolling_max(baseline, window, t_axis, GPU_AVAILABLE=True)

            dFF_chunk = (chunk_gpu - baseline) / baseline
            dFF_chunk = dFF_chunk.get()

            # trim padding
            # first calculate how much padding was actually added
            left_trim  = start - chunk_start
            right_trim = chunk_end - (start + chunk_size)

            # define output slice
            slicer_out = [slice(None)] * dFF_chunk.ndim
            time_len = dFF_chunk.shape[t_axis]

            start_idx = left_trim
            end_idx = time_len - right_trim if right_trim > 0 else time_len

            # actual output slice
            slicer_out[t_axis] = slice(start_idx, end_idx)

            slices.append(dFF_chunk[tuple(slicer_out)])

        return np.concatenate(slices, axis=t_axis)

    else:
        for start in tqdm(range(0, T, chunk_size),
                          desc='chunked dFF calculation on CPU...'):
            chunk_start = max(0, start - pad)
            chunk_end = min(T, start + chunk_size + pad)

            slicer = [slice(None)] * F_array.ndim
            slicer[t_axis] = slice(chunk_start, chunk_end)
            chunk = F_array[tuple(slicer)].astype(np.float32, copy=False)

            baseline = convolve_gaussian(chunk, sigma, t_axis, GPU_AVAILABLE=False)
            baseline = rolling_min(baseline, window, t_axis, GPU_AVAILABLE=False)
            baseline = rolling_max(baseline, window, t_axis, GPU_AVAILABLE=False)

            dFF_chunk = (chunk - baseline) / baseline

            # trim padding
            # first calculate how much padding was actually added
            left_trim  = start - chunk_start
            right_trim = chunk_end - (start + chunk_size)

            # define output slice
            slicer_out = [slice(None)] * dFF_chunk.ndim
            time_len = dFF_chunk.shape[t_axis]

            start_idx = left_trim
            end_idx = time_len - right_trim if right_trim > 0 else time_len

            # actual output slice
            slicer_out[t_axis] = slice(start_idx, end_idx)

            slices.append(dFF_chunk[tuple(slicer_out)])

        return np.concatenate(slices, axis=t_axis)

def calculate_dFF_percentile(
        F_array,
        t_axis,
        window_size,
        GPU_AVAILABLE=False,
        CHUNK=False,
        chunk_size=2000,
        pct=20,
        return_baseline=False,
        device_name='',
        progress_desc=None,
        ):
    '''
    calculate dF/F using rolling-percentile baseline
    '''
    # originally used a window_size = 1800 frames; now requires input
    T = F_array.shape[t_axis]

    print(f'calculating dF/F using rolling {pct}th-percentile {window_size}-frame baselines')

    # ------------------
    # non-chunked
    # ------------------
    if not CHUNK:
        if GPU_AVAILABLE:
            import cupy as cp
            if device_name != '':
                print(f'device: {device_name}')
            if progress_desc:
                print(f'{progress_desc}: transferring fluorescence array to GPU...', flush=True)
            F = cp.asarray(F_array, dtype=cp.float32)
            cp.cuda.Stream.null.synchronize()
            if progress_desc:
                print(f'{progress_desc}: transfer complete.', flush=True)
        else:
            if progress_desc:
                print(f'{progress_desc}: preparing fluorescence array on CPU...', flush=True)
            F = F_array.astype(np.float32, copy=False)

        baseline = rolling_percentile(
            F,
            window_size,
            pct,
            t_axis,
            GPU_AVAILABLE,
            progress_desc=progress_desc or 'rolling percentile baseline',
            )
        dFF = (F - baseline) / baseline

        if GPU_AVAILABLE:
            dFF_out = dFF.get()
            baseline_out = baseline.get() if return_baseline else None
            del F, baseline, dFF
            cp.cuda.Stream.null.synchronize()
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
            return (dFF_out, baseline_out) if return_baseline else dFF_out

        return (dFF, baseline) if return_baseline else dFF

    # ------------------
    # chunked
    # ------------------
    pad = window_size // 2
    dff_slices = []
    baseline_slices = [] if return_baseline else None

    if GPU_AVAILABLE:
        import cupy as cp
        xp = cp
        device = 'GPU'
    else:
        xp = np
        device = 'CPU'

    for start in tqdm(range(0, T, chunk_size),
                      desc=f'Chunked dFF calculation on {device}...'):

        chunk_start = max(0, start - pad)
        chunk_end   = min(T, start + chunk_size + pad)

        slicer = [slice(None)] * F_array.ndim
        slicer[t_axis] = slice(chunk_start, chunk_end)

        chunk = F_array[tuple(slicer)].astype(np.float32, copy=False)
        if GPU_AVAILABLE:
            chunk = xp.asarray(chunk)

        baseline = rolling_percentile(chunk, window_size, pct, t_axis, GPU_AVAILABLE)
        dFF_chunk = (chunk - baseline) / baseline

        if GPU_AVAILABLE:
            dFF_chunk = dFF_chunk.get()
            baseline = baseline.get()
            del chunk
            cp.cuda.Stream.null.synchronize()
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()

        # trim padding
        # first calculate how much padding was actually added
        left_trim  = start - chunk_start
        right_trim = chunk_end - (start + chunk_size)

        # define output slice
        slicer_out = [slice(None)] * dFF_chunk.ndim
        time_len = dFF_chunk.shape[t_axis]

        start_idx = left_trim
        end_idx = time_len - right_trim if right_trim > 0 else time_len

        # actual output slice
        slicer_out[t_axis] = slice(start_idx, end_idx)

        dff_slices.append(dFF_chunk[tuple(slicer_out)])

        if return_baseline:
            baseline_slices.append(baseline[tuple(slicer_out)])

    dFF_out = np.concatenate(dff_slices, axis=t_axis)

    if return_baseline:
        baseline_out = np.concatenate(baseline_slices, axis=t_axis)
        return dFF_out, baseline_out

    return dFF_out
# ------------------------------
# dFF calculation functions end
# ------------------------------

def spatial_gaussian_filter(
        movie,
        sigma_spatial=1,
        GPU_AVAILABLE=False,
        CHUNK=False,
        chunk_size=2000,
        inplace=False,
        show_progress=True
        ):
    '''
    apply spatial gaussian filtering (only over y and x) to a 3D movie in memory-efficient chunks.

    parameters:
    - movie: np.ndarray
        3D array (frames x height x width).
    - sigma_spatial: float, default=1
        spatial std for gaussian filter over height and width.
    - GPU_AVAILABLE: bool, default=False
        use CuPy GPU filtering if available.
    - CHUNK: bool, default=False
        whether to filter in smaller chunks along the frame axis.
    - chunk_size: int, default=2000
        number of frames per chunk.
    - inplace: bool, default=False
        when CHUNK is True, write filtered chunks back into movie to avoid
        allocating a second full-size movie.
    - show_progress: bool, default=True
        show tqdm progress bars during chunked filtering.

    returns:
    - filtered_movie: np.ndarray
        array of same shape as movie, spatially filtered.
    '''
    T, H, W = movie.shape

    if not CHUNK:
        # full-array processing
        if GPU_AVAILABLE:
            import cupy as cp
            from cupyx.scipy.ndimage import gaussian_filter
            movie_gpu = cp.array(movie, dtype=cp.float32)
            filtered = gaussian_filter(movie_gpu, sigma=(0, sigma_spatial, sigma_spatial))
            return cp.asnumpy(filtered)
        else:
            from scipy.ndimage import gaussian_filter
            return gaussian_filter(movie, sigma=(0, sigma_spatial, sigma_spatial))

    # chunked processing
    filtered = movie if inplace else np.empty_like(movie)

    if GPU_AVAILABLE:
        import cupy as cp
        from cupyx.scipy.ndimage import gaussian_filter

        for start in tqdm(range(0, T, chunk_size),
                          desc='chunked spatial-filtering frames on GPU...',
                          disable=not show_progress):
            end = min(start + chunk_size, T)
            chunk_gpu = cp.asarray(movie[start:end], dtype=cp.float32)
            filtered_chunk = gaussian_filter(chunk_gpu, sigma=(0, sigma_spatial, sigma_spatial))
            filtered[start:end] = cp.asnumpy(filtered_chunk)
            del chunk_gpu, filtered_chunk
            cp.cuda.Stream.null.synchronize()
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()

    else:
        from scipy.ndimage import gaussian_filter

        for start in tqdm(range(0, T, chunk_size),
                          desc='chunked spatial-filtering frames on CPU...',
                          disable=not show_progress):
            end = min(start + chunk_size, T)
            filtered[start:end] = gaussian_filter(movie[start:end], sigma=(0, sigma_spatial, sigma_spatial))

    return filtered


#%% reference images
def plot_reference(mov, outpath=r'', recname='',
        grids=-1, stride=-1, dim=512, channel=1, frames=1000,
        GPU_AVAILABLE=False, save_figure=True, figure_outpath=None):
    '''
    plot a reference image (mean Z-projection) with optional grid annotations

    parameters:
    mov : np.ndarray or cupy.ndarray
        imaging data as a 3D array (frames x height x width)
    grids : list[int] or int, optional
        grid line positions for annotation; set to -1 to disable grid processing
    stride : int, optional
        stride length between grid lines; required if grids is not -1
    dim : int, optional
        dimension of the imaging data for plotting (default: 512)
    channel : int, optional
        channel number to annotate in the title (default: 1)
    outpath : str, optional
        output path to save the reference image and array (default: '')
    GPU_AVAILABLE : bool, optional
        flag to enable GPU processing using cupy (default: False)

    returns:
    ref_im : np.ndarray
        processed reference image (mean Z-projection)
    '''
    proc_path = Path(outpath) / 'processed_data'
    proc_path.mkdir(parents=True, exist_ok=True)

    # choose how many frames to use
    n_frames = min(frames, mov.shape[0])

    # compute mean z-projection over subset of frames
    if GPU_AVAILABLE:
        import cupy as cp
        mov_gpu = cp.array(mov[:n_frames])  # move only the needed frames to vram
        ref_im_gpu = cp.mean(mov_gpu, axis=0)
        ref_im = ref_im_gpu.get()
    else:
        ref_im = np.mean(mov[:n_frames], axis=0)

    ref_im = post_processing_suite2p_gui(ref_im)

    if save_figure:
        figure_outpath = outpath if figure_outpath is None else figure_outpath
        figure_outpath = Path(figure_outpath)
        figure_outpath.mkdir(parents=True, exist_ok=True)

        if grids != -1:
            boundary_low = grids[0]
            boundary_high = grids[-1] + stride

        fig, ax = plt.subplots(figsize=(4, 4))
        ax.imshow(
            ref_im,
            aspect='auto',
            cmap='gist_gray',
            interpolation='none',
            extent=[0, dim, dim, 0],
            )

        if grids != -1:
            for i in range(len(grids)):
                ax.plot(
                    [grids[i], grids[i]],
                    [boundary_low, boundary_high],
                    color='grey',
                    linewidth=1,
                    alpha=.5,
                    )
                ax.plot(
                    [boundary_low, boundary_high],
                    [grids[i], grids[i]],
                    color='grey',
                    linewidth=1,
                    alpha=.5,
                    )

            ax.plot(
                [grids[-1] + stride, grids[-1] + stride],
                [boundary_low, boundary_high],
                color='grey',
                linewidth=1,
                alpha=.5,
                )
            ax.plot(
                [boundary_low, boundary_high],
                [grids[-1] + stride, grids[-1] + stride],
                color='grey',
                linewidth=1,
                alpha=.5,
                )

        ax.set(xlim=(0, dim), ylim=(0, dim))
        fig.suptitle(f'ref ch{channel}')
        fig.tight_layout()

        if grids != -1:
            out_stem = figure_outpath / f'{recname}_ref_ch{channel}_{stride}'
        else:
            out_stem = figure_outpath / f'{recname}_ref_ch{channel}'

        for ext in ['.png', '.pdf']:
            fig.savefig(f'{out_stem}{ext}', dpi=300, bbox_inches='tight')
        plt.close(fig)

    np.save(proc_path / f'ref_mat_ch{channel}.npy', ref_im)

    if GPU_AVAILABLE:
        del mov_gpu, ref_im_gpu
        gc.collect()
        cp.get_default_memory_pool().free_all_blocks()

    return ref_im

def post_processing_suite2p_gui(img_orig):
    '''
    apply percentile-based contrast normalisation and rescale to 8-bit for GUI display.

    parameters:
    - img_orig: np.ndarray
        input 2D image.

    returns:
    - img_proc: np.ndarray
        normalised image (uint8, range 0–255).
    '''
    # normalise to 1st and 99th percentile
    perc_low, perc_high = np.percentile(img_orig, [1, 99])
    img_proc = (img_orig - perc_low) / (perc_high - perc_low)
    img_proc = np.maximum(0, np.minimum(1, img_proc))

    # convert to uint8
    img_proc *= 255
    img_proc = img_proc.astype(np.uint8)

    return img_proc

#%% behaviour file processing
def process_txt_nobeh(txtfile):
    '''
    parse the reduced imaging-side txt format used when no wheel-based behaviour structure is needed.

    parameters:
    - txtfile: str
        path to the log file.

    returns:
    - logfile: dict
    '''
    curr_logfile = {}
    file = open(txtfile, 'r')

    line = get_next_line(file)

    pulse_times       = []
    pulse_parameters  = []
    frame_times       = []
    shutter_ON_times  = []
    shutter_OFF_times = []
    buzzer_times      = []
    wheel_dummy = 0

    while line[0].find('$') == 0:
        if line[0] == '$WE':
            wheel_dummy += 1
        if line[0] == '$PC':
            pulse_times.append(float(line[1]))
        if line[0] == '$PP':
            pulse_parameters.append([s for s in line[1:]])
        if line[0] == '$FM' and line[2] == '0':
            frame_times.append(float(line[1]))
        if line[0] == '$SC':
            shutter_ON_times.append(float(line[1]))
        if line[0] == '$SO':
            shutter_OFF_times.append(float(line[1]))
        if line[0] == '$BZ':
            buzzer_times.append(float(line[1]))
        line = get_next_line(file)

    curr_logfile['pulse_times']       = pulse_times
    curr_logfile['pulse_parameters']  = pulse_parameters
    curr_logfile['frame_times']       = frame_times
    curr_logfile['shutter_on_times']  = shutter_ON_times
    curr_logfile['shutter_off_times'] = shutter_OFF_times
    curr_logfile['buzzer_times']      = buzzer_times

    # flag for 'well, actually this one HAS behaviour', 24 June 2025
    curr_logfile['behaviour'] = False
    if wheel_dummy > 5:
        curr_logfile['behaviour'] = True

    return curr_logfile

def extract_opto_pulse_metadata(logfile, samp_freq, gap_ms=1000):
    '''
    extract shared pulse-train metadata for opto-imaging preprocessing scripts.

    returns a dictionary with pulse widths, taper info, split pulse trains,
    and derived train-duration/frame-interval quantities used downstream.
    '''
    pulse_parameters = logfile['pulse_parameters']
    pulse_params = pulse_parameters[-1]

    pulse_width_on = float(pulse_params[2]) / 1e6
    pulse_width = float(pulse_params[3]) / 1e6
    pulse_number = int(pulse_params[4])
    taper_enabled = int(pulse_params[7])
    taper_raw = int(pulse_params[8])
    taper_duration = taper_raw if taper_raw > 1000 else 0
    duty_cycle = f'{int(round(100 * pulse_width_on / pulse_width))}%'

    pulse_times = np.asarray(logfile['pulse_times'], dtype=float)
    split_idx = np.where(np.diff(pulse_times) >= gap_ms)[0] + 1
    pulse_trains = np.split(pulse_times, split_idx)
    total_pulses = int(len(pulse_times) / pulse_number)

    example_train = pulse_trains[0]
    train_rel = example_train - example_train[0]
    total_train_duration_ms = round(train_rel[-1] + pulse_width * 1000)

    total_train_duration_s = total_train_duration_ms / 1000
    total_train_duration_frames = int(total_train_duration_s * samp_freq)
    min_interval_frames = 2 * total_train_duration_frames

    return {
        'pulse_times': pulse_times,
        'pulse_params': pulse_params,
        'pulse_trains': pulse_trains,
        'pulse_width_on': pulse_width_on,
        'pulse_width': pulse_width,
        'pulse_number': pulse_number,
        'taper_enabled': taper_enabled,
        'taper_duration': taper_duration,
        'duty_cycle': duty_cycle,
        'total_pulses': total_pulses,
        'total_train_duration_ms': total_train_duration_ms,
        'total_train_duration_s': total_train_duration_s,
        'total_train_duration_frames': total_train_duration_frames,
        'min_interval_frames': min_interval_frames,
    }

def detect_step_pairs(
        trace, zthr=100, min_interval_frames=30,
        plot_debug=False, debug_savepath=None):
    '''
    detect shutter/PMT step onset & offset pairs.
    works even if step directions flip across recordings.

    returns:
        onsets  : list of frame indices
        offsets : list of frame indices
    '''
    d = np.diff(trace)

    # MAD noise estimate
    m = np.nanmedian(d)
    mad = np.nanmedian(np.abs(d - m))
    z = (d - m) / mad

    if plot_debug:
        fig, ax = plt.subplots(figsize=(3,3))
        ax.plot(z)
        ax.set(ylim=(-100,100))
        if debug_savepath is None:
            plt.show()
        else:
            fig.savefig(debug_savepath, dpi=150, bbox_inches='tight')
            plt.close(fig)

    # detect all large steps (ignore sign)
    cand = np.where(np.abs(z) > zthr)[0] + 1

    if len(cand) == 0:
        return [], []

    # find edges
    gaps = np.where(np.diff(cand) > min_interval_frames)[0] + 1

    # split into separate trains
    trains = np.split(cand, gaps)

    onsets = []
    offsets = []

    for t in trains:
        if len(t) < 2:
            continue
        onsets.append(t[0])
        offsets.append(t[-1])

    return onsets, offsets
