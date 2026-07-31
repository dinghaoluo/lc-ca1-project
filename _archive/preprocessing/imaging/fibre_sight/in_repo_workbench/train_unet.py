'''
Created on 10 April 2026
Modified on 19 May 2026
Modified on 23 June 2026
Modified on 23 July 2026 to keep the loss experiments beside the training loop
train the small U-Net and keep its run record

@author: Dinghao Luo
'''

#%% imports
from contextlib import nullcontext
from datetime import datetime
from pathlib import Path
import argparse
import csv
import sys

import numpy as np
import torch
import torch.nn.functional as F

if __package__ in {None, ''}:
    repo_root = Path(__file__).resolve().parents[3]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from preprocessing.imaging.fibre_sight._repo import add_repo_paths, default_figure_root
from preprocessing.imaging.fibre_sight.config import get_section, load_config, resolve_path, save_config
from preprocessing.imaging.fibre_sight.dataset import AxonROIDataset
from preprocessing.imaging.fibre_sight.model import build_model

REPO_ROOT = add_repo_paths()


#%% loss
# BCE was the first baseline; Dice helped once the foreground occupied only a small part of each crop.
def soft_dice_loss(logits, targets, eps=1e-6):
    probs = torch.sigmoid(logits)
    dims = tuple(range(1, probs.ndim))

    intersection = torch.sum(probs * targets, dim=dims)
    denominator = torch.sum(probs, dim=dims) + torch.sum(targets, dim=dims)
    dice = (2 * intersection + eps) / (denominator + eps)
    return 1 - dice.mean()


def bce_dice_loss(logits, targets, bce_weight=1.0, dice_weight=1.0):
    bce = F.binary_cross_entropy_with_logits(logits, targets)
    dice = soft_dice_loss(logits, targets)
    loss = bce_weight * bce + dice_weight * dice

    return loss, {
        'loss': float(loss.detach().cpu()),
        'bce': float(bce.detach().cpu()),
        'dice_loss': float(dice.detach().cpu()),
        }


def soft_tversky_loss(logits, targets, alpha=0.3, beta=0.7, eps=1e-6):
    probs = torch.sigmoid(logits)
    dims = tuple(range(2, probs.ndim))

    tp = torch.sum(probs * targets, dim=dims)
    fp = torch.sum(probs * (1 - targets), dim=dims)
    fn = torch.sum((1 - probs) * targets, dim=dims)

    score = (tp + eps) / (tp + alpha * fp + beta * fn + eps)
    return 1 - score.mean()


def segmentation_loss(logits, targets, config):
    mode = config.get('mode', 'bce_dice')
    channel_weights = config.get('channel_weights', None)
    pos_weight = config.get('pos_weight', None)

    if channel_weights is not None:
        weights = torch.as_tensor(channel_weights, device=logits.device, dtype=logits.dtype)
        weights = weights.view(1, -1, 1, 1)
    else:
        weights = 1.0

    if pos_weight is not None:
        pos_weight = torch.as_tensor(pos_weight, device=logits.device, dtype=logits.dtype)
        pos_weight = pos_weight.view(1, -1, 1, 1)

    bce_raw = F.binary_cross_entropy_with_logits(
        logits,
        targets,
        pos_weight=pos_weight,
        reduction='none',
        )
    bce = torch.mean(bce_raw * weights)

    if mode == 'bce_dice':
        seg = soft_dice_loss(logits * weights, targets)
        seg_name = 'dice_loss'
    elif mode == 'bce_tversky':
        seg = soft_tversky_loss(
            logits,
            targets,
            alpha=config.get('tversky_alpha', 0.3),
            beta=config.get('tversky_beta', 0.7),
            )
        seg_name = 'tversky_loss'
    else:
        raise ValueError(f'unknown loss mode: {mode}')

    loss = config.get('bce_weight', 1.0) * bce + config.get('seg_weight', 1.0) * seg
    return loss, {
        'loss': float(loss.detach().cpu()),
        'bce': float(bce.detach().cpu()),
        seg_name: float(seg.detach().cpu()),
        }


#%% cli
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--config',
        type=Path,
        default=Path(__file__).resolve().parent / 'configs' / 'dlight_hpc_lc_opto_unet.yaml',
        )
    return parser.parse_args()


#%% setup
def get_device(device_name):
    import torch

    if device_name == 'auto':
        return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    return torch.device(device_name)


def make_dataloaders(config):
    from torch.utils.data import DataLoader

    data_cfg = get_section(config, 'data')
    train_cfg = get_section(config, 'train')

    manifest_path = resolve_path(data_cfg['manifest'], REPO_ROOT)
    normalise_percentiles = data_cfg.get('normalise_percentiles', [1, 99.7])

    train_dataset = AxonROIDataset(
        manifest_path,
        split=data_cfg.get('train_split', 'train'),
        patch_size=data_cfg.get('patch_size', 256),
        patches_per_image=data_cfg.get('patches_per_image', 24),
        foreground_fraction=data_cfg.get('foreground_fraction', 0.75),
        normalise_percentiles=normalise_percentiles,
        augment=True,
        cache_images=data_cfg.get('cache_images', False),
        target_mode=data_cfg.get('target_mode', 'foreground'),
        support_radius=data_cfg.get('support_radius', 3),
        seed=train_cfg.get('seed', 7),
        )
    val_dataset = AxonROIDataset(
        manifest_path,
        split=data_cfg.get('val_split', 'val'),
        patch_size=data_cfg.get('patch_size', 256),
        patches_per_image=max(1, data_cfg.get('val_patches_per_image', 8)),
        foreground_fraction=data_cfg.get('foreground_fraction', 0.75),
        normalise_percentiles=normalise_percentiles,
        augment=False,
        cache_images=data_cfg.get('cache_images', False),
        target_mode=data_cfg.get('target_mode', 'foreground'),
        support_radius=data_cfg.get('support_radius', 3),
        seed=train_cfg.get('seed', 7) + 1000,
        )

    loader_args = {
        'batch_size': train_cfg.get('batch_size', 8),
        'num_workers': train_cfg.get('num_workers', 0),
        'pin_memory': train_cfg.get('pin_memory', True),
        }
    train_loader = DataLoader(train_dataset, shuffle=True, **loader_args)
    # validation keeps the image unchanged; only crop selection moves with the epoch
    val_loader = DataLoader(val_dataset, shuffle=False, **loader_args)

    return train_dataset, val_dataset, train_loader, val_loader


def prepare_run_dir(config):
    train_cfg = get_section(config, 'train')
    out_root = resolve_path(train_cfg.get('out_dir', 'data/imaging/fibre_sight/runs'), REPO_ROOT)
    run_name = train_cfg.get('run_name', '')

    if not run_name:
        run_name = datetime.now().strftime('%Y%m%d_%H%M%S')

    run_dir = out_root / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    save_config(config, run_dir / 'config.yaml')
    return run_dir


def prepare_figure_run_dir(run_name):
    run_dir = default_figure_root() / 'runs' / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


#%% training
def train_one_epoch(model, loader, optimiser, scaler, device, config, epoch):
    from tqdm import tqdm

    model.train()
    loader.dataset.set_epoch(epoch)
    rows = []
    amp = use_amp(config, device)

    for batch in tqdm(loader, desc=f'train {epoch}', leave=False):
        images = batch['image'].to(device, non_blocking=True)
        masks = batch['mask'].to(device, non_blocking=True)

        optimiser.zero_grad(set_to_none=True)
        with autocast_context(amp):
            logits = model(images)
            loss, loss_parts = loss_from_config(logits, masks, config)

        if amp:
            scaler.scale(loss).backward()
            scaler.step(optimiser)
            scaler.update()
        else:
            loss.backward()
            optimiser.step()

        rows.append(loss_parts)

    return average_rows(rows)


def validate(model, loader, device, config, epoch):
    import torch
    from tqdm import tqdm

    model.eval()
    loader.dataset.set_epoch(epoch)
    rows = []
    amp = use_amp(config, device)

    with torch.no_grad():
        for batch in tqdm(loader, desc=f'val {epoch}', leave=False):
            images = batch['image'].to(device, non_blocking=True)
            masks = batch['mask'].to(device, non_blocking=True)

            with autocast_context(amp):
                logits = model(images)
                _, loss_parts = loss_from_config(logits, masks, config)

            loss_parts['pixel_dice'] = batch_dice(logits, masks)
            rows.append(loss_parts)

    return average_rows(rows)


def loss_from_config(logits, masks, config):
    loss_cfg = get_section(config, 'loss', {'bce_weight': 1.0, 'dice_weight': 1.0})
    return segmentation_loss(logits, masks, loss_cfg)


def batch_dice(logits, masks, threshold=0.5, eps=1e-8):
    import torch

    pred = torch.sigmoid(logits[:, :1]) >= threshold
    target = masks[:, :1] >= 0.5
    dims = tuple(range(1, pred.ndim))
    tp = torch.sum(pred & target, dim=dims)
    fp = torch.sum(pred & ~target, dim=dims)
    fn = torch.sum(~pred & target, dim=dims)
    dice = (2 * tp.float() + eps) / (2 * tp.float() + fp.float() + fn.float() + eps)
    return float(dice.mean().detach().cpu())


def use_amp(config, device):
    train_cfg = get_section(config, 'train')
    return bool(train_cfg.get('amp', True) and device.type == 'cuda')


def autocast_context(enabled):
    if not enabled:
        return nullcontext()

    import torch
    return torch.amp.autocast('cuda')


#%% output
def save_checkpoint(path, model, optimiser, epoch, best_score, config):
    import torch

    checkpoint = {
        'epoch': epoch,
        'best_score': best_score,
        'model_state': model.state_dict(),
        'optimiser_state': optimiser.state_dict(),
        'model_config': get_section(config, 'model'),
        'data_config': get_section(config, 'data'),
        'postprocess_config': get_section(config, 'postprocess'),
        'config': config,
        }
    torch.save(checkpoint, path)


def write_history(history, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = sorted({key for row in history for key in row.keys()})

    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(history)


def plot_history(history, path):
    import matplotlib.pyplot as plt

    add_repo_paths()
    from common_functions import mpl_formatting
    mpl_formatting()

    epochs = [row['epoch'] for row in history]
    fig, axes = plt.subplots(1, 2, figsize=(7, 3), constrained_layout=True)

    axes[0].plot(epochs, [row['train_loss'] for row in history], label='train')
    axes[0].plot(epochs, [row['val_loss'] for row in history], label='val')
    axes[0].set_xlabel('epoch')
    axes[0].set_ylabel('loss')
    axes[0].legend(frameon=False)

    axes[1].plot(epochs, [row.get('val_pixel_dice', np.nan) for row in history])
    axes[1].set_xlabel('epoch')
    axes[1].set_ylabel('validation Dice')

    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)


def average_rows(rows):
    if len(rows) == 0:
        return {}

    keys = rows[0].keys()
    return {
        key: float(np.mean([row[key] for row in rows]))
        for key in keys
        }


#%% main
def main():
    import torch

    args = parse_args()
    config = load_config(args.config)
    run_dir = prepare_run_dir(config)
    figure_run_dir = prepare_figure_run_dir(run_dir.name)

    train_cfg = get_section(config, 'train')
    device = get_device(train_cfg.get('device', 'auto'))
    print(f'training on {device}')
    print(f'run directory: {run_dir}')

    train_dataset, val_dataset, train_loader, val_loader = make_dataloaders(config)
    print(f'train sessions: {len(train_dataset.rows)}')
    print(f'validation sessions: {len(val_dataset.rows)}')

    model = build_model(get_section(config, 'model')).to(device)
    optimiser = torch.optim.AdamW(
        model.parameters(),
        lr=train_cfg.get('learning_rate', 3e-4),
        weight_decay=train_cfg.get('weight_decay', 1e-5),
        )
    scaler = torch.amp.GradScaler('cuda', enabled=use_amp(config, device))

    history = []
    best_score = -np.inf
    epochs = int(train_cfg.get('epochs', 80))

    for epoch in range(1, epochs + 1):
        train_stats = train_one_epoch(model, train_loader, optimiser, scaler, device, config, epoch)
        val_stats = validate(model, val_loader, device, config, epoch)

        row = {'epoch': epoch}
        row.update({f'train_{key}': value for key, value in train_stats.items()})
        row.update({f'val_{key}': value for key, value in val_stats.items()})
        history.append(row)

        val_score = row.get('val_pixel_dice', -np.inf)
        print(
            f'epoch {epoch:03d} | '
            f'train loss {row["train_loss"]:.4f} | '
            f'val loss {row["val_loss"]:.4f} | '
            f'val Dice {val_score:.4f}'
            )

        # Keep the checkpoint with the best validation Dice; this is the mask overlap I inspect.
        if val_score > best_score:
            best_score = val_score
            save_checkpoint(run_dir / 'best.pt', model, optimiser, epoch, best_score, config)
        # best.pt is used for prediction; latest.pt records where each long run ended
        save_checkpoint(run_dir / 'latest.pt', model, optimiser, epoch, best_score, config)

        write_history(history, run_dir / 'history.csv')
        plot_history(history, figure_run_dir / 'history.png')


if __name__ == '__main__':
    main()
