'''
Created on 8 April 2026
Archived on 23 July 2026
losses used while training the small U-Net

@author: Dinghao Luo
'''

#%% imports
import torch
import torch.nn.functional as F


#%% losses
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
