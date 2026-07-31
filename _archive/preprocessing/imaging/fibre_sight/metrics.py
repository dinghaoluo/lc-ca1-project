'''
Created on 20 April 2026
Archived on 23 July 2026
pixel and connected-component scores for held-out labels

@author: Dinghao Luo
'''

#%% imports
import numpy as np


#%% masks
def mask_scores(pred_mask, target_mask, eps=1e-8):
    pred_mask = np.asarray(pred_mask).astype(bool)
    target_mask = np.asarray(target_mask).astype(bool)

    tp = np.sum(pred_mask & target_mask)
    fp = np.sum(pred_mask & ~target_mask)
    fn = np.sum(~pred_mask & target_mask)

    precision = tp / (tp + fp + eps)
    recall = tp / (tp + fn + eps)
    dice = 2 * tp / (2 * tp + fp + fn + eps)
    iou = tp / (tp + fp + fn + eps)
    f2 = 5 * tp / (5 * tp + 4 * fn + fp + eps)

    return {
        'dice': float(dice),
        'f2': float(f2),
        'iou': float(iou),
        'precision': float(precision),
        'recall': float(recall),
        'true_positive_pixels': int(tp),
        'false_positive_pixels': int(fp),
        'false_negative_pixels': int(fn),
        }


def component_counts(pred_labelled, target_labelled):
    pred_ids = set(np.unique(pred_labelled)) - {0}
    target_ids = set(np.unique(target_labelled)) - {0}

    return {
        'predicted_components': len(pred_ids),
        'target_components': len(target_ids),
        }
