'''
Created on 4 April 2026
small U-Net used for channel-2 axon masks

@author: Dinghao Luo
'''

#%% imports
import torch
import torch.nn as nn


#%% blocks
class DoubleConv(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_channels, out_channels, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(out_channels),
            nn.ReLU(inplace=True),
            )

    def forward(self, x):
        return self.net(x)


class Down(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.net = nn.Sequential(
            nn.MaxPool2d(2),
            DoubleConv(in_channels, out_channels),
            )

    def forward(self, x):
        return self.net(x)


class Up(nn.Module):
    def __init__(self, in_channels, skip_channels, out_channels):
        super().__init__()
        self.up = nn.ConvTranspose2d(in_channels, out_channels, kernel_size=2, stride=2)
        self.conv = DoubleConv(out_channels + skip_channels, out_channels)

    def forward(self, x, skip):
        x = self.up(x)
        x = pad_to_match(x, skip)
        x = torch.cat([skip, x], dim=1)
        return self.conv(x)


class OutConv(nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.conv = nn.Conv2d(in_channels, out_channels, kernel_size=1)

    def forward(self, x):
        return self.conv(x)


#%% model
class SmallUNet(nn.Module):
    def __init__(self, in_channels=1, out_channels=1, base_channels=24, depth=4):
        super().__init__()
        if depth < 2:
            raise ValueError('depth must be at least 2')

        # Keep the first working model small here; the 512-pixel trial widens it in its YAML recipe.
        channels = [base_channels * (2 ** idx) for idx in range(depth)]
        self.in_conv = DoubleConv(in_channels, channels[0])
        self.downs = nn.ModuleList([
            Down(channels[idx - 1], channels[idx])
            for idx in range(1, depth)
            ])
        self.ups = nn.ModuleList([
            Up(channels[idx], channels[idx - 1], channels[idx - 1])
            for idx in range(depth - 1, 0, -1)
            ])
        self.out_conv = OutConv(channels[0], out_channels)

    def forward(self, x):
        features = [self.in_conv(x)]
        for down in self.downs:
            features.append(down(features[-1]))

        x = features[-1]
        for up, skip in zip(self.ups, reversed(features[:-1])):
            x = up(x, skip)

        return self.out_conv(x)


def build_model(config):
    return SmallUNet(
        in_channels=config.get('in_channels', 1),
        out_channels=config.get('out_channels', 1),
        base_channels=config.get('base_channels', 24),
        depth=config.get('depth', 4),
        )


def pad_to_match(x, reference):
    diff_y = reference.size(2) - x.size(2)
    diff_x = reference.size(3) - x.size(3)

    if diff_y == 0 and diff_x == 0:
        return x

    return nn.functional.pad(
        x,
        [diff_x // 2, diff_x - diff_x // 2,
         diff_y // 2, diff_y - diff_y // 2],
        )
