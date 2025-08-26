# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "April 14, 2025"
__description__ = "This script make giff file"

import imageio.v2 as imageio
import os

# Define input/output paths
img_dir = '/leonardo/home/userexternal/mdasilva/leonardo_work/WORKSHOP2025/figs'
output_gif = '/leonardo/home/userexternal/mdasilva/leonardo_work/WORKSHOP2025/figs/animated_psl_RegCM5_2023Oct19-27th.gif'

# Get all image file names sorted by timestamp
images = sorted([img for img in os.listdir(img_dir) if img.endswith('.png')])

# Create gif
frames = []
for img_name in images:
    img_path = os.path.join(img_dir, img_name)
    image = imageio.imread(img_path)
    frames.append(image)

# Save as gif (duration in seconds per frame)
imageio.mimsave(output_gif, frames, duration=1, loop=0)

print(f"GIF saved to: {output_gif}")
