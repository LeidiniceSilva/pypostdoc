import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap

# Create a modified jet colormap starting with white
jet = cm.get_cmap('jet', 256)
new_colors = jet(np.linspace(0, 1, 256))
new_colors[0] = np.array([1, 1, 1, 1])  # Replace the first color with white (RGBA)

# Create a new colormap
new_cmap = ListedColormap(new_colors)

# Test the new colormap with a plot
data = np.random.rand(10, 10)

plt.imshow(data, cmap=new_cmap)
plt.colorbar()
plt.show()
exit()
