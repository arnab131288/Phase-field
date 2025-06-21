"""
Plot and save 2D fields using Matplotlib.
"""
import os
import matplotlib.pyplot as plt
from config import params

OUT = params['OUT_DIR']

def plot_field(field, title, prefix=''):
    """Render a field and save figure to output directory."""
    fig, ax = plt.subplots()
    cax = ax.imshow(
        field.T, origin='lower',
        extent=[0, params['NX']*params['DELTA_X'], 0, params['NY']*params['DELTA_Y']]
    )
    fig.colorbar(cax, label=title)
    ax.set(title=title, xlabel='x', ylabel='y')
    # Sanitize filename
    fname = f"{prefix}{title.lower().replace(' ', '_').replace(':','')}.png"
    path = os.path.join(OUT, fname)
    plt.savefig(path)
    plt.close(fig)
