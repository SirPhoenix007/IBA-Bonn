#!/usr/bin/env python3
# --------------------------------------------------------------
# second_xaxis_simple_line.py
#
# Minimal example:
#   * Primary x‑axis (top) → seconds
#   * Secondary x‑axis (bottom) → minutes (seconds / 60)
#   * Both axes share the same y‑range – no extra y‑axes appear.
#
# Requirements
#   pip install matplotlib numpy
#
# --------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

def sec_to_min(sec):
    """Conversion used for the secondary axis."""
    return sec / 60.0

def min_to_sec(minute):
    """Inverse conversion (Matplotlib needs the inverse)."""
    return minute * 60.0

def main():
    # ------------------------------------------------------------------
    # 1️⃣  Create some test data (seconds → sin wave)
    # ------------------------------------------------------------------
    x_seconds = np.linspace(0, 600, 500)                # 0 … 600 seconds
    y = np.sin(2 * np.pi * x_seconds / 120)            # 5‑period sine wave

    # ------------------------------------------------------------------
    # 2️⃣  Build the figure + primary axes
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 5))

    # Plot the data
    ax.plot(x_seconds, y, color='tab:blue', lw=2,
            label=r'$\sin\!\bigl(2\pi t/120\bigr)$')
    ax.set_xlabel('Time (seconds)', fontsize=12, color='tab:blue')
    ax.set_ylabel('Amplitude', fontsize=12, color='tab:blue')
    ax.tick_params(axis='x', colors='tab:blue')
    ax.tick_params(axis='y', colors='tab:blue')
    ax.grid(True, ls='--', alpha=0.5)

    # ------------------------------------------------------------------
    # 3️⃣  Add the **secondary** x‑axis *below* the primary one
    # ------------------------------------------------------------------
    #   - `secondary_xaxis('bottom')` creates a new axis at the bottom.
    #   - The two functions (`sec_to_min` and `min_to_sec`) define the
    #     forward and inverse mapping between the two scales.
    #   - The result is just a thin horizontal line with tick marks.
    #
    sec_ax = ax.secondary_xaxis('bottom',
                                functions=(sec_to_min, min_to_sec))

    sec_ax.set_xlabel('Time (minutes)', fontsize=12,
                      color='dimgray', labelpad=10)
    sec_ax.tick_params(axis='x', colors='dimgray')
    sec_ax.xaxis.set_tick_params(which='major', length=6)
    sec_ax.xaxis.set_tick_params(which='minor', length=3)

    # Optional: make the bottom line a little thinner / dimmer
    for spine in sec_ax.spines.values():
        spine.set_color('dimgray')
        spine.set_linewidth(0.8)

    # ------------------------------------------------------------------
    # 4️⃣  Finish and show
    # ------------------------------------------------------------------
    plt.title('Sin wave – secondary X‑axis below (sec → min)', fontsize=14,
              pad=15)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()