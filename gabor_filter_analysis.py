import numpy as np
import matplotlib.pyplot as plt
from skimage import io, color
import os,sys
sys.path.append('../physion/src')
from physion.visual_stim.preprocess_NI import img_after_hist_normalization
from skimage.filters import gabor

#%%

img_path = os.path.join(os.path.expanduser('~'), 'work', 'lab-notebook', 'physion', 'src', 'physion', 'visual_stim', 'NI_bank', '4.jpg')


img = color.rgb2gray(io.imread(img_path)).astype(np.float32)
img_norm = img_after_hist_normalization(img).astype(np.float32)
plt.imshow(img_norm, cmap='gray')
plt.title("Input image")
plt.axis('off')
#%%
frequencies = [0.05, 0.1, 0.2, 0.3]
orientations = np.linspace(0, np.pi, 8, endpoint=False)

def compute_gabor_magnitudes(img, frequencies, orientations):
    mags = []

    for theta in orientations:
        for freq in frequencies:
            real, imag = gabor(img, frequency=freq, theta=theta)
            mag = np.sqrt(real**2 + imag**2)
            mags.append(mag)

    return np.stack(mags)

mags = compute_gabor_magnitudes(img_norm, frequencies, orientations)
print("Shape:", mags.shape)

#%%
# strongest response per pixel
mags_max = mags.max(axis=0)

plt.imshow(mags_max, cmap='inferno')
plt.title("Max Gabor response")
plt.colorbar()

#%% 
from scipy.ndimage import binary_dilation

# ----------------------------------------
# Edge tracing from Gabor energy
# ----------------------------------------

def trace_edges(m):

    # robust thresholds
    T_low  = np.percentile(m, 60)
    T_high = np.percentile(m, 90)

    # strong edges
    strong = m >= T_high

    # weak / low-contrast edges
    weak = (m >= T_low) & (m < T_high)

    # connect weak edges to strong ones
    connected = weak & binary_dilation(strong, iterations=2)

    return strong, connected

strong_edges, low_contrast_edges = trace_edges(mags_max)

#%%
fig, ax = plt.subplots(1, 3, figsize=(15,5))

ax[0].imshow(img_norm, cmap='gray')
ax[0].set_title("Input")
ax[0].axis('off')

ax[1].imshow(strong_edges, cmap='gray')
ax[1].set_title("High-contrast edges")
ax[1].axis('off')

ax[2].imshow(low_contrast_edges, cmap='gray')
ax[2].set_title("Low-contrast traced edges")
ax[2].axis('off')

plt.tight_layout()
# 

# 
# #%%
def estimate_noise_floor(m):
    return np.percentile(m, 20)

def strong_threshold(m):
    return np.percentile(m, 80)

def low_contrast_score(m):
    T_low = estimate_noise_floor(m)
    T_high = strong_threshold(m)

    norm = (m - T_low) / (T_high - T_low + 1e-8)
    norm = np.clip(norm, 0, 1)

    weights = norm * (1 - norm)
    return weights.mean()

lc_score = low_contrast_score(mags_max)
print("Low-contrast edge score:", lc_score)

#%%
overlay = np.zeros((*img_norm.shape, 3), dtype=np.float32)

# grayscale background
base = (img_norm - img_norm.min()) / (img_norm.max() - img_norm.min())
overlay[..., 0] = base
overlay[..., 1] = base
overlay[..., 2] = base

# strong edges in red
overlay[strong_edges] = [1, 0, 0]

# low-contrast edges in cyan
overlay[low_contrast_edges] = [0, 1, 1]

plt.figure(figsize=(8,8))
plt.imshow(overlay)
plt.title("Edge tracing overlay")
plt.axis('off')
# 
# #%%
T_low = estimate_noise_floor(mags_max)
T_high = strong_threshold(mags_max)

low_contrast_mask = (mags_max > T_low) & (mags_max < T_high)

plt.imshow(low_contrast_mask, cmap='gray')
plt.title("Low-contrast edge regions")
plt.axis('off')

#%%
def compute_cnr(img):
    # simple approximation:
    # signal = bright regions
    # background = dark regions

    p10 = np.percentile(img, 10)
    p90 = np.percentile(img, 90)

    background = img[img <= p10]
    signal = img[img >= p90]

    mu_signal = signal.mean()
    mu_background = background.mean()
    sigma_noise = img.std()

    return (mu_signal - mu_background) / (sigma_noise + 1e-8)

cnr = compute_cnr(img_norm)
print("CNR:", cnr)

#%%
print(f"""
Summary:
--------
Low-contrast edge score : {lc_score:.4f}
CNR                     : {cnr:.4f}
""")

#%%
def interpret(lc, cnr):
    if cnr > 2 and lc < 0.1:
        return "Clean but flat"
    elif cnr < 1 and lc > 0.2:
        return "Textured but noisy"
    elif cnr > 2 and lc > 0.2:
        return "Rich and sharp"
    else:
        return "Low detail or blurry"

print("Interpretation:", interpret(lc_score, cnr))
# %%
