import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import factorial
import matplotlib as mp

# CONFIG
plt.style.use("default")
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

# Function to calculate Zernike radial polynomial
def zernike_radial(m, n, rho):
    if (n - m) % 2 != 0:
        return np.zeros_like(rho)
    
    radial = np.zeros_like(rho)
    for k in range((n - m) // 2 + 1):
        radial += rho**(n - 2 * k) * ((-1)**k * factorial(n - k) /
                                     (factorial(k) * factorial((n + m) // 2 - k) * factorial((n - m) // 2 - k)))
    return radial

# Function to calculate Zernike polynomial (both radial and angular parts)
def zernike(m, n, rho, theta):
    if m >= 0:
        return zernike_radial(m, n, rho) * np.cos(m * theta)
    else:
        return zernike_radial(-m, n, rho) * np.sin(-m * theta)

# Zernike indices (m, n) for the first 10 polynomials
zernike_indices = [
    (0, 0),  # Z0: Piston
    (1, 1),  # Z1: Tilt X
    (-1, 1), # Z2: Tilt Y
    (0, 2),  # Z3: Defocus
    (-2, 2), # Z4: Astigmatism 45 degrees
    (2, 2),  # Z5: Astigmatism 0 degrees
    (1, 3),  # Z6: Coma X
    (-1, 3), # Z7: Coma Y
    (3, 3),  # Z8: Primary Spherical Aberration
    (-3, 3), # Z9: Secondary Astigmatism 45 degrees
]

titles = [r"$Z _ {0} ^{0}$", 
          r"$Z _ {1} ^{1}$",
          r"$Z _ {1} ^{-1}$", 
          r"$Z _ {2} ^{0}$",
          r"$Z _ {2} ^{-2}$",
          r"$Z _ {2} ^{2}$",
          r"$Z _ {3} ^{1}$",
          r"$Z _ {3} ^{-1}$",
          r"$Z _ {3} ^{3}$",
          r"$Z _ {3} ^{-3}$"]

# Generate polar grid
num_points = 300
rho = np.linspace(0, 1, num_points)
theta = np.linspace(0, 2 * np.pi, num_points)
rho, theta = np.meshgrid(rho, theta)

# Convert polar coordinates to Cartesian for plotting
x = rho * np.cos(theta)
y = rho * np.sin(theta)

# Set wavelength (in nanometers)
wavelength = 500  # Example wavelength in nm

# Desired RMS wavefront error (WFE)
rms_wavefront_error = wavelength / 10  # Example RMS error

# Coefficients for the first 10 Zernike polynomials (scaled to simulate desired RMS WFE)
# These coefficients are chosen to ensure the RMS WFE is λ/10
coefficients = np.array([
    0.00,  # Adjusted for illustrative purposes
    0.03,
    -0.03,
    0.01,
    -0.01,
    0.02,
    -0.01,
    0.02,
    0.02,
    0.01
])

# Compute the wavefront as a sum of the first 10 Zernike polynomials
wavefront = np.zeros_like(rho)
for i, (m, n) in enumerate(zernike_indices):
    wavefront += coefficients[i] * zernike(m, n, rho, theta)

# Normalize wavefront data to WFE/λ
wavefront_normalized = wavefront #/ rms_wavefront_error

# Create figure and subplots
fig = plt.figure(figsize=(15, 5))

# 3D plot with modified background color
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(x, y, wavefront_normalized, cmap='plasma', edgecolor='none')
ax1.set_title('3D Projection of Wavefront')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_zticklabels([])

# Set 3D plot background color
ax1.w_xaxis.pane.fill = True
ax1.w_xaxis.pane.set_facecolor('ghostwhite')
ax1.w_yaxis.pane.fill = True
ax1.w_yaxis.pane.set_facecolor('ghostwhite')
ax1.w_zaxis.pane.fill = True
ax1.w_zaxis.pane.set_facecolor('ghostwhite')

# Polar plot with only angular labels
ax2 = fig.add_subplot(122, projection='polar')
c = ax2.pcolormesh(theta, rho, wavefront_normalized, shading='auto', cmap='plasma')
ax2.set_title('Polar cut of WFE')

# Remove radial labels but keep angular labels
ax2.set_yticklabels([])  # Remove radial tick labels
ax2.set_yticks([])       # Remove radial ticks
ax2.set_xticks(np.linspace(0, 2 * np.pi, 12, endpoint=False))  # Add angular ticks for better visibility
# Add colorbar for reference
cbar = fig.colorbar(c, ax=ax2,  fraction=0.046, pad=0.04,label='WFE [Arb. units]')

plt.tight_layout()
#plt.show()
plt.savefig("Wavefront_reconstruction.pdf", bbox_inches = "tight")
#plt.show()
fig, axs = plt.subplots(2, 5, figsize=(10.5, 5), subplot_kw={'projection': '3d'})
axes = axs.flatten()

row = 0
col = 0


cmaps = ["magma", "hot", "hot", "viridis", "viridis","plasma","plasma","plasma","plasma","plasma" ]

for i in range(10):
    (m, n) = zernike_indices[i]

    zz = coefficients[i] * zernike(m, n, rho, theta)
    axes[i].set_title(titles[i])
    axes[i].plot_surface(x, y, zz, cmap="plasma", edgecolor='none')
    axes[i].set_axis_off()

plt.tight_layout()
plt.show()
plt.savefig("zernikes.pdf", bbox_inches = "tight")