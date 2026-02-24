import numpy as np, csv, os
# Create a small synthetic HSI cube
np.random.seed(0)
cube = np.random.rand(10, 10, 15).astype(np.float32)
# Inject a synthetic fungus pattern (bright spots)
cube[3:6, 3:6, :] *= 1.5
# Save as .npy (segmentation script can read .npy)
np.save('test_cube.npy', cube)
# Create reference CSV
wavelengths = np.arange(400, 400 + 15*10, 10)  # 400,410,...,550 nm
reflectance = np.linspace(0.1, 0.3, 15)
with open('test_fungus_ref.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['wavelength', 'fungus_reflection'])
    writer.writerows(zip(wavelengths, reflectance))
print('Test data files written: test_cube.npy, test_fungus_ref.csv')