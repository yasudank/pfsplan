# Moon Brightness Model
# Based on Krisciunas & Schaefer (1991) PASP 103, 1033

import numpy as np

class MoonBrightnessModel:
    def __init__(self):
        self.k = {
            "g": 0.15, "r": 0.10, "i": 0.09, "z": 0.08, "y": 0.1
        }
        self.Q = {
            "g": 0.129, "r": 0.152, "i": 0.114, "z": 0.048, "y": 0.038
        }
        self.mu_sky = {
            "g": 22.25, "r": 21.18, "i": 20.32, "z": 19.59, "y": 18.28
        }
        self.lam_eff = {
            "g": 478.0, "r": 617.0, "i": 766.0, "z": 888.0, "y": 974.0
        }
        self.Msun = {
            "g": -26.520, "r": -26.922, "i": -27.042, "z": -27.054, "y": -27.059, "V": -26.756
        }

    # Airmass for zenith distance z (deg)
    def X(self, z):
        z = np.radians(z)
        return 1.0 / np.sqrt(1.0 - 0.96 * np.sin(z)**2)

    # Rayleigh scattering
    def tR(self, band, X):
        p = 608.0  # Pressure at Mauna Kea (hPa)
        H = 4.2  # Height of Mauna Kea (km)
        lam = self.lam_eff[band] * 1.0E-03
        tauR = p / 1013.25 * (0.00864 + 6.5E-06 * H) * lam**(-(3.916 + 0.074 * lam + 0.050 / lam))
        return np.exp(-tauR * X)

    # Mie scattering
    def tM(self, band, X):
        lam = self.lam_eff[band] * 1.0E-03
        alpha = -1.38
        kM = np.where(lam < 0.4, 0.050, 0.013 * lam**alpha)
        return 10.0**(-0.4 * kM * X)

    def Bmoon(self, band, alpha, z_moon, X_sky, rho):
        # band: filter band
        # alpha: lunar phase angle (deg)
        # z_moon: lunar zenith distance (deg)
        # X_sky: airmass for the sky
        # rho: angular separation between the Moon and the field (deg)
        XV = self.Msun[band] - self.Msun["V"]
        phi = 180 - alpha
        Istar = 10.0**(-0.4 * (3.84 + 0.026 * np.abs(phi) + 4.0E-09 * phi**4)) * 10.0**(-0.4 * XV)
        X_moon = self.X(z_moon)
        rho = np.radians(rho)
        fR = 10.0**0.92 * (1.06 + np.cos(rho)**2)
        fM = 10.0**(2.44 - np.degrees(rho) / 40.0)
        BmoonR = fR * Istar * 10.0**(-0.4 * self.k[band] * X_moon) * (1.0 - self.tR(band, X_sky))
        BmoonM = fM * Istar * 10.0**(-0.4 * self.k[band] * X_moon) * (1.0 - self.tM(band, X_sky))
        return BmoonR + BmoonM

    def deltaMag(self, band, alpha, z_moon, z_sky, rho):
        if np.any(z_moon < 90.0):
            X_sky = self.X(z_sky)
            B0 = self.Q[band] * 5.48E+06 * 10.0**(-0.4 * self.mu_sky[band]) * X_sky
            Bm = self.Bmoon(band, alpha, z_moon, X_sky, rho)
            return -2.5 * np.log10((Bm + B0) / B0)
        else:
            return np.zeros_like(z_moon)

if __name__ == "__main__":
    # Example usage
    moon = MoonBrightnessModel()
    # alpha = 60 deg, z_moon = 60 deg, z_sky = 30 deg, rho = 30 deg
    alpha = np.array([60])
    z_moon = np.array([60])
    z_sky = np.array([30])
    rho = np.array([30])
    print(moon.deltaMag("g", alpha, z_moon, z_sky, rho))
