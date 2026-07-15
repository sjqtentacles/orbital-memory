"""Real physical units for the co-orbital memory cell.

The integrators run in canonical restricted-three-body units — G = 1, total
mass 1, primary separation 1, mean motion n = 1, orbital period 2*pi. In those
units the dynamics depend ONLY on the mass ratio ``mu``; the physical size,
speed, and timescale are pure rescalings applied afterward. A ``System`` holds
the three scale factors implied by a real pair of bodies:

    length_scale = a                      (metres per nondimensional length)
    time_scale   = 1 / n = sqrt(a^3 / (G * M_total))   (seconds per nondim time)
    vel_scale    = a * n = length_scale / time_scale   (m/s per nondim speed)

and converts simulation output into years / AU / m·s⁻¹ / kg. Nothing here
touches the engine, so every canonical test stays exactly as it was; these
factors just let us say the simulated bit librates at 148 real years — the
period actually observed for Jupiter's Trojans.
"""

from dataclasses import dataclass
from math import pi, sqrt

# Physical constants (SI).
G_SI = 6.674e-11              # m^3 kg^-1 s^-2
AU_M = 1.495978707e11         # metres in an astronomical unit
SECONDS_PER_YEAR = 365.25 * 86400.0


@dataclass(frozen=True)
class System:
    """A real two-body pair, carrying the scales that dimensionalize a run.

    ``a_m`` is the primary separation in metres; masses in kg. The secondary
    is the lighter body (the planet or moon whose L4/L5 points hold the bit).
    """

    name: str
    m_primary_kg: float
    m_secondary_kg: float
    a_m: float
    G: float = G_SI

    # --- derived dimensionless / scale quantities -----------------------
    @property
    def M_tot_kg(self) -> float:
        return self.m_primary_kg + self.m_secondary_kg

    @property
    def mu(self) -> float:
        """Secondary mass fraction — the one number the engine actually uses."""
        return self.m_secondary_kg / self.M_tot_kg

    @property
    def length_scale(self) -> float:
        """Metres per nondimensional length unit."""
        return self.a_m

    @property
    def time_scale(self) -> float:
        """Seconds per nondimensional time unit (= 1/n = sqrt(a^3/GM))."""
        return sqrt(self.a_m ** 3 / (self.G * self.M_tot_kg))

    @property
    def vel_scale(self) -> float:
        """m/s per nondimensional speed unit (= a*n)."""
        return self.length_scale / self.time_scale

    @property
    def period_seconds(self) -> float:
        return 2 * pi * self.time_scale

    @property
    def period_years(self) -> float:
        return self.period_seconds / SECONDS_PER_YEAR

    # --- nondimensional -> real -----------------------------------------
    def years(self, t_nd: float) -> float:
        return t_nd * self.time_scale / SECONDS_PER_YEAR

    def au(self, len_nd: float) -> float:
        return len_nd * self.length_scale / AU_M

    def mps(self, v_nd: float) -> float:
        return v_nd * self.vel_scale

    def kmps(self, v_nd: float) -> float:
        return self.mps(v_nd) / 1000.0

    def kg(self, m_nd: float) -> float:
        return m_nd * self.M_tot_kg

    # --- real -> nondimensional -----------------------------------------
    def nd_time(self, t_years: float) -> float:
        return t_years * SECONDS_PER_YEAR / self.time_scale

    def nd_length(self, len_au: float) -> float:
        return len_au * AU_M / self.length_scale

    def nd_vel(self, v_mps: float) -> float:
        return v_mps / self.vel_scale


# --- real systems (masses in kg, separations in metres) -----------------
# Sun-Jupiter: the classic Trojan swarm — the validation anchor.
SUN_JUPITER = System("Sun-Jupiter", 1.989e30, 1.898e27, 5.2044 * AU_M)

# Sun-Earth: home of 2010 TK7, the (large-amplitude) Earth Trojan.
SUN_EARTH = System("Sun-Earth", 1.989e30, 5.972e24, 1.0 * AU_M)

# Saturn's co-orbital moons: Telesto/Calypso ride Tethys, Helene/Polydeuces
# ride Dione — real hardware that has held its bit for the age of the system.
SATURN_TETHYS = System("Saturn-Tethys", 5.6834e26, 6.174e20, 2.947e8)
SATURN_DIONE = System("Saturn-Dione", 5.6834e26, 1.0955e21, 3.774e8)

SYSTEMS = {s.name: s for s in (SUN_JUPITER, SUN_EARTH, SATURN_TETHYS, SATURN_DIONE)}
