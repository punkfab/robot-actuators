"""
First-order sizing of a roller-screw linear actuator in the layout Tesla showed for
Optimus: a frameless hollow brushless motor whose rotor carries the nut of a planetary
roller screw, the screw as the output rod, a thrust bearing on the rotor, a position
sensor on the rotor and a force sensor in the load path. Sized here for the 8 kN class
driving a knee through a short crank.

Nothing here is Tesla data beyond the layout and the reported force classes
(500 / 3,900 / 8,000 N, secondary sources, see README). The requirements come from a
knee duty cycle on a 57 kg robot pushed through the linkage, not from a datasheet.

What it models
  screw    planetary roller screw kinematics (d_s = (n - 2) d_r, rollers single start,
           screw and nut n starts), how many rollers fit, Hertz contact at the
           screw/roller and roller/nut flanks, efficiency forward and back, dynamic load
           rating from a power-law fit to catalogue values, L10 life over the duty cycle
  motor    inner-rotor surface-magnet machine around the nut: torque from shear stress,
           copper loss, winding temperature against a lumped resistance to the limb,
           voltage limit at 48 V, rotor inertia
  linkage  actuator on the thigh, rod on a crank on the shank: stroke and moment arm
           against knee angle
  impact   the rotor's inertia reflected to the rod, which a short lead makes very large

    python rollerscrew/design.py          (or: make rollerscrew)
"""

import json
from dataclasses import asdict, dataclass, field, replace
from math import asin, atan, cos, degrees, floor, hypot, pi, radians, sin, sqrt, tan, tau
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"

# ---- materials, SI ------------------------------------------------------------
E_STEEL, NU_STEEL = 206e9, 0.3
RHO_STEEL, RHO_FE, RHO_CU, RHO_MAG, RHO_AL = 7850, 7650, 8900, 7500, 2700
RHO_E_CU20, ALPHA_CU = 1.72e-8, 0.00393
C_CU, C_FE = 385, 460                       # J/kg K
P_STATIC = 4.0e9                            # contact pressure for permanent set ~1e-4 of roller dia (ISO 76 uses 4.2 GPa for ball bearings)

# Catalogue dynamic load ratings of planetary roller screws, (d mm, lead mm, Ca kN).
# Quoted by humanoid.guide (?p=20588) from the Ewellix catalogue; not checked against it.
CATALOGUE = [(15, 5, 26.0), (15, 8, 27.4), (21, 5, 50.6), (21, 8, 57.0), (21, 10, 59.2)]


def fit_ca():
    """Ca = c d^a lead^b, least squares in logs over CATALOGUE."""
    A = np.array([[1, np.log(d), np.log(l)] for d, l, _ in CATALOGUE])
    y = np.log([c for *_, c in CATALOGUE])
    (lc, a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(np.exp(lc)), float(a), float(b)


CA_FIT = fit_ca()


# ---- planetary roller screw ----------------------------------------------------
@dataclass
class Screw:
    d_s: float = 14.0          # screw pitch diameter, mm
    lead: float = 4.0          # mm per nut revolution
    starts: int = 5            # screw and nut starts; rollers single start
    flank: float = 45.0        # flank angle to the radial plane, deg
    nut_len_k: float = 1.6     # threaded nut length / d_s (assumed proportion)
    mu: float = 0.019          # effective friction coefficient in the efficiency formula (gives ~0.85 at 15 x 5)
    k_share: float = 1.3       # peak/mean thread load (elastic load sharing along the nut, assumed)
    drag: float = 0.05         # N m: nut seal, wiper and bearing preload drag (assumed)

    @property
    def pitch(self): return self.lead / self.starts
    @property
    def d_r(self): return self.d_s / (self.starts - 2)        # same helix angle as the nut thread
    @property
    def d_n(self): return self.d_s + 2 * self.d_r
    @property
    def depth(self): return 0.6 * self.pitch                  # thread depth, each member
    @property
    def nut_len(self): return self.nut_len_k * self.d_s
    @property
    def roller_env(self): return self.d_r + 2 * self.depth    # roller gear tip dia ~ thread major
    @property
    def nut_od(self): return self.d_n + 2 * self.depth + 2 * 2.5    # 2.5 mm wall carries the ring gear

    @property
    def n_rollers(self):
        orbit = (self.d_s + self.d_r) / 2
        return min(12, floor(pi / asin(min(1.0, (self.roller_env + 0.4) / (2 * orbit)))))

    @property
    def contacts(self):
        """threads engaged per roller, per side"""
        return floor(self.nut_len / self.pitch)

    @property
    def helix(self): return atan(self.lead / (pi * self.d_s))

    def eta(self):
        lam, phi = self.helix, atan(self.mu)
        fwd = tan(lam) / tan(lam + phi)
        back = tan(lam - phi) / tan(lam) if lam > phi else 0.0
        return fwd, back

    def ca(self):
        c, a, b = CA_FIT
        return c * self.d_s ** a * self.lead ** b * 1e3        # N

    @classmethod
    def for_lead(cls, d_s, lead, min_pitch=0.5):
        """5 starts if the pitch stays >= min_pitch, else 4 (rollers d_s/2), else None"""
        for n in (5, 4):
            if lead / n >= min_pitch - 1e-9:
                return cls(d_s=d_s, lead=lead, starts=n)
        return None

    def roller_phase(self):
        """Axial thread phase each roller must carry relative to its gear teeth, as a
        fraction of the pitch: rollers are index-matched to their positions."""
        return [(self.starts * k / self.n_rollers) % 1 for k in range(self.n_rollers)]

    # Hertz point contact, Hamrock-Brewe approximations
    def _hertz(self, q, rx_curv, ry_curv):
        e_eff = E_STEEL / (1 - NU_STEEL ** 2)                  # E' for two steel bodies
        rx, ry = 1 / rx_curv, 1 / ry_curv
        if rx > ry:
            rx, ry = ry, rx
        alpha = ry / rx
        k = 1.0339 * alpha ** 0.636
        ell = 1.0003 + 0.5968 / alpha
        r = 1 / (1 / rx + 1 / ry)
        a = (6 * k * k * ell * q * r / (pi * e_eff)) ** (1 / 3)
        b = (6 * ell * q * r / (pi * k * e_eff)) ** (1 / 3)
        return 3 * q / (2 * pi * a * b), a, b

    def contact(self, force):
        """Peak normal load per flank and peak pressure at the screw and nut contacts."""
        al = radians(self.flank)
        q = force / (self.n_rollers * self.contacts * cos(al) * cos(self.helix)) * self.k_share
        rs, rr, rn = (x / 2 * 1e-3 for x in (self.d_s, self.d_r, self.d_n))
        prof = sin(al) / rr                                    # roller flank: sphere centred on the roller axis
        p_s, *_ = self._hertz(q, sin(al) / rs + sin(al) / rr, prof)
        p_n, *_ = self._hertz(q, sin(al) / rr - sin(al) / rn, prof)
        return dict(q=q, p_screw=p_s, p_nut=p_n)

    def mass(self, stroke):
        """screw (threaded over stroke + nut), rollers, nut; kg"""
        l_screw = stroke + self.nut_len + 20
        m_s = RHO_STEEL * pi / 4 * (self.d_s * 1e-3) ** 2 * l_screw * 1e-3
        m_r = self.n_rollers * RHO_STEEL * pi / 4 * (self.d_r * 1e-3) ** 2 * (self.nut_len + 4) * 1e-3
        m_n = RHO_STEEL * pi / 4 * ((self.nut_od * 1e-3) ** 2 - (self.d_n * 1e-3) ** 2) * (self.nut_len + 8) * 1e-3
        return m_s, m_r, m_n


# ---- motor around the nut --------------------------------------------------------
@dataclass
class Motor:
    bore: float                 # rotor inner diameter = nut OD, mm
    stack: float = 30.0         # mm
    slots: int = 12
    poles: int = 14
    kw: float = 0.933           # 12s/14p concentrated winding
    rotor_yoke: float = 2.5     # mm
    magnet: float = 2.5
    gap: float = 0.5
    tooth_tip: float = 1.0
    slot_depth: float = 8.0
    br: float = 1.3             # N48H-class, T
    k_flux: float = 0.80        # leakage and slotting on the fundamental (assumed)
    fill: float = 0.45
    j_peak: float = 20e6        # A/m2 rms: peak torque, a few seconds
    t_hot: float = 140.0        # winding hot-spot limit, C (class H with margin)
    t_amb: float = 40.0
    h_eff: float = 40.0         # W/m2K over the housing: conduction into the limb structure + air (assumed)
    r_int: float = 0.4          # K/W winding to housing (assumed)
    v_bus: float = 48.0
    turns: int = 0              # per coil; 0 = choose in size()

    @property
    def r_gap(self): return self.bore / 2 + self.rotor_yoke + self.magnet          # magnet outer, mm
    @property
    def r_si(self): return self.r_gap + self.gap
    @property
    def b1(self):
        bg = self.br * self.magnet / (self.magnet + 1.05 * self.gap)
        return 4 / pi * bg * sin(0.85 * pi / 2) * self.k_flux                      # fundamental peak, T

    @property
    def tooth_w(self):
        """sized so the tooth carries the gap flux of one slot pitch at 1.6 T"""
        return 2 * pi * self.r_si / self.slots * self.b1 * 2 / pi / 1.6

    @property
    def stator_yoke(self):
        phi_half = self.b1 * 2 / pi * (2 * pi * self.r_si / self.poles) / 2       # T*mm per mm stack
        return max(2.0, phi_half / 1.5)

    @property
    def r_so(self): return self.r_si + self.tooth_tip + self.slot_depth + self.stator_yoke
    @property
    def od(self): return 2 * (self.r_so + 2.0)                                    # 2 mm housing wall

    @property
    def slot_area(self):
        ra, rb = self.r_si + self.tooth_tip, self.r_si + self.tooth_tip + self.slot_depth
        return (pi * (rb * rb - ra * ra) - self.slots * self.tooth_w * self.slot_depth) / self.slots * 1e-6   # m2

    @property
    def coil_len(self):
        """mean turn length, m: two stack runs plus two end arcs round the tooth"""
        return (2 * self.stack + pi * (self.tooth_w + self.slot_depth / 2)) * 1e-3

    sigma_max: float = 40e3     # Pa: gap shear cap standing in for saturation (surface-magnet servo motors peak ~30-50 kPa)

    def torque(self, j):
        lin = self.r_gap * 1e-3 * self.stack * 1e-3 * self.kw * self.b1 * self.slots * self.slot_area * self.fill * j / sqrt(2)
        return min(lin, tau * (self.r_gap * 1e-3) ** 2 * self.stack * 1e-3 * self.sigma_max)

    def torque_lin(self, j):
        return self.r_gap * 1e-3 * self.stack * 1e-3 * self.kw * self.b1 * self.slots * self.slot_area * self.fill * j / sqrt(2)

    def p_cu(self, j, temp=None):
        rho = RHO_E_CU20 * (1 + ALPHA_CU * ((temp or self.t_hot) - 20))
        return rho * j * j * self.slots * self.slot_area * self.fill / 2 * self.coil_len    # each coil fills half of two slots

    @property
    def r_th(self):
        length = self.stack + 2 * 6 + 30                     # winding overhang, caps, bearings
        area = pi * self.od * length * 1e-6 + 2 * pi / 4 * (self.od * 1e-3) ** 2
        return self.r_int + 1 / (self.h_eff * area)

    @property
    def j_cont(self):
        return sqrt((self.t_hot - self.t_amb) / self.r_th / self.p_cu(1.0))

    def heat_capacity(self):
        """copper plus the teeth, the mass the hot spot heats first; J/K"""
        m_cu = RHO_CU * self.slots * self.slot_area * self.fill * self.coil_len / 2
        m_t = RHO_FE * self.slots * self.tooth_w * self.slot_depth * self.stack * 1e-9
        return m_cu * C_CU + m_t * C_FE

    def phase(self, turns=None):
        n = turns or self.turns
        a_cond = self.slot_area * self.fill / (2 * n)
        r = RHO_E_CU20 * (1 + ALPHA_CU * (self.t_hot - 20)) * n * self.coil_len * (self.slots // 3) / a_cond
        i_per_j = self.slot_area * self.fill / (2 * n)         # phase rms A per (A/m2)
        kt = self.torque_lin(1.0) / i_per_j                        # N m per A rms
        return dict(r=r, kt=kt, i_per_j=i_per_j)

    def speed(self, j, turns=None):
        """rotor rad/s reachable at current density j (resistive drop only, no field weakening)"""
        ph = self.phase(turns)
        v_ph = 0.95 * self.v_bus / sqrt(6)                     # phase rms, SVPWM linear range
        return max(0.0, 3 * (v_ph - ph["r"] * ph["i_per_j"] * j) / ph["kt"])

    def masses(self):
        r = lambda a, b: pi * ((b * 1e-3) ** 2 - (a * 1e-3) ** 2) * self.stack * 1e-3
        m_ry = RHO_FE * r(self.bore / 2, self.bore / 2 + self.rotor_yoke)
        m_mag = RHO_MAG * r(self.r_gap - self.magnet, self.r_gap) * 0.85
        m_st = RHO_FE * (r(self.r_si, self.r_so) - self.slots * self.slot_area * self.stack * 1e-3)
        m_cu = RHO_CU * self.slots * self.slot_area * self.fill * self.coil_len / 2
        return dict(rotor_yoke=m_ry, magnets=m_mag, stator=m_st, copper=m_cu)

    def rotor_inertia(self, nut_mass, nut_r):
        """rotor yoke, magnets and the nut, kg m2 (hub and rollers added by the caller)"""
        r = lambda a, b: (a * a + b * b) / 2 * 1e-6
        m = self.masses()
        return (m["rotor_yoke"] * r(self.bore / 2, self.bore / 2 + self.rotor_yoke)
                + m["magnets"] * r(self.r_gap - self.magnet, self.r_gap)
                + nut_mass * r(nut_r[0], nut_r[1]))


# ---- knee linkage -----------------------------------------------------------------
@dataclass
class Linkage:
    """Knee at the origin, thigh up +y. Actuator body pinned on the thigh at `base`;
    rod pinned on a crank of radius `crank` on the shank, at `crank_ang` from the
    shank's axis (shank hangs along -y at knee angle 0, flexion rotates it +). mm."""
    base: tuple = (-40.0, 180.0)   # defaults are fit_linkage()'s result
    crank: float = 35.0
    crank_ang: float = -20.0       # deg from the shank axis
    flex_max: float = 125.0

    def pin(self, flex):
        a = radians(-90 + self.crank_ang - flex)            # shank axis at -90 deg, flexion turns it backward (-)
        return self.crank * cos(a), self.crank * sin(a)

    def length(self, flex):
        x, y = self.pin(flex)
        return hypot(x - self.base[0], y - self.base[1])

    def arm(self, flex, h=0.1):
        """moment arm, mm: d length / d flexion (rad); positive = rod extends as the knee bends"""
        return (self.length(flex + h) - self.length(flex - h)) / radians(2 * h)

    def stroke(self):
        ls = [self.length(f) for f in np.linspace(0, self.flex_max, 200)]
        return max(ls) - min(ls), min(ls), max(ls)


def fit_linkage(max_stroke=65.0, min_len=0.0):
    """Grid search: the largest smallest moment arm over 0..flex_max, the arm keeping one
    sign (the rod only pushes to extend the knee), stroke under max_stroke."""
    best = None
    for bx in np.arange(-60, 61, 10):
        for by in (180, 220, 260, 300, 340):
            for c in (25, 30, 35, 40, 45):
                for ca in np.arange(-180, 180, 10):
                    lk = Linkage(base=(float(bx), float(by)), crank=float(c), crank_ang=float(ca))
                    arms = np.array([lk.arm(f) for f in np.linspace(0, lk.flex_max, 40)])
                    st, lo, _ = lk.stroke()
                    if not np.all(arms < 0) or st > max_stroke or lo < min_len:
                        continue
                    if best is None or np.abs(arms).min() > best[0]:
                        best = (np.abs(arms).min(), lk)
    return best[1]



# ---- duty: the knee over a day -------------------------------------------------------
ROBOT_MASS = 57.0          # kg, the planned Optimus mass (Wikipedia, citing Motor Authority 2021)
PAYLOAD = 20.0
STRIDE_S, STRIDE_M = 1.1, 1.3
# One day's work, per knee (assumed): 10 km walking, 200 squats to 100 deg holding the
# payload, 30 floors up and down (16 steps a floor, one knee leads every other step).
DAY = dict(walk=10e3 / STRIDE_M, squat=200, stairs_up=240, stairs_down=240)


def _pchip(keys, vals, s):
    from scipy.interpolate import PchipInterpolator
    return PchipInterpolator(keys, vals)(s)


def gait(n=400):
    """Knee flexion (deg) and extensor moment (N m) over one stride. Synthetic: shaped
    like normal human walking (loading-response flexion ~18 deg with an extensor peak of
    0.5 N m/kg, swing flexion ~60 deg); keypoints by hand, not digitised from a data set."""
    s = np.linspace(0, 1, n)
    pk = [0, .06, .15, .30, .40, .50, .62, .72, .85, .95, 1.0]
    ang = [3, 10, 18, 9, 5, 12, 40, 60, 30, 5, 3]
    mom = [0.0, 0.25, 0.50, 0.15, -0.15, -0.10, -0.05, -0.02, -0.12, -0.15, 0.0]
    return s * STRIDE_S, _pchip(pk, ang, s), _pchip(pk, mom, s) * ROBOT_MASS


def squat(n=300, period=3.0, depth=100.0):
    """Down to `depth` and up again holding the payload. Knee moment taken as the static
    one for half the weight on a 0.12 m lever at 90 deg, scaled by sin(flexion), x1.5 on
    the way up (assumed)."""
    s = np.linspace(0, 1, n)
    flex = depth * np.sin(pi * s) ** 2
    m90 = (ROBOT_MASS + PAYLOAD) * 9.81 * 0.5 * 0.12
    mom = m90 * np.sin(np.radians(flex)) * np.where(s > 0.5, 1.5, 1.0)
    return s * period, flex, mom


def stairs(up, n=300):
    """One stance on a stair: up, the knee extends from 60 deg under a 1.0 N m/kg extensor
    peak; down, it flexes to 85 deg against one of 1.1 N m/kg. Synthetic, by hand."""
    s = np.linspace(0, 1, n)
    if up:
        flex, mom = _pchip([0, .3, .7, 1], [60, 45, 12, 5], s), _pchip([0, .25, .6, 1], [0.3, 1.0, 0.4, 0.0], s)
    else:
        flex, mom = _pchip([0, .3, .7, 1], [5, 20, 70, 85], s), _pchip([0, .3, .7, 1], [0.0, 0.5, 1.1, 0.2], s)
    return s * 1.2, flex, mom * ROBOT_MASS


def rod_history(link, t, flex, moment):
    """rod force (N, + = pushing) and speed (m/s) from knee angle and moment. With
    arm = dl/dθ negative here (the rod shortens as the knee bends), F = -M / arm."""
    arm = np.array([link.arm(f) for f in flex]) * 1e-3
    length = np.array([link.length(f) for f in flex]) * 1e-3
    return -moment / arm, np.gradient(length, t), length


def duties(link):
    """each duty's rod force and speed history, and how many a day"""
    out = {}
    for name, (t, fl, mo) in dict(walk=gait(), squat=squat(), stairs_up=stairs(True), stairs_down=stairs(False)).items():
        f, v, _ = rod_history(link, t, fl, mo)
        out[name] = dict(t=t, flex=fl, moment=mo, force=f, vel=v, count=DAY[name])
    return out


# ---- a whole actuator --------------------------------------------------------------
B_W = 7.0                  # bearing width, mm (thin section, stand-in)
B_SECTION = 6.0            # bearing radial section, mm
CAP = 4.0
ENCODER = 6.0              # rotor magnet + sensor board
CLEVIS = 12.0              # rear pin to the load cell
LOAD_CELL = 8.0
ROD_OUT = 6.0              # rod past the front cap, fully retracted
EYE = 10.0                 # rod end face to the rod-eye pin
NUT_END = 4.0              # each end of the nut: ring gear 2.5 + carrier ring 1.5


@dataclass
class Actuator:
    screw: Screw = field(default_factory=Screw)
    stack: float = 30.0
    slot_depth: float = 8.0
    turns: int = 0
    link: Linkage = field(default_factory=Linkage)
    spacer: float = 0.0          # aluminium sleeve between the nut and the rotor yoke, mm radial: a larger air gap

    @property
    def motor(self):
        return Motor(bore=self.screw.nut_od + 2 * self.spacer, stack=self.stack, slot_depth=self.slot_depth, turns=self.turns)

    def spacer_mass(self):
        r0 = self.screw.nut_od / 2
        return RHO_AL * pi * ((r0 + self.spacer) ** 2 - r0 ** 2) * self.stack * 1e-9

    @property
    def overhang(self):
        return self.slot_depth / 2 + 2.0                 # end winding past the stack

    def stations(self):
        """axial layout, mm from the rear pin, rod retracted; the same numbers drive cad.py.
        The nut sits at the front of the rotor, between the thrust bearings. Retracting, the
        screw's tail runs back through the rotor bore and, if the bore is too short, on into
        a tail tube behind the rear cap. Extending, the threads leave through a wiper in the
        front cap under a boot."""
        stroke = self.link.stroke()[0]
        b = {}                                           # body, from the rear cap's outer face
        b["encoder"] = CAP
        b["bearing_rear"] = b["encoder"] + ENCODER
        b["stack"] = b["bearing_rear"] + B_W + self.overhang
        b["bearing_front"] = b["stack"] + self.stack + self.overhang
        b["cap_front"] = b["bearing_front"] + 2 * B_W
        b["front"] = b["cap_front"] + CAP
        b["nut"] = b["cap_front"] - self.screw.nut_len - 2 * NUT_END     # whole nut incl. ring gears and carriers
        z = {"clevis": 0.0, "load_cell": CLEVIS, "tail": CLEVIS + LOAD_CELL}
        z["tail_len"] = max(0.0, stroke + 2.0 - (b["nut"] - b["encoder"]))
        z["cap_rear"] = z["tail"] + z["tail_len"]
        z.update({k: z["cap_rear"] + v for k, v in b.items()})
        z["nut_thread"] = z["nut"] + NUT_END
        z["screw_rear"] = z["nut_thread"] - stroke - 1.0
        z["rod_end"] = z["front"] + ROD_OUT
        z["pin_min"] = z["rod_end"] + EYE
        z["stroke"] = stroke
        return z

    def force(self, j):
        """rod force at current density j, driving (motor torque through η forward), minus drag"""
        eta, _ = self.screw.eta()
        return max(0.0, (self.motor.torque(j) - self.screw.drag) * tau / (self.screw.lead * 1e-3) * eta)

    def rod_speed(self, j):
        return self.motor.speed(j) * self.screw.lead * 1e-3 / tau

    def j_for(self, force):
        eta, _ = self.screw.eta()
        t = abs(force) * self.screw.lead * 1e-3 / tau / eta + self.screw.drag
        return t / self.motor.torque_lin(1.0)

    def mass(self):
        m = self.motor.masses()
        z = self.stations()
        stroke = self.link.stroke()[0]
        ms, mr, mn = self.screw.mass(stroke)
        l_screw = max(self.link.stroke()[1] - EYE, z["rod_end"]) - z["screw_rear"]      # the rod reaches the linkage's pin
        ms = RHO_STEEL * pi / 4 * (self.screw.d_s * 1e-3) ** 2 * l_screw * 1e-3
        od = self.motor.od
        body = z["front"] - z["cap_rear"]
        housing = RHO_AL * (pi * od * 2.0 * body + 2 * pi / 4 * od * od * CAP) * 1e-9
        b_out = self.screw.nut_od / 2 + 2 + B_SECTION                  # bearing OD radius (hub + section)
        webs = RHO_AL * pi * (self.motor.r_so ** 2 - b_out ** 2) * 3 * B_W * 1e-9    # bearing seats, found by cad.py's mass check
        housing += webs
        tail = RHO_AL * pi * (self.screw.d_s + 4) * 1.5 * z["tail_len"] * 1e-9
        bo = self.screw.nut_od + 4
        bearings = 3 * RHO_STEEL * pi / 4 * ((bo + 2 * B_SECTION) ** 2 - bo ** 2) * 1e-6 * B_W * 1e-3 * 0.6
        rotor_hub = RHO_STEEL * pi / 4 * ((self.screw.nut_od + 4) ** 2 - self.screw.nut_od ** 2) * 1e-6 * (body - self.stack) * 1e-3
        extras = 0.10              # encoder, load cell, rod eye, clevis, seals, fasteners, leads (assumed)
        parts = dict(**m, screw=ms, rollers=mr, nut=mn, rotor_hub=rotor_hub + self.spacer_mass(), housing=housing + tail,
                     bearings=bearings, extras=extras)
        return sum(parts.values()), parts

    def reflected_mass(self):
        """rotor + nut + rollers seen at the rod, kg. Screw does not turn; the roller carrier
        turns at d_n/(d_s+d_n) of the nut speed."""
        s = self.screw
        _, mr, mn = s.mass(self.link.stroke()[0])
        j = self.motor.rotor_inertia(mn, (s.d_n / 2, s.nut_od / 2))
        j += mr * (((s.d_s + s.d_r) / 2) * 1e-3) ** 2 * (s.d_n / (s.d_s + s.d_n)) ** 2
        r_hub = (s.nut_od / 2, s.nut_od / 2 + 2)
        _, parts = self.mass()
        j += (parts["rotor_hub"] - self.spacer_mass()) * (r_hub[0] ** 2 + r_hub[1] ** 2) / 2 * 1e-6
        j += self.spacer_mass() * (r_hub[0] ** 2 + (s.nut_od / 2 + self.spacer) ** 2) / 2 * 1e-6
        return j * (tau / (s.lead * 1e-3)) ** 2, j

    def sized(self, demand):
        """largest turns per coil that still reaches every (force, speed) demand at 48 V"""
        best = None
        for n in range(1, 60):
            a = replace(self, turns=n)
            if all(a.rod_speed(a.j_for(f)) >= abs(v) for f, v in demand):
                best = a
        return best


def life(act, ds):
    """L10 of the screw, in days of DAY: cubic mean load weighted by nut revolutions."""
    num = den = 0.0
    for d in ds.values():
        w = np.abs(d["vel"]) / (act.screw.lead * 1e-3) * (d["t"][1] - d["t"][0])
        num += d["count"] * np.sum(np.abs(d["force"]) ** 3 * w)
        den += d["count"] * np.sum(w)
    f_m = (num / den) ** (1 / 3)
    l10 = (act.screw.ca() / f_m) ** 3 * 1e6
    return l10 / den, f_m, den


# ---- study ----------------------------------------------------------------------------
F_PEAK = 8000.0           # the reported class
SHOCK = 1.5               # landing / stumble overload on the screw, x F_PEAK (assumed)
LIFE_DAYS = 3 * 300       # L10 target: three working years of DAY
LEADS = (2.0, 2.5, 3.0, 4.0, 5.0)
IMPACT_V, IMPACT_K = 0.15, 10e6     # rod speed change at a hard landing (m/s), stiffness in series at the rod (N/m); assumed


def landing(a):
    """peak screw force when the rod is suddenly driven at IMPACT_V: the rotor, seen through the
    lead as a mass at the rod, is accelerated through the series stiffness IMPACT_K"""
    return IMPACT_V * sqrt(IMPACT_K * (a.reflected_mass()[0] + a.mass()[1]["screw"]))
DIAMETERS = (10, 12, 14, 16, 18)


def requirements(link):
    ds = duties(link)
    allf = np.concatenate([d["force"] for d in ds.values()])
    w = ds["walk"]
    # thermal: rms over a day's active time, all duties
    t_act = {k: d["t"][-1] * d["count"] for k, d in ds.items()}
    f2 = sum(np.mean(d["force"] ** 2) * t_act[k] for k, d in ds.items()) / sum(t_act.values())
    demand = []
    for d in ds.values():
        demand += list(zip(d["force"], d["vel"]))[::3]
    return dict(ds=ds, demand=demand, f_rms=float(sqrt(f2)), f_rms_walk=float(np.sqrt(np.mean(w["force"] ** 2))),
                f_max=float(np.max(np.abs(allf))), v_max=float(max(np.max(np.abs(d["vel"])) for d in ds.values())),
                f_squat=float(np.max(np.abs(ds["squat"]["force"]))))


def sweep(req, link):
    rows = []
    for lead in LEADS:
        for d_s in DIAMETERS:
            s = Screw.for_lead(d_s, lead)
            if s is None:
                continue
            stat = s.contact(F_PEAK * SHOCK)
            if max(stat["p_screw"], stat["p_nut"]) > P_STATIC:
                continue
            for depth, spacer in ((d, sp) for d in (6.0, 8.0, 10.0) for sp in (0.0, 3.0, 6.0, 9.0)):
                for stack in np.arange(10, 80, 2.0):
                    a = Actuator(screw=s, stack=float(stack), slot_depth=depth, link=link, spacer=spacer)
                    if a.force(a.motor.j_peak) < F_PEAK or a.force(a.motor.j_cont) < 1.2 * req["f_rms"]:
                        continue
                    a = a.sized(req["demand"])
                    if a is None:
                        continue
                    days, *_ = life(a, req["ds"])
                    if days < LIFE_DAYS:
                        break
                    mass, _ = a.mass()
                    rows.append(dict(lead=lead, d_s=d_s, depth=depth, spacer=spacer, stack=float(stack), turns=a.turns, mass=mass,
                                     impact=landing(a),
                                     od=a.motor.od, pin_min=a.stations()["pin_min"], days=days,
                                     m_ref=a.reflected_mass()[0], eta=s.eta()[0], eta_back=s.eta()[1], act=a))
                    break               # shortest stack that works for this screw and slot depth
    return rows


def describe(a, req):
    s, m = a.screw, a.motor
    ph = m.phase()
    mass, parts = a.mass()
    m_ref, j_rot = a.reflected_mass()
    stroke, l_min, l_max = a.link.stroke()
    eta_f, eta_b = s.eta()
    st, pk = s.contact(F_PEAK * SHOCK), s.contact(F_PEAK)
    days, f_m, rev_day = life(a, req["ds"])
    p_peak = m.p_cu(m.j_peak)
    return dict(
        screw=dict(d_s=s.d_s, d_r=round(s.d_r, 3), d_n=round(s.d_n, 3), lead=s.lead, pitch=s.pitch, starts=s.starts,
                   rollers=s.n_rollers, nut_len=round(s.nut_len, 1), nut_od=round(s.nut_od, 2), threads_engaged=s.contacts,
                   helix_deg=round(degrees(s.helix), 2), eta_fwd=eta_f, eta_back=eta_b, ca=s.ca(),
                   p_at_peak=max(pk["p_screw"], pk["p_nut"]), p_at_shock=max(st["p_screw"], st["p_nut"]),
                   flank_load_at_shock=st["q"], roller_phase=s.roller_phase()),
        motor=dict(od=m.od, stack=m.stack, slot_depth=m.slot_depth, r_gap=m.r_gap, b1=m.b1, tooth_w=m.tooth_w,
                   stator_yoke=m.stator_yoke, turns=m.turns, kt=ph["kt"], r_phase=ph["r"], slot_area_mm2=m.slot_area * 1e6,
                   j_peak=m.j_peak, j_cont=m.j_cont, t_peak=m.torque(m.j_peak), t_cont=m.torque(m.j_cont),
                   i_peak=ph["i_per_j"] * m.j_peak, i_cont=ph["i_per_j"] * m.j_cont, p_cu_peak=p_peak,
                   p_cu_cont=m.p_cu(m.j_cont), r_th=m.r_th,
                   burst_s_from_cold=m.heat_capacity() * (m.t_hot - m.t_amb) / p_peak,
                   burst_s_from_40K_below=m.heat_capacity() * 40 / p_peak),
        actuator=dict(f_peak=a.force(m.j_peak), f_cont=a.force(m.j_cont), v_noload=a.rod_speed(0.0),
                      v_at_cont=a.rod_speed(m.j_cont), v_at_peak=a.rod_speed(m.j_peak),
                      mass=mass, parts=parts, od=m.od, stroke=stroke, pin_min_needed=a.stations()["pin_min"],
                      pin_min=l_min, pin_max=l_max, reflected_mass=m_ref, rotor_inertia=j_rot,
                      backdrive_force=s.drag * tau / (s.lead * 1e-3) / max(eta_b, 1e-6),
                      impact_force=landing(a)),
        duty=dict(f_rms_day=req["f_rms"], f_rms_walk=req["f_rms_walk"], f_max=req["f_max"], v_max=req["v_max"],
                  f_squat=req["f_squat"], cubic_mean=f_m, nut_rev_per_day=rev_day, l10_days=days, l10_years=days / 300,
                  arm_0=a.link.arm(0), arm_90=a.link.arm(90),
                  arm_min=min(abs(a.link.arm(f)) for f in np.linspace(0, a.link.flex_max, 60))),
    )


def ca_check():
    """Peak contact pressure each catalogue screw would see at its own Ca in this model:
    if the Hertz model and the proportions are sane these come out roughly equal."""
    out = []
    for d, l, c in CATALOGUE:
        s = Screw(d_s=d, lead=l)
        p = s.contact(c * 1e3)
        out.append(dict(d=d, lead=l, ca_kn=c, fit_kn=s.ca() / 1e3, rollers=s.n_rollers,
                        p_at_ca=max(p["p_screw"], p["p_nut"]) / 1e9))
    return out


def plots(best, rows, req):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.8))
    a, m, link = best, best.motor, best.link
    col = dict(walk="#c87a2f", squat="#b03a2e", stairs_up="#6a51a3", stairs_down="#2a7f3f")

    flex = np.linspace(0, link.flex_max, 200)
    p = ax[0, 0]
    p.plot(flex, [abs(link.arm(f)) for f in flex], color="#1f4e79", lw=1.6)
    p.set(xlabel="knee flexion (deg)", ylabel="moment arm (mm)", ylim=(0, None),
          title=f"crank linkage: {link.stroke()[0]:.0f} mm stroke over 0–{link.flex_max:.0f}°")
    q = p.twinx()
    q.spines["right"].set_visible(True)
    q.plot(flex, [abs(link.arm(f)) * 1e-3 * a.force(m.j_peak) for f in flex], color="#b03a2e", lw=1.2, label="at peak force")
    q.plot(flex, [abs(link.arm(f)) * 1e-3 * a.force(m.j_cont) for f in flex], "--", color="#b03a2e", lw=1.0, label="at continuous force")
    q.set_ylabel("knee torque (N·m)", color="#b03a2e")
    q.set_ylim(0, None)
    q.legend(loc="lower center", frameon=False, fontsize=8)

    p = ax[0, 1]
    js = np.linspace(0, m.j_peak, 80)
    p.plot([a.rod_speed(j) * 1e3 for j in js], [a.force(j) / 1e3 for j in js], color="#1f4e79", lw=1.6, label="peak")
    jc = np.linspace(0, m.j_cont, 60)
    p.plot([a.rod_speed(j) * 1e3 for j in jc], [a.force(j) / 1e3 for j in jc], color="0.35", lw=1.3, ls="--", label="continuous")
    for k, d in req["ds"].items():
        p.plot(np.abs(d["vel"]) * 1e3, np.abs(d["force"]) / 1e3, ".", ms=2, color=col[k], label=k.replace("_", " "))
    p.set(xlabel="rod speed (mm/s)", ylabel="rod force (kN)", title="force against speed at 48 V, and the duty")
    p.legend(frameon=False, fontsize=8)

    p = ax[1, 0]
    for d in sorted({r["d_s"] for r in rows}):
        rr = {}
        for r in rows:
            if r["d_s"] == d and (r["lead"] not in rr or r["mass"] < rr[r["lead"]]["mass"]):
                rr[r["lead"]] = r
        rr = [rr[k] for k in sorted(rr)]
        p.plot([r["lead"] for r in rr], [r["mass"] for r in rr], "o-", ms=3, lw=1.2, label=f"screw Ø{d}")
    p.plot([best.screw.lead], [best.mass()[0]], "*", ms=12, color="k")
    p.set(xlabel="lead (mm)", ylabel="actuator mass (kg)", title="lightest meeting the duty and the landing, by screw")
    p.legend(frameon=False, fontsize=8)

    p = ax[1, 1]
    leads = np.array([2, 2.5, 3, 4, 5, 6, 8])
    p.plot(leads, [Screw(d_s=best.screw.d_s, lead=l).eta()[0] for l in leads], "--", color="0.4", lw=1.2, label="forward η")
    p.plot(leads, [Screw(d_s=best.screw.d_s, lead=l).eta()[1] for l in leads], color="0.1", lw=1.4, label="back-driving η")
    p.set(xlabel="lead (mm)", ylabel="efficiency", ylim=(0, 1), title=f"lead trade, Ø{best.screw.d_s} screw, chosen motor")
    q = p.twinx()
    q.spines["right"].set_visible(True)
    sub = [replace(best, screw=replace(best.screw, lead=float(l))) for l in leads]
    q.semilogy(leads, [x.reflected_mass()[0] for x in sub], color="#b03a2e", lw=1.4, label="rotor mass seen at the rod")
    q.semilogy(leads, [landing(x) for x in sub], color="#b03a2e", lw=1.0, ls=":",
               label=f"landing force, {IMPACT_V} m/s")
    q.axhline(F_PEAK * SHOCK, color="#b03a2e", lw=0.6, alpha=0.5)
    q.set_ylabel("kg  |  N", color="#b03a2e")
    p.legend(frameon=False, fontsize=8, loc="center right")
    q.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(0.0, 0.42))
    fig.tight_layout()
    fig.savefig(OUT / "design.png", dpi=150)


def main():
    OUT.mkdir(exist_ok=True)
    c, ea, eb = CA_FIT
    print(f"catalogue fit: Ca = {c:.3f} d^{ea:.2f} lead^{eb:.2f} kN; contact pressure each catalogue screw sees at its Ca here:")
    for r in ca_check():
        print(f"  {r['d']}x{r['lead']}: Ca {r['ca_kn']:.1f} (fit {r['fit_kn']:.1f}) kN, {r['rollers']} rollers, {r['p_at_ca']:.2f} GPa")

    # the actuator's own retracted length constrains the linkage: size, refit, size again
    link = fit_linkage()
    for it in range(3):
        req = requirements(link)
        rows = sweep(req, link)
        if not rows:
            raise SystemExit("nothing meets the duty")
        ok = [r for r in rows if r["impact"] <= F_PEAK * SHOCK]
        lightest = min(rows, key=lambda r: r["mass"])
        best = min(ok, key=lambda r: r["mass"])["act"]
        need = max(r["act"].stations()["pin_min"] for r in (lightest, min(ok, key=lambda r: r["mass"])))
        print(f"pass {it}: linkage pin length {link.stroke()[1]:.0f}-{link.stroke()[2]:.0f} mm, actuator needs >= {need:.0f} mm")
        if link.stroke()[1] >= need:
            break
        link = fit_linkage(min_len=need + 5)

    stroke, l_min, l_max = link.stroke()
    print(f"\nlinkage: base {link.base}, crank {link.crank} mm at {link.crank_ang} deg; stroke {stroke:.1f} mm, "
          f"moment arm {abs(link.arm(0)):.1f} mm straight, {abs(link.arm(90)):.1f} at 90 deg")
    print("duty per day: " + ", ".join(f"{k} x{v:.0f}" for k, v in DAY.items()))
    print(f"  rod force rms over the active day {req['f_rms']:.0f} N (walking alone {req['f_rms_walk']:.0f}), "
          f"max {req['f_max']:.0f} N, max rod speed {req['v_max']*1e3:.0f} mm/s")

    print(f"\n{'lead':>5} {'d_s':>4} {'slot':>5} {'spcr':>5} {'stack':>6} {'turns':>5} {'mass kg':>8} {'OD':>5} {'pins':>5} {'L10 yr':>7} "
          f"{'m_ref kg':>9} {'η fwd':>6} {'η back':>7}")
    for r in sorted(rows, key=lambda r: r["mass"])[:12]:
        print(f"{r['lead']:5.1f} {r['d_s']:4d} {r['depth']:5.0f} {r['spacer']:5.0f} {r['stack']:6.0f} {r['turns']:5d} {r['mass']:8.3f} {r['od']:5.1f} "
              f"{r['pin_min']:5.0f} {r['days']/300:7.1f} {r['m_ref']:9.0f} {r['eta']:6.2f} {r['eta_back']:7.2f}")

    print("\nlightest surviving the landing, for each lead (landing: rod speed change %.2f m/s into %.0f N/um at the rod; screw static limit %.0f kN):"
          % (IMPACT_V, IMPACT_K / 1e6, F_PEAK * SHOCK / 1e3))
    by_lead = []
    for lead in LEADS:
        rr = [r for r in rows if r["lead"] == lead and r["impact"] <= F_PEAK * SHOCK]
        if not rr:
            print(f"  {lead:3.1f} mm: nothing survives the landing")
            continue
        r = min(rr, key=lambda r: r["mass"])
        x = r["act"]
        f_imp = r["impact"]
        bd = x.screw.drag * tau / (lead * 1e-3) / r["eta_back"]
        by_lead.append(dict(lead=lead, d_s=r["d_s"], spacer=r["spacer"], stack=r["stack"], od=r["od"], mass=r["mass"],
                            m_ref=r["m_ref"], impact=f_imp, eta=r["eta"], eta_back=r["eta_back"], backdrive=bd,
                            l10_years=r["days"] / 300))
        print(f"  {lead:3.1f} mm: Ø{r['d_s']} screw, Ø{r['od']:.0f} x {r['stack']:.0f} motor, {r['mass']:.2f} kg; "
              f"rotor seen at the rod {r['m_ref']:.0f} kg, landing {f_imp/1e3:.1f} kN; η {r['eta']:.2f} / back {r['eta_back']:.2f}; "
              f"drag to back-drive {bd:.0f} N")

    d = describe(best, req)
    d["by_lead"] = by_lead
    d["lightest_ignoring_landing"] = dict(lead=lightest["lead"], d_s=lightest["d_s"], spacer=lightest["spacer"],
                                          od=lightest["od"], mass=lightest["mass"], impact=lightest["impact"])
    print(f"\nlightest of all: {lightest['lead']} mm lead, Ø{lightest['od']:.0f} motor, {lightest['mass']:.3f} kg, "
          f"but landing {lightest['impact']/1e3:.1f} kN > {F_PEAK*SHOCK/1e3:.0f} kN")
    print("chosen (lightest that survives the landing):")
    for k, v in d.items():
        if isinstance(v, dict):
            print(f"  {k}: " + ", ".join(f"{kk}={vv:.4g}" if isinstance(vv, float) else f"{kk}={vv}" for kk, vv in v.items()
                                         if kk not in ("roller_phase", "parts")))
    print("  mass (g): " + ", ".join(f"{k} {v*1e3:.0f}" for k, v in d["actuator"]["parts"].items()))
    print("  roller thread phase by position (fraction of pitch): " + ", ".join(f"{x:.2f}" for x in d["screw"]["roller_phase"]))
    plots(best, [r for r in rows if r["impact"] <= F_PEAK * SHOCK], req)
    spec = dict(params=dict(screw=asdict(best.screw), stack=best.stack, slot_depth=best.slot_depth, turns=best.turns,
                            spacer=best.spacer, link=asdict(best.link)), stations=best.stations(), **d)
    (OUT / "design.json").write_text(json.dumps(spec, indent=1, default=float))
    print(f"\nwrote {OUT/'design.json'} and {OUT/'design.png'}")


def chosen():
    """the design main() picked, rebuilt from out/design.json (for cad.py)"""
    p = json.loads((OUT / "design.json").read_text())["params"]
    return Actuator(screw=Screw(**p["screw"]), stack=p["stack"], slot_depth=p["slot_depth"], turns=p["turns"], spacer=p["spacer"],
                    link=Linkage(**{**p["link"], "base": tuple(p["link"]["base"])}))


if __name__ == "__main__":
    main()
