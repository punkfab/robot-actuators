"""
RobStride MIT-protocol frames, with each model's own field ranges.

From the manuals' "6.5 Command 3: MIT Dynamic Parameters" (11-bit ID = motor CAN ID,
8 data bytes, big-endian bit fields):

    p   16 bit  ±12.57 rad (all models)
    v   12 bit  ±v_max   (RS00 33, RS01/02 44, RS03 20, RS04 15, RS05/06 50, RS10P 13)
    kp  12 bit  0–500 (RS00/01/02/05/10P)   0–5000 (RS03/04/06)
    kd  12 bit  0–5                         0–100
    τ   12 bit  ±τ_peak

The RANGES differ per model and the firmware decodes with its own. lerobot's
RobstrideMotorsBus (motors/robstride, 2026-06) encodes kp over 0–500 and kd over
0–5 for every model, and uses v ±33 for RS02/03/04. On an RS03/04/06 that would make a
commanded kp arrive 10× higher and kd 20× higher, if the manuals are right; not
verified on hardware here. Use encode() below, or check with a
motor on a stand at low gain first.

(RobStride's private protocol, type-1 extended frames, carries the same fields as
16-bit values; servo.py quantizes as the MIT frame does.)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from specs import SPECS, Spec  # noqa: E402

CMD_ENABLE = bytes([0xFF] * 7 + [0xFC])
CMD_DISABLE = bytes([0xFF] * 7 + [0xFD])
CMD_ZERO = bytes([0xFF] * 7 + [0xFE])
CMD_CLEAR = bytes([0xFF] * 7 + [0xFB])     # clear fault; the reply carries the state


def _u(x, lo, hi, bits):
    return int(np.round((np.clip(x, lo, hi) - lo) / (hi - lo) * ((1 << bits) - 1)))


def _f(u, lo, hi, bits):
    return lo + u * (hi - lo) / ((1 << bits) - 1)


def ranges(spec: Spec):
    return dict(p=(-spec.p_max, spec.p_max), v=(-spec.v_max, spec.v_max),
                kp=(0.0, spec.kp_max), kd=(0.0, spec.kd_max), t=(-spec.t_max, spec.t_max))


def encode(spec: Spec | str, p, v, kp, kd, t) -> bytes:
    s = SPECS[spec] if isinstance(spec, str) else spec
    r = ranges(s)
    pu, vu = _u(p, *r["p"], 16), _u(v, *r["v"], 12)
    ku, du, tu = _u(kp, *r["kp"], 12), _u(kd, *r["kd"], 12), _u(t, *r["t"], 12)
    return bytes([pu >> 8, pu & 0xFF, vu >> 4, ((vu & 0xF) << 4) | (ku >> 8), ku & 0xFF,
                  du >> 4, ((du & 0xF) << 4) | (tu >> 8), tu & 0xFF])


def decode_command(spec: Spec | str, data: bytes):
    """Inverse of encode(): what the firmware reads out of a frame."""
    s = SPECS[spec] if isinstance(spec, str) else spec
    r = ranges(s)
    pu = (data[0] << 8) | data[1]
    vu = (data[2] << 4) | (data[3] >> 4)
    ku = ((data[3] & 0xF) << 8) | data[4]
    du = (data[5] << 4) | (data[6] >> 4)
    tu = ((data[6] & 0xF) << 8) | data[7]
    return (_f(pu, *r["p"], 16), _f(vu, *r["v"], 12), _f(ku, *r["kp"], 12),
            _f(du, *r["kd"], 12), _f(tu, *r["t"], 12))


def decode_feedback(spec: Spec | str, data: bytes):
    """Reply frame: id, p (16), v (12), τ (12), MOS temperature ×10 (16)."""
    s = SPECS[spec] if isinstance(spec, str) else spec
    r = ranges(s)
    pu = (data[1] << 8) | data[2]
    vu = (data[3] << 4) | (data[4] >> 4)
    tu = ((data[4] & 0xF) << 8) | data[5]
    return dict(id=data[0], p=_f(pu, *r["p"], 16), v=_f(vu, *r["v"], 12),
                tau=_f(tu, *r["t"], 12), temp_mos=((data[6] << 8) | data[7]) / 10)


def lerobot_as_received(kind: str, kp, kd, v=0.0):
    """What an RS model would decode from lerobot's frame for (kp, kd, v), if the
    manual's ranges are what the firmware uses (lerobot: kp 0–500, kd 0–5, its v)."""
    lerobot_v = {"RS00": 33, "RS01": 44, "RS02": 33, "RS03": 33, "RS04": 33,
                 "RS05": 50, "RS06": 50}[kind]
    s = SPECS[kind]
    ku, du, vu = _u(kp, 0, 500, 12), _u(kd, 0, 5, 12), _u(v, -lerobot_v, lerobot_v, 12)
    return (_f(ku, 0, s.kp_max, 12), _f(du, 0, s.kd_max, 12), _f(vu, -s.v_max, s.v_max, 12))


if __name__ == "__main__":
    for kind in ("RS00", "RS02", "RS03", "RS04", "RS06"):
        cmd = (0.5, 1.0, 40.0, 2.0, 3.0)
        back = decode_command(kind, encode(kind, *cmd))
        err = max(abs(a - b) for a, b in zip(cmd, back))
        kp, kd, v = lerobot_as_received(kind, 40.0, 2.0, 1.0)
        print(f"{kind}: round-trip max error {err:.4f}   lerobot kp 40 kd 2 v 1.0 -> "
              f"kp {kp:6.1f} kd {kd:5.1f} v {v:5.2f}")
