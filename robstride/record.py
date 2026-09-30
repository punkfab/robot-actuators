"""
Record a RobStride on the pendulum bench, for bench.py to fit.

NOT YET RUN ON HARDWARE. Written from the manuals' MIT-protocol section (frame
layout in can_mit.py); start with --check, then --zero, then a low load.

Setup
  - The motor in MIT protocol mode (RobStride's MotorStudio, comm mode 2), CAN at 1 Mbps,
    Linux socketcan:  sudo ip link set can0 up type can bitrate 1000000
  - Output axis horizontal, the bench arm (bench.R = 0.25 m to the mass centre) with
    the mass bench.py asks for (it prints kg per load), hanging straight down
  - pip install python-can (the robot-actuators venv or any other)

    python robstride/record.py --kind RS03 --id 1 --check          # read state only
    python robstride/record.py --kind RS03 --id 1 --zero           # arm hanging: set zero
    python robstride/record.py --kind RS03 --id 1 --load 0.25 --traj coast
    ...                                                            # the bench.PLAN set
    python robstride/bench.py RS03 logs/*.npz

Each run: enable, ease from the current pose to the trajectory's start (min-jerk,
not logged), run the trajectory at bench.HOST_HZ with bench.command(), logging the
frame sent and the reply (q, q̇, τ = driver's Kt·iq estimate, MOS temperature), then
disable. The same command() drives the sim, so bench.py replays exactly what was sent.
Safety: |τ| clamped to --max-torque (default: the load's gravity torque × 2.5), the
soft limit in bench.command(), and any missed reply or CTRL-C disables the motor.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bench  # noqa: E402
from can_mit import CMD_CLEAR, CMD_DISABLE, CMD_ENABLE, CMD_ZERO, decode_feedback, encode  # noqa: E402
from specs import SPECS  # noqa: E402


class Motor:
    def __init__(self, kind, can_id, channel="can0"):
        import can
        self.can = can
        self.kind, self.id = kind, can_id
        self.bus = can.Bus(interface="socketcan", channel=channel, bitrate=1_000_000)

    def xfer(self, data, timeout=0.005):
        self.bus.send(self.can.Message(arbitration_id=self.id, data=data, is_extended_id=False))
        t_end = time.perf_counter() + timeout
        while time.perf_counter() < t_end:
            msg = self.bus.recv(timeout=max(0.0, t_end - time.perf_counter()))
            if msg is not None and len(msg.data) >= 8 and msg.data[0] == self.id:
                return decode_feedback(self.kind, msg.data)
        return None

    def mit(self, p, v, kp, kd, t):
        return self.xfer(encode(self.kind, p, v, kp, kd, t))

    def close(self):
        try:
            self.xfer(CMD_DISABLE)
        finally:
            self.bus.shutdown()


def ease_to(motor, q_start, q_goal, kind, load, seconds=3.0):
    s = SPECS[kind]
    mgr = bench.bench_mass(kind, load) * 9.81 * bench.R
    kp = min(4 * mgr, s.kp_max)
    kd = min(2 * 0.7 * np.sqrt(kp * bench.bench_mass(kind, load) * bench.R ** 2), s.kd_max)
    n = int(seconds * bench.HOST_HZ)
    for k in range(n):
        u = k / (n - 1)
        p = q_start + (q_goal - q_start) * (10 * u**3 - 15 * u**4 + 6 * u**5)
        fb = motor.mit(p, 0.0, kp, kd, mgr * np.sin(p))
        if fb is None:
            raise RuntimeError("no reply while easing to the start pose")
        time.sleep(1 / bench.HOST_HZ)
    return fb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", required=True, choices=sorted(SPECS))
    ap.add_argument("--id", type=int, required=True)
    ap.add_argument("--channel", default="can0")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--zero", action="store_true")
    ap.add_argument("--traj", choices=bench.TRAJS)
    ap.add_argument("--load", type=float, default=bench.LOAD_FRACTION)
    ap.add_argument("--max-torque", type=float)
    ap.add_argument("--out", default=str(Path.home() / "robstride_logs"))
    a = ap.parse_args()

    motor = Motor(a.kind, a.id, a.channel)
    try:
        if a.check:
            fb = motor.xfer(CMD_CLEAR)
            print(fb if fb else "no reply: check wiring, ID, bitrate and MIT protocol mode")
            return
        if a.zero:
            print(motor.xfer(CMD_ZERO))
            return
        if not a.traj:
            ap.error("--traj, --check or --zero")
        mgr = bench.bench_mass(a.kind, a.load) * 9.81 * bench.R
        tmax = a.max_torque or 2.5 * mgr
        print(f"{a.kind} id {a.id}: bench mass {bench.bench_mass(a.kind, a.load):.2f} kg at "
              f"{bench.R} m, |τ| ≤ {tmax:.1f} N·m. Enter to start, CTRL-C to abort.")
        input()
        fb = motor.xfer(CMD_ENABLE)
        if fb is None or abs(fb["p"]) > np.radians(15):
            raise RuntimeError(f"not hanging near zero ({fb}); run --zero with the arm hanging")
        ease_to(motor, fb["p"], bench.start_angle(a.traj), a.kind, a.load)
        n = int(bench.DURATION[a.traj] * bench.HOST_HZ)
        log = {k: [] for k in ("t", "cmd", "q", "dq", "tau", "temp", "wall")}
        q, dq, t0 = fb["p"], fb["v"], time.perf_counter()
        for k in range(n):
            t_next = t0 + k / bench.HOST_HZ
            while time.perf_counter() < t_next:
                pass
            c = list(bench.command(a.traj, k / bench.HOST_HZ, q, dq, a.kind, a.load))
            c[4] = float(np.clip(c[4], -tmax, tmax))
            fb = motor.mit(*c)
            if fb is None:
                raise RuntimeError(f"missed reply at step {k}")
            q, dq = fb["p"], fb["v"]
            for key, v in zip(log, (k / bench.HOST_HZ, c, q, dq, fb["tau"], fb["temp_mos"],
                                    time.perf_counter() - t0)):
                log[key].append(v)
        out = Path(a.out)
        out.mkdir(exist_ok=True)
        path = out / f"{a.kind}_{a.traj}_{a.load:.2f}_{time.strftime('%Y%m%d-%H%M%S')}.npz"
        np.savez(path, kind=a.kind, traj=a.traj, load=a.load,
                 **{k: np.array(v) for k, v in log.items()})
        print(f"wrote {path}  (late frames: {int(np.sum(np.diff(log['wall']) > 1.5 / bench.HOST_HZ))})")
    finally:
        motor.close()


if __name__ == "__main__":
    main()
