"""
Record real SO-101 follower trajectories for real2sim (runs on the HARDWARE).

One joint moves through BAM-style excitation trajectories while the others hold a
base pose; every sample logs the goal and all present positions (+ current, load and
supply voltage for the moving joint). replay.py plays the same goals in MuJoCo with
each actuator model and scores them; fit.py tunes the model to the logs.

Run it from so101-lab's venv (it has lerobot + the Feetech SDK and your port/calibration
config), with this folder on the path:

    cd ~/sandbox/dnewcome/so101-lab
    uv run python ~/sandbox/punkfab/robot-actuators/sts3215/record.py \\
        --joint elbow_flex --traj all --payload-g 0 --out ~/so101_logs/

THIS MOVES THE ARM. Clear the workspace; keep a hand on the power. Like lerobot, it
disables torque on exit, so support the arm when the run ends. Amplitudes are
modest (≤ 20°) around a mid pose and every goal goes through lerobot's
max_relative_target clamp. The torque-off "drop" trajectory from the BAM paper is
deliberately NOT included: on an arm the link would fall.

Joint angles are lerobot degrees (use_degrees=True), which so101-lab's sim backend
already treats as the MJCF new-calibration angles. Check with --check before trusting
a replay: it reads the pose so you can compare with the sim's.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / "sandbox/dnewcome/so101-lab"))

JOINTS = ("shoulder_pan", "shoulder_lift", "elbow_flex", "wrist_flex", "wrist_roll")
BASE_POSE_DEG = dict(shoulder_pan=0.0, shoulder_lift=-35.0, elbow_flex=35.0,
                     wrist_flex=30.0, wrist_roll=0.0)       # a mid pose, all joints loaded
RATE_HZ = 100.0
KP_FIRMWARE = 16


def trajectories(name: str):
    """Goal offset (deg) as a function of t, and duration (s)."""
    if name == "chirp":          # accelerated oscillations: 0.3 -> 2 Hz, ±15°
        T = 12.0
        f0, f1 = 0.3, 2.0
        return (lambda t: 15.0 * np.sin(2 * np.pi * (f0 * t + (f1 - f0) * t * t / (2 * T)))), T
    if name == "ramp":           # raise and lower slowly: stiction / backdrive plateaus
        return (lambda t: 20.0 * np.interp(t, [0, 1, 7, 8, 14, 15], [0, 0, 1, 1, 0, 0]) - 10.0), 15.0
    if name == "steps":          # settling error and approach-direction hysteresis
        seq = [0, 5, 0, -5, 0, 10, 0, -10, 0, 20, 0, -20, 0]
        return (lambda t: float(seq[min(int(t / 1.5), len(seq) - 1)])), 1.5 * len(seq)
    if name == "sub":            # slow oscillation with fast small ripple
        return (lambda t: 12.0 * np.sin(2 * np.pi * 0.2 * t) + 3.0 * np.sin(2 * np.pi * 3.0 * t)), 12.0
    raise ValueError(name)


def connect():
    from lerobot.robots.so_follower.config_so_follower import SOFollowerRobotConfig
    from lerobot.robots.so_follower.so_follower import SOFollower
    from so101_config import ARMS
    cfg = SOFollowerRobotConfig(port=ARMS["follower"]["port"], id=ARMS["follower"]["id"],
                                max_relative_target=20.0, use_degrees=True)
    robot = SOFollower(cfg)
    robot.connect(calibrate=False)
    for m in robot.bus.motors:                 # make sure the gain matches the sim
        if m != "gripper":
            robot.bus.write("P_Coefficient", m, KP_FIRMWARE)
    return robot


def read_extra(robot, joint):
    out = {}
    for reg in ("Present_Current", "Present_Load", "Present_Voltage"):
        try:
            out[reg] = float(robot.bus.read(reg, joint, normalize=False))
        except Exception:
            out[reg] = np.nan
    return out


def move_to(robot, pose_deg, seconds=3.0):
    obs = robot.get_observation()
    start = {j: obs[f"{j}.pos"] for j in JOINTS}
    n = int(seconds * RATE_HZ)
    for k in range(1, n + 1):
        a = k / n
        robot.send_action({f"{j}.pos": start[j] + a * (pose_deg[j] - start[j]) for j in JOINTS})
        time.sleep(1 / RATE_HZ)


def record(robot, joint, traj, payload_g, out_dir: Path):
    f, T = trajectories(traj)
    move_to(robot, BASE_POSE_DEG)
    time.sleep(1.0)
    rows, extra = [], []
    t0 = time.perf_counter()
    while (t := time.perf_counter() - t0) < T:
        goal = dict(BASE_POSE_DEG)
        goal[joint] = BASE_POSE_DEG[joint] + f(t)
        robot.send_action({f"{j}.pos": goal[j] for j in JOINTS})
        obs = robot.get_observation()
        rows.append([t, goal[joint]] + [obs[f"{j}.pos"] for j in JOINTS])
        extra.append(list(read_extra(robot, joint).values()))
        time.sleep(max(0.0, (len(rows) / RATE_HZ) - (time.perf_counter() - t0)))
    rows, extra = np.array(rows), np.array(extra)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    path = out_dir / f"{joint}_{traj}_{payload_g:.0f}g_{stamp}.npz"
    meta = dict(joint=joint, traj=traj, payload_g=payload_g, kp=KP_FIRMWARE,
                base_pose_deg=BASE_POSE_DEG, joints=JOINTS, rate_hz=RATE_HZ,
                columns=["t", "goal_deg"] + [f"{j}_deg" for j in JOINTS],
                extra_columns=["current_raw", "load_raw", "voltage_raw_0.1V"])
    np.savez(path, data=rows, extra=extra, meta=json.dumps(meta))
    dt = np.diff(rows[:, 0])
    print(f"  {path.name}: {len(rows)} samples, {1/np.median(dt):.0f} Hz median, "
          f"supply {np.nanmedian(extra[:, 2])/10:.1f} V")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--joint", default="elbow_flex", choices=JOINTS)
    ap.add_argument("--traj", default="all", help="chirp | ramp | steps | sub | all")
    ap.add_argument("--payload-g", type=float, default=0.0, help="mass held at the gripper")
    ap.add_argument("--out", default=str(Path.home() / "so101_logs"))
    ap.add_argument("--check", action="store_true", help="just read and print the pose")
    args = ap.parse_args()
    robot = connect()
    try:
        if args.check:
            obs = robot.get_observation()
            print({j: round(obs[f"{j}.pos"], 1) for j in JOINTS})
            print("voltage", read_extra(robot, "shoulder_lift")["Present_Voltage"] / 10, "V")
            return
        out = Path(args.out).expanduser()
        out.mkdir(parents=True, exist_ok=True)
        names = ["chirp", "ramp", "steps", "sub"] if args.traj == "all" else [args.traj]
        for name in names:
            record(robot, args.joint, name, args.payload_g, out)
        move_to(robot, BASE_POSE_DEG)
    finally:
        robot.bus.disable_torque()
        robot.disconnect()


if __name__ == "__main__":
    main()
