"""
Put an extended friction model (models.py) on a MuJoCo hinge, per step.

MuJoCo already implements the budget-and-clip that BAM formalizes: dof_frictionloss
is a friction constraint whose force is clipped to ±frictionloss. So each step we
only rewrite dof_frictionloss (static part) and dof_damping (Kv) from the model,
the way bam/mujoco.py does:

  τm = qfrc_actuator[dof]                    ideal motor torque at the output
  τe = −qfrc_bias + qfrc_spring + qfrc_applied + (qfrc_constraint − own friction row)

qfrc_bias is gravity + Coriolis, so inertial reaction is NOT load (as in BAM).
Equality constraints (gear couplings, tendons) do count as external load. τe comes
from the previous step (the friction and the constraint forces are solved together);
MuJoCo's soft constraints keep it from jumping, which is BAM's argument too.

Creep. MuJoCo's friction constraint is soft: inside the budget the joint still slips
at a speed ∝ (1 − solimp). With the default solimp (0.9, 0.95) a 0.25 N·m load held by
0.08 N·m of friction creeps 8e-3 rad/s — 0.5°/s, enough to erase the static window
that the extended models exist to capture. We stiffen this dof's solimp to
(0.999, 0.9999): creep drops 100× (≈8e-5 rad/s), stable with implicitfast.

Usage:
    fr = FrictionUpdater(model, "joint_out", Friction.from_efficiency(0.83, 0.08))
    ...
    data.ctrl[:] = ...
    fr.step(data)          # update friction, then mj_step
"""

import mujoco
import numpy as np

from models import Friction


class FrictionUpdater:
    def __init__(self, model: mujoco.MjModel, joint: str, friction: Friction,
                 solimp=(0.999, 0.9999)):
        self.m = model
        self.f = friction
        jid = model.joint(joint).id
        if model.jnt_type[jid] not in (mujoco.mjtJoint.mjJNT_HINGE, mujoco.mjtJoint.mjJNT_SLIDE):
            raise ValueError(f"{joint}: friction models need a 1-dof hinge/slide joint")
        self.dof = model.jnt_dofadr[jid]
        model.dof_damping[self.dof] = friction.kv
        if solimp is not None:
            model.dof_solimp[self.dof, :2] = solimp
        self.tau_m = self.tau_e = 0.0

    def external_torque(self, d: mujoco.MjData) -> float:
        i = self.dof
        own = (d.efc_type == mujoco.mjtConstraint.mjCNSTR_FRICTION_DOF) & (d.efc_id == i)
        constraint = d.qfrc_constraint[i] - float(np.sum(d.efc_force[own]))
        return (-d.qfrc_bias[i] + d.qfrc_spring[i] + d.qfrc_applied[i] + constraint)

    def update(self, d: mujoco.MjData):
        """Refresh frictionloss from the current ctrl and the last solved state."""
        mujoco.mj_fwdActuation(self.m, d)          # qfrc_actuator for the NEW ctrl
        self.tau_m = float(d.qfrc_actuator[self.dof])
        self.tau_e = self.external_torque(d)
        self.m.dof_frictionloss[self.dof] = float(
            self.f.budget(self.tau_m, self.tau_e, d.qvel[self.dof]))

    def step(self, d: mujoco.MjData, n: int = 1):
        for _ in range(n):
            self.update(d)
            mujoco.mj_step(self.m, d)

    def prime(self, d: mujoco.MjData):
        """Call once after creating MjData so τe starts from a solved state."""
        mujoco.mj_forward(self.m, d)
        self.update(d)
        mujoco.mj_forward(self.m, d)
