"""User-editable contact shear-strength callback, returns magnitude in Pa.

Example only: this reproduces a supplied constant, without material fitting.
Replace the body to implement a researched pressure/temperature/speed relation.
"""
import numpy as np


def shear_strength(state,parameters):
    # Available arrays, on real-contact cells only:
    # state.x_m, state.relative_velocity_m_s (signed), state.temperature_k,
    # state.nominal_contact_pressure_pa, state.real_contact_pressure_pa (=p_c/a),
    # state.relative_contact_area, state.mean_gap_m.
    # Use abs(state.relative_velocity_m_s) for a magnitude speed dependence.
    tau=float(parameters['shear_strength_pa'])
    return np.full_like(state.x_m,tau)
