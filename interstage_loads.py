"""
interstage_loads.py  --  Load and boundary-condition derivation for the
hot-staging lattice-truss interstage (Section 3 and Appendix B).

Requirements (versions used): Python 3.11, numpy 2.x, scipy 1.x,
ambiance 1.3.1 (standard atmosphere),
rocketcea 1.2.3 (Python wrapper of CEA2, the second version of the NASA Chemical
Equilibrium with Applications (CEA) Fortran code), cantera 3.2.0 (validation only).
Run:  python interstage_loads.py

Names used in the code and output:
  NASA  National Aeronautics and Space Administration
  NACA  National Advisory Committee for Aeronautics
  USSA  U.S. Standard Atmosphere 1976
  tf    tonne-force;  MECO  main engine cutoff;  RVac  Raptor Vacuum engine
  mq / maxq  max-q (maximum dynamic pressure);  pm / premeco  pre-MECO
  qs    quasi-static;  _d  after the dynamic amplification factor (DAF)
  CNa   normal-force slope C_N_alpha;  cdc  crossflow drag coefficient
  Pc    chamber pressure;  MR  oxidizer-to-fuel mass ratio (O/F);  eps  area ratio
  Tc, Te, pe, Me, Ve, ae  chamber temperature, exit temperature, exit pressure,
        exit Mach number, exit velocity, exit speed of sound
  gam   ratio of specific heats;  MW  molecular weight (kg/kmol)
  hp / tp  CEA problem with assigned enthalpy and pressure / temperature and pressure
  stag  stagnation state;  FE  finite element model;  2-D  two-dimensional
  y / u yield (x 1.25) / ultimate (x 1.4) design load;  m/s2  m/s^2;  deg  degrees
"""
import numpy as np
from scipy.optimize import brentq
from rocketcea.cea_obj_w_units import CEA_Obj

g0 = 9.80665          # standard gravity used to convert tf -> N
g = 9.81              # gravitational acceleration used in the load factors
R = 4.5               # interstage outer radius (m), Section 2
S = np.pi * R**2      # reference area (m^2)
R_m = 4.45            # radius of the ring wall midline, on which the members are centred (m)

# --------------------------------------------------------------------------
# Vehicle data (Section 2, Table 2.1)
# --------------------------------------------------------------------------
m_up = 1800e3                       # upper stage wet + payload (kg)
m_tot = {"liftoff": 5700e3, "maxq": 4300e3, "premeco": 2500e3}
T = 8240 * 1e3 * g0                 # Super Heavy thrust, 8,240 tf -> N
print(f"Thrust T = {T/1e6:.2f} MN")

# --------------------------------------------------------------------------
# 3.2  Axial compression
# --------------------------------------------------------------------------
rho_10 = 0.41271     # USSA, H = 10 km (kg/m^3)
T_10 = 223.15        # USSA, H = 10 km (K)
V_mq = 416.0         # max-q velocity (m/s)
q_mq = 0.5 * rho_10 * V_mq**2
CD = 0.125
D = q_mq * S * CD
print(f"\n[3.2] q_maxq = {q_mq/1e3:.2f} kPa, D = {D/1e6:.3f} MN, D/T = {100*D/T:.2f} %")

def axial(case, drag=0.0, D_up=0.0):
    W = m_tot[case] * g
    nx = (T - drag) / W
    N = m_up * nx * g + D_up
    return nx, N, N / (2 * np.pi * R_m)

ax = {}
ax["liftoff"] = axial("liftoff")
ax["maxq"] = axial("maxq", drag=D, D_up=D)        # all drag assigned to upper stage
ax["maxq_noDup"] = axial("maxq", drag=D, D_up=0.0)
ax["premeco"] = axial("premeco")
for k, (nx, N, nN) in ax.items():
    print(f"  {k:11s} n_x = {nx:.3f}  N = {N/1e6:6.2f} MN  n_N = {nN/1e6:.3f} MN/m")
N_pad = m_up * g                                   # 1 g weight on the pad
N_hs = 1614 * 1e3 * g0                             # full upper-stage thrust, 1,614 tf
print(f"  pad        N = {N_pad/1e6:6.2f} MN  n_N = {N_pad/(2*np.pi*R_m)/1e6:.3f} MN/m")
print(f"  hot-stage  N = {N_hs/1e6:6.2f} MN  n_N = {N_hs/(2*np.pi*R_m)/1e6:.3f} MN/m "
      f"({100*N_hs/ax['premeco'][1]:.0f} % of pre-MECO)")

# Roughness-limited drag regime (Niskanen 2009, Eq. 3.79)
beta_s, S_suth = 1.458e-6, 110.4                  # USSA 1976 Eq. (51)
mu_10 = beta_s * T_10**1.5 / (T_10 + S_suth)
nu_10 = mu_10 / rho_10
L_veh, Rs = 124.0, 60e-6
Re = V_mq * L_veh / nu_10
Re_crit = 51 * (Rs / L_veh) ** -1.039
print(f"  mu_10 = {mu_10:.4e} Pa s, Re = {Re:.3e}, Re_crit = {Re_crit:.3e}, Re/Re_crit = {Re/Re_crit:.1f}")

# --------------------------------------------------------------------------
# 3.3  Shear and bending
# --------------------------------------------------------------------------
V_wind = 50.0
alpha = V_wind / V_mq
# slender-body + viscous crossflow (Allen & Perkins, NACA Rep. 1048)
eta, cdc, Ap = 0.692, 1.2, 52.0 * 9.0
CN = np.sin(2 * alpha) * np.cos(alpha / 2) + eta * cdc * (Ap / S) * np.sin(alpha) ** 2
CNa_AP = CN / alpha
CNa = 3.0
F_aero = q_mq * S * CNa * alpha
# crossflow force on the booster (constant-diameter cylinder: potential term is zero)
Ap_b = 72.0 * 9.0
F_boost = q_mq * eta * cdc * Ap_b * np.sin(alpha) ** 2
# center of pressure (Allen & Perkins, NACA Rep. 1048, Eq. 14): the potential force acts
# at L - Q/S_b from the nose, the crossflow force at the planform centroid
L_up, L_n = 52.0, 12.0                             # upper-stage length, nose length (m)
def h_cp_of(Ln):
    rho_o = (R**2 + Ln**2) / (2 * R)                # tangent-ogive radius
    xs = np.linspace(0.0, Ln, 20001)
    ys = np.sqrt(rho_o**2 - (Ln - xs)**2) + R - rho_o
    Qn = np.trapezoid(np.pi * ys**2, xs)            # nose volume
    h_pot = L_up - (Ln - Qn / S)                    # height above the interstage
    C_pot = np.sin(2 * alpha) * np.cos(alpha / 2)
    C_cf = eta * cdc * (Ap / S) * np.sin(alpha) ** 2
    return h_pot, (C_pot * h_pot + C_cf * L_up / 2) / (C_pot + C_cf), C_pot / (C_pot + C_cf)
h_pot, h_cp, w_pot = h_cp_of(L_n)
k_trim = 1.25
a_lat_mq = (k_trim * F_aero + F_boost) / m_tot["maxq"]
F_in = m_up * a_lat_mq
h_cg = 25.0
V_qs = F_aero - F_in
M_qs = F_aero * h_cp - F_in * h_cg
DAF = 2.0
V_mq_d, M_mq_d = DAF * V_qs, DAF * M_qs
print(f"\n[3.3] alpha = {alpha:.4f} rad ({np.degrees(alpha):.2f} deg), C_Na(Allen-Perkins) = {CNa_AP:.3f} /rad")
print(f"  potential share = {w_pot:.3f}, h_pot = {h_pot:.2f} m, h_cp = {h_cp:.2f} m "
      f"(L_n 10 m: {h_cp_of(10.0)[1]:.2f} m, L_n 14 m: {h_cp_of(14.0)[1]:.2f} m)")
print(f"  F_boost = {F_boost/1e6:.3f} MN, total lateral force = {(k_trim*F_aero+F_boost)/1e6:.3f} MN")
print(f"  F_aero = {F_aero/1e6:.3f} MN, a_lat = {a_lat_mq:.4f} m/s2 ({a_lat_mq/g:.4f} g), F_in = {F_in/1e6:.3f} MN")
print(f"  V_qs = {V_qs/1e6:.3f} MN, M_qs = {M_qs/1e6:.2f} MN m -> V = {V_mq_d/1e6:.3f} MN, M = {M_mq_d/1e6:.2f} MN m")
t_rise = 30.0 / V_mq
print(f"  gust rise time = {t_rise:.3f} s")

delta = np.radians(1.0)
F_side = T * np.sin(delta)
a_lat_pm = F_side / m_tot["premeco"]
V_pm = m_up * a_lat_pm
M_pm = V_pm * h_cg
print(f"  pre-MECO: F_side = {F_side/1e6:.3f} MN, a_lat = {a_lat_pm:.4f} m/s2 ({a_lat_pm/g:.4f} g), V = {V_pm/1e6:.3f} MN, M = {M_pm/1e6:.2f} MN m")

def line(Vs, Ms):
    return Vs / (np.pi * R_m), Ms / (np.pi * R_m**2)
qmq, nMmq = line(V_mq_d, M_mq_d)
qpm, nMpm = line(V_pm, M_pm)
print(f"  max-q:    q_max = {qmq/1e6:.4f} MN/m, n_M = {nMmq/1e6:.4f} MN/m")
print(f"  pre-MECO: q_max = {qpm/1e6:.4f} MN/m, n_M = {nMpm/1e6:.4f} MN/m")

# --------------------------------------------------------------------------
# 3.4  Pre-staging thermal state
# --------------------------------------------------------------------------
M_mq = V_mq / np.sqrt(1.4 * 8314.32 / 28.9644 * T_10)      # USSA Eq. (50)
r = 0.89
Tr_air = T_10 * (1 + r * 0.2 * M_mq**2)
print(f"\n[3.4] M_maxq = {M_mq:.3f}, T_r(air, max-q) = {Tr_air:.1f} K")

# Ascent heating from max-q to separation (flight 13 webcast telemetry, booster)
from ambiance import Atmosphere                    # standard atmosphere (checked below)
a10, a55 = Atmosphere(10e3), Atmosphere(55e3)      # USSA Table I, geometric altitude
print(f"  ambiance check: rho(10 km) = {a10.density[0]:.5f} (USSA 0.41351), "
      f"rho(55 km) = {a55.density[0]:.2e} (USSA 5.7e-4), p(55 km) = {a55.pressure[0]:.1f} Pa (USSA 42.5)")
t_tl = np.array([58.0, 75.0, 90.0, 105.0, 120.0, 138.0, 141.0])          # s after liftoff
z_tl = np.array([10.0, 14.3, 20.7, 29.1, 39.6, 53.0, 55.0]) * 1e3        # m
v_tl = np.array([1500, 1960, 2783, 3800, 4943, 5900, 5900]) / 3.6        # m/s
gam_a, R_a = 1.4, 8314.32 / 28.9644                # USSA Eq. (50): gamma = 1.40
cp_a = gam_a * R_a / (gam_a - 1)
x_tube = np.pi * 0.100                             # running length, as in Section 3.5
rho_ti, c_ti, k_ti = 4420.0, 580.0, 6.6            # Ti-6Al-4V at 20 C (TIMET datasheet)

def ascent_h(tt, Tw):
    a_ = Atmosphere(np.interp(tt, t_tl, z_tl)); V_ = np.interp(tt, t_tl, v_tl)
    Ta, pa = a_.temperature[0], a_.pressure[0]
    M_ = V_ / np.sqrt(gam_a * R_a * Ta)
    Tr_ = Ta * (1 + r * (gam_a - 1) / 2 * M_**2)
    Ts_ = Ta + 0.5 * (Tw - Ta) + 0.22 * (Tr_ - Ta)                     # Eckert
    rs_ = pa / (R_a * Ts_)
    mu_ = 1.458e-6 * Ts_**1.5 / (Ts_ + 110.4)                           # USSA Eq. (51)
    k_ = 2.64638e-3 * Ts_**1.5 / (Ts_ + 245.4 * 10 ** (-12 / Ts_))      # USSA Eq. (53)
    Pr_ = mu_ * cp_a / k_
    Re_ = rs_ * V_ * x_tube / mu_
    h_ = 0.0296 * Re_**-0.2 * Pr_**(-2 / 3) * rs_ * V_ * cp_a          # turbulent flat plate
    return h_, Tr_, M_, Re_

def ascent_wall(t_w, dt=0.01):
    """lumped tube wall heated on its outer face from T+0:58 (300 K) to separation"""
    Tw, tt, qmax, Q = 300.0, t_tl[0], (0.0, 0.0), 0.0
    while tt < t_tl[-1] - 1e-9:
        h_, Tr_, _, _ = ascent_h(tt, Tw)
        q_ = h_ * (Tr_ - Tw)
        if q_ > qmax[0]:
            qmax = (q_, tt)
        Q += q_ * dt
        Tw += q_ * dt / (rho_ti * c_ti * t_w)
        tt += dt
    return Tw, q_, qmax, Q

print("  ascent heating (turbulent, x = 0.314 m):")
for tt in t_tl:
    h_, Tr_, M_, Re_ = ascent_h(tt, 300.0)
    print(f"    t = {tt:5.0f} s  M = {M_:.2f}  T_r = {Tr_:6.1f} K  h = {h_:6.1f} W/m2K  "
          f"Re_x = {Re_:.2e}  q(300 K) = {h_*(Tr_-300)/1e3:5.2f} kW/m2")
for t_w in (0.002, 0.003, 0.004, 0.005, 0.006, 0.008, 0.010, 0.012):
    Tw, q_end, qmax, Q = ascent_wall(t_w)
    print(f"    tube wall {t_w*1e3:4.0f} mm: T_w(separation) = {Tw:6.1f} K, "
          f"q_max = {qmax[0]/1e3:.1f} kW/m2 at T+{qmax[1]:.0f} s, Q = {Q/1e3:.0f} kJ/m2, "
          f"through-wall dT at separation = {q_end*t_w/(2*k_ti):.2f} K")
# 100 mm ring: 1-D conduction, outer face heated, inner face adiabatic
n_r, L_r = 101, 0.100
dx_r = L_r / (n_r - 1); al_r = k_ti / (rho_ti * c_ti); dt_r = 0.4 * dx_r**2 / al_r
T_ring, tt = np.full(n_r, 300.0), t_tl[0]
while tt < t_tl[-1] - 1e-9:
    h_, Tr_, _, _ = ascent_h(tt, T_ring[0])
    Tn = T_ring.copy()
    Tn[1:-1] += al_r * dt_r / dx_r**2 * (T_ring[2:] - 2 * T_ring[1:-1] + T_ring[:-2])
    Tn[0] += 2 * al_r * dt_r / dx_r**2 * (T_ring[1] - T_ring[0]) \
        + 2 * dt_r * h_ * (Tr_ - T_ring[0]) / (rho_ti * c_ti * dx_r)
    Tn[-1] += 2 * al_r * dt_r / dx_r**2 * (T_ring[-2] - T_ring[-1])
    T_ring, tt = Tn, tt + dt_r
print(f"    ring 100 mm: surface {T_ring[0]:.1f} K, at 10 mm depth {T_ring[10]:.1f} K, "
      f"mean {T_ring.mean():.1f} K")

# --------------------------------------------------------------------------
# 3.5  Plume impingement chain
# --------------------------------------------------------------------------
cea = CEA_Obj(oxName="LOX", fuelName="CH4", pressure_units="bar",
              temperature_units="K", cstar_units="m/s",
              sonic_velocity_units="m/s", isp_units="sec",
              enthalpy_units="J/kg", density_units="kg/m^3",
              viscosity_units="millipoise")
Pc, MR = 350.0, 3.6
print("\n[3.5] CEA (RocketCEA) exit states")
ex = {}
for name, eps in (("Raptor", 34.34), ("RVac", 80.0), ("eps945", 945.0)):
    Tc, Tt, Te = cea.get_Temperatures(Pc=Pc, MR=MR, eps=eps)
    MWe, ge = cea.get_exit_MolWt_gamma(Pc=Pc, MR=MR, eps=eps)
    pr = cea.get_PcOvPe(Pc=Pc, MR=MR, eps=eps)
    Me = cea.get_MachNumber(Pc=Pc, MR=MR, eps=eps)
    ae = cea.get_SonicVelocities(Pc=Pc, MR=MR, eps=eps)[2]
    ex[name] = dict(Tc=Tc, Te=Te, pe=Pc / pr, Me=Me, g=ge, MW=MWe, ae=ae)
    print(f"  {name:6s} Tc={Tc:.1f} Te={Te:.1f} K pe={Pc/pr:.4f} bar Me={Me:.4f} "
          f"gam={ge:.4f} MW={MWe:.3f} Ve={Me*ae:.1f} m/s")

gam = ex["eps945"]["g"]
Me, pe, Te = ex["RVac"]["Me"], ex["RVac"]["pe"], ex["RVac"]["Te"]
theta = np.arctan(0.30 / 1.35)

def nu(M, g_):
    k_ = np.sqrt((g_ + 1) / (g_ - 1))
    return k_ * np.arctan(np.sqrt((M**2 - 1) / k_**2)) - np.arctan(np.sqrt(M**2 - 1))

def chain(g_, theta_=theta, Te_=Te, pe_=pe, Me_=Me):
    M1 = brentq(lambda M: nu(M, g_) - (nu(Me_, g_) + theta_), Me_, 50.0)
    f = (1 + (g_ - 1) / 2 * Me_**2) / (1 + (g_ - 1) / 2 * M1**2)
    p1 = pe_ * f ** (g_ / (g_ - 1))
    T1 = Te_ * f
    po2 = p1 * ((g_ + 1)**2 * M1**2 / (4 * g_ * M1**2 - 2 * (g_ - 1)))**(g_ / (g_ - 1)) \
        * (1 - g_ + 2 * g_ * M1**2) / (g_ + 1)
    return M1, p1, T1, po2

M1, p1, T1, po2 = chain(gam)
T0 = cea.get_Tcomb(Pc=po2, MR=MR)          # hp equilibrium at the pitot pressure
Tr = T1 + r * (T0 - T1)
print(f"  theta = {np.degrees(theta):.8f} deg, nu(Me) = {np.degrees(nu(Me,gam)):.2f} deg")
print(f"  M1 = {M1:.4f}, p1 = {p1:.5f} bar, T1 = {T1:.1f} K, p_o2 = {po2:.4f} bar")
print(f"  T0 = {T0:.1f} K, T_r = {Tr:.1f} K")

# Stagnation-point heating: Fay-Riddell form (Hamilton, NASA Technical Memorandum 84580, Eqs. 4a, 6)
# edge (stagnation) state: CEA equilibrium at (p_o2, h0)  -> RocketCEA chamber at Pc = p_o2
_, mu_s_mP, _, Pr_s = cea.get_Chamber_Transport(Pc=po2, MR=MR, frozen=0)
mu_s = mu_s_mP * 1e-4                                   # millipoise -> Pa s
rho_s = cea.get_Chamber_Density(Pc=po2, MR=MR)
H0 = cea.get_Chamber_H(Pc=po2, MR=MR)                # J/kg
# wall state: CEA tp equilibrium at (300 K, p_o2) -- CEARUN run "TPprops" (Appendix A.3)
rho_w, mu_w, h_w = 1.1675, 0.15048e-4, -12028.5e3
Tw, Rn = 300.0, 0.100
dudx = (1 / Rn) * np.sqrt(2 * (po2 - p1) * 1e5 / rho_s)
q_s = 0.76 / Pr_s**0.6 * (rho_w * mu_w)**0.1 * (rho_s * mu_s)**0.4 * (H0 - h_w) * np.sqrt(dudx)
h_s = q_s / (T0 - Tw)
print(f"  stag: rho_s={rho_s:.4f} mu_s={mu_s:.4e} Pr_s={Pr_s:.4f} H0={H0/1e3:.2f} kJ/kg")
print(f"  du/dx = {dudx:.0f} 1/s, q_s = {q_s/1e6:.3f} MW/m2, h_s = {h_s:.0f} W/m2K "
      f"(2-D line: {h_s/np.sqrt(2):.0f})")

# Tangential surfaces: Eckert reference temperature + turbulent flat plate
Tstar = T1 + 0.5 * (Tw - T1) + 0.22 * (Tr - T1)
# CEA tp equilibrium at (T*, p1) -- CEARUN run "TPprops" (Appendix A.3)
rho_st, mu_st, cp_st, Pr_st = 2.6337e-2, 0.52230e-4, 1958.4, 0.7437
R_gas = 8314.46 / ex["eps945"]["MW"]
u1 = M1 * np.sqrt(gam * R_gas * T1)
x = np.pi * Rn                                   # half circumference of the tube
Re_x = rho_st * u1 * x / mu_st
St = 0.0296 * Re_x**-0.2 * Pr_st**(-2 / 3)
h_t = St * rho_st * u1 * cp_st
print(f"  T* = {Tstar:.1f} K, u1 = {u1:.0f} m/s, Re_x = {Re_x:.3e}, St = {St:.5f}, h_t = {h_t:.0f} W/m2K")

q_dyn1 = 0.5 * gam * p1 * 1e5 * M1**2
print(f"  plume dynamic pressure at state 1: {q_dyn1/1e3:.1f} kPa")
f_s, T_ref_FE = 0.10, 3000.0
q_blend = f_s * h_s * (T0 - Tw) + (1 - f_s) * h_t * (Tr - Tw)
h_blend = q_blend / (T_ref_FE - Tw)
h_FE = 750.0
print(f"  q_blend = {q_blend:.0f} W/m2 -> h_blend = {h_blend:.0f} W/m2K; adopted {h_FE:.0f} "
      f"(q_FE(300 K) = {h_FE*(T_ref_FE-Tw):.0f} W/m2, +{100*(h_FE/h_blend-1):.1f} %)")

# --------------------------------------------------------------------------
# 3.6  Limit and ultimate line loads
# --------------------------------------------------------------------------
FS_u, FS_y = 1.4, 1.25
print("\n[3.6] line loads (MN/m): limit / yield (x1.25) / ultimate (x1.4)")
for nm, nN, nM, qv in (("max-q", ax["maxq"][2], nMmq, qmq), ("pre-MECO", ax["premeco"][2], nMpm, qpm)):
    peak = nN + nM
    print(f"  {nm:8s} n_N={nN/1e6:.3f} n_M={nM/1e6:.3f} q={qv/1e6:.4f} peak={peak/1e6:.3f} | "
          f"y: {FS_y*peak/1e6:.3f} | u: {FS_u*peak/1e6:.3f}  (n_M/n_N = {nM/nN:.2f})")
print("  resultants (FE input): limit / yield / ultimate")
for nm, (Nr, Vr, Mr) in (("max-q", (ax["maxq"][1], V_mq_d, M_mq_d)), ("pre-MECO", (ax["premeco"][1], V_pm, M_pm))):
    for lab, val, u in (("N", Nr, "MN"), ("V", Vr, "MN"), ("M", Mr, "MN m")):
        print(f"  {nm:8s} {lab} = {val/1e6:.3f} / {FS_y*val/1e6:.3f} / {FS_u*val/1e6:.3f} {u}")

# --------------------------------------------------------------------------
# 3.7  Validation and sensitivity
# --------------------------------------------------------------------------
print("\n[3.7] validation")
# energy check at the RVac exit with the CEARUN enthalpies (Appendix A.1, Table A.1)
H_c_run, H_e_run, Ve_run = -1526.63e3, -8018.22e3, 3603.2
print(f"  V_e from enthalpy = {np.sqrt(2*(H_c_run-H_e_run)):.1f} m/s vs CEA exit velocity {Ve_run} m/s (CEARUN, eps = 80)")
e_p1 = cea.get_eps_at_PcOvPe(Pc=Pc, MR=MR, PcOvPe=Pc / p1)
T_dir = cea.get_Temperatures(Pc=Pc, MR=MR, eps=e_p1)[2]
M_dir = cea.get_MachNumber(Pc=Pc, MR=MR, eps=e_p1)
print(f"  direct CEA expansion to p1: eps = {e_p1:.1f}, T = {T_dir:.1f} K, M = {M_dir:.3f} "
      f"(chain: {T1:.1f} K, {M1:.3f})")
try:
    import cantera as ct
    gas = ct.Solution("gri30.yaml")
    X0 = {"CH4": 1.0, "O2": MR * 16.04246 / 31.9988}
    H0_c = H0
    def heq(Tg):
        gas.TPX = Tg, po2 * 1e5, X0
        gas.equilibrate("TP")
        return gas.enthalpy_mass
    T0_ct = brentq(lambda Tg: heq(Tg) - H0_c, 2000, 4500, xtol=1e-4)
    print(f"  Cantera {ct.__version__} (GRI-Mech 3.0) T0 = {T0_ct:.1f} K vs CEA {T0:.1f} K "
          f"({100*(T0_ct/T0-1):+.2f} %)")
except ImportError:
    print("  Cantera not installed; skipped")

print("\n  sensitivity (one input changed at a time)")
base_pm = ax["premeco"][1]
def N_pm(mu=m_up, mt=m_tot["premeco"], TT=T):
    return mu * TT / mt
print(f"  pre-MECO N: m_up +100 t -> {100*(N_pm(mu=m_up+100e3)/base_pm-1):+.1f} %; "
      f"booster dry -50 t -> {100*(N_pm(mt=m_tot['premeco']-50e3)/base_pm-1):+.1f} %")
for lab, kw in (("C_Na 2.71", dict(c=CNa_AP)), ("h_cp = h_pot", dict(hac=h_pot)),
                ("V_wind 40 m/s", dict(vw=40.0)), ("no trim", dict(kt=1.0)),
                ("no booster force", dict(fb=0.0))):
    c = kw.get("c", CNa); hac = kw.get("hac", h_cp); vw = kw.get("vw", V_wind); kt = kw.get("kt", k_trim)
    al = vw / V_mq
    Fa = q_mq * S * c * al
    Fb = kw.get("fb", q_mq * eta * cdc * Ap_b * np.sin(al) ** 2)
    Fi = m_up * (kt * Fa + Fb) / m_tot["maxq"]
    Mx = DAF * (Fa * hac - Fi * h_cg)
    print(f"  max-q M with {lab:12s}: {Mx/1e6:6.2f} MN m ({100*(Mx/M_mq_d-1):+.1f} %), "
          f"peak line load {(ax['maxq'][2] + Mx/(np.pi*R_m**2))/1e6:.3f} MN/m")
for lab, gg in (("gamma 1.1944 (exit)", ex["RVac"]["g"]), ("gamma 1.2033 (at p1)", 1.2033)):
    m_, p_, t_, po_ = chain(gg)
    T0_ = cea.get_Tcomb(Pc=po_, MR=MR)
    print(f"  {lab:22s}: T1={t_:.1f} K  T0={T0_:.1f} K  T_r={t_+r*(T0_-t_):.1f} K")
r2 = 0.7437 ** (1 / 3)
print(f"  r = Pr*^(1/3) = {r2:.3f}: T_r = {T1 + r2*(T0-T1):.1f} K")
for th in (np.radians(10.0), np.radians(15.0)):
    m_, p_, t_, po_ = chain(gam, theta_=th)
    T0_ = cea.get_Tcomb(Pc=po_, MR=MR)
    print(f"  theta = {np.degrees(th):.0f} deg: M1={m_:.3f} p1={p_:.4f} T1={t_:.1f} T_r={t_+r*(T0_-t_):.1f} K")
for fs in (0.05, 0.20):
    hb = (fs * h_s * (T0 - Tw) + (1 - fs) * h_t * (Tr - Tw)) / (T_ref_FE - Tw)
    print(f"  f_s = {fs:.2f}: h_blend = {hb:.0f} W/m2K")
print(f"  2-D stagnation line (k=0): h_blend = "
      f"{(f_s*h_s/np.sqrt(2)*(T0-Tw)+(1-f_s)*h_t*(Tr-Tw))/(T_ref_FE-Tw):.0f} W/m2K")
