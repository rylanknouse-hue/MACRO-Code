"""
Proj2_FermModel_TeamXX.py

ENGR 16100 - Project 2: Ethanol Production Facility Computational Model
Mass & Energy Balance, Fluid Mechanics, Capital Cost, and EROI Calculator.

This Python model supports engineering design teams in optimizing localized 
biofuel production from agricultural corn waste.

Key Features & Engineering Modeling:
1. Universal Accounting Equation (UAE) Mass Balances:
   - Fermenter: Conversions, ethanol generation (0.51 wt sugar), CO2 release (0.49 wt sugar).
   - Filter: Insoluble fiber removal efficiency (eta_filt).
   - Distiller: Concentrates ethanol exit stream (eta_dist).
   - Dehydrator: Excess water removal efficiency (eta_dehy).
2. Fluid Mechanics & Energy Losses:
   - Weighted average fluid density calculations based on component mass fractions.
   - Darcy-Weisbach friction head loss (H_DW) and fitting losses (bends, valves).
   - Mechanical pump work and daily electrical energy consumption.
3. Financial & Economic Analysis:
   - Complete capital cost breakdown including Site acquisition, taxes, mold remediation,
     unit operation equipment costs, pumps, piping, valves, and fittings.
   - Energy Return on Investment (EROI = Usable Energy Out / Electrical Energy In).
"""

import math

# --- PHYSICAL CONSTANTS & CONVERSIONS ---
G = 9.81  # m/s^2
SECONDS_PER_DAY = 86400.0
GALLONS_TO_M3 = 0.00378541  # 1 gallon = 0.00378541 m^3
M3_TO_GALLONS = 1.0 / GALLONS_TO_M3
KWH_TO_MJ = 3.6  # 1 kWh = 3.6 MJ
ENERGY_DENSITY_E100_MJ_PER_GAL = 80.1  # MJ / gallon of E100 ethanol

# Component Densities (kg/m^3)
RHO_SUGAR = 1599.0
RHO_ETHANOL = 789.0
RHO_FIBER = 1311.0
RHO_WATER = 997.0

# --- CATALOG DATA TABLES ---

# Fermenters: efficiency (% sugar converted), electrical (kWh/day), capital cost ($ per m^3/s flow rate)
FERMENTER_CATALOG = {
    "Scrap": {"eta": 0.50, "kwh_day": 46600.0, "cost_per_m3s": 320000.0},
    "Average": {"eta": 0.75, "kwh_day": 47200.0, "cost_per_m3s": 380000.0},
    "Premium": {"eta": 0.90, "kwh_day": 47500.0, "cost_per_m3s": 460000.0},
    "World-class": {"eta": 0.95, "kwh_day": 48000.0, "cost_per_m3s": 1100000.0},
}

# Distillers: efficiency (% alcohol exit purity by weight), electrical (kWh/day), capital cost ($ per m^3/s flow rate)
DISTILLER_CATALOG = {
    "Scrap": {"eta": 0.81, "kwh_day": 47004.0, "cost_per_m3s": 390000.0},
    "Average": {"eta": 0.90, "kwh_day": 47812.0, "cost_per_m3s": 460000.0},
    "Premium": {"eta": 0.915, "kwh_day": 48200.0, "cost_per_m3s": 560000.0},
    "World-class": {"eta": 0.98, "kwh_day": 49500.0, "cost_per_m3s": 1370000.0},
}

# Material Removal (Filters & Dehydrators): efficiency (% contaminant/water removed), electrical (kWh/day), capital cost ($ per m^3/s flow rate)
MATERIAL_REMOVAL_CATALOG = {
    "Scrap": {"eta": 0.50, "kwh_day": 48800.0, "cost_per_m3s": 200000.0},
    "Average": {"eta": 0.75, "kwh_day": 49536.0, "cost_per_m3s": 240000.0},
    "Premium": {"eta": 0.90, "kwh_day": 50350.0, "cost_per_m3s": 280000.0},
    "World-class": {"eta": 0.98, "kwh_day": 51000.0, "cost_per_m3s": 480000.0},
}

# Pipe Friction Factors (f) and unit prices per meter for D = 0.10 m
PIPE_CATALOG_D010 = {
    "Salvage": {"f": 0.05, "cost_per_m": 1000.0},
    "Questionable": {"f": 0.03, "cost_per_m": 1200.0},
    "Better": {"f": 0.02, "cost_per_m": 1440.0},
    "Nice": {"f": 0.01, "cost_per_m": 2160.0},
    "Outstanding": {"f": 0.005, "cost_per_m": 2700.0},
    "Glorious": {"f": 0.002, "cost_per_m": 2970.0},
}

# Site Survey Specs
SITE_CATALOG = {
    "Site 1": {
        "purchase_cost": 1500000.0,
        "tax_rate": 0.0,
        "remediation_cost": 0.0,
        "description": "Two 2-story buildings (2400 sqft total). Demolishes 10 low-income housing units.",
    },
    "Site 2": {
        "purchase_cost": 1300000.0,
        "tax_rate": 0.0,
        "remediation_cost": 0.0,
        "description": "L-shaped warehouse (3300 sqft total). Requires 4-mile access road near state park.",
    },
    "Site 3": {
        "purchase_cost": 1600000.0,
        "tax_rate": 0.10,  # 10% tax upfront
        "remediation_cost": 12000.0,  # mold remediation cost
        "description": "3-story tower (2300 sqft total). Vertical gravity flow. Job creation in low-income area.",
    },
}


def calculate_density(m_sugar, m_ethanol, m_fiber, m_water):
    """Calculates weighted fluid density (kg/m^3) based on mass fractions of components."""
    m_total = m_sugar + m_ethanol + m_fiber + m_water
    if m_total <= 0:
        return RHO_WATER
    w_sugar = m_sugar / m_total
    w_ethanol = m_ethanol / m_total
    w_fiber = m_fiber / m_total
    w_water = m_water / m_total
    return (w_sugar * RHO_SUGAR) + (w_ethanol * RHO_ETHANOL) + (w_fiber * RHO_FIBER) + (w_water * RHO_WATER)


def run_mass_balance(slurry_vol_gal_day, ferm_line, filt_line, dist_line, dehy_line):
    """
    Simulates UAE mass balance through the 4 sequential unit operations.
    """
    # Step 0: Raw Slurry Inputs (20% wt sugar, 20% wt fiber, 60% wt water)
    slurry_vol_m3_day = slurry_vol_gal_day * GALLONS_TO_M3
    rho_slurry = calculate_density(0.20, 0.0, 0.20, 0.60)
    slurry_mass_kg_day = slurry_vol_m3_day * rho_slurry

    m_sugar_in = 0.20 * slurry_mass_kg_day
    m_fiber_in = 0.20 * slurry_mass_kg_day
    m_water_in = 0.60 * slurry_mass_kg_day

    # Step 1: Fermenter
    eta_ferm = FERMENTER_CATALOG[ferm_line]["eta"]
    m_sugar_converted = m_sugar_in * eta_ferm
    m_ethanol_ferm_out = 0.51 * m_sugar_converted
    m_co2_waste = 0.49 * m_sugar_converted
    m_sugar_ferm_out = m_sugar_in * (1.0 - eta_ferm)
    m_water_ferm_out = m_water_in
    m_fiber_ferm_out = m_fiber_in

    # Step 2: Filter (removes fiber)
    eta_filt = MATERIAL_REMOVAL_CATALOG[filt_line]["eta"]
    m_fiber_filt_out = m_fiber_ferm_out * (1.0 - eta_filt)
    m_fiber_filt_waste = m_fiber_ferm_out * eta_filt
    m_sugar_filt_out = m_sugar_ferm_out
    m_water_filt_out = m_water_ferm_out
    m_ethanol_filt_out = m_ethanol_ferm_out

    # Step 3: Distiller (concentrates ethanol)
    eta_dist = DISTILLER_CATALOG[dist_line]["eta"]
    m_ethanol_dist_out = m_ethanol_filt_out

    denom = m_water_filt_out + m_sugar_filt_out + m_fiber_filt_out
    if denom > 0 and eta_dist > 0:
        ratio_factor = (m_ethanol_filt_out * (1.0 / eta_dist - 1.0)) / denom
    else:
        ratio_factor = 1.0

    m_water_dist_out = m_water_filt_out * ratio_factor
    m_sugar_dist_out = m_sugar_filt_out * ratio_factor
    m_fiber_dist_out = m_fiber_filt_out * ratio_factor

    m_water_dist_waste = m_water_filt_out - m_water_dist_out
    m_sugar_dist_waste = m_sugar_filt_out - m_sugar_dist_out
    m_fiber_dist_waste = m_fiber_filt_out - m_fiber_dist_out

    # Step 4: Dehydrator (removes excess water)
    eta_dehy = MATERIAL_REMOVAL_CATALOG[dehy_line]["eta"]
    m_water_dehy_out = m_water_dist_out * (1.0 - eta_dehy)
    m_water_dehy_waste = m_water_dist_out * eta_dehy
    m_ethanol_dehy_out = m_ethanol_dist_out
    m_sugar_dehy_out = m_sugar_dist_out
    m_fiber_dehy_out = m_fiber_dist_out

    # Final Output Summary
    m_total_final = m_ethanol_dehy_out + m_water_dehy_out + m_sugar_dehy_out + m_fiber_dehy_out
    purity_ethanol = (m_ethanol_dehy_out / m_total_final) if m_total_final > 0 else 0.0

    rho_final = calculate_density(m_sugar_dehy_out, m_ethanol_dehy_out, m_fiber_dehy_out, m_water_dehy_out)
    vol_final_m3_day = m_total_final / rho_final
    vol_final_gal_day = vol_final_m3_day * M3_TO_GALLONS

    # Volumetric flow rates into each unit (m^3/s)
    q_in_slurry_m3s = slurry_vol_m3_day / SECONDS_PER_DAY

    rho_ferm_out = calculate_density(m_sugar_ferm_out, m_ethanol_ferm_out, m_fiber_ferm_out, m_water_ferm_out)
    m_ferm_out_total = m_sugar_ferm_out + m_ethanol_ferm_out + m_fiber_ferm_out + m_water_ferm_out
    q_ferm_out_m3s = (m_ferm_out_total / rho_ferm_out) / SECONDS_PER_DAY

    rho_filt_out = calculate_density(m_sugar_filt_out, m_ethanol_filt_out, m_fiber_filt_out, m_water_filt_out)
    m_filt_out_total = m_sugar_filt_out + m_ethanol_filt_out + m_fiber_filt_out + m_water_filt_out
    q_filt_out_m3s = (m_filt_out_total / rho_filt_out) / SECONDS_PER_DAY

    rho_dist_out = calculate_density(m_sugar_dist_out, m_ethanol_dist_out, m_fiber_dist_out, m_water_dist_out)
    m_dist_out_total = m_sugar_dist_out + m_ethanol_dist_out + m_fiber_dist_out + m_water_dist_out
    q_dist_out_m3s = (m_dist_out_total / rho_dist_out) / SECONDS_PER_DAY

    return {
        "slurry_mass_kg_day": slurry_mass_kg_day,
        "slurry_vol_gal_day": slurry_vol_gal_day,
        "q_in_slurry_m3s": q_in_slurry_m3s,
        "q_ferm_out_m3s": q_ferm_out_m3s,
        "q_filt_out_m3s": q_filt_out_m3s,
        "q_dist_out_m3s": q_dist_out_m3s,
        "m_ethanol_out_kg_day": m_ethanol_dehy_out,
        "m_total_final_kg_day": m_total_final,
        "purity_ethanol": purity_ethanol,
        "vol_final_gal_day": vol_final_gal_day,
        "m_sugar_dist_waste": m_sugar_dist_waste,
        "m_co2_waste": m_co2_waste,
        "m_fiber_filt_waste": m_fiber_filt_waste,
    }


def run_full_model(
    slurry_vol_gal_day=690000.0,
    site_name="Site 3",
    ferm_line="World-class",
    filt_line="Premium",
    dist_line="World-class",
    dehy_line="World-class",
    pump_line="Standard",
    pump_eff=0.86,
    pump_lift_m=18.3,  # 60 ft vertical elevation (top floor fermenter)
    pipe_tier="Nice",
    pipe_d=0.10,
    pipe_len_m=18.3,
    num_bends=4,
    bend_xi=0.3,  # 90 degree bend
    num_valves=8,
    valve_k=800.0,
):
    """
    Executes complete system modeling (Mass Balance, Hydraulics, Electrical Energy, EROI, and Capital Cost).
    """
    # 1. Mass Balance Calculation
    mb = run_mass_balance(slurry_vol_gal_day, ferm_line, filt_line, dist_line, dehy_line)

    # 2. Hydraulic Loss Calculation
    q_m3s = mb["q_in_slurry_m3s"]
    rho_slurry = calculate_density(0.20, 0.0, 0.20, 0.60)
    area_pipe = math.pi * (pipe_d / 2.0) ** 2
    v_pipe = q_m3s / area_pipe

    f_friction = PIPE_CATALOG_D010[pipe_tier]["f"]

    # Head losses (m)
    h_dw = f_friction * (pipe_len_m / pipe_d) * (v_pipe**2 / (2.0 * G))
    h_bends = num_bends * (bend_xi * (v_pipe**2 / (2.0 * G)))
    h_valves = num_valves * (valve_k * (v_pipe**2 / (2.0 * G)))
    h_total = pump_lift_m + h_dw + h_bends + h_valves

    # Mechanical Work done by pump per day (Joules)
    e_pump_work_j_day = mb["slurry_mass_kg_day"] * G * h_total
    e_pump_work_kwh_day = e_pump_work_j_day / (KWH_TO_MJ * 1e6)

    # Electrical Energy required for pump
    e_pump_elec_kwh_day = e_pump_work_kwh_day / pump_eff

    # Electrical Energy required for Unit Operations
    e_ferm_elec_kwh_day = FERMENTER_CATALOG[ferm_line]["kwh_day"]
    e_filt_elec_kwh_day = MATERIAL_REMOVAL_CATALOG[filt_line]["kwh_day"]
    e_dist_elec_kwh_day = DISTILLER_CATALOG[dist_line]["kwh_day"]
    e_dehy_elec_kwh_day = MATERIAL_REMOVAL_CATALOG[dehy_line]["kwh_day"]

    e_total_elec_kwh_day = (
        e_pump_elec_kwh_day
        + e_ferm_elec_kwh_day
        + e_filt_elec_kwh_day
        + e_dist_elec_kwh_day
        + e_dehy_elec_kwh_day
    )
    e_total_elec_mj_day = e_total_elec_kwh_day * KWH_TO_MJ

    # Energy Output (MJ/day)
    e_out_mj_day = mb["vol_final_gal_day"] * ENERGY_DENSITY_E100_MJ_PER_GAL

    # EROI
    eroi = e_out_mj_day / e_total_elec_mj_day if e_total_elec_mj_day > 0 else 0.0

    # 3. Capital Cost Calculation ($)
    site_info = SITE_CATALOG[site_name]
    cost_site = site_info["purchase_cost"] * (1.0 + site_info["tax_rate"]) + site_info["remediation_cost"]

    cost_ferm = FERMENTER_CATALOG[ferm_line]["cost_per_m3s"] * mb["q_in_slurry_m3s"]
    cost_filt = MATERIAL_REMOVAL_CATALOG[filt_line]["cost_per_m3s"] * mb["q_ferm_out_m3s"]
    cost_dist = DISTILLER_CATALOG[dist_line]["cost_per_m3s"] * mb["q_filt_out_m3s"]
    cost_dehy = MATERIAL_REMOVAL_CATALOG[dehy_line]["cost_per_m3s"] * mb["q_dist_out_m3s"]

    cost_pipe = PIPE_CATALOG_D010[pipe_tier]["cost_per_m"] * pipe_len_m
    cost_bends = num_bends * 128.0  # 90 degree bend for 0.10m pipe
    cost_valves = num_valves * 100.0  # Salvage valve for K=800

    cost_pump = 460000.0 * mb["q_in_slurry_m3s"]

    capital_cost_total = (
        cost_site + cost_ferm + cost_filt + cost_dist + cost_dehy + cost_pipe + cost_bends + cost_valves + cost_pump
    )

    return {
        "site_name": site_name,
        "slurry_input_gal_day": slurry_vol_gal_day,
        "ethanol_output_gal_day": mb["vol_final_gal_day"],
        "ethanol_purity_percent": mb["purity_ethanol"] * 100.0,
        "fluid_velocity_m_s": v_pipe,
        "head_loss_dw_m": h_dw,
        "head_loss_total_m": h_total,
        "total_elec_kwh_day": e_total_elec_kwh_day,
        "total_elec_mj_day": e_total_elec_mj_day,
        "usable_energy_out_mj_day": e_out_mj_day,
        "eroi": eroi,
        "capital_cost_total": capital_cost_total,
        "cost_site": cost_site,
        "cost_equipment": cost_ferm + cost_filt + cost_dist + cost_dehy + cost_pump,
        "cost_piping_fittings": cost_pipe + cost_bends + cost_valves,
        "mass_balance": mb,
    }


def find_target_slurry_volume(target_ethanol_gal_day=100000.0, **model_kwargs):
    """
    Uses linear scaling to find the exact required corn slurry input volume (gal/day)
    to achieve the target ethanol output volume (e.g. 100,000 gal/day).
    """
    base_results = run_full_model(slurry_vol_gal_day=500000.0, **model_kwargs)
    base_out = base_results["ethanol_output_gal_day"]
    scale_factor = target_ethanol_gal_day / base_out
    required_slurry = 500000.0 * scale_factor
    return run_full_model(slurry_vol_gal_day=required_slurry, **model_kwargs)


def print_summary_report(results):
    """Prints a formatted execution summary of the model results."""
    print("=" * 70)
    print("      ENGR 16100 PROJECT 2: ETHANOL PRODUCTION MODEL RESULTS")
    print("=" * 70)
    print(f"Selected Facility Site     : {results['site_name']}")
    print(f"Corn Slurry Input Volume  : {results['slurry_input_gal_day']:,.0f} gal/day")
    print(f"Ethanol Output Volume     : {results['ethanol_output_gal_day']:,.2f} gal/day")
    print(f"Ethanol Purity             : {results['ethanol_purity_percent']:.2f}% (Requirement: >= 98.0%)")
    print("-" * 70)
    print("HYDRAULIC & FLUID METRICS:")
    print(f"  Fluid Velocity in Pipe   : {results['fluid_velocity_m_s']:.2f} m/s")
    print(f"  Friction Head Loss (H_DW): {results['head_loss_dw_m']:.2f} m")
    print(f"  Total Effective Head     : {results['head_loss_total_m']:.2f} m")
    print("-" * 70)
    print("ENERGY & EROI METRICS:")
    print(f"  Total Electrical In      : {results['total_elec_kwh_day']:,.2f} kWh/day ({results['total_elec_mj_day']:,.2f} MJ/day)")
    print(f"  Usable Ethanol Energy Out: {results['usable_energy_out_mj_day']:,.2f} MJ/day")
    print(f"  Energy ROI (EROI)        : {results['eroi']:.2f}x")
    print("-" * 70)
    print("FINANCIAL & CAPITAL COSTS:")
    print(f"  Site Acquisition & Tax   : ${results['cost_site']:,.2f}")
    print(f"  Equipment Capital Cost   : ${results['cost_equipment']:,.2f}")
    print(f"  Piping & Fittings Cost   : ${results['cost_piping_fittings']:,.2f}")
    print(f"  TOTAL CAPITAL COST       : ${results['capital_cost_total']:,.2f}")
    print("=" * 70)


if __name__ == "__main__":
    # Test Case 1: Model Validation Benchmark Case (from Model Validation Document)
    print("RUNNING MODEL VALIDATION BENCHMARK CASE:")
    val_results = run_full_model(
        slurry_vol_gal_day=835200.0,
        site_name="Site 3",
        ferm_line="Scrap",        # 50%
        filt_line="Scrap",        # 50%
        dist_line="Scrap",        # 81%
        dehy_line="Scrap",        # 50%
        pump_line="Cheap",        # 80%
        pump_eff=0.80,
    )
    print(f"Validation Ethanol Output : {val_results['ethanol_output_gal_day']/1000:.2f}k gal/day (Target Benchmark: ~68.70k gal/day)")
    print(f"Validation Purity         : {val_results['ethanol_purity_percent']:.1f}% (Target Benchmark: 87%)")
    print(f"Validation EROI           : {val_results['eroi']:.2f} (Target Benchmark: ~5.55)")
    print(f"Validation Capital Cost   : ${val_results['capital_cost_total']:,.0f} (Target Benchmark: ~$1,867,000)")
    print("\n" + "#"*70 + "\n")

    # Test Case 2: Optimal Recommended Configuration (Exact 100,000 gal/day @ >=98% Purity)
    print("RUNNING RECOMMENDED OPTIMAL DESIGN CONFIGURATION (Target: 100,000 gal/day E100 @ >=98% Purity):")
    opt_results = find_target_slurry_volume(
        target_ethanol_gal_day=100000.0,
        site_name="Site 3",
        ferm_line="World-class",      # 95% conversion
        filt_line="Premium",          # 90% filtration
        dist_line="World-class",      # 98% distillation purity
        dehy_line="World-class",      # 98% water removal
        pump_line="High-Grade",       # 89% pump efficiency
        pump_eff=0.89,
    )
    print_summary_report(opt_results)

