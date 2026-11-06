"""
Proj2_FermModel_TeamXX_v4.py

ENGR 16100 - Project 2: Ethanol Production Facility Computational Model (Version 4)
Mass & Energy Balance, Fluid Mechanics, Site Selection, Capital Cost, EROI, and Financial ROI Calculator.

Key Enhancements in Version 4:
1. Low-Resistance Glorious Valve Integration:
   - Configured with 'Glorious' low-resistance valves (K = 500, $297/valve) as default.
   - Reduces valve head loss from 4,826.8 m to 3,016.8 m (-37.4% total effective head reduction).
   - Reduces daily electrical energy draw from 244,544 kWh/day to 227,471 kWh/day.
   - Boosts Energy Return on Investment (EROI) from 9.10x to 9.78x.
2. Complete 4-Unit Operation Mass Balance: Fermenter, Filter, Distiller, and Dehydrator.
3. Interactive User Input (CLI prompts for Target Volume, Site, Fermenter, Filter, Distiller, Dehydrator, Pump, and Valve Tiers).
4. Detailed Equipment-by-Equipment Capital Cost & Efficiency Breakdown.
5. Financial ROI & Payback Period Analysis across all three candidate sites.
"""

import math
import sys

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

# Fermenters
FERMENTER_CATALOG = {
    "Scrap": {"eta": 0.50, "kwh_day": 46600.0, "cost_per_m3s": 320000.0},
    "Average": {"eta": 0.75, "kwh_day": 47200.0, "cost_per_m3s": 380000.0},
    "Premium": {"eta": 0.90, "kwh_day": 47500.0, "cost_per_m3s": 460000.0},
    "World-class": {"eta": 0.95, "kwh_day": 48000.0, "cost_per_m3s": 1100000.0},
}

# Distillers
DISTILLER_CATALOG = {
    "Scrap": {"eta": 0.81, "kwh_day": 47004.0, "cost_per_m3s": 390000.0},
    "Average": {"eta": 0.90, "kwh_day": 47812.0, "cost_per_m3s": 460000.0},
    "Premium": {"eta": 0.915, "kwh_day": 48200.0, "cost_per_m3s": 560000.0},
    "World-class": {"eta": 0.98, "kwh_day": 49500.0, "cost_per_m3s": 1370000.0},
}

# Material Removal (Filters & Dehydrators)
MATERIAL_REMOVAL_CATALOG = {
    "Scrap": {"eta": 0.50, "kwh_day": 48800.0, "cost_per_m3s": 200000.0},
    "Average": {"eta": 0.75, "kwh_day": 49536.0, "cost_per_m3s": 240000.0},
    "Premium": {"eta": 0.90, "kwh_day": 50350.0, "cost_per_m3s": 280000.0},
    "World-class": {"eta": 0.98, "kwh_day": 51000.0, "cost_per_m3s": 480000.0},
}

# Pumps Catalog
PUMP_CATALOG = {
    "Cheap": {"eff": 0.80, "cost_per_m3s": 260000.0},
    "Value": {"eff": 0.83, "cost_per_m3s": 310000.0},
    "Standard": {"eff": 0.86, "cost_per_m3s": 380000.0},
    "High-Grade": {"eff": 0.89, "cost_per_m3s": 460000.0},
    "Premium": {"eff": 0.92, "cost_per_m3s": 552000.0},
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

# Valves Catalog (Loss coefficient K and unit cost)
VALVE_CATALOG = {
    "Salvage": {"k": 800.0, "cost_per_unit": 100.0},
    "Glorious": {"k": 500.0, "cost_per_unit": 297.0},
}

# Site Survey Specifications & Impact Metrics
SITE_CATALOG = {
    "Site 1": {
        "purchase_cost": 1500000.0,
        "tax_rate": 0.0,
        "remediation_cost": 0.0,
        "infra_add_on_cost": 0.0,
        "pump_lift_m": 12.2,       # 40 ft lift (2-story buildings)
        "pipe_len_m": 30.5,        # 100 ft piping between 2 separate buildings
        "num_bends": 6,            # 6 pipe bends
        "layout_type": "Two 2-story buildings (2,400 sqft total)",
        "social_impact": "Demolishes 10 low-income housing units (housing replacement required).",
        "environmental_impact": "Low environmental disturbance; modern HVAC/electrical infrastructure.",
    },
    "Site 2": {
        "purchase_cost": 1300000.0,
        "tax_rate": 0.0,
        "remediation_cost": 0.0,
        "infra_add_on_cost": 150000.0, # Estimated 4-mile access road construction cost
        "pump_lift_m": 7.62,       # 25 ft lift (ledge height)
        "pipe_len_m": 36.6,        # 120 ft piping around L-shaped warehouse
        "num_bends": 6,            # 6 pipe bends
        "layout_type": "L-shaped warehouse w/ ledge (3,300 sqft total)",
        "social_impact": "Close to corn agribusiness supplier.",
        "environmental_impact": "4-mile access road impacts state park wildlife.",
    },
    "Site 3": {
        "purchase_cost": 1600000.0,
        "tax_rate": 0.10,          # 10% real estate tax upfront
        "remediation_cost": 12000.0, # Upper floor mold remediation
        "infra_add_on_cost": 0.0,
        "pump_lift_m": 15.24,      # 50 ft initial lift to Floor 3 fermenter (gravity cascade after)
        "pipe_len_m": 18.3,        # 60 ft total piping
        "num_bends": 4,            # 4 pipe bends
        "layout_type": "3-story tower (2,300 sqft total) with vertical gravity cascade",
        "social_impact": "Creates dozens of stable jobs in a low-income community.",
        "environmental_impact": "Minimal land disturbance; urban infrastructure revitalization.",
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
    Simulates UAE mass balance through all 4 sequential unit operations:
    Fermenter -> Filter -> Distiller -> Dehydrator.
    """
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

    # Step 2: Filter
    eta_filt = MATERIAL_REMOVAL_CATALOG[filt_line]["eta"]
    m_fiber_filt_out = m_fiber_ferm_out * (1.0 - eta_filt)
    m_fiber_filt_waste = m_fiber_ferm_out * eta_filt
    m_sugar_filt_out = m_sugar_ferm_out
    m_water_filt_out = m_water_ferm_out
    m_ethanol_filt_out = m_ethanol_ferm_out

    # Step 3: Distiller
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

    # Step 4: Dehydrator
    eta_dehy = MATERIAL_REMOVAL_CATALOG[dehy_line]["eta"]
    m_water_dehy_out = m_water_dist_out * (1.0 - eta_dehy)
    m_water_dehy_waste = m_water_dist_out * eta_dehy
    m_ethanol_dehy_out = m_ethanol_dist_out
    m_sugar_dehy_out = m_sugar_dist_out
    m_fiber_dehy_out = m_fiber_dist_out

    m_total_final = m_ethanol_dehy_out + m_water_dehy_out + m_sugar_dehy_out + m_fiber_dehy_out
    purity_ethanol = (m_ethanol_dehy_out / m_total_final) if m_total_final > 0 else 0.0

    rho_final = calculate_density(m_sugar_dehy_out, m_ethanol_dehy_out, m_fiber_dehy_out, m_water_dehy_out)
    vol_final_m3_day = m_total_final / rho_final
    vol_final_gal_day = vol_final_m3_day * M3_TO_GALLONS

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
    slurry_vol_gal_day=689575.0,
    site_name="Site 3",
    ferm_line="World-class",
    filt_line="Premium",
    dist_line="World-class",
    dehy_line="World-class",
    pump_line="High-Grade",
    pipe_tier="Nice",
    pipe_d=0.10,
    bend_xi=0.3,
    num_valves=8,
    valve_tier="Glorious",
    price_per_gal_ethanol=2.50,   # Wholesale ethanol selling price ($/gal)
    elec_cost_per_kwh=0.10,       # Electricity cost ($/kWh)
):
    """
    Executes complete system modeling (Mass Balance, Site Hydraulics, Energy/EROI, Capital Cost & Financial ROI).
    """
    # 1. Mass Balance
    mb = run_mass_balance(slurry_vol_gal_day, ferm_line, filt_line, dist_line, dehy_line)

    # 2. Site Hydraulics & Energy Loss
    site_info = SITE_CATALOG[site_name]
    pump_lift_m = site_info["pump_lift_m"]
    pipe_len_m = site_info["pipe_len_m"]
    num_bends = site_info["num_bends"]

    q_m3s = mb["q_in_slurry_m3s"]
    area_pipe = math.pi * (pipe_d / 2.0) ** 2
    v_pipe = q_m3s / area_pipe
    f_friction = PIPE_CATALOG_D010[pipe_tier]["f"]

    # Valve specs
    valve_info = VALVE_CATALOG.get(valve_tier, VALVE_CATALOG["Glorious"])
    valve_k = valve_info["k"]

    # Head losses (m)
    h_dw = f_friction * (pipe_len_m / pipe_d) * (v_pipe**2 / (2.0 * G))
    h_bends = num_bends * (bend_xi * (v_pipe**2 / (2.0 * G)))
    h_valves = num_valves * (valve_k * (v_pipe**2 / (2.0 * G)))
    h_total = pump_lift_m + h_dw + h_bends + h_valves

    pump_eff = PUMP_CATALOG[pump_line]["eff"]

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

    # Energy Return on Investment (EROI)
    eroi = e_out_mj_day / e_total_elec_mj_day if e_total_elec_mj_day > 0 else 0.0

    # 3. Capital Costs ($)
    cost_site = (
        site_info["purchase_cost"] * (1.0 + site_info["tax_rate"])
        + site_info["remediation_cost"]
        + site_info["infra_add_on_cost"]
    )

    cost_ferm = FERMENTER_CATALOG[ferm_line]["cost_per_m3s"] * mb["q_in_slurry_m3s"]
    cost_filt = MATERIAL_REMOVAL_CATALOG[filt_line]["cost_per_m3s"] * mb["q_ferm_out_m3s"]
    cost_dist = DISTILLER_CATALOG[dist_line]["cost_per_m3s"] * mb["q_filt_out_m3s"]
    cost_dehy = MATERIAL_REMOVAL_CATALOG[dehy_line]["cost_per_m3s"] * mb["q_dist_out_m3s"]
    cost_pump = PUMP_CATALOG[pump_line]["cost_per_m3s"] * mb["q_in_slurry_m3s"]

    cost_pipe = PIPE_CATALOG_D010[pipe_tier]["cost_per_m"] * pipe_len_m
    cost_bends = num_bends * 128.0
    cost_valves = num_valves * valve_info["cost_per_unit"]

    capital_cost_total = (
        cost_site + cost_ferm + cost_filt + cost_dist + cost_dehy + cost_pipe + cost_bends + cost_valves + cost_pump
    )

    # 4. Financial ROI & Operational Economics
    daily_revenue = mb["vol_final_gal_day"] * price_per_gal_ethanol
    annual_revenue = daily_revenue * 365.0

    daily_opex_elec = e_total_elec_kwh_day * elec_cost_per_kwh
    annual_opex_elec = daily_opex_elec * 365.0

    # Assume maintenance/labor/overhead is ~15% of annual revenue
    annual_other_opex = 0.15 * annual_revenue
    annual_net_profit = annual_revenue - annual_opex_elec - annual_other_opex

    financial_roi_percent = (annual_net_profit / capital_cost_total * 100.0) if capital_cost_total > 0 else 0.0
    payback_years = (capital_cost_total / annual_net_profit) if annual_net_profit > 0 else float("inf")

    return {
        "site_name": site_name,
        "site_info": site_info,
        "ferm_line": ferm_line,
        "filt_line": filt_line,
        "dist_line": dist_line,
        "dehy_line": dehy_line,
        "pump_line": pump_line,
        "valve_tier": valve_tier,
        "valve_k": valve_k,
        "slurry_input_gal_day": slurry_vol_gal_day,
        "ethanol_output_gal_day": mb["vol_final_gal_day"],
        "ethanol_purity_percent": mb["purity_ethanol"] * 100.0,
        "fluid_velocity_m_s": v_pipe,
        "head_loss_dw_m": h_dw,
        "head_loss_valves_m": h_valves,
        "head_loss_total_m": h_total,
        "total_elec_kwh_day": e_total_elec_kwh_day,
        "total_elec_mj_day": e_total_elec_mj_day,
        "usable_energy_out_mj_day": e_out_mj_day,
        "eroi": eroi,
        "capital_cost_total": capital_cost_total,
        "cost_site": cost_site,
        "cost_ferm": cost_ferm,
        "cost_filt": cost_filt,
        "cost_dist": cost_dist,
        "cost_dehy": cost_dehy,
        "cost_pump": cost_pump,
        "cost_equipment": cost_ferm + cost_filt + cost_dist + cost_dehy + cost_pump,
        "cost_piping_fittings": cost_pipe + cost_bends + cost_valves,
        "daily_revenue": daily_revenue,
        "annual_revenue": annual_revenue,
        "daily_opex_elec": daily_opex_elec,
        "annual_opex_elec": annual_opex_elec,
        "annual_net_profit": annual_net_profit,
        "financial_roi_percent": financial_roi_percent,
        "payback_years": payback_years,
        "mass_balance": mb,
    }


def find_target_slurry_volume(target_ethanol_gal_day=100000.0, **model_kwargs):
    """
    Solves for the required corn slurry input volume (gal/day) to achieve the target ethanol output.
    """
    base_results = run_full_model(slurry_vol_gal_day=500000.0, **model_kwargs)
    base_out = base_results["ethanol_output_gal_day"]
    scale_factor = target_ethanol_gal_day / base_out
    required_slurry = 500000.0 * scale_factor
    return run_full_model(slurry_vol_gal_day=required_slurry, **model_kwargs)


def print_summary_report(results):
    """Prints a detailed formatted execution summary of the model results."""
    mb = results["mass_balance"]
    print("\n" + "=" * 78)
    print("      ENGR 16100 PROJECT 2: ETHANOL PRODUCTION MODEL RESULTS (v4)")
    print("=" * 78)
    print(f"Selected Facility Site     : {results['site_name']} ({results['site_info']['layout_type']})")
    print(f"Corn Slurry Input Volume  : {results['slurry_input_gal_day']:,.0f} gal/day")
    print(f"Ethanol Output Volume     : {results['ethanol_output_gal_day']:,.2f} gal/day")
    print(f"Ethanol Exit Purity        : {results['ethanol_purity_percent']:.2f}% (Requirement: >= 98.0%)")
    print("-" * 78)
    print("UNIT OPERATIONS & EFFICIENCIES:")
    print(f"  1. Fermenter ({results['ferm_line']:<11}): {FERMENTER_CATALOG[results['ferm_line']]['eta']*100:.1f}% sugar conversion -> ${results['cost_ferm']:,.2f}")
    print(f"  2. Filter    ({results['filt_line']:<11}): {MATERIAL_REMOVAL_CATALOG[results['filt_line']]['eta']*100:.1f}% fiber removal    -> ${results['cost_filt']:,.2f}")
    print(f"     [Filter Fiber Removal: {mb['m_fiber_filt_waste']:,.1f} kg/day of dry waste fiber]")
    print(f"  3. Distiller ({results['dist_line']:<11}): {DISTILLER_CATALOG[results['dist_line']]['eta']*100:.1f}% alcohol exit purity -> ${results['cost_dist']:,.2f}")
    print(f"  4. Dehydrator({results['dehy_line']:<11}): {MATERIAL_REMOVAL_CATALOG[results['dehy_line']]['eta']*100:.1f}% water removal     -> ${results['cost_dehy']:,.2f}")
    print(f"  5. Pump      ({results['pump_line']:<11}): {PUMP_CATALOG[results['pump_line']]['eff']*100:.1f}% mechanical eff    -> ${results['cost_pump']:,.2f}")
    print(f"  6. Valves    ({results['valve_tier']:<11}): K = {results['valve_k']:.0f} loss coefficient")
    print("-" * 78)
    print("SITE & HYDRAULIC METRICS:")
    print(f"  Pumping Height Lift      : {results['site_info']['pump_lift_m']:.2f} m ({results['site_info']['pump_lift_m']*3.28084:.0f} ft)")
    print(f"  Pipe Length & Bends      : {results['site_info']['pipe_len_m']:.1f} m, {results['site_info']['num_bends']} bends")
    print(f"  Fluid Velocity in Pipe   : {results['fluid_velocity_m_s']:.2f} m/s")
    print(f"  Friction Head Loss (H_DW): {results['head_loss_dw_m']:.2f} m")
    print(f"  Valves Head Loss (H_valv): {results['head_loss_valves_m']:.2f} m")
    print(f"  Total Effective Head     : {results['head_loss_total_m']:.2f} m")
    print("-" * 78)
    print("ENERGY & EROI METRICS:")
    print(f"  Total Electrical In      : {results['total_elec_kwh_day']:,.2f} kWh/day ({results['total_elec_mj_day']:,.2f} MJ/day)")
    print(f"  Usable Ethanol Energy Out: {results['usable_energy_out_mj_day']:,.2f} MJ/day")
    print(f"  Energy ROI (EROI)        : {results['eroi']:.2f}x")
    print("-" * 78)
    print("FINANCIAL & CAPITAL COSTS (CapEx):")
    print(f"  Site Acquisition & Tax   : ${results['cost_site']:,.2f}")
    print(f"  Equipment Capital Cost   : ${results['cost_equipment']:,.2f}")
    print(f"  Piping & Fittings Cost   : ${results['cost_piping_fittings']:,.2f}")
    print(f"  TOTAL CAPITAL COST       : ${results['capital_cost_total']:,.2f}")
    print("-" * 78)
    print("FINANCIAL ROI & OPERATIONAL ECONOMICS:")
    print(f"  Annual Revenue ($2.50/gal): ${results['annual_revenue']:,.2f}/yr")
    print(f"  Annual Elec OpEx ($0.10/kwh): ${results['annual_opex_elec']:,.2f}/yr")
    print(f"  Estimated Annual Net Profit : ${results['annual_net_profit']:,.2f}/yr")
    print(f"  Simple Financial ROI (%)    : {results['financial_roi_percent']:.2f}% per year")
    print(f"  Capital Payback Period      : {results['payback_years']*365.0:.2f} days ({results['payback_years']:.4f} years)")
    print("-" * 78)
    print("COMMUNITY & ENVIRONMENTAL IMPACTS:")
    print(f"  Social Impact            : {results['site_info']['social_impact']}")
    print(f"  Environmental Impact     : {results['site_info']['environmental_impact']}")
    print("=" * 78 + "\n")


def compare_all_sites(target_ethanol_gal_day=100000.0, valve_tier="Glorious"):
    """Generates a comparative side-by-side summary table of Site 1, Site 2, and Site 3."""
    print("\n" + "#" * 84)
    print(f"    COMPARATIVE ANALYSIS OF FACILITY SITES ({valve_tier} Valves, Target: 100k gal/day)")
    print("#" * 84)
    header = f"{'Metric / Criteria':<28} | {'Site 1':<16} | {'Site 2':<16} | {'Site 3 (Recommended)':<18}"
    print(header)
    print("-" * 84)

    res1 = find_target_slurry_volume(target_ethanol_gal_day, site_name="Site 1", valve_tier=valve_tier)
    res2 = find_target_slurry_volume(target_ethanol_gal_day, site_name="Site 2", valve_tier=valve_tier)
    res3 = find_target_slurry_volume(target_ethanol_gal_day, site_name="Site 3", valve_tier=valve_tier)

    s1_cost = f"${res1['cost_site']:,.0f}"
    s2_cost = f"${res2['cost_site']:,.0f}"
    s3_cost = f"${res3['cost_site']:,.0f}"
    print(f"{'Site Acquisition Cost':<28} | {s1_cost:<16} | {s2_cost:<16} | {s3_cost:<18}")

    print(f"{'Layout / Architecture':<28} | {'2-Story Bldgs':<16} | {'L-Warehouse':<16} | {'3-Story Tower':<18}")

    s1_lift = f"{res1['site_info']['pump_lift_m']:.1f} m"
    s2_lift = f"{res2['site_info']['pump_lift_m']:.1f} m"
    s3_lift = f"{res3['site_info']['pump_lift_m']:.1f} m"
    print(f"{'Pumping Lift Height':<28} | {s1_lift:<16} | {s2_lift:<16} | {s3_lift:<18}")

    s1_head = f"{res1['head_loss_total_m']:.2f} m"
    s2_head = f"{res2['head_loss_total_m']:.2f} m"
    s3_head = f"{res3['head_loss_total_m']:.2f} m"
    print(f"{'Total Effective Head':<28} | {s1_head:<16} | {s2_head:<16} | {s3_head:<18}")

    s1_elec = f"{res1['total_elec_kwh_day']:,.0f} kWh"
    s2_elec = f"{res2['total_elec_kwh_day']:,.0f} kWh"
    s3_elec = f"{res3['total_elec_kwh_day']:,.0f} kWh"
    print(f"{'Daily Elec Consumption':<28} | {s1_elec:<16} | {s2_elec:<16} | {s3_elec:<18}")

    s1_eroi = f"{res1['eroi']:.2f}x"
    s2_eroi = f"{res2['eroi']:.2f}x"
    s3_eroi = f"{res3['eroi']:.2f}x"
    print(f"{'Energy ROI (EROI)':<28} | {s1_eroi:<16} | {s2_eroi:<16} | {s3_eroi:<18}")

    s1_cap = f"${res1['capital_cost_total']:,.0f}"
    s2_cap = f"${res2['capital_cost_total']:,.0f}"
    s3_cap = f"${res3['capital_cost_total']:,.0f}"
    print(f"{'Total Capital Cost (CapEx)':<28} | {s1_cap:<16} | {s2_cap:<16} | {s3_cap:<18}")

    s1_froi = f"{res1['financial_roi_percent']:.1f}%"
    s2_froi = f"{res2['financial_roi_percent']:.1f}%"
    s3_froi = f"{res3['financial_roi_percent']:.1f}%"
    print(f"{'Financial ROI (%/yr)':<28} | {s1_froi:<16} | {s2_froi:<16} | {s3_froi:<18}")

    s1_pay = f"{res1['payback_years']*365.0:.1f} days"
    s2_pay = f"{res2['payback_years']*365.0:.1f} days"
    s3_pay = f"{res3['payback_years']*365.0:.1f} days"
    print(f"{'Capital Payback Period':<28} | {s1_pay:<16} | {s2_pay:<16} | {s3_pay:<18}")

    print("#" * 84 + "\n")


def interactive_user_input():
    """Interactive CLI interface allowing custom parameter selection."""
    print("\n" + "=" * 60)
    print("   ETHANOL PRODUCTION MODEL: INTERACTIVE DESIGN SESSION (v4)")
    print("=" * 60)

    try:
        vol_str = input("Enter target ethanol output (gal/day) [Default: 100000]: ").strip()
        target_vol = float(vol_str) if vol_str else 100000.0
    except ValueError:
        target_vol = 100000.0

    print("\nSelect Facility Site:")
    print("  1) Site 1 (Two 2-story buildings, $1.5M)")
    print("  2) Site 2 (L-shaped warehouse, $1.3M + $150k road)")
    print("  3) Site 3 (3-story tower, $1.6M + 10% tax + $12k mold) [Recommended]")
    site_choice = input("Choice (1, 2, or 3) [Default: 3]: ").strip()
    site_map = {"1": "Site 1", "2": "Site 2", "3": "Site 3"}
    site_name = site_map.get(site_choice, "Site 3")

    tier_map = {"1": "Scrap", "2": "Average", "3": "Premium", "4": "World-class"}

    print("\nSelect Fermenter Line:")
    print("  1) Scrap (50%) 2) Average (75%) 3) Premium (90%) 4) World-class (95%) [Default]")
    f_choice = input("Choice (1-4) [Default: 4]: ").strip()
    ferm_line = tier_map.get(f_choice, "World-class")

    print("\nSelect Filter Line:")
    print("  1) Scrap (50%) 2) Average (75%) 3) Premium (90%) [Default] 4) World-class (98%)")
    filt_choice = input("Choice (1-4) [Default: 3]: ").strip()
    filt_line = tier_map.get(filt_choice, "Premium")

    print("\nSelect Distiller Line:")
    print("  1) Scrap (81%) 2) Average (90%) 3) Premium (91.5%) 4) World-class (98%) [Default]")
    d_choice = input("Choice (1-4) [Default: 4]: ").strip()
    dist_line = tier_map.get(d_choice, "World-class")

    print("\nSelect Dehydrator Line:")
    print("  1) Scrap (50%) 2) Average (75%) 3) Premium (90%) 4) World-class (98%) [Default]")
    dehy_choice = input("Choice (1-4) [Default: 4]: ").strip()
    dehy_line = tier_map.get(dehy_choice, "World-class")

    print("\nSelect Valve Tier:")
    print("  1) Salvage (K = 800, $100/valve)")
    print("  2) Glorious (K = 500, $297/valve) [Recommended default]")
    v_choice = input("Choice (1 or 2) [Default: 2]: ").strip()
    valve_tier = "Salvage" if v_choice == "1" else "Glorious"

    print("\nRunning simulation with your custom parameters...")
    results = find_target_slurry_volume(
        target_ethanol_gal_day=target_vol,
        site_name=site_name,
        ferm_line=ferm_line,
        filt_line=filt_line,
        dist_line=dist_line,
        dehy_line=dehy_line,
        valve_tier=valve_tier,
    )
    print_summary_report(results)


if __name__ == "__main__":
    if sys.stdin.isatty() and "--batch" not in sys.argv:
        interactive_user_input()
    else:
        print("RUNNING MODEL VALIDATION BENCHMARK CASE (Salvage Valves K=800):")
        val_results = run_full_model(
            slurry_vol_gal_day=835200.0,
            site_name="Site 3",
            ferm_line="Scrap",
            filt_line="Scrap",
            dist_line="Scrap",
            dehy_line="Scrap",
            pump_line="Cheap",
            valve_tier="Salvage",
        )
        print(f"Validation Ethanol Output : {val_results['ethanol_output_gal_day']/1000:.2f}k gal/day (Target Benchmark: ~68.70k gal/day)")
        print(f"Validation Purity         : {val_results['ethanol_purity_percent']:.1f}% (Target Benchmark: 87%)")
        print(f"Validation EROI           : {val_results['eroi']:.2f} (Target Benchmark: ~5.55)")
        print(f"Validation Capital Cost   : ${val_results['capital_cost_total']:,.0f} (Target Benchmark: ~$1,867,000)")

        compare_all_sites(target_ethanol_gal_day=100000.0, valve_tier="Glorious")

        opt_results = find_target_slurry_volume(target_ethanol_gal_day=100000.0, site_name="Site 3", valve_tier="Glorious")
        print_summary_report(opt_results)
