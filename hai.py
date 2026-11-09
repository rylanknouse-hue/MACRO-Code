"""ENGR 161 Project 2 - ethanol plant model + exhaustive design optimizer (see report)"""
import numpy as np, itertools, math

G = 9.80665; SEC = 86400.0; FT = 0.3048
RHO = dict(s=1599.0, e=789.0, f=1311.0, w=997.0)
E_DENS = 21160.1772            # MJ/m^3 ethanol
V_OUT_TARGET = 378.54125342579829   # m^3/day (100,000 gal)
PURITY_MIN = 0.98
SLURRY_DENS = 0.2*1599 + 0.2*1311 + 0.6*997   # 1180.2

# ---------------- catalog data ----------------
F_EFF  = [0.5, 0.75, 0.9, 0.95];  F_KWH  = [46600, 47200, 47500, 48000];  F_COST  = [320000, 380000, 460000, 1100000]
FIL_EFF= [0.5, 0.75, 0.9, 0.98];  FIL_KWH= [48800, 49536, 50350, 51000];  FIL_COST= [200000, 240000, 280000, 480000]
D_EFF  = [0.81, 0.9, 0.915, 0.98];D_KWH  = [47004, 47812, 48200, 49500];  D_COST  = [390000, 460000, 560000, 1370000]
H_EFF, H_KWH, H_COST = FIL_EFF, FIL_KWH, FIL_COST
LINE4 = ['Scrap', 'Average', 'Premium', 'World-class']

PUMP_EFF = [0.80, 0.83, 0.86, 0.89, 0.92]; PUMP_COLS = ['Cheap','Value','Standard','High-Grade','Premium']
PUMP_ROWS = [6,9,12,15,18,21,24,27,30,33,36]
PUMP_PRICE = np.array([
 [200000,240000,280000,340000,415000],[220000,260000,310000,380000,456000],
 [240000,290000,340000,410000,502000],[260000,310000,380000,460000,552000],
 [290000,350000,420000,500000,607000],[320000,380000,460000,550000,668000],
 [350000,420000,510000,610000,735000],[390000,460000,560000,670000,808000],
 [420000,510000,610000,740000,889000],[470000,560000,670000,810000,978000],
 [510000,620000,740000,890000,1076000]], float)

DIAMS = [0.10,0.11,0.12,0.13,0.14,0.15]
PIPE_GRADES = ['Salvage','Questionable','Better','Nice','Outstanding','Glorious']
PIPE_F = [0.05,0.03,0.02,0.01,0.005,0.002]
PIPE_COST = np.array([[1000,1200,1440,2160,2700,2970],[1200,1440,1720,2580,3230,3550],
  [2570,3080,3700,5550,6940,7640],[6300,7560,9070,14000,17000,19000],
  [14000,16000,20000,29000,37000,40000],[26000,31000,37000,55000,69000,76000]], float)
VALVE_GRADES = ['Salvage','Questionable','Outstanding','Glorious']
VALVE_K = [800,700,600,500]
VALVE_COST = np.array([[100,120,270,297],[120,144,323,355],[257,308,694,764],
  [630,756,1700,1900],[1400,1600,3700,4000],[2600,3100,6900,7600]], float)
BEND_ANGLES = [20,30,45,60,75,90]; BEND_K = [0.1,0.15,0.2,0.22,0.27,0.3]
BEND_COST = np.array([[100,105,110,116,122,128],[149,157,164,173,181,190],
  [493,517,543,570,599,700],[1400,1500,1600,1600,1700,1800],
  [3200,3400,3600,3800,3900,4100],[6200,6500,6900,7200,7600,8000]], float)
DUCT = {1.0:228.0, 1.25:414.0, 1.5:700.0}

SITE_COST = {1:1_500_000.0, 2:1_300_000.0, 3:1_600_000.0*1.10}
MOLD_COST = 12_000.0   # $20/ft^2 x 600 ft^2, only if the moldy top floor is used

# ---------------- process (mass balance) ----------------
def dens(s,e,f,w):
    t = s+e+f+w
    return (s*RHO['s']+e*RHO['e']+f*RHO['f']+w*RHO['w'])/t

def process(m0, ef, efil, ed, eh):
    """m0 = slurry mass flow kg/s. Returns dict of mass flows (kg/s), densities, flows (m3/s)."""
    s0, f0, w0 = 0.2*m0, 0.2*m0, 0.6*m0
    conv = ef*s0; eth = 0.51*conv; co2 = 0.49*conv
    s1, f1, w1 = s0-conv, f0, w0
    m1 = s1+eth+f1+w1
    frem = efil*f1; f2 = f1-frem; s2, w2 = s1, w1
    m2 = s2+eth+f2+w2
    oth_in = s2+f2+w2; oth_out = eth*(1-ed)/ed
    fr = min(1.0, oth_out/oth_in)
    s3, f3, w3 = s2*fr, f2*fr, w2*fr
    m3 = eth+s3+f3+w3
    wrem = eh*w3; w4 = w3-wrem
    m4 = eth+s3+f3+w4
    m = [m0,m1,m2,m3,m4]
    rho = [dens(s0,0,f0,w0), dens(s1,eth,f1,w1), dens(s2,eth,f2,w2), dens(s3,eth,f3,w3), dens(s3,eth,f3,w4)]
    Q = [m[i]/rho[i] for i in range(5)]
    waste = [co2, frem, oth_in-(s3+f3+w3), wrem]       # F, Fil, D, H waste (kg/s)
    return dict(m=m, rho=rho, Q=Q, waste=waste, eth=eth, purity=eth/m4, m4=m4)

def vout_per_vin(ef, efil, ed, eh):
    """output stream volume (m3/day) per m3/day slurry input - linear scaling."""
    p = process(SLURRY_DENS/SEC, ef, efil, ed, eh)   # V_in = 1 m3/day
    return p['Q'][4]*SEC, p['purity']

# ---------------- layout definitions ----------------
# Each site: list of slot groups (name, building, level_ft, capacity), building roof heights (ft),
# and whether different buildings are physically separate (need transfer piping).
SITES = {
 1: dict(slots=[('A-L0','A',0,1),('A-L1','A',20,1),('B-L0','B',0,1),('B-L1','B',20,1)],
         roof={'A':40,'B':40}, separate=True),
 2: dict(slots=[('Ground','W',0,4),('Ledge','W',25,1)], roof={'W':50}, separate=False),
 3: dict(slots=[('Tall-L0','T',0,1),('Small-L0','S',0,1),('Tall-L1','T',20,1),('Tall-L2','T',40,1)],
        roof={'T':60,'S':20}, separate=False),
}
CROSS_LINK_FT = 30.0   # 20 ft gap between Site-1 buildings + 10 ft rise to inlet
SAME_FLOOR_FT = 10.0   # rise to the top inlet of the next unit on the same floor

def layouts(site):
    """All feasible assignments of (F,Fil,D,H) to slots; downhill (non-ascending) cascade."""
    S = SITES[site]; slots = S['slots']; out = []
    for combo in itertools.product(range(len(slots)), repeat=4):
        cnt = {}
        for c in combo: cnt[c] = cnt.get(c,0)+1
        if any(cnt[c] > slots[c][3] for c in cnt): continue
        lv = [slots[c][2] for c in combo]
        if any(lv[i] < lv[i+1] for i in range(3)): continue      # must be non-ascending
        out.append(combo)
    return out

def layout_geometry(site, combo):
    S = SITES[site]; sl = [S['slots'][c] for c in combo]
    z_ft = [s[2] for s in sl]; bl = [s[1] for s in sl]
    links = []   # for streams 1..3: (length_m, n_bends)
    for i in range(3):
        if S['separate'] and bl[i] != bl[i+1]:
            links.append((CROSS_LINK_FT*FT, 2))
        elif z_ft[i] == z_ft[i+1]:
            links.append((SAME_FLOOR_FT*FT, 2))
        else:
            links.append((0.0, 0))
    L0 = (z_ft[0]+10)*FT
    duct = max(0.0, (S['roof'][bl[0]] - (z_ft[0]+10)))*FT
    mold = MOLD_COST if site == 3 else 0.0  
    return dict(z=[z*FT for z in z_ft], z_ft=z_ft, bl=bl, links=links, L0=L0, duct=duct, mold=mold,
                names=[s[0] for s in sl])

def pump_price(eff_idx, lift_m):
    for r, rating in enumerate(PUMP_ROWS):
        if rating >= lift_m: return PUMP_PRICE[r, eff_idx], rating
    raise ValueError('lift too high')

# ---------------- evaluator ----------------

def evaluate(site, combo, grades, pump_i, Vin, diam, pipe_g, valve_g, bend_angle_idx=5, duct_d=1.0,
             legacy_flow_cost=False, legacy_bend_price=None, layout_override=None):
    """
    grades=(f,fil,d,h) indices 0..3 ; diam=[D0..D4] per stream ; pipe_g=[g0,g1,g2,g3] (stream0 pipe + links)
    valve_g = list of 8 valve grade idx: [F_in, F_out, Fil_in, Fil_out, D_in, D_out, H_in, H_out]
    Returns full dict of results (energies in kJ/day).
    """
    geo = layout_override or layout_geometry(site, combo)
    gf, gfil, gd, gh = grades
    m0 = Vin*SLURRY_DENS/SEC
    P = process(m0, F_EFF[gf], FIL_EFF[gfil], D_EFF[gd], H_EFF[gh])
    Q, rho, m = P['Q'], P['rho'], P['m']
    A = [math.pi*d*d/4 for d in diam]
    v = [Q[k]/A[k] for k in range(5)]
    dyn = [0.5*rho[k]*v[k]**2 for k in range(5)]                  # Pa
    kJ = lambda watts: watts*SEC/1000.0
    # valve -> stream map
    vstream = [0,1,1,2,2,3,3,4]
    E = {}
    E['valves'] = [kJ(VALVE_K[valve_g[i]]*dyn[vstream[i]]*Q[vstream[i]]) for i in range(8)]
    # pipes
    fric = []; Lpipes=[(0, geo['L0'])]+[(k+1, geo['links'][k][0]) for k in range(3)]
    for gi,(k,L) in zip(range(4),Lpipes):
        fric.append(kJ(PIPE_F[pipe_g[gi]]*(L/diam[k])*dyn[k]*Q[k]) if L>0 else 0.0)
    E['friction'] = fric
    nb = [2]+[geo['links'][k][1] for k in range(3)]
    E['bends'] = [kJ(nb[i]*BEND_K[bend_angle_idx]*dyn[i]*Q[i]) for i in range(4)]
    # waste
    wm = P['waste']
    E['waste_gpe'] = [kJ(wm[i]*G*geo['z'][i]) for i in range(4)]
    E['waste_ke']  = [kJ(0.5*wm[i]*v[i+1]**2) for i in range(4)]
    h_out = max(3.0, geo['z'][3])
    E['out_ke'] = kJ(0.5*m[4]*v[4]**2); E['out_pe'] = kJ(m[4]*G*h_out)
    useful = (sum(E['valves'])+sum(E['friction'])+sum(E['bends'])+sum(E['waste_gpe'])+sum(E['waste_ke'])
              +E['out_ke']+E['out_pe'])
    eta = PUMP_EFF[pump_i]
    E_pump_in = useful/eta
    units_kWh = [F_KWH[gf], FIL_KWH[gfil], D_KWH[gd], H_KWH[gh]]
    E_units = sum(units_kWh)*3600.0
    E_in = E_pump_in + E_units
    Vout = Q[4]*SEC
    E_out = Vout*P['purity']*E_DENS*1000.0
    # ---- capital cost ----
    Q0 = Q[0]
    pp, rating = pump_price(pump_i, geo['L0'])
    if legacy_flow_cost:
        qf = [Q0, m[0]/rho[1], m[1]/rho[2], m[2]/rho[3]]
    else:
        qf = [Q[0], Q[1], Q[2], Q[3]]
    C = {}
    C['site'] = SITE_COST[site]
    C['mold'] = MOLD_COST if site == 3 else 0.0
    C['pump'] = pp*Q0
    C['fermenter'] = F_COST[gf]*qf[0]; C['filter'] = FIL_COST[gfil]*qf[1]
    C['distiller'] = D_COST[gd]*qf[2];  C['dehydrator'] = H_COST[gh]*qf[3]
    di = lambda d: DIAMS.index(round(d,2))
    C['pipes'] = sum(PIPE_COST[di(diam[k]), pipe_g[i]]*L for i,(k,L) in enumerate(Lpipes))
    C['duct'] = DUCT[duct_d]*geo['duct']
    vc = [VALVE_COST[di(diam[vstream[i]]), valve_g[i]] for i in range(8)]
    C['valves'] = sum(vc)
    bprice = (lambda d: legacy_bend_price) if legacy_bend_price else (lambda d: BEND_COST[di(d), bend_angle_idx])
    C['bends'] = sum(nb[i]*bprice(diam[i]) for i in range(4))
    cap = sum(C.values())
    return dict(P=P, Vin=Vin, useful=useful, E=E, E_pump_in=E_pump_in, E_units=E_units, E_in=E_in, E_out=E_out,
                ROI=E_out/E_in, purity=P['purity'], Vout=Vout, C=C, capital=cap, rating=rating, geo=geo, v=v, Q=Q,
                eta=eta, units_kWh=units_kWh)


# ---------------- search ----------------
import sys, time

def prune(cost, en):
    """indices of Pareto-optimal (min cost, min energy) points"""
    idx = np.lexsort((en, cost)); keep=[]; best=np.inf
    for i in idx:
        if en[i] < best - 1e-9:
            keep.append(i); best = en[i]
    return np.array(keep)

def seg_options(k, Q, rho, L, nb, waste_kg, m4, legacy_bend=None):
    """enumerate (cost, energy kJ/day) for stream k over diameter, pipe grade, valve grades"""
    nv = [1,2,2,2,1][k]
    has_pipe = L > 0
    out_c=[]; out_e=[]; desc=[]
    for di,D in enumerate(DIAMS):
        A = math.pi*D*D/4; v = Q/A; dyn = 0.5*rho*v*v
        base_kJ = dyn*Q*SEC/1000.0
        ke = 0.0
        if k>=1: ke += 0.5*waste_kg[k-1]*v*v*SEC/1000.0
        if k==4: ke += 0.5*m4*v*v*SEC/1000.0
        bend_e = nb*BEND_K[5]*base_kJ
        bend_c = nb*(legacy_bend if legacy_bend else BEND_COST[di,5])
        pgs = range(6) if has_pipe else [None]
        for pg in pgs:
            fr_e = PIPE_F[pg]*(L/D)*base_kJ if has_pipe else 0.0
            pc = PIPE_COST[di,pg]*L if has_pipe else 0.0
            for vg in itertools.product(range(4), repeat=nv):
                ve = sum(VALVE_K[g] for g in vg)*base_kJ
                vc = sum(VALVE_COST[di,g] for g in vg)
                out_c.append(pc+vc+bend_c); out_e.append(ke+bend_e+fr_e+ve); desc.append((di,pg,vg))
    return np.array(out_c), np.array(out_e), desc

def combine(front_a, front_b):
    ca, ea, ia = front_a; cb, eb, ib = front_b
    C = (ca[:,None]+cb[None,:]).ravel(); E = (ea[:,None]+eb[None,:]).ravel()
    k = prune(C,E)
    na, nb_ = len(ca), len(cb)
    ai = k//nb_; bi = k%nb_
    idx = np.concatenate([ia[ai], ib[bi].reshape(len(bi),-1)], axis=1)
    return C[k], E[k], idx

def search_site(site, legacy=False, verbose=True, purity_min=PURITY_MIN):
    t0=time.time(); cands=[]
    LAY = layouts(site)
    nlay = len(LAY)
    for gf,gfil,gd,gh in itertools.product(range(4),repeat=4):
        vpv, pur = vout_per_vin(F_EFF[gf],FIL_EFF[gfil],D_EFF[gd],H_EFF[gh])
        if pur < purity_min: continue
        Vin = V_OUT_TARGET/vpv
        P = process(Vin*SLURRY_DENS/SEC, F_EFF[gf],FIL_EFF[gfil],D_EFF[gd],H_EFF[gh])
        Q,rho,m = P['Q'],P['rho'],P['m']
        units_kWh = F_KWH[gf]+FIL_KWH[gfil]+D_KWH[gd]+H_KWH[gh]
        E_units = units_kWh*3600.0
        E_out = (Q[4]*SEC)*P['purity']*E_DENS*1000.0
        if legacy: qf=[Q[0], m[0]/rho[1], m[1]/rho[2], m[2]/rho[3]]
        else: qf=[Q[0],Q[1],Q[2],Q[3]]
        unit_cost = F_COST[gf]*qf[0]+FIL_COST[gfil]*qf[1]+D_COST[gd]*qf[2]+H_COST[gh]*qf[3]
        for li,combo in enumerate(LAY):
            geo = layout_geometry(site, combo)
            const_E = 0.0
            for i in range(4): const_E += WG(P['waste'][i],geo['z'][i])
            const_E += m[4]*G*max(3.0,geo['z'][3])*SEC/1000.0
            fixed = SITE_COST[site]+geo['mold']+unit_cost+DUCT[1.0]*geo['duct']
            Ls=[geo['L0']]+[l[0] for l in geo['links']]+[0.0]; Nb=[2]+[l[1] for l in geo['links']]+[0]
            fronts=[]
            for k in range(5):
                c,e,d = seg_options(k,Q[k],rho[k],Ls[k],Nb[k],P['waste'],m[4],100.0 if legacy else None)
                kp = prune(c,e); fronts.append((c[kp],e[kp],kp[:,None].copy()))
            cur = (fronts[0][0],fronts[0][1],np.arange(len(fronts[0][0]))[:,None])
            for k in range(1,5):
                cur = combine(cur,(fronts[k][0],fronts[k][1],np.arange(len(fronts[k][0]))[:,None]))
            Cs,Es,idx = cur
            capl=[];roil=[];pil=[];jl=[]
            for pi in range(5):
                pp,rating = pump_price(pi, geo['L0'])
                cap = fixed + pp*Q[0] + Cs
                Ein = (const_E+Es)/PUMP_EFF[pi] + E_units
                capl.append(cap); roil.append(E_out/Ein); pil.append(np.full(len(cap),pi)); jl.append(np.arange(len(cap)))
            cap=np.concatenate(capl); roi=np.concatenate(roil); pa=np.concatenate(pil); ja=np.concatenate(jl)
            order=np.lexsort((-roi,cap)); best=-1
            for o in order:
                if roi[o]>best+1e-12:
                    best=roi[o]; cands.append((cap[o],roi[o],site,li,(gf,gfil,gd,gh),int(pa[o]),tuple(idx[ja[o]]),Vin))
    if verbose: print(f'site {site}: {nlay} layouts, {len(cands)} pareto-candidates, {time.time()-t0:.1f}s')
    return cands, LAY

def WG(mass, z): return mass*G*z*SEC/1000.0
cands_front_store = {}

def global_front(cands):
    arr = np.array([(c[0],c[1]) for c in cands])
    idx = np.lexsort((-arr[:,1], arr[:,0])); keep=[]; best=-1
    for i in idx:
        if arr[i,1] > best+1e-12: keep.append(i); best=arr[i,1]
    return [cands[i] for i in keep]


# ---------------- decode ----------------

def decode(c, legacy=False):
    cap,roi,site,li,grades,pi,idx,Vin = c
    LAY = layouts(site); combo = LAY[li]; geo = layout_geometry(site, combo)
    gf,gfil,gd,gh = grades
    P = process(Vin*SLURRY_DENS/SEC, F_EFF[gf],FIL_EFF[gfil],D_EFF[gd],H_EFF[gh])
    Ls=[geo['L0']]+[l[0] for l in geo['links']]+[0.0]; Nb=[2]+[l[1] for l in geo['links']]+[0]
    diam=[];pipe_g=[0,0,0,0];valve_g=[0]*8
    vmap=[[0],[1,2],[3,4],[5,6],[7]]
    for k in range(5):
        c_,e_,d_ = seg_options(k,P['Q'][k],P['rho'][k],Ls[k],Nb[k],P['waste'],P['m'][4],100.0 if legacy else None)
        kp = prune(c_,e_); di,pg,vg = d_[kp[idx[k]]]
        diam.append(DIAMS[di])
        if k<4 and pg is not None: pipe_g[k]=pg
        for a,g in zip(vmap[k],vg): valve_g[a]=g
    r = evaluate(site, combo, grades, pi, Vin, diam, pipe_g, valve_g, legacy_flow_cost=legacy,
                 legacy_bend_price=100.0 if legacy else None)
    return r, dict(diam=diam,pipe_g=pipe_g,valve_g=valve_g,pi=pi,grades=grades,site=site,combo=combo,Vin=Vin)


# =====================================================================================
#  MAIN : (1) validate against course Validation Model, (2) exhaustive search, (3) report
# =====================================================================================
def validate():
    geo = dict(z=[12.192, 6.096, 0.0, 0.0], z_ft=[40,20,0,0], bl=['T']*4, links=[(0,0),(0,0),(10*FT,2)],
               L0=50*FT, duct=10*FT, mold=12000.0, names=[])
    r = evaluate(3, None, (0,0,0,0), 0, 3161.58, [0.1]*5, [3,3,3,3], [0]*8, legacy_flow_cost=True,
                 legacy_bend_price=100.0, layout_override=geo)
    print('VALIDATION vs course sheet (Site 3, scrap equipment, 835,200 gal/day input)')
    print(f"  capital cost {r['capital']:>14,.2f}   sheet 1,866,773.74")
    print(f"  energy ROI   {r['ROI']:>14.4f}   sheet 5.5459")
    print(f"  purity       {r['purity']:>14.5f}   sheet 0.87214")
    print(f"  output m3/d  {r['Vout']:>14.3f}   sheet 260.059")

def run(sites=(1,2,3)):
    results = {}
    for s in sites:
        cands, LAY = search_site(s)
        results[s] = global_front(cands)
    return results

def describe(r, d):
    print(f"  capital ${r['capital']:,.0f} | ROI {r['ROI']:.3f} | purity {r['purity']:.2%} | output {r['Vout']:.2f} m3/d")
    print(f"  slurry input {r['Vin']:,.1f} m3/d = {r['Vin']*264.172:,.0f} gal/d = {r['Vin']*SLURRY_DENS/1000:,.0f} t/d")
    print('  units placed at:', r['geo']['names'])
    print('  equipment:', dict(zip(['fermenter','filter','distiller','dehydrator'],[LINE4[g] for g in d['grades']])),
          '| pump', PUMP_COLS[d['pi']], f"({r['eta']}) rating {r['rating']} m")
    print('  stream diameters (m):', d['diam'], '| pipe grades', [PIPE_GRADES[g] for g in d['pipe_g']])
    print('  valve grades:', [VALVE_GRADES[g] for g in d['valve_g']])
    print('  cost breakdown:', {k: round(v) for k, v in r['C'].items()})


# =====================================================================================
#  INTERACTIVE WHAT-IF MODE : type inputs (Enter = default), get results
# =====================================================================================
def ask(prompt, default, cast=str, choices=None):
    while True:
        raw = input(f"{prompt} [{default}]: ").strip()
        if raw == '': raw = str(default)
        try:
            v = cast(raw)
            if choices is not None and v not in choices:
                print('   choose one of:', choices); continue
            return v
        except ValueError:
            print('   invalid entry, try again')

def grade_menu(title, names, default):
    print(f"{title}: " + ', '.join(f"{i+1}={n}" for i, n in enumerate(names)))
    return ask('   choice', default + 1, int, list(range(1, len(names)+1))) - 1

def interactive():
    print('\n=== WHAT-IF DESIGN EVALUATOR (press Enter to accept the default shown) ===')
    print('Defaults = recommended design (Site 3, ROI/$-optimal).\n')
    site = ask('Site (1, 2 or 3)', 3, int, [1, 2, 3])
    LAY = layouts(site)
    print('Feasible layouts (Fermenter, Filter, Distiller, Dehydrator):')
    default_li = 0
    for i, c in enumerate(LAY):
        g = layout_geometry(site, c)
        print(f"  {i+1:2d}. {g['names']}")
        if site == 3 and g['names'] == ['Tall-L1','Tall-L0','Tall-L0','Small-L0']: default_li = i
        if site == 2 and g['names'] == ['Ledge','Ground','Ground','Ground']: default_li = i
        if site == 1 and g['names'] == ['A-L1','B-L1','B-L0','A-L0']: default_li = i
    li = ask('Layout number', default_li + 1, int, list(range(1, len(LAY)+1))) - 1
    combo = LAY[li]

    print('\n-- Unit operations (1=Scrap, 2=Average, 3=Premium, 4=World-class) --')
    gf   = grade_menu('Fermenter ', LINE4, 2)
    gfil = grade_menu('Filter    ', LINE4, 2)
    gd   = grade_menu('Distiller ', LINE4, 1)
    gh   = grade_menu('Dehydrator', LINE4, 3)
    pi   = grade_menu('Pump grade', PUMP_COLS, 3 if site != 2 else 2)

    print('\n-- Slurry input --')
    vpv, pur = vout_per_vin(F_EFF[gf], FIL_EFF[gfil], D_EFF[gd], H_EFF[gh])
    auto = V_OUT_TARGET / vpv * 264.172052
    gal = ask('Slurry input, gal/day (type 0 to auto-size for exactly 100,000 gal/day output)', 0, float)
    Vin = (auto/264.172052) if gal <= 0 else gal/264.172052
    if gal <= 0: print(f'   auto-sized input: {auto:,.0f} gal/day')

    print('\n-- Piping (diameter m: 0.10-0.15) for the 5 stream segments --')
    names = ['Pump -> Fermenter', 'Fermenter -> Filter', 'Filter -> Distiller', 'Distiller -> Dehydrator', 'Dehydrator -> Output']
    ddef = [0.11, 0.15, 0.12, 0.10, 0.11]; pdef = [0, 0, 0, 0]; vdef = [3, 1, 3, 0, 0]
    geo = layout_geometry(site, combo)
    diam = []; pipe_g = []; valve_g = []
    vcount = [1, 2, 2, 2, 1]
    for k in range(5):
        print(f' Segment {k+1}: {names[k]}')
        diam.append(ask('   diameter (m)', ddef[k], float, DIAMS))
        if k < 4:
            L = [geo['L0']] + [l[0] for l in geo['links']]
            if L[k] > 0:
                pg = grade_menu(f'   pipe grade ({L[k]:.2f} m long)', PIPE_GRADES, pdef[k])
            else:
                pg = 0; print('   (no pipe in this segment for the chosen layout)')
            pipe_g.append(pg)
        vg = grade_menu(f'   valve grade (applies to {vcount[k]} valve(s))', VALVE_GRADES, vdef[k])
        valve_g += [vg] * vcount[k]
    duct_d = ask('\nCO2 duct diameter (m: 1.0, 1.25, 1.5)', 1.0, float, [1.0, 1.25, 1.5])
    legacy = ask('Price like the course validation sheet? (y/n)', 'n', str, ['y', 'n']) == 'y'

    r = evaluate(site, combo, (gf, gfil, gd, gh), pi, Vin, diam, pipe_g, valve_g, duct_d=duct_d,
                 legacy_flow_cost=legacy, legacy_bend_price=100.0 if legacy else None)
    report(r, site)
    return r

def report(r, site):
    P = r['P']; ok = lambda b: 'PASS' if b else 'FAIL'
    print('\n' + '=' * 70 + '\nRESULTS\n' + '=' * 70)
    print(f"Layout: {r['geo']['names']}   (pump height rating {r['rating']} m)")
    print(f"\nAGRICULTURAL SLURRY NEEDED: {r['Vin']*264.172052:,.0f} gal/day = {r['Vin']:,.1f} m3/day = {r['Vin']*SLURRY_DENS/1000:,.0f} t/day")
    print(f"                            ~{r['Vin']*264.172052/11600:,.1f} tanker loads/day (11,600 gal each)")
    print('\nMass / flow balance:')
    print('  stream                 kg/s    rho kg/m3   m3/s      m/s')
    for n, m, rho, q, v in zip(['slurry (pump)', 'fermenter out', 'filter out', 'distiller out', 'product out'],
                               P['m'], P['rho'], r['Q'], r['v']):
        print(f'  {n:<20s}{m:9.3f}{rho:11.1f}{q:11.5f}{v:9.2f}')
    print('  waste kg/s: CO2 %.3f | fiber %.3f | distiller %.3f | water %.3f' % tuple(P['waste']))
    print('\nEnergy (kJ/day):')
    print(f"  pump input            {r['E_pump_in']:>16,.0f}   (hydraulic losses {r['useful']:,.0f} / eff {r['eta']})")
    print(f"  unit operations       {r['E_units']:>16,.0f}")
    print(f"  TOTAL IN              {r['E_in']:>16,.0f}")
    print(f"  ethanol energy OUT    {r['E_out']:>16,.0f}")
    print('\nCapital cost:')
    for k, v in r['C'].items(): print(f'  {k:<12s}${v:>13,.2f}')
    print(f"  {'TOTAL':<12s}${r['capital']:>13,.2f}")
    print('\nSUMMARY vs requirements')
    print(f"  Output volume : {r['Vout']*264.172052:>10,.0f} gal/day  (need >= 100,000)  {ok(r['Vout'] >= V_OUT_TARGET*0.9999)}")
    print(f"  Purity        : {r['purity']:>10.2%}      (need >= 98%)       {ok(r['purity'] >= 0.98)}")
    print(f"  Energy ROI    : {r['ROI']:>10.3f}      (target 10)         {ok(r['ROI'] >= 10)}")
    print(f"  Capital cost  : ${r['capital']:,.0f}")

def menu():
    while True:
        print('\n==== ETHANOL PLANT MODEL ====')
        print(' 1) Evaluate a design (what-if: enter your own inputs)')
        print(' 2) Run the exhaustive optimizer (about 1 minute)')
        print(' 3) Validate against the course validation sheet')
        print(' 4) Quit')
        c = input('Choice [1]: ').strip() or '1'
        if c == '1': interactive()
        elif c == '2':
            res = run()
            for s in (1, 2, 3):
                front = res[s]; k = int(np.argmax([f[1]/f[0] for f in front]))
                r, d = decode(front[k]); print(f'\nSITE {s}: best ROI-per-dollar design'); describe(r, d)
        elif c == '3': validate()
        elif c == '4': break

if __name__ == '__main__':
    menu()