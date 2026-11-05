// engine.js — thermodynamic state and process engine built on CoolProp (WebAssembly)
// All values inside this module are SI: K, Pa, J/kg, J/kg/K, kg/m3, kg/s, W.
// Johann Cruz · engcalc

export const FLASH_PAIRS = {
  TP:  { label: 'T – P',            q: ['T', 'P'] },
  PH:  { label: 'P – H',            q: ['P', 'H'] },
  PS:  { label: 'P – S',            q: ['P', 'S'] },
  HS:  { label: 'H – S',            q: ['H', 'S'] },
  TS:  { label: 'T – S',            q: ['T', 'S'] },
  PQ:  { label: 'P – Q (saturated / two-phase)', q: ['P', 'Q'] },
  TQ:  { label: 'T – Q (saturated / two-phase)', q: ['T', 'Q'] },
  PdT: { label: 'P – superheat / subcool ΔT',    q: ['P', 'dT'] },
  DT:  { label: 'D – T',            q: ['D', 'T'] },
  DP:  { label: 'D – P',            q: ['D', 'P'] },
};

export const PROCESS_TYPES = {
  comp:  { label: 'Compressor',            needs: ['P2', 'eta'] },
  pump:  { label: 'Pump',                  needs: ['P2', 'eta'] },
  turb:  { label: 'Turbine / expander',    needs: ['P2', 'eta'] },
  valve: { label: 'Valve / choke (isenthalpic)', needs: ['P2'] },
  hx:    { label: 'Heater / cooler',       needs: ['spec', 'dP'] },
};

export const HX_SPECS = {
  T:    'Outlet temperature',
  Q:    'Outlet quality Q',
  dT:   'Outlet superheat (+) / subcool (−) ΔT',
  duty: 'Duty (+ heating, − cooling)',
};

const PHASE_NAMES = {
  0: 'Liquid', 1: 'Supercritical', 2: 'Supercritical gas', 3: 'Supercritical liquid',
  4: 'Critical point', 5: 'Vapour', 6: 'Two-phase', 7: 'Unknown', 8: 'Unknown',
};
export const PHASE_CLASS = {
  0: 'liq', 1: 'sc', 2: 'vap', 3: 'liq', 4: 'sc', 5: 'vap', 6: 'tp', 7: 'unk', 8: 'unk',
};

export class Engine {
  constructor(cp) {
    this.cp = cp;
    this.ip = cp.input_pairs;
    this.par = cp.parameters;
    this.cache = {};
  }

  fluids() {
    return this.cp.get_global_param_string('FluidsList').split(',').sort((a, b) => a.localeCompare(b));
  }

  // One AbstractState for flashes and one for side calculations (saturation), per fluid.
  _states(fluid) {
    let c = this.cache[fluid];
    if (!c) {
      const main = this.cp.factory('HEOS', fluid);
      const side = this.cp.factory('HEOS', fluid);
      const k = (key) => { try { return main.keyed_output(this.par[key]); } catch (e) { return NaN; } };
      let pure = true;
      try { pure = this.cp.get_fluid_param_string(fluid, 'pure') === 'true'; } catch (e) {}
      c = {
        main, side, pure,
        Tc: k('iT_critical'), Pc: k('iP_critical'), Dc: k('irhomass_critical'),
        Tt: k('iT_triple'), Pt: k('iP_triple'),
        Tmin: k('iT_min'), Tmax: k('iT_max'), Pmax: k('iP_max'),
        M: k('imolar_mass'),
      };
      this.cache[fluid] = c;
    }
    return c;
  }

  info(fluid) {
    const c = this._states(fluid);
    return { Tc: c.Tc, Pc: c.Pc, Tt: c.Tt, Pt: c.Pt, Tmin: c.Tmin, Tmax: c.Tmax, Pmax: c.Pmax, M: c.M, pure: c.pure };
  }

  _resetMain(fluid) {
    const c = this.cache[fluid];
    try { c.main.delete(); } catch (e) {}
    c.main = this.cp.factory('HEOS', fluid);
  }

  // Saturation temperature at pressure P (NaN if not below the critical pressure)
  Tsat(fluid, P) {
    const c = this._states(fluid);
    if (!(P > 0) || !(P < c.Pc)) return NaN;
    try { c.side.update(this.ip.PQ_INPUTS, P, 0); return c.side.T(); }
    catch (e) { return NaN; }
  }

  Psat(fluid, T) {
    const c = this._states(fluid);
    if (!(T < c.Tc) || T < c.Tmin) return NaN;
    try { c.side.update(this.ip.QT_INPUTS, 0, T); return c.side.p(); }
    catch (e) { return NaN; }
  }

  // Throws a FlashError with a plain-language message.
  _check(fluid, q, v, fmt) {
    const c = this._states(fluid);
    if (!Number.isFinite(v)) throw new FlashError(`${qName(q)} is missing or not a number.`);
    if (q === 'T') {
      if (v <= 0) throw new FlashError('Temperature must be above absolute zero.');
      if (v < c.Tmin) throw new FlashError(`Temperature ${fmt('T', v)} is below the lowest temperature the ${fluid} equation of state covers (${fmt('T', c.Tmin)}).`);
      if (v > c.Tmax * 1.5) throw new FlashError(`Temperature ${fmt('T', v)} is well above the upper limit of the ${fluid} equation of state (${fmt('T', c.Tmax)}).`);
    }
    if (q === 'P') {
      if (v <= 0) throw new FlashError('Absolute pressure must be above zero. If you use a gauge unit, check that the value is above −1 atm.');
      if (v > c.Pmax) throw new FlashError(`Pressure ${fmt('P', v)} is above the upper limit of the ${fluid} equation of state (${fmt('P', c.Pmax)}).`);
    }
    if (q === 'Q' && (v < 0 || v > 1)) throw new FlashError('Quality Q must be between 0 (saturated liquid) and 1 (saturated vapour).');
    if (q === 'D' && v <= 0) throw new FlashError('Density must be above zero.');
  }

  // Flash a state from a pair key and two SI values. Returns a props object.
  flash(fluid, pair, a, b, fmt = defaultFmt) {
    const c = this._states(fluid);
    const [qa, qb] = FLASH_PAIRS[pair].q;
    this._check(fluid, qa, a, fmt);
    if (qb !== 'dT') this._check(fluid, qb, b, fmt);
    const ip = this.ip;
    let call;
    switch (pair) {
      case 'TP': {
        const Ts = this.Tsat(fluid, b);
        if (Number.isFinite(Ts) && Math.abs(a - Ts) < Math.max(1e-6 * Ts, 1e-4)) {
          throw new FlashError(`T and P are on the saturation line (Tsat at this pressure is ${fmt('T', Ts)}). T and P alone do not fix a two-phase state. Use P – Q or T – Q, or move T slightly off saturation.`);
        }
        call = [ip.PT_INPUTS, b, a]; break;
      }
      case 'PH': call = [ip.HmassP_INPUTS, b, a]; break;
      case 'PS': call = [ip.PSmass_INPUTS, a, b]; break;
      case 'HS': call = [ip.HmassSmass_INPUTS, a, b]; break;
      case 'TS': call = [ip.SmassT_INPUTS, b, a]; break;
      case 'PQ':
        if (a >= c.Pc) throw new FlashError(`Pressure ${fmt('P', a)} is at or above the critical pressure (${fmt('P', c.Pc)}). There is no liquid–vapour split above the critical pressure, so quality cannot be set.`);
        if (a < c.Pt * 0.999) throw new FlashError(`Pressure ${fmt('P', a)} is below the triple-point pressure (${fmt('P', c.Pt)}). No liquid–vapour equilibrium exists there.`);
        call = [ip.PQ_INPUTS, a, b]; break;
      case 'TQ':
        if (b >= c.Tc) throw new FlashError(`Temperature ${fmt('T', b)} is at or above the critical temperature (${fmt('T', c.Tc)}). Quality cannot be set above the critical temperature.`);
        call = [ip.QT_INPUTS, b, a]; break;
      case 'PdT': {
        if (a >= c.Pc) throw new FlashError(`Pressure ${fmt('P', a)} is at or above the critical pressure (${fmt('P', c.Pc)}), so there is no saturation temperature to measure superheat or subcooling from.`);
        if (!Number.isFinite(b)) throw new FlashError('ΔT is missing or not a number.');
        if (b === 0) throw new FlashError('ΔT = 0 is exactly on the saturation line, so the state is not fixed. Use P – Q with Q = 0 (saturated liquid) or Q = 1 (saturated vapour).');
        const Ts = this.Tsat(fluid, a);
        if (!Number.isFinite(Ts)) throw new FlashError('Could not find the saturation temperature at this pressure.');
        this._check(fluid, 'T', Ts + b, fmt);
        call = [ip.PT_INPUTS, a, Ts + b]; break;
      }
      case 'DT': call = [ip.DmassT_INPUTS, a, b]; break;
      case 'DP': call = [ip.DmassP_INPUTS, a, b]; break;
      default: throw new FlashError('Unknown input pair.');
    }
    return this._update(fluid, call, fmt);
  }

  _update(fluid, call, fmt) {
    const c = this._states(fluid);
    try {
      c.main.update(call[0], call[1], call[2]);
    } catch (e) {
      this._resetMain(fluid);
      throw new FlashError(translateError(e, fmt));
    }
    return this._read(fluid);
  }

  flashPH(fluid, P, H, fmt = defaultFmt) {
    this._check(fluid, 'P', P, fmt);
    return this._update(fluid, [this.ip.HmassP_INPUTS, H, P], fmt);
  }

  flashPS(fluid, P, S, fmt = defaultFmt) {
    this._check(fluid, 'P', P, fmt);
    return this._update(fluid, [this.ip.PSmass_INPUTS, P, S], fmt);
  }

  _read(fluid) {
    const c = this._states(fluid);
    const s = c.main;
    const g = (f) => { try { const v = f(); return Number.isFinite(v) ? v : NaN; } catch (e) { return NaN; } };
    let ph = 7;
    try { ph = s.phase().value; } catch (e) {}
    const T = s.T(), P = s.p();
    const two = ph === 6;
    const p = {
      T, P,
      D: g(() => s.rhomass()),
      H: g(() => s.hmass()),
      S: g(() => s.smass()),
      U: g(() => s.umass()),
      Q: two ? g(() => s.Q()) : NaN,
      phase: ph,
      Cp: two ? NaN : g(() => s.cpmass()),
      Cv: two ? NaN : g(() => s.cvmass()),
      Z: g(() => s.compressibility_factor()),
      mu: two ? NaN : g(() => s.viscosity()),
      k: two ? NaN : g(() => s.conductivity()),
      a: two ? NaN : g(() => s.speed_sound()),
      sigma: two ? g(() => s.surface_tension()) : NaN,
    };
    p.gamma = p.Cp / p.Cv;
    p.Tsat = this.Tsat(fluid, P);
    p.margin = two ? 0 : (Number.isFinite(p.Tsat) ? T - p.Tsat : NaN);
    // Refine single-phase names below the critical pressure
    p.phaseName = PHASE_NAMES[ph] || 'Unknown';
    if (ph === 0 && Number.isFinite(p.Tsat)) p.phaseName = 'Subcooled liquid';
    if (ph === 5 && Number.isFinite(p.Tsat)) p.phaseName = 'Superheated vapour';
    if (ph === 2) p.phaseName = 'Gas (T > Tc)';
    if (ph === 1) p.phaseName = 'Supercritical (T > Tc, P > Pc)';
    if (ph === 3) p.phaseName = 'Compressed liquid (P > Pc)';
    if (ph === 6) p.phaseName = p.Q < 1e-9 ? 'Saturated liquid' : (p.Q > 1 - 1e-9 ? 'Saturated vapour' : 'Two-phase');
    p.phaseClass = PHASE_CLASS[ph] || 'unk';
    return p;
  }

  // Process step from inlet props. Returns outlet props plus step data.
  process(fluid, inlet, proc, mdot, fmt = defaultFmt) {
    const t = proc.type;
    const warn = [];
    let out, h2s = NaN;
    if (t === 'comp' || t === 'pump' || t === 'turb') {
      const P2 = proc.P2, eta = proc.eta;
      if (!(eta > 0 && eta <= 1)) throw new FlashError('Isentropic efficiency must be greater than 0 and not more than 1 (enter it as a fraction, e.g. 0.75).');
      if (t !== 'turb' && !(P2 > inlet.P)) throw new FlashError(`Outlet pressure must be higher than the inlet pressure (${fmt('P', inlet.P)}) for a ${t === 'comp' ? 'compressor' : 'pump'}.`);
      if (t === 'turb' && !(P2 < inlet.P)) throw new FlashError(`Outlet pressure must be lower than the inlet pressure (${fmt('P', inlet.P)}) for a turbine or expander.`);
      const iso = this.flashPS(fluid, P2, inlet.S, fmt);
      h2s = iso.H;
      const h2 = t === 'turb' ? inlet.H - eta * (inlet.H - h2s) : inlet.H + (h2s - inlet.H) / eta;
      out = this.flashPH(fluid, P2, h2, fmt);
      out.T2s = iso.T;
      if (t === 'comp' && (inlet.phase === 0 || inlet.phase === 3 || (inlet.phase === 6 && inlet.Q < 1 - 1e-9))) warn.push('Liquid at compressor inlet.');
      if (t === 'pump' && (inlet.phase === 5 || inlet.phase === 2 || (inlet.phase === 6 && inlet.Q > 1e-9))) warn.push('Vapour at pump inlet.');
      if (t === 'turb' && out.phase === 6 && out.Q < 0.88) warn.push(`Outlet quality ${out.Q.toFixed(3)} is below about 0.88; blade erosion risk for steam turbines.`);
    } else if (t === 'valve') {
      if (!(proc.P2 <= inlet.P)) throw new FlashError(`Outlet pressure must not be higher than the inlet pressure (${fmt('P', inlet.P)}) across a valve.`);
      out = this.flashPH(fluid, proc.P2, inlet.H, fmt);
      if (out.phase === 6 && inlet.phase !== 6) warn.push('Fluid flashes to two-phase across the valve.');
    } else if (t === 'hx') {
      const P2 = inlet.P - (proc.dP || 0);
      if (!(P2 > 0)) throw new FlashError('Pressure drop is larger than the inlet pressure.');
      const sp = proc.spec;
      if (sp === 'T') out = this.flash(fluid, 'TP', proc.specVal, P2, fmt);
      else if (sp === 'Q') out = this.flash(fluid, 'PQ', P2, proc.specVal, fmt);
      else if (sp === 'dT') out = this.flash(fluid, 'PdT', P2, proc.specVal, fmt);
      else if (sp === 'duty') {
        if (!(mdot > 0)) throw new FlashError('A duty specification needs a mass flow above zero.');
        out = this.flashPH(fluid, P2, inlet.H + proc.specVal / mdot, fmt);
      } else throw new FlashError('Unknown heater / cooler specification.');
    } else {
      throw new FlashError('Unknown equipment type.');
    }
    const dh = out.H - inlet.H;
    return { out, dh, h2s, warn };
  }

  // Points along a process between two states, as a straight line in (log P, h).
  path(fluid, a, b, n = 40) {
    const pts = [];
    const lpa = Math.log(a.P), lpb = Math.log(b.P);
    for (let i = 0; i <= n; i++) {
      const f = i / n;
      const P = Math.exp(lpa + (lpb - lpa) * f);
      const H = a.H + (b.H - a.H) * f;
      if (i === 0) { pts.push({ T: a.T, P: a.P, H: a.H, S: a.S }); continue; }
      if (i === n) { pts.push({ T: b.T, P: b.P, H: b.H, S: b.S }); continue; }
      try {
        const c = this._states(fluid);
        c.side.update(this.ip.HmassP_INPUTS, H, P);
        pts.push({ T: c.side.T(), P, H, S: c.side.smass() });
      } catch (e) {
        this._resetSide(fluid);
      }
    }
    return pts;
  }

  _resetSide(fluid) {
    const c = this.cache[fluid];
    try { c.side.delete(); } catch (e) {}
    c.side = this.cp.factory('HEOS', fluid);
  }

  // Saturation dome: liquid and vapour lines from the triple point (or Tmin) to near the critical point.
  dome(fluid, n = 90) {
    const c = this._states(fluid);
    const T0 = Math.max(c.Tt, c.Tmin), Tc = c.Tc;
    const liq = [], vap = [];
    if (!(Tc > T0)) return { liq, vap, crit: null };
    for (let i = 0; i <= n; i++) {
      // cluster points near the critical point
      const x = i / n;
      const T = Tc - (Tc - T0) * Math.pow(1 - x, 2.2) - (i === n ? (Tc - T0) * 2e-5 : 0);
      for (const [Q, arr] of [[0, liq], [1, vap]]) {
        try {
          c.side.update(this.ip.QT_INPUTS, Q, T);
          arr.push({ T, P: c.side.p(), H: c.side.hmass(), S: c.side.smass() });
        } catch (e) { this._resetSide(fluid); }
      }
    }
    let crit = null;
    try {
      c.side.update(this.ip.DmassT_INPUTS, c.Dc, c.Tc);
      crit = { T: c.Tc, P: c.side.p(), H: c.side.hmass(), S: c.side.smass() };
    } catch (e) { this._resetSide(fluid); }
    return { liq, vap, crit };
  }
}

export class FlashError extends Error {}

function qName(q) {
  return { T: 'Temperature', P: 'Pressure', H: 'Enthalpy', S: 'Entropy', D: 'Density', Q: 'Quality', dT: 'ΔT' }[q] || q;
}

function defaultFmt(q, v) { return `${+v.toPrecision(6)} (SI)`; }

// Turn a CoolProp exception into a short plain-language message.
export function translateError(e, fmt) {
  let m = (e && e.message) ? e.message : String(e);
  const i = m.indexOf(']: ');
  if (i >= 0) m = m.slice(i + 3);
  m = m.replace(/\s*\(Operation:.*$/s, '').replace(/\s*\[Context:.*$/s, '').trim();
  if (/within 1e-4 % of given p/i.test(m) || /saturation pressure.*corresponding/i.test(m)) {
    return 'T and P are on the saturation line. T and P alone do not fix a two-phase state. Use P – Q or T – Q.';
  }
  if (/below Tmelt/i.test(m)) return 'Temperature is below the melting line at this pressure (solid region). CoolProp covers liquid, vapour and supercritical states only.';
  if (/HSU_P_flash|PHSU_D_flash|solver|did not converge|not able to|failed|Brent|Secant/i.test(m)) {
    return `CoolProp could not converge on this state. The inputs may be outside the fluid's valid range, or very close to the critical point. Detail: ${m}`;
  }
  return m;
}
