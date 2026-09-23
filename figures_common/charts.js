// Shared SVG chart helpers for the ExoRT figure pages (inlined by
// tools/build_figure_pages.py). No dependencies. Colors come from CSS tokens.
const Charts = (() => {
  const NS = "http://www.w3.org/2000/svg";
  const W = 380, H = 250, M = { l: 54, r: 12, t: 14, b: 40 };
  const PW = W - M.l - M.r, PH = H - M.t - M.b;

  function el(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function num(v) {
    const a = Math.abs(v);
    if (a === 0) return "0";
    if (a >= 100) return v.toFixed(1);
    if (a >= 1) return v.toFixed(2);
    if (a >= 0.01) return v.toFixed(3);
    return v.toExponential(2);
  }
  function signed(v, d = 2) { return (v >= 0 ? "+" : "−") + Math.abs(v).toFixed(d); }
  function niceTicks(lo, hi, n = 5) {
    const span = hi - lo || Math.abs(hi) || 1, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag;
    const out = [];
    for (let v = Math.floor(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toFixed(10));
    if (out[out.length - 1] < hi) out.push(+(out[out.length - 1] + step).toFixed(10));  // scale covers the data
    return out;
  }
  function card(parent) {
    const c = document.createElement("div"); c.className = "card chart"; parent.appendChild(c); return c;
  }
  function frame(c, title, sub) {
    const h = document.createElement("h3"); h.textContent = title; c.appendChild(h);
    if (sub) { const s = document.createElement("div"); s.className = "sub"; s.textContent = sub; c.appendChild(s); }
    const svg = el("svg", { class: "plot", viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": title });
    c.appendChild(svg);
    const tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true; c.appendChild(tip);
    return { svg, tip };
  }
  function swatch(s, w = 14) {
    return `<svg width="${w}" height="8" viewBox="0 0 ${w} 8"><line x1="0" y1="4" x2="${w}" y2="4" stroke="${s.color}" stroke-width="2.2" ${s.dash ? `stroke-dasharray="${w > 20 ? "5 4" : "3 2"}"` : ""}/>${s.marker ? `<circle cx="${w / 2}" cy="4" r="3" fill="${s.color}"/>` : ""}</svg>`;
  }
  function showTip(c, tip, svg, px, py, title, rows) {
    tip.innerHTML = `<div class="t">${title}</div>` + rows.map(r =>
      `<div class="r">${swatch(r)}<span>${r.name}</span><span class="v">${r.value}</span></div>`).join("");
    tip.hidden = false;
    const box = svg.getBoundingClientRect(), cb = c.getBoundingClientRect();
    const x = box.left - cb.left + px * (box.width / W), y = box.top - cb.top + py * (box.height / H);
    const tw = tip.offsetWidth;
    tip.style.left = (x + 14 + tw > cb.width ? Math.max(4, x - tw - 14) : x + 14) + "px";
    tip.style.top = Math.max(4, y - 20) + "px";
  }
  function axes(svg, xlabel, ylabel) {
    el("text", { class: "axlabel", x: M.l + PW / 2, y: H - 4, "text-anchor": "middle" }, svg).textContent = xlabel;
    el("text", { class: "axlabel", x: 12, y: M.t + PH / 2, "text-anchor": "middle", transform: `rotate(-90 12 ${M.t + PH / 2})` }, svg).textContent = ylabel;
  }
  function yLinear(svg, lo, hi) {
    const yt = niceTicks(lo, hi, 5), y0 = yt[0], y1 = yt[yt.length - 1];
    const fy = v => M.t + PH * (1 - (v - y0) / (y1 - y0));
    for (const v of yt) {
      el("line", { class: "gridl", x1: M.l, x2: M.l + PW, y1: fy(v), y2: fy(v) }, svg);
      el("text", { x: M.l - 6, y: fy(v) + 3.5, "text-anchor": "end" }, svg).textContent = +v.toPrecision(3);
    }
    if (y0 < 0 && y1 > 0) el("line", { class: "axis", x1: M.l, x2: M.l + PW, y1: fy(0), y2: fy(0) }, svg);
    return fy;
  }
  function yLog(svg, lo, hi) {
    const d0 = Math.floor(Math.log10(lo)), d1 = Math.ceil(Math.log10(hi));
    const fy = v => M.t + PH * (1 - (Math.log10(v) - d0) / (d1 - d0));
    for (let d = d0; d <= d1; d++) {
      el("line", { class: "gridl", x1: M.l, x2: M.l + PW, y1: fy(10 ** d), y2: fy(10 ** d) }, svg);
      el("text", { x: M.l - 6, y: fy(10 ** d) + 3.5, "text-anchor": "end" }, svg).textContent = `1e${d}`;
    }
    return fy;
  }
  function order(cfg) { return cfg.drawOrder || cfg.series.map((_, i) => i); }
  function path(svg, d, s, width = 1.8) {
    el("path", { d, fill: "none", stroke: s.color, "stroke-width": width, "stroke-dasharray": s.dash ? "5 4" : "none", "stroke-linejoin": "round" }, svg);
  }

  // Per-band step chart. series[i].v is indexed by band; edges may run either way.
  function spectrum(c, cfg) {
    const { svg, tip } = frame(c, cfg.title, cfg.sub);
    const [x0, x1] = cfg.xdom, log = !!cfg.log;
    const fx = log ? v => M.l + PW * (Math.log(v) - Math.log(x0)) / (Math.log(x1) - Math.log(x0))
                   : v => M.l + PW * (v - x0) / (x1 - x0);
    const bands = [];
    for (let b = 0; b < cfg.edges.length - 1; b++) {
      const a = Math.min(cfg.edges[b], cfg.edges[b + 1]), z = Math.max(cfg.edges[b], cfg.edges[b + 1]);
      if (z > x0 && a < x1) bands.push({ b, a: Math.max(a, x0), z: Math.min(z, x1), ra: a, rz: z });
    }
    bands.sort((p, q) => p.a - q.a);
    let lo = 0, hi = 0;
    for (const s of cfg.series) for (const bd of bands) { lo = Math.min(lo, s.v[bd.b]); hi = Math.max(hi, s.v[bd.b]); }
    const fy = yLinear(svg, lo, hi * (hi > 0 ? 1.04 : 1));
    for (const v of cfg.xticks) {
      el("line", { class: "gridl", x1: fx(v), x2: fx(v), y1: M.t, y2: M.t + PH }, svg);
      el("text", { x: fx(v), y: M.t + PH + 15, "text-anchor": "middle" }, svg).textContent = v;
    }
    el("line", { class: "axis", x1: M.l, x2: M.l + PW, y1: M.t + PH, y2: M.t + PH }, svg);
    axes(svg, cfg.xlabel, cfg.ylabel);
    if (cfg.annot) {
      const [v, lab] = cfg.annot, x = fx(v);
      el("line", { class: "annot", x1: x, x2: x, y1: M.t + 14, y2: M.t + PH }, svg);
      el("text", { class: "annot-t", x: x + 4, y: M.t + 10, "text-anchor": x > M.l + PW * 0.7 ? "end" : "start" }, svg).textContent = lab;
    }
    const hl = el("rect", { class: "hl", y: M.t, height: PH, width: 0, x: 0 }, svg);
    for (const i of order(cfg)) {
      const s = cfg.series[i];
      let d = "";
      bands.forEach((bd, k) => {
        const y = fy(s.v[bd.b]);
        d += (k ? "L" : "M") + `${fx(bd.a).toFixed(2)},${y.toFixed(2)}L${fx(bd.z).toFixed(2)},${y.toFixed(2)}`;
      });
      path(svg, d, s);
    }
    const hit = el("rect", { x: M.l, y: M.t, width: PW, height: PH, fill: "transparent" }, svg);
    const inv = log ? px => Math.exp(Math.log(x0) + (px - M.l) / PW * (Math.log(x1) - Math.log(x0)))
                    : px => x0 + (px - M.l) / PW * (x1 - x0);
    hit.addEventListener("pointermove", ev => {
      const r = svg.getBoundingClientRect(), px = (ev.clientX - r.left) * W / r.width, xv = inv(px);
      const bd = bands.find(q => xv >= q.a && xv <= q.z);
      if (!bd) { hl.setAttribute("width", 0); tip.hidden = true; return; }
      hl.setAttribute("x", fx(bd.a)); hl.setAttribute("width", Math.max(1, fx(bd.z) - fx(bd.a)));
      showTip(c, tip, svg, px, M.t + 20, cfg.bandLabel(bd),
        cfg.series.map(s => ({ ...s, value: (cfg.fmtV || num)(s.v[bd.b]) })));
    });
    hit.addEventListener("pointerleave", () => { hl.setAttribute("width", 0); tip.hidden = true; });
  }

  // Vertical profile against log pressure (top of atmosphere at the top).
  function profile(c, cfg) {
    const { svg, tip } = frame(c, cfg.title, cfg.sub);
    const p = cfg.p, pmin = Math.min(...p), pmax = Math.max(...p);
    const d0 = Math.floor(Math.log10(pmin)), d1 = Math.ceil(Math.log10(pmax));
    const fy = v => M.t + PH * (Math.log10(v) - d0) / (d1 - d0);
    let lo = 0, hi = 0;
    for (const s of cfg.series) for (const v of s.v) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    const pad = (hi - lo) * 0.04 || 0.01, xt = niceTicks(lo - pad, hi + pad, 5);
    const xa = xt[0], xb = xt[xt.length - 1];
    const fx = v => M.l + PW * (v - xa) / (xb - xa);
    for (const v of xt) {
      el("line", { class: "gridl", x1: fx(v), x2: fx(v), y1: M.t, y2: M.t + PH }, svg);
      el("text", { x: fx(v), y: M.t + PH + 15, "text-anchor": "middle" }, svg).textContent = +v.toPrecision(3);
    }
    for (let d = d0; d <= d1; d++) {
      el("line", { class: "gridl", x1: M.l, x2: M.l + PW, y1: fy(10 ** d), y2: fy(10 ** d) }, svg);
      el("text", { x: M.l - 6, y: fy(10 ** d) + 3.5, "text-anchor": "end" }, svg).textContent = +(10 ** d).toPrecision(3);
    }
    if (xa < 0 && xb > 0) el("line", { class: "axis", x1: fx(0), x2: fx(0), y1: M.t, y2: M.t + PH }, svg);
    el("line", { class: "axis", x1: M.l, x2: M.l + PW, y1: M.t + PH, y2: M.t + PH }, svg);
    axes(svg, cfg.xlabel || "heating rate [K day⁻¹]", "pressure [hPa]");
    const cross = el("line", { class: "cross", x1: M.l, x2: M.l + PW, y1: 0, y2: 0, visibility: "hidden" }, svg);
    for (const i of order(cfg)) {
      const s = cfg.series[i];
      path(svg, s.v.map((v, k) => `${k ? "L" : "M"}${fx(v).toFixed(2)},${fy(p[k]).toFixed(2)}`).join(""), s);
    }
    const hit = el("rect", { x: M.l, y: M.t, width: PW, height: PH, fill: "transparent" }, svg);
    hit.addEventListener("pointermove", ev => {
      const r = svg.getBoundingClientRect(), py = (ev.clientY - r.top) * H / r.height, px = (ev.clientX - r.left) * W / r.width;
      const lp = d0 + (py - M.t) / PH * (d1 - d0);
      let k = 0, best = Infinity;
      for (let j = 0; j < p.length; j++) { const dd = Math.abs(Math.log10(p[j]) - lp); if (dd < best) { best = dd; k = j; } }
      cross.setAttribute("y1", fy(p[k])); cross.setAttribute("y2", fy(p[k])); cross.setAttribute("visibility", "visible");
      showTip(c, tip, svg, px, fy(p[k]), `${num(p[k])} hPa`,
        cfg.series.map(s => ({ ...s, value: num(s.v[k]) + " K/d" })));
    });
    hit.addEventListener("pointerleave", () => { cross.setAttribute("visibility", "hidden"); tip.hidden = true; });
  }

  // Line chart with markers over a shared x array; hover snaps to the nearest x.
  function xy(c, cfg) {
    const { svg, tip } = frame(c, cfg.title, cfg.sub);
    const xs = cfg.x, x0 = cfg.xdom ? cfg.xdom[0] : Math.min(...xs), x1 = cfg.xdom ? cfg.xdom[1] : Math.max(...xs);
    const fx = v => M.l + PW * (v - x0) / (x1 - x0);
    const all = cfg.series.flatMap(s => s.v.filter(v => v !== null));
    let fy;
    if (cfg.ylog) fy = yLog(svg, Math.min(...all), Math.max(...all));
    else {
      let lo = Math.min(...all), hi = Math.max(...all);
      if (cfg.zero) { lo = Math.min(lo, 0); hi = Math.max(hi, 0); }
      const pad = (hi - lo) * 0.06 || 1;
      fy = yLinear(svg, lo - (lo < 0 || !cfg.zero ? pad : 0), hi + pad);
    }
    for (const v of cfg.xticks || xs) {
      el("line", { class: "gridl", x1: fx(v), x2: fx(v), y1: M.t, y2: M.t + PH }, svg);
      el("text", { x: fx(v), y: M.t + PH + 15, "text-anchor": "middle" }, svg).textContent = v;
    }
    el("line", { class: "axis", x1: M.l, x2: M.l + PW, y1: M.t + PH, y2: M.t + PH }, svg);
    axes(svg, cfg.xlabel, cfg.ylabel);
    const cross = el("line", { class: "cross", x1: 0, x2: 0, y1: M.t, y2: M.t + PH, visibility: "hidden" }, svg);
    for (const i of order(cfg)) {
      const s = cfg.series[i];
      const pts = s.v.map((v, k) => v === null ? null : [fx(xs[k]), fy(v)]).filter(Boolean);
      path(svg, pts.map((q, k) => `${k ? "L" : "M"}${q[0].toFixed(2)},${q[1].toFixed(2)}`).join(""), s, 2);
      if (s.marker !== false) for (const q of pts)
        el("circle", { cx: q[0], cy: q[1], r: 4, fill: s.color, stroke: "var(--surface)", "stroke-width": 2 }, svg);
    }
    const hit = el("rect", { x: M.l - 8, y: M.t, width: PW + 16, height: PH, fill: "transparent" }, svg);
    hit.addEventListener("pointermove", ev => {
      const r = svg.getBoundingClientRect(), px = (ev.clientX - r.left) * W / r.width;
      let k = 0, best = Infinity;
      xs.forEach((x, j) => { const dd = Math.abs(fx(x) - px); if (dd < best) { best = dd; k = j; } });
      cross.setAttribute("x1", fx(xs[k])); cross.setAttribute("x2", fx(xs[k])); cross.setAttribute("visibility", "visible");
      showTip(c, tip, svg, fx(xs[k]), M.t + 20, (cfg.fmtX || (v => v))(xs[k]),
        cfg.series.filter(s => s.v[k] !== null).map(s => ({ ...s, value: (cfg.fmtY || num)(s.v[k]) })));
    });
    hit.addEventListener("pointerleave", () => { cross.setAttribute("visibility", "hidden"); tip.hidden = true; });
  }

  function legend(target, series) {
    target.innerHTML = series.map(s => `<span>${swatch(s, 26)}${s.name}</span>`).join("");
  }

  return { card, spectrum, profile, xy, legend, num, signed };
})();
