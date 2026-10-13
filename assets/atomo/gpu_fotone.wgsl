// IL FOTONE SULLA GPU: misure, palline di energia e disegno, senza rileggere il campo sulla CPU
// (studio motore_grafico_gpu, docs/MOTORE_GRAFICO_GPU.md). La regola è in regola_luce.wgsl.
//
// Ogni kernel rifà, con le stesse operazioni f32 nello stesso ordine, il codice CPU del Laboratorio:
//   misura   -> Fotone3D::misura + Fotone3D::elicita + somma di E x B (fotone3d.rs, gpu_fotone.rs)
//   avanza   -> Traccianti::avanza (fotdisegno.rs): passo di Eulero con v = 2 E x B / (E^2 + B^2), |v| <= 1
//   rinasci  -> la pallina persa prende il posto di una viva a caso (+-1 cella): DIVERSO dalla CPU (dichiarato)
//   spin     -> somme per lo spin letto dalle palline
//   istanze  -> fotdisegno::palline (colore = direzione di E attorno al volo, allungamento lungo E, scia)
// Le somme di un gruppo di lavoro sono in f32; la CPU le somma tra i gruppi in f64.

struct Par {
    n: u32, tick: u32, n_palline: u32, ns: u32,
    verso: f32, scia_op: u32, allunga: u32, n_sonde: u32,
    rif: vec4<f32>,   // centro atteso (celle dal centro della scatola)
    e1: vec4<f32>,    // base trasversa e1; w = allungamento massimo (celle)
    e2: vec4<f32>,
    dir: vec4<f32>,
};

@group(0) @binding(0)  var<uniform>             par : Par;
@group(0) @binding(1)  var<storage, read>       fx : array<vec2<f32>>;
@group(0) @binding(2)  var<storage, read>       fy : array<vec2<f32>>;
@group(0) @binding(3)  var<storage, read>       fz : array<vec2<f32>>;
@group(0) @binding(4)  var<storage, read_write> ex : array<f32>;     // E del tick prima (per l'elicità)
@group(0) @binding(5)  var<storage, read_write> ey : array<f32>;
@group(0) @binding(6)  var<storage, read_write> ez : array<f32>;
@group(0) @binding(7)  var<storage, read_write> parz : array<f32>;   // 20 numeri per gruppo di lavoro
@group(0) @binding(8)  var<storage, read_write> stato : array<f32>;  // [0] = max |F|^2, [1] = max |E|^2
@group(0) @binding(9)  var<storage, read_write> pos : array<vec4<f32>>;
@group(0) @binding(10) var<storage, read_write> vel : array<vec4<f32>>;
@group(0) @binding(11) var<storage, read_write> scia : array<vec4<f32>>; // 4 per pallina, la più recente prima
@group(0) @binding(12) var<storage, read_write> eta : array<u32>;
@group(0) @binding(13) var<storage, read_write> ist : array<Istanza>;
@group(0) @binding(14) var<storage, read>       sonde_in : array<vec4<f32>>;
@group(0) @binding(15) var<storage, read_write> sonde_out : array<vec4<f32>>;
@group(0) @binding(16) var<storage, read_write> parzp : array<f32>;  // 8 numeri per gruppo (spin dalle palline)

struct Istanza { pos_size: vec4<f32>, allunga: vec3<f32>, colore: u32 };

const MORTA: u32 = 0xffffffffu;
const SCIA_MAX: u32 = 4u;
const NP: u32 = 20u;
const TAU: f32 = 6.28318530717958647692;

// ───────────────────────── misure dal campo ─────────────────────────
var<workgroup> rid: array<f32, 1280>; // 64 righe x 20 numeri

fn avvolgi(v: f32, nf: f32) -> f32 { return v - nf * round(v / nf); }

@compute @workgroup_size(64)
fn misura(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let n = par.n;
    let r = gid.x;           // una riga (i, j) di n celle lungo k
    let l = lid.x;
    var a: array<f32, 20>;
    for (var q = 0u; q < NP; q++) { a[q] = 0.0; }
    if (r < n * n) {
        let i = r / n;
        let j = r % n;
        let nf = f32(n);
        let h = f32(n / 2u);
        let ddx = avvolgi(f32(i) - h - par.rif.x, nf);
        let ddy = avvolgi(f32(j) - h - par.rif.y, nf);
        let e1 = par.e1.xyz;
        let e2 = par.e2.xyz;
        for (var k = 0u; k < n; k++) {
            let c = r * n + k;
            let px = fx[c]; let py = fy[c]; let pz = fz[c];
            let e = vec3<f32>(px.x, py.x, pz.x);
            let b = vec3<f32>(px.y, py.y, pz.y);
            let w = (px.x * px.x + px.y * px.y) + (py.x * py.x + py.y * py.y) + (pz.x * pz.x + pz.y * pz.y);
            // elicità: E prima x E ora, proiettati sul piano trasverso
            let ep = vec3<f32>(ex[c], ey[c], ez[c]);
            let a1 = dot(ep, e1); let a2 = dot(ep, e2);
            let b1 = dot(e, e1);  let b2 = dot(e, e2);
            a[10] += a1 * b2 - a2 * b1;
            a[11] += 0.5 * (a1 * a1 + a2 * a2 + b1 * b1 + b2 * b2);
            ex[c] = e.x; ey[c] = e.y; ez[c] = e.z;
            // Poynting
            let s = cross(e, b);
            a[12] += s.x; a[13] += s.y; a[14] += s.z;
            a[15] = max(a[15], w);
            a[16] = max(a[16], dot(e, e));
            if (w == 0.0) { continue; }
            let ddz = avvolgi(f32(k) - h - par.rif.z, nf);
            a[0] += w;
            a[1] += w * ddx; a[2] += w * ddy; a[3] += w * ddz;
            a[4] += w * ddx * ddx; a[5] += w * ddy * ddy; a[6] += w * ddz * ddz;
            a[7] += w * ddx * ddy; a[8] += w * ddx * ddz; a[9] += w * ddy * ddz;
        }
    }
    for (var q = 0u; q < NP; q++) { rid[l * NP + q] = a[q]; }
    workgroupBarrier();
    for (var s = 32u; s > 0u; s = s >> 1u) {
        if (l < s) {
            for (var q = 0u; q < 15u; q++) { rid[l * NP + q] += rid[(l + s) * NP + q]; }
            rid[l * NP + 15u] = max(rid[l * NP + 15u], rid[(l + s) * NP + 15u]);
            rid[l * NP + 16u] = max(rid[l * NP + 16u], rid[(l + s) * NP + 16u]);
        }
        workgroupBarrier();
    }
    if (l == 0u) { for (var q = 0u; q < NP; q++) { parz[wid.x * NP + q] = rid[q]; } }
}

// i massimi di tutti i gruppi -> stato (per le palline, senza passare dalla CPU)
var<workgroup> mx: array<vec2<f32>, 256>;

@compute @workgroup_size(256)
fn finalizza(@builtin(local_invocation_id) lid: vec3<u32>) {
    let ng = (par.n * par.n + 63u) / 64u;
    var m = vec2<f32>(0.0, 0.0);
    for (var g = lid.x; g < ng; g += 256u) { m = max(m, vec2<f32>(parz[g * NP + 15u], parz[g * NP + 16u])); }
    mx[lid.x] = m;
    workgroupBarrier();
    for (var s = 128u; s > 0u; s = s >> 1u) {
        if (lid.x < s) { mx[lid.x] = max(mx[lid.x], mx[lid.x + s]); }
        workgroupBarrier();
    }
    if (lid.x == 0u) { stato[0] = mx[0].x; stato[1] = mx[0].y; }
}

// ───────────────────────── campo in un punto (trilineare, periodico) ─────────────────────────
fn campo_in(p: vec3<f32>) -> array<vec2<f32>, 3> {
    let n = par.n;
    let ni = i32(n);
    let h = f32(n / 2u);
    var i0: array<u32, 3>;
    var w: array<f32, 3>;
    for (var c = 0; c < 3; c++) {
        let x = p[c] + h;
        let fl = floor(x);
        w[c] = x - fl;
        i0[c] = u32(((i32(fl) % ni) + ni) % ni);
    }
    var o: array<vec2<f32>, 3>;
    o[0] = vec2<f32>(0.0); o[1] = vec2<f32>(0.0); o[2] = vec2<f32>(0.0);
    for (var dx = 0u; dx < 2u; dx++) {
        for (var dy = 0u; dy < 2u; dy++) {
            for (var dz = 0u; dz < 2u; dz++) {
                let wx = select(1.0 - w[0], w[0], dx == 1u);
                let wy = select(1.0 - w[1], w[1], dy == 1u);
                let wz = select(1.0 - w[2], w[2], dz == 1u);
                let wt = wx * wy * wz;
                if (wt == 0.0) { continue; }
                let c = (((i0[0] + dx) % n) * n + (i0[1] + dy) % n) * n + (i0[2] + dz) % n;
                let qx = fx[c]; let qy = fy[c]; let qz = fz[c];
                o[0].x += wt * qx.x; o[0].y += wt * qx.y;
                o[1].x += wt * qy.x; o[1].y += wt * qy.y;
                o[2].x += wt * qz.x; o[2].y += wt * qz.y;
            }
        }
    }
    return o;
}

// velocità dell'energia v = 2 E x B / (E^2 + B^2) (|v| <= 1) e densità u
fn vel_energia(p: vec3<f32>) -> vec4<f32> {
    let q = campo_in(p);
    let e = vec3<f32>(q[0].x, q[1].x, q[2].x);
    let b = vec3<f32>(q[0].y, q[1].y, q[2].y);
    let u = (e.x * e.x + e.y * e.y + e.z * e.z) + (b.x * b.x + b.y * b.y + b.z * b.z);
    if (u <= 0.0) { return vec4<f32>(0.0); }
    let s = vec3<f32>(e.y * b.z - e.z * b.y, e.z * b.x - e.x * b.z, e.x * b.y - e.y * b.x);
    var v = vec3<f32>(2.0 * s.x / u, 2.0 * s.y / u, 2.0 * s.z / u);
    let m = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (m > 1.0) { v = v / m; }
    return vec4<f32>(v, u);
}

// ───────────────────────── palline di energia ─────────────────────────
@compute @workgroup_size(256)
fn avanza(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= par.n_palline) { return; }
    let p = pos[i].xyz;
    let b = i * SCIA_MAX;
    // la scia: avanti, ogni 2 tick la storia scorre; indietro si azzera
    if (par.scia_op == 1u) {
        for (var s = SCIA_MAX - 1u; s > 0u; s--) { scia[b + s] = scia[b + s - 1u]; }
        scia[b] = vec4<f32>(p, 0.0);
    } else if (par.scia_op == 2u) {
        for (var s = 0u; s < SCIA_MAX; s++) { scia[b + s] = vec4<f32>(p, 0.0); }
    }
    let imax = max(stato[0], 1e-30);
    let vu = vel_energia(p);
    if (vu.w < 1e-4 * imax) { eta[i] = MORTA; return; }
    pos[i] = vec4<f32>(p.x + par.verso * vu.x, p.y + par.verso * vu.y, p.z + par.verso * vu.z, 0.0);
    vel[i] = vec4<f32>(vu.xyz, 0.0);
    eta[i] = min(eta[i] + 1u, 65535u);
}

fn hash(x: u32) -> u32 {
    var h = x * 747796405u + 2891336453u;
    h = ((h >> ((h >> 28u) + 4u)) ^ h) * 277803737u;
    return (h >> 22u) ^ h;
}
fn hash01(x: u32) -> f32 { return f32(hash(x) >> 8u) * (1.0 / 16777216.0); }

@compute @workgroup_size(256)
fn rinasci(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= par.n_palline) { return; }
    if (eta[i] != MORTA) { return; }
    var p = pos[i].xyz;
    let base = hash(i ^ hash(par.tick * 9781u + 17u));
    for (var t = 0u; t < 16u; t++) {
        let j = hash(base + t * 2654435761u) % par.n_palline;
        let ej = eta[j];
        // una pallina viva che non sta rinascendo adesso (eta 0 = rinata in questo passo)
        if (ej != MORTA && ej != 0u) {
            let q = pos[j].xyz;
            let s = hash(base ^ (t + 101u));
            p = q + vec3<f32>(2.0 * hash01(s) - 1.0, 2.0 * hash01(s + 1u) - 1.0, 2.0 * hash01(s + 2u) - 1.0);
            break;
        }
    }
    pos[i] = vec4<f32>(p, 0.0);
    vel[i] = vec4<f32>(0.0);
    for (var s = 0u; s < SCIA_MAX; s++) { scia[i * SCIA_MAX + s] = vec4<f32>(p, 0.0); }
    eta[i] = 0u;
}

// i surfisti (pochi, gestiti dalla CPU): velocità e densità nei loro punti
@compute @workgroup_size(8)
fn sonde(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x >= par.n_sonde) { return; }
    sonde_out[gid.x] = vel_energia(sonde_in[gid.x].xyz);
}

// spin dalle palline: somme di (p - rif) x v e di v (la CPU corregge per il centro misurato)
var<workgroup> rs: array<f32, 2048>; // 256 x 8

@compute @workgroup_size(256)
fn spin(@builtin(global_invocation_id) gid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>, @builtin(num_workgroups) nw: vec3<u32>) {
    var a: array<f32, 8>;
    for (var q = 0u; q < 8u; q++) { a[q] = 0.0; }
    let passo = nw.x * 256u;
    for (var i = gid.x; i < par.n_palline; i += passo) {
        let r = pos[i].xyz - par.rif.xyz;
        let v = vel[i].xyz;
        let l = cross(r, v);
        a[0] += l.x; a[1] += l.y; a[2] += l.z;
        a[3] += v.x; a[4] += v.y; a[5] += v.z;
    }
    let t = lid.x;
    for (var q = 0u; q < 8u; q++) { rs[t * 8u + q] = a[q]; }
    workgroupBarrier();
    for (var s = 128u; s > 0u; s = s >> 1u) {
        if (t < s) { for (var q = 0u; q < 8u; q++) { rs[t * 8u + q] += rs[(t + s) * 8u + q]; } }
        workgroupBarrier();
    }
    if (t == 0u) { for (var q = 0u; q < 8u; q++) { parzp[wid.x * 8u + q] = rs[q]; } }
}

// ───────────────────────── disegno: istanze direttamente dal buffer delle palline ─────────────────────────
fn colore_fase(ang: f32, alfa: f32) -> vec4<f32> {
    var st = array<vec3<f32>, 5>(vec3<f32>(0.20, 0.95, 1.00), vec3<f32>(0.25, 0.45, 1.00), vec3<f32>(0.70, 0.30, 1.00),
                                 vec3<f32>(1.00, 0.35, 0.70), vec3<f32>(1.00, 0.85, 0.30));
    let a = ang - TAU * floor(ang / TAU);
    let x = (a / TAU) * 5.0;
    let i = u32(floor(x)) % 5u;
    let s = x - floor(x);
    let c = st[i] + (st[(i + 1u) % 5u] - st[i]) * s;
    return vec4<f32>(c, alfa);
}

// coordinate della scena -> schermo (render.rs::schermo: (x, z, -y))
fn schermo(p: vec3<f32>) -> vec3<f32> { return vec3<f32>(p.x, p.z, -p.y); }

@compute @workgroup_size(256)
fn istanze(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= par.n_palline) { return; }
    let imax = max(stato[0], 1e-30);
    let emax = max(sqrt(stato[1]), 1e-30);
    let p = pos[i].xyz;
    let q = campo_in(p);
    let e = vec3<f32>(q[0].x, q[1].x, q[2].x);
    let b = vec3<f32>(q[0].y, q[1].y, q[2].y);
    let u = (dot(e, e) + dot(b, b)) / imax;
    let ang = atan2(dot(e, par.e2.xyz), dot(e, par.e1.xyz));
    let fade = min(f32(eta[i]) / 6.0, 1.0);
    let a = fade * min(0.12 + 0.6 * sqrt(u), 0.72);
    let c = colore_fase(ang, a);
    var st = vec3<f32>(0.0);
    if (par.allunga != 0u) { st = schermo(e * (par.e1.w / emax)); }
    let ns = par.ns;
    let base = i * (ns + 1u);
    // prima la scia (disegnata dietro la testa): la stessa pallina λ/3, 2λ/3, λ, 4λ/3 fa
    for (var s = 0u; s < ns; s++) {
        let j = ns - 1u - s;
        let w = 0.5 * (1.0 - f32(j) / f32(ns));
        let sj = scia[i * SCIA_MAX + j].xyz;
        // allungata a metà strada verso il campione più nuovo (la testa per j = 0), come fotdisegno::palline
        var nuovo = p;
        if (j > 0u) { nuovo = scia[i * SCIA_MAX + j - 1u].xyz; }
        ist[base + s] = Istanza(vec4<f32>(schermo(sj), 0.85 - 0.08 * f32(j)), schermo(0.5 * (nuovo - sj)),
                                pack4x8unorm(vec4<f32>(c.rgb, w * a)));
    }
    ist[base + ns] = Istanza(vec4<f32>(schermo(p), 1.0), st, pack4x8unorm(c));
}
