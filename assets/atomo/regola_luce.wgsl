// LA REGOLA ZIG-ZAG DELLA LUCE SULLA GPU (porta di src/atomo/fotone3d.rs, Fotone3D::shift).
//
// Campo F = E + icB: 3 numeri complessi per cella, salvati come 3 vec2<f32> consecutivi
// (indice = 3*cella + componente, cella = (i*n + j)*n + k, i = x, j = y, k = z). Cubo periodico n^3.
// Una chiamata = una traslazione lungo l'asse `ax`: la parte +1 di S_ax va avanti di una cella, la parte -1
// indietro, la parte 0 resta. Con a = F_{ax+1}, b = F_{ax+2}, m = cella dietro, p = cella davanti:
//   a' = 1/2 (a_m + a_p) + 1/2 i (b_p - b_m),  b' = 1/2 (b_m + b_p) + 1/2 i (a_m - a_p),  F_ax' = F_ax.
// verso = -1 scambia m e p (l'inverso esatto). Un tick = assi x,y,z (tick pari) o z,y,x (tick dispari):
// l'ordine lo decide il lato Rust, una dispatch per asse, due buffer a ping-pong.

struct Parametri {
    n:     u32,
    ax:    u32,
    verso: i32,
    _pad:  u32,
}

@group(0) @binding(0) var<uniform>             par : Parametri;
@group(0) @binding(1) var<storage, read>       src : array<vec2<f32>>;
@group(0) @binding(2) var<storage, read_write> dst : array<vec2<f32>>;

@compute @workgroup_size(8, 8, 4)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let n = par.n;
    let k = gid.x;
    let j = gid.y;
    let i = gid.z;
    if (i >= n || j >= n || k >= n) { return; }
    let cella = (i * n + j) * n + k;

    // vicini lungo l'asse (periodico)
    var c = vec3<u32>(i, j, k);
    var cm = c;
    var cp = c;
    let a = par.ax;
    cm[a] = select(c[a] - 1u, n - 1u, c[a] == 0u);
    cp[a] = select(c[a] + 1u, 0u, c[a] + 1u == n);
    var im = (cm.x * n + cm.y) * n + cm.z;
    var ip = (cp.x * n + cp.y) * n + cp.z;
    if (par.verso < 0) { let t = im; im = ip; ip = t; }

    let j1 = (a + 1u) % 3u;
    let j2 = (a + 2u) % 3u;
    let am = src[3u * im + j1];
    let bm = src[3u * im + j2];
    let ap = src[3u * ip + j1];
    let bp = src[3u * ip + j2];

    // stesse espressioni del lato CPU (f32), nello stesso ordine
    let a_new = vec2<f32>(0.5 * (am.x + ap.x) - 0.5 * (bp.y - bm.y),
                          0.5 * (am.y + ap.y) + 0.5 * (bp.x - bm.x));
    let b_new = vec2<f32>(0.5 * (bm.x + bp.x) - 0.5 * (am.y - ap.y),
                          0.5 * (bm.y + bp.y) + 0.5 * (am.x - ap.x));

    dst[3u * cella + a]  = src[3u * cella + a];
    dst[3u * cella + j1] = a_new;
    dst[3u * cella + j2] = b_new;
}
