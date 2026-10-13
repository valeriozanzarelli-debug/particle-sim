## R1 — regola sulla GPU

**Cosa si prova.** La regola zig-zag della luce (`Fotone3D::tick` in `src/atomo/fotone3d.rs`) riscritta
come compute shader WGSL (`assets/atomo/regola_luce.wgsl`, lato Rust `src/atomo/gpu_regola.rs`). Stessa
regola, niente di nuovo: per ogni asse (x, y, z nei tick pari; z, y, x nei dispari) F si divide sugli
autovettori di S_j, la parte +1 va avanti di una cella, la −1 indietro, la 0 resta. Una dispatch per asse,
due buffer che si passano il campo a ping-pong (3 numeri complessi f32 per cella), cubo periodico.

**Esame pre-registrato** (`lab/out/mg_prereg.json`): stesso impulso del Laboratorio (parametri di default:
8 celle per λ, n = 96, w0 = 1 λ, circolare +, lungo x), 96 tick su CPU e su GPU,
max|F_gpu − F_cpu| / max|F_cpu| ≤ 1e-4.

**Risultato (cloud, Vulkan software llvmpipe): PASSA.**

| misura | valore |
|---|---|
| max\|F_gpu − F_cpu\| / max\|F_cpu\| dopo 96 tick | **0** (identici bit per bit) |
| soglia | 1e-4 |
| quanto si è mosso il campo in 96 tick (max\|F₉₆ − F₀\| / max\|F\|) | 1.29 (il confronto non è banale) |
| 96 tick avanti + 96 indietro sulla GPU, errore rispetto al lancio | 7.5e-7 (arrotondamenti f32) |
| ms per tick, CPU (rayon, 4 core) | 8.8 |
| ms per tick, GPU llvmpipe (cioè Vulkan che gira sulla CPU) | 41 |

**Come leggerlo.** La differenza è zero perché lo shader fa le stesse operazioni f32 nello stesso ordine
e llvmpipe non le fonde (niente FMA). Su una scheda vera il compilatore può fondere moltiplicazione e
somma: lì ci aspettiamo differenze dell'ordine di 1e-7, sempre molto sotto la soglia — da verificare sul
portatile via ponte. I tempi del cloud **non dicono niente** sulla GPU: llvmpipe è un emulatore su CPU,
più lento di rayon. È una validazione (la GPU rifà esattamente la regola), non una scoperta.

**Manopole:** zero. **Debiti nuovi:** nessuno. **Limite noto:** un buffer di campo pesa 24·n³ byte
(n = 192: 170 MB; n = 384: 1.36 GB, due buffer); serve che il device abbia
`max_storage_buffer_binding_size` abbastanza grande (il default di wgpu è 128 MB).
