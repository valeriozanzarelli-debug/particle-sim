# Atomo supercritico (Z oltre 137): l'1s che si tuffa nel mare di Dirac

Studio `supercrit_` — file: `lab/supercrit_*.py`, `lab/out/supercrit_*`. Stato: **pre-registrato, corse in arrivo**.

## 1. Il problema (cosa è noto, cosa è aperto)

- Con un nucleo **puntiforme** l'equazione di Dirac dà E_1s = mc² √(1 − (Zα)²): a Zα = 1 (Z ≈ 137) l'energia va
  a zero e oltre diventa immaginaria: la formula si rompe ("caduta sul centro").
- Con un nucleo **esteso** (raggio di qualche fm) l'1s continua a scendere e tocca −mc² (il bordo del "mare di
  Dirac") a una carica critica Z_cr ≈ 170–173: Pomeranchuk–Smorodinsky (J. Phys. USSR 9, 97, 1945, stima troppo
  alta); Gershtein–Zeldovich (Lett. Nuovo Cimento 1, 835, 1969); Pieper–Greiner (Z. Phys. 218, 327, 1969);
  Müller–Rafelski–Greiner (Z. Phys. 257, 62 e 183, 1972); Zeldovich–Popov (Sov. Phys. Usp. 14, 673, 1972).
  Valori moderni: Z_cr ≈ 173 (Maltsev et al., PRL 123, 113401, 2019, con modello di nucleo realistico);
  ≈ 170 con R = 1,2 (2,5 Z)^{1/3} fm (Popov et al., arXiv:2209.00125).
- Oltre Z_cr lo stato non è più legato: diventa una **risonanza** dentro il continuo negativo. Se il guscio K è
  vuoto, un elettrone del mare lo riempie e il buco che resta scappa: **emissione spontanea di un positrone**;
  il vuoto diventa "carico". Le larghezze Γ della risonanza sono piccole (eV–keV): esempio, tabella esatta di
  Godunov, Machet, Vysotsky (EPJC 77, 782, 2017) per un nucleo "a guscio" di raggio R = 0,031 ħ/mc (12 fm):
  Z_cr ≈ 177,5, Γ(184) = 5,2·10⁻⁴ mc² = 265 eV, Γ(190) = 6,1·10⁻³ mc², Γ(232) = 0,19 mc².
  Vicino alla soglia Γ è soppressa dalla barriera di Coulomb vista dal positrone: Γ ∝ exp(−2πη),
  η = Zα E/p (Zeldovich–Popov).
- **Mai osservato.** Nessun nucleo ha Z > 118; si usano collisioni lente U+U (Z = 184 per ~10⁻²¹ s). I picchi di
  positroni visti a GSI negli anni '80 (EPOS, ORANGE) non sono stati confermati (APEX, anni '90). Oggi si
  propongono collisioni U–U a energie sotto la barriera con traiettorie scelte (Maltsev et al. 2019; Popov et
  al.), esperimenti previsti a GSI/FAIR e HIAF. È quindi una previsione standard della QED in campo forte
  mai verificata: il problema aperto è **sperimentale**, non teorico.
- **Analogo in materiali**: il "collasso atomico" in grafene (Shytov–Katsnelson–Levitov 2007; Pereira–Nilsson–
  Castro Neto 2007) è stato osservato (Wang et al., Science 340, 734, 2013). Lì il reticolo di carbonio fa da
  taglio vicino all'impurità, proprio come farà il nostro reticolo.
- **Quantum walk / automi di Dirac**: ci sono walk con campo elettrico e con Coulomb (Di Molfetta–Debbasch–
  Brachet 2014; walk di Dirac con Coulomb per la ricerca spaziale, arXiv:2210.13920), automi di Dirac 1-D con
  scattering da potenziale (Bisio–D'Ariano–Tosini 2013), simulazioni su ioni intrappolati (Huerta Alderete et al.
  2020). **Non abbiamo trovato** un lavoro che faccia con un walk/automa di Dirac il tuffo dell'1s oltre Z_cr e il
  decadimento della risonanza nel tempo. Le simulazioni tempo-dipendenti della risonanza esistono, ma con
  l'equazione di Dirac ordinaria (Müller, Reinhardt, Greiner anni '70–'80; Maltsev, Shabaev, Tupitsyn 2015–2019).
  Quindi: se funziona, è **nuovo come metodo, non come fisica**.

(Letteratura controllata online il 13/10/2026: arXiv:2008.05005, arXiv:2209.00125, arXiv:1707.07497, tesi di
R.V. Popov 2024 — U+U: R_cr = 32,7 fm, larghezza naturale < 2 keV a R = 16 fm. Gli altri riferimenti storici sono
citati a memoria.)

## 2. Previsione prima dei dati

Vedi `lab/out/supercrit_prereg.json` (generato da `lab/supercrit_prereg.py`). In breve:

1. Nucleo puntiforme, Z ≤ 137: la regola ritrova Dirac, ma con errore che si chiude come a^{2γ} (γ = √(1−(Zα)²)):
   a Z = 137, 2γ = 0,05, cioè **in pratica non converge** (limite atteso).
2. Nucleo puntiforme, Zα > 1: il reticolo non ha punti dentro r = a/2, quindi fa da **nucleo finto grande ~a**.
   Previsione: l'1s tocca −mc² a un Z_cr(a) finito, che scende verso 137 come 1 + π²/(2 ln²(1/(c·a))) con c fisso.
   È un artefatto (dichiarato): non è fisica.
3. Nucleo esteso (DEBITO: il raggio): Z_cr = 177,46 per il guscio R = 0,031 (letteratura), 170–173 per la sfera.
4. Oltre Z_cr: la risonanza si vede sia dalla fase (salto di π) sia dal tempo (lo stato si svuota con legge
   esponenziale e un'onda di carica positiva esce). Γ entro il 10 % della tabella di Godunov et al.

## 3–8. Risultati, esami, debiti

(in arrivo)
