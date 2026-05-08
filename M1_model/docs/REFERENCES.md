# Quellen und Zuordnung

## Aktiver HOLO/M1PN-Pfad

- **M1-Closure / Eddington-Faktor**
  - Implementiert in `src/calc_psi2.m`.
  - Verwendet die 1D-M1-Closure ueber den reduzierten Fluss `f = |j| / rho`.
  - Die konkrete Formel war bereits im Projekt vorhanden; wenn du dafuer die exakte Literaturstelle mitgeben willst, kann ich sie hier noch sauber ergaenzen.

- **Fail-safe flux limiting**
  - Implementiert im aktiven gekoppelten Pfad in `src/failsafe_limit_density.m`.
  - Orientierung: Abschnitt **2.5.6.2 "Fail-safe flux limiting"** aus deinem Buch *Property-Preserving Numerical Schemes for Conservation Laws*.
  - In der aktuellen Version wird nur die Bedingung `rho >= 0` erzwungen, nicht das volle M1-Realisierbarkeitsgebiet.

- **HOLO-Idee fuer die Kopplung**
  - Implementiert in `src/heuns_method.m` und `src/real_HO_flux.m`.
  - M1 bleibt das eigentliche Modell; PN liefert nur den High-Order-Kandidatenfluss.
  - Diese Struktur folgt deiner beschriebenen HOLO-Kopplung:
    `M1-LO + PN-HO + fail-safe limiting auf den Korrekturfluss`.

## PN-Transport und Benchmarks

- **PN-Transport in charakteristischen Variablen**
  - Historischer Ausgangspunkt des PN-FV-Teils in den Legacy-Dateien:
    `doi:10.1137/15M1052871`
  - Der aktive Code nutzt jetzt die gemeinsame Flussroutine `src/pn_interface_flux.m`, die denselben Grundgedanken in einer zentralen Schnittstelle abbildet.

- **Plane-source Benchmark**
  - Implementiert in `scripts/m1_plane_source_benchmark.m` und
    `scripts/pn_m1pn_negativity_benchmark.m`.
  - Orientierung: `Seminarquelle 2`, Abschnitt **6.1.1 Plane source**.

- **Positive PN Closures, Figure 3**
  - Implementiert in `src/simulate_pn_siu.m` und
    `scripts/positive_pn_closures_figure3.m`.
  - Orientierung: C. Hauck and R. McClarren,
    *Positive PN Closures*, SIAM J. Sci. Comput. 32(5), 2010,
    DOI `10.1137/090764918`.
  - Verwendete Stellen:
    - Abschnitt **3.2 One-dimensional setting**
    - diskrete Ordinatenform **(13)**
    - SIU-Schema **(15)-(20)**
    - Diagnosegroesse `F*` aus **(21)**
    - Plotziel: **Figure 3**
  - Im aktiven Skript wird der Figure-3-Testfall jetzt auf das aktuelle
    `M1PN`-Modell angewendet; der SIU-PN-Lauf dient dort nur noch als
    Referenz auf demselben Testfall.

## Offene / von dir nachlieferbare Referenzen

- **Paper mit Paul zur realizierbaren Closure**
  - Darauf verweisen wir inhaltlich, aber ich habe aktuell keine exakte bibliographische Angabe.
  - Wenn du mir PDF, DOI, BibTeX oder einen Screenshot der Titelseite gibst, trage ich die Referenz hier sauber ein und kommentiere die betroffenen Stellen genauer.
