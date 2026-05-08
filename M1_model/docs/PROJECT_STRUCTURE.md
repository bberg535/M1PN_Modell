# M1\_model Projektstruktur

## Ordner
- `src/`: aktive MATLAB-Funktionen des M1PN-Modells und der Referenzloeser
- `scripts/`: ausfuehrbare MATLAB-Skripte fuer Hauptlauf und Benchmarks
- `docs/`: kurze Projektdokumentation
- `results/`: Ablage fuer Plots, Exporte und numerische Ausgaben
- `legacy/`: alte oder nicht mehr benoetigte Skripte/Hilfsdateien

## Aktive Einstiegspunkte
- `scripts/m1pn_model.m`: Hauptskript fuer das gekoppelte M1PN-Modell
- `scripts/m1_plane_source_benchmark.m`: Plane-Source-Benchmark gegen eine S_N-Referenz
- `scripts/pn_m1pn_negativity_benchmark.m`: Plane-Source-Stresstest fuer PN-vs-M1PN mit Fokus auf negativen PN-Dichten

## Hinweise
- Die MATLAB-Pfade werden ueber `setup_project_paths.m` initialisiert.
- Alte `*.asv`-Backupdateien werden nicht mitgefuehrt.
