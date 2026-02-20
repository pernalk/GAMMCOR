# CPU Profiling z perf stat — instrukcja dla agenta kodującego

## Cel

Zebrać hardware counter metrics z programu Fortran (obliczenia quantum chemistry RPA/ERPA) używając `perf stat`. Wynik: plik tekstowy z metrykami do analizy bottlenecków i oceny potencjału GPU acceleration.

## Wymagania systemowe

### Sprawdź dostępność narzędzi

```bash
# perf musi być zainstalowany
perf --version

# Jeśli brak:
# Ubuntu/Debian:
sudo apt install linux-tools-common linux-tools-$(uname -r)
# Fedora:
sudo dnf install perf
# RHEL/CentOS:
sudo yum install perf
```

### Sprawdź uprawnienia

```bash
cat /proc/sys/kernel/perf_event_paranoid
```

Wymagana wartość: `1` lub niżej. Jeśli wyższa:

```bash
sudo sysctl kernel.perf_event_paranoid=1
# Lub na stałe:
echo 'kernel.perf_event_paranoid=1' | sudo tee -a /etc/sysctl.conf
```

### Sprawdź dostępne eventy hardware

```bash
perf list hw cache
```

Zanotuj które eventy są dostępne — nie każdy CPU wspiera wszystkie countery.

## Kompilacja programu

Skompiluj z optymalizacją + debug symbols + frame pointer:

```bash
gfortran -O2 -g -fno-omit-frame-pointer -march=native \
    -o mycode [LISTA_PLIKÓW_F90] -lblas -llapack
```

Kluczowe flagi:

| Flaga | Cel |
|---|---|
| `-O2` | Profilujemy zoptymalizowany kod, NIE debug build |
| `-g` | Symbole debug — nazwy funkcji zamiast adresów w raporcie |
| `-fno-omit-frame-pointer` | Umożliwia perf chodzenie po call stacku |
| `-march=native` | Optymalizacje pod konkretny CPU |

Jeśli kompilator to nvfortran:

```bash
nvfortran -O2 -g -Mframe -o mycode [LISTA_PLIKÓW_F90] -lblas -llapack
```

Jeśli kompilator to ifort/ifx:

```bash
ifx -O2 -g -fno-omit-frame-pointer -xHost -o mycode [LISTA_PLIKÓW_F90] -lblas -llapack
```

## Przygotowanie środowiska

Przed profilowaniem wyłącz źródła szumu:

```bash
# Wyłącz frequency scaling (governor na performance)
sudo cpupower frequency-set -g performance 2>/dev/null || true

# Wyłącz turbo boost (opcjonalnie, dla powtarzalności)
# Intel:
echo 1 | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || true
# AMD:
echo 0 | sudo tee /sys/devices/system/cpu/cpufreq/boost 2>/dev/null || true
```

## Uruchomienie perf stat

### Krok 1: Podstawowe metryki (obowiązkowy)

```bash
perf stat -e cycles,instructions,cache-references,cache-misses,\
branches,branch-misses,L1-dcache-loads,L1-dcache-load-misses,\
LLC-loads,LLC-load-misses,task-clock,context-switches,cpu-migrations,page-faults \
    ./mycode < input.dat 2> results/perf_stat_basic.txt
```

Jeśli program nie czyta stdin, pomiń `< input.dat`. Dostosuj wywołanie do sposobu uruchamiania programu (argumenty, pliki konfiguracyjne, itp.).

### Krok 2: Metryki FP (opcjonalny, jeśli dostępne)

Sprawdź najpierw czy eventy istnieją:

```bash
perf list | grep fp_arith
```

Jeśli dostępne (Intel Skylake+):

```bash
perf stat -e fp_arith_inst_retired.scalar_double,\
fp_arith_inst_retired.128b_packed_double,\
fp_arith_inst_retired.256b_packed_double,\
fp_arith_inst_retired.512b_packed_double \
    ./mycode < input.dat 2> results/perf_stat_fp.txt
```

Jeśli AMD (Zen3+), szukaj:

```bash
perf list | grep -i 'fp\|flop\|sse\|avx'
```

I użyj dostępnych eventów. Jeśli żadne FP eventy nie są dostępne, pomiń ten krok.

### Krok 3: Metryki pamięci (opcjonalny)

```bash
perf stat -e dTLB-loads,dTLB-load-misses,\
iTLB-loads,iTLB-load-misses,\
bus-cycles,ref-cycles \
    ./mycode < input.dat 2> results/perf_stat_memory.txt
```

### Krok 4: Powtarzalność (obowiązkowy)

Uruchom 3-5 razy z obliczeniem statystyki:

```bash
perf stat -r 5 -e cycles,instructions,cache-references,cache-misses,\
branches,branch-misses,LLC-loads,LLC-load-misses,task-clock \
    ./mycode < input.dat 2> results/perf_stat_repeated.txt
```

W wyniku pojawią się wartości `( +- X.XX% )` — jeśli variancja jest powyżej 5%, zanotuj to, oznacza niedeterminizm.

## Struktura wyników

Utwórz katalog i zbierz wyniki:

```bash
mkdir -p results
```

Po zakończeniu w `results/` powinny być:

```
results/
├── perf_stat_basic.txt      # główne metryki (obowiązkowy)
├── perf_stat_fp.txt         # metryki FP (jeśli dostępne)
├── perf_stat_memory.txt     # metryki TLB/bus (opcjonalny)
├── perf_stat_repeated.txt   # statystyka powtarzalności (obowiązkowy)
└── system_info.txt          # info o systemie (patrz niżej)
```

## Zbierz informacje o systemie

```bash
{
  echo "=== CPU ==="
  lscpu | grep -E 'Model name|Socket|Core|Thread|MHz|Cache|Architecture'
  echo ""
  echo "=== Memory ==="
  free -h
  echo ""
  echo "=== Kernel ==="
  uname -r
  echo ""
  echo "=== Compiler ==="
  gfortran --version 2>/dev/null || nvfortran --version 2>/dev/null || ifx --version 2>/dev/null
  echo ""
  echo "=== BLAS/LAPACK ==="
  ldconfig -p | grep -E 'blas|lapack' || echo "statically linked or not found"
  echo ""
  echo "=== NUMA topology ==="
  numactl --hardware 2>/dev/null || echo "numactl not available"
  echo ""
  echo "=== GPU (jeśli obecny) ==="
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>/dev/null || echo "no NVIDIA GPU"
} > results/system_info.txt
```

## Troubleshooting

| Problem | Rozwiązanie |
|---|---|
| `perf_event_paranoid` error | `sudo sysctl kernel.perf_event_paranoid=1` |
| `<not supported>` przy evencie | Dany CPU nie wspiera tego countera — pomiń go |
| `<not counted>` przy evencie | Za dużo eventów naraz — rozdziel na osobne uruchomienia (max ~4-8 eventów hardware naraz, reszta jest multipleksowana) |
| Brak nazw funkcji (same adresy) | Brakuje `-g` przy kompilacji |
| Wynik `0` dla cache-misses | Program zbyt krótki — użyj większego inputu |
| Duża variancja (>10%) | Wyłącz turbo boost, ustaw governor na performance, sprawdź czy nie ma innych procesów |

## Uwagi

- `perf stat` nie mówi KTÓRE funkcje są wolne — mówi JAKI CHARAKTER ma bottleneck (memory/compute/branch). Do identyfikacji konkretnych funkcji potrzebny będzie `perf record` + flame graph (następny krok).
- Wszystkie pliki z `results/` są potrzebne do dalszej analizy — nie usuwaj ich.
- Jeśli program wykonuje się krócej niż 1 sekundę, użyj większego inputu. Perf stat potrzebuje wystarczająco dużo sampli żeby metryki były statystycznie istotne.
- Jeśli `perf stat -r 5` trwa zbyt długo (bo program liczy godzinami), użyj `-r 3` albo pomiń i uruchom jednokrotnie.