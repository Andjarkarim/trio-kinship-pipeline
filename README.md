# Verwandtschaftsanalyse aus Ganzgenomdaten (CEPH 1463)

Snakemake-Pipeline, die aus Ganzgenomdaten von vier Personen ableitet, wer mit wem verwandt ist, ohne dass Rollen vorgegeben werden. Grundlage sind Mendel-Fehler und der KING-robust-Schätzer. Ausgewertet werden die Chromosomen 20 bis 22 von NA12877, NA12878, NA12891 und NA12892 (Platinum Genomes).

## Ergebnis

- NA12878 ist das Kind von NA12891 und NA12892. Die Mendel-Fehlerrate liegt bei 0,015 % und 0,096 % pro Paar, bei allen anderen Paaren zwischen 5,4 und 5,9 %.
- Von den zwölf möglichen Trios hat nur NA12891 + NA12892 → NA12878 eine niedrige Fehlerrate (0,654 %).
- KING bestätigt das über IBS0: 0,00016 und 0,0010 bei den Eltern-Kind-Paaren, 0,058 bis 0,064 bei den übrigen.
- NA12877 ist mit den anderen drei nicht verwandt.
- Grundlage sind 329.422 gefilterte SNPs.

## Pipeline

FASTQ (ENA) → `bwa mem` + `samtools sort` → GATK `MarkDuplicates` → `HaplotypeCaller` (GVCF) → `CombineGVCFs` → `GenotypeGVCFs` → `bcftools` (Genotypfilter, SNP-Filter) → `concat` → Mendel-Fehler (`scripts/mendel_check.py`) und KING (PLINK 2).

Für NA12877 und NA12878 werden fertige BAM-Dateien verwendet, NA12891 und NA12892 werden aus dem ENA geladen und selbst abgeglichen.

## Inhalt

| Datei / Ordner | Inhalt |
|---|---|
| `Snakefile` | alle Regeln des Workflows |
| `config.example.yaml` | Vorlage für die Konfiguration |
| `scripts/mendel_check.py` | Mendel-Fehler für alle Paare und Trios |
| `scripts/plot_results.py` | vier Abbildungen (PNG und PDF) |
| `results/` | `mendel.tsv`, `king.kin0`, Abbildungen |
| `bericht.tex`, `bericht.pdf` | Bericht |
| `dokumentation.tex`, `dokumentation.pdf` | ausführliche Dokumentation |
| `literatur.bib` | Literatur |

## Benutzung

Umgebung mit conda anlegen (getestet mit Snakemake 9.27.0, GATK 4.6.2.0, bcftools 1.24, PLINK 2.0.0a.7.8):

```bash
conda create -n trio_test -c conda-forge -c bioconda \
  snakemake bwa samtools gatk4 bcftools plink2 matplotlib python=3.11
conda activate trio_test
```

Konfiguration anpassen: `config.example.yaml` nach `config.yaml` kopieren und die Pfade zum Referenzgenom (GRCh38, Ensembl Primary Assembly, mit BWA-, samtools- und GATK-Index) und zu den BAM-Dateien eintragen.

```bash
cp config.example.yaml config.yaml
mkdir -p data work
snakemake -n
snakemake -j 8 --resources mem_gb=16 --rerun-triggers mtime
python scripts/plot_results.py
```

## Grenzen

- Nur drei Chromosomen und vier Personen. Die Aussage gilt für diese Daten.
- Vater und Mutter lassen sich mit autosomalen Daten nicht unterscheiden, dafür wäre Chromosom X nötig.
- Die Genotypfilter (Tiefe ≥ 15, GQ ≥ 30) sind nicht kalibriert, und eine Basenqualitäts-Rekalibrierung (BQSR) fehlt.
- Die absoluten KINSHIP-Werte sind in dieser kleinen Gruppe nach oben verschoben. Belastbar sind IBS0 und die Rangfolge.

## Quellen

Eberle et al. 2017 (Platinum Genomes), Manichaikul et al. 2010 (KING), Chang et al. 2015 (PLINK 2), Mölder et al. 2021 (Snakemake), McKenna et al. 2010 (GATK), Li 2013 (BWA-MEM). Vollständige Angaben in `literatur.bib`.
