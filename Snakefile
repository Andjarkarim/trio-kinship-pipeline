configfile: "config.yaml"

SAMPLES = list(config["samples"])
CHROMS = [str(c) for c in config["chroms"]] or [str(i) for i in range(1, 23)]
REF = config["reference"]
ALIGNER = config.get("aligner", "bwa")


def bam_of(w):
    # fertige BAM aus der Config, sonst die selbst erzeugte
    return config["samples"][w.sample].get("bam") or f"work/{w.sample}.md.bam"


wildcard_constraints:
    sample="|".join(SAMPLES),
    chrom="|".join(CHROMS),


rule all:
    input: "results/king.kin0", "results/mendel.tsv"


def ena_url(acc):
    num = acc[3:]
    sub = {6: "", 7: f"/00{num[-1]}", 8: f"/0{num[-2:]}", 9: f"/{num[-3:]}"}[len(num)]
    return f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{acc[:6]}{sub}/{acc}"


rule download_fastq:
    output:
        temp("data/{sample}_{read}.fastq.gz"),
    wildcard_constraints:
        read="1|2",
    params:
        url=lambda w: ena_url(config["samples"][w.sample]["accession"]).replace("ftp://", "https://"),
        acc=lambda w: config["samples"][w.sample]["accession"],
    shell:
        """
        n=0
        until curl -fL -C - --connect-timeout 60 --speed-limit 50000 --speed-time 120 -o {output}.part {params.url}/{params.acc}_{wildcards.read}.fastq.gz; do
            n=$((n+1))
            if [ $n -ge 50 ]; then echo "Download nach 50 Versuchen aufgegeben" >&2; exit 1; fi
            echo "Versuch $n fehlgeschlagen, warte 30 s und setze fort ..." >&2
            sleep 30
        done
        gzip -t < {output}.part
        mv {output}.part {output}
        """


rule align:
    input:
        r1="data/{sample}_1.fastq.gz",
        r2="data/{sample}_2.fastq.gz",
    output: bam=temp("work/{sample}.sorted.bam")
    threads: config["threads"]
    resources: mem_gb=config.get("align_mem_gb", 8)
    params:
        ref=REF, aligner=ALIGNER,
        rg=lambda w: f"@RG\\tID:{w.sample}\\tSM:{w.sample}\\tPL:ILLUMINA\\tLB:{w.sample}",
    shell:
        "{params.aligner} mem -t {threads} -R '{params.rg}' {params.ref} {input.r1} {input.r2} "
        "| samtools sort -@ 4 -m 2G -o {output.bam} -"


rule markdup:
    input: "work/{sample}.sorted.bam"
    output:
        bam=temp("work/{sample}.md.bam"),
        bai=temp("work/{sample}.md.bam.bai"),
        metrics="qc/{sample}.dup_metrics.txt",
    resources: mem_gb=8
    shell:
        "gatk --java-options '-Xmx6g' MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} "
        "&& samtools index {output.bam}"


rule haplotypecaller:
    input:
        bam=bam_of,
        bai=lambda w: bam_of(w) + ".bai",
    output: "gvcf/{sample}.{chrom}.g.vcf.gz"
    threads: 2
    resources: mem_gb=8
    params: ref=REF
    shell:
        "gatk --java-options '-Xmx6g' HaplotypeCaller -R {params.ref} -I {input.bam} "
        "-L {wildcards.chrom} -O {output} -ERC GVCF --native-pair-hmm-threads {threads}"


rule combine_gvcfs:
    input: expand("gvcf/{sample}.{{chrom}}.g.vcf.gz", sample=SAMPLES)
    output: "joint/{chrom}.g.vcf.gz"
    resources: mem_gb=8
    params:
        ref=REF,
        v=lambda w, input: " ".join(f"-V {f}" for f in input),
    shell:
        "gatk --java-options '-Xmx6g' CombineGVCFs -R {params.ref} {params.v} "
        "-L {wildcards.chrom} -O {output}"


rule genotype:
    input: "joint/{chrom}.g.vcf.gz"
    output: "joint/{chrom}.vcf.gz"
    resources: mem_gb=8
    params: ref=REF
    shell:
        "gatk --java-options '-Xmx6g' GenotypeGVCFs -R {params.ref} -V {input} "
        "-L {wildcards.chrom} -O {output}"


rule filter_snps:
    input: "joint/{chrom}.vcf.gz"
    output: "filtered/{chrom}.snps.vcf.gz"
    params: dp=config["min_dp"], gq=config["min_gq"]
    shell:
        "bcftools +setGT {input} -- -t q -n . -i 'FMT/DP<{params.dp} || FMT/GQ<{params.gq}' "
        "| bcftools view -m2 -M2 -v snps -i 'QUAL>30 && F_MISSING=0' -Oz -o {output} "
        "&& bcftools index -t {output}"


rule concat:
    input: expand("filtered/{chrom}.snps.vcf.gz", chrom=CHROMS)
    output: "results/snps.vcf.gz"
    shell: "bcftools concat -Oz -o {output} {input} && bcftools index -t {output}"


rule king:
    input: "results/snps.vcf.gz"
    output: "results/king.kin0"
    params: out="results/king"
    shell:
        """
        plink2 --vcf {input} --double-id --vcf-half-call m \
               --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 50 \
               --make-bed --out work/cohort
        plink2 --bfile work/cohort --mac 1 --make-king-table --out {params.out}
        """


rule mendel:
    input: "results/snps.vcf.gz"
    output: "results/mendel.tsv"
    shell: "python scripts/mendel_check.py --vcf {input} > {output}"
