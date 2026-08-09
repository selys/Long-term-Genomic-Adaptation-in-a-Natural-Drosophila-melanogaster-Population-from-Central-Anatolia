#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH -J gene_length_control
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/gene_length_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/gene_length_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

OUT_DIR="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL"
mkdir -p $OUT_DIR

# -----------------------------------------------------------------------
# 1. Get gene coordinates/lengths for GCF_000001215.4_Release_6_plus_ISO1_MT
#    NOTE: SnpEff's BDGP6.32.105 database only stores a compiled binary
#    predictor (snpEffPredictor.bin) -- the source GTF/GFF is NOT retained
#    after the database is built, so we cannot get gene models from there.
#    Instead we use the annotation that matches the exact reference genome
#    used throughout this pipeline (same RefSeq assembly accession).
# -----------------------------------------------------------------------
echo ">>> Locating a gene GFF/GTF for GCF_000001215.4_Release_6_plus_ISO1_MT..."

BASE_DIR="/arf/scratch/ssenkal/selin/sekanslar"
REF_DIR="$BASE_DIR/RAW_DIR/dmel_reference_genome"

GENE_ANNOT=""

# a) Look right next to the .fna reference genome first (fast, scoped)
if [ -d "$REF_DIR" ]; then
    GENE_ANNOT=$(find "$REF_DIR" -maxdepth 2 \( -iname "*.gff*" -o -iname "*.gtf*" \) 2>/dev/null | head -n 1)
fi

# b) Broaden slightly within sekanslar/ if not found next to the .fna
if [ -z "$GENE_ANNOT" ]; then
    echo "Not found in $REF_DIR, checking within $BASE_DIR ..."
    GENE_ANNOT=$(find "$BASE_DIR" -maxdepth 3 \( -iname "*GCF_000001215*gff*" -o -iname "*GCF_000001215*gtf*" \) 2>/dev/null | head -n 1)
fi

# c) Last resort: download the matching NCBI RefSeq GFF3 for the SAME
#    assembly accession used to build the BAMs/VCF, so gene coordinates
#    are guaranteed consistent with the rest of the pipeline.
if [ -z "$GENE_ANNOT" ]; then
    echo "No local GFF/GTF found. Attempting to download the matching NCBI RefSeq GFF3..."
    mkdir -p "$REF_DIR"
    NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz"
    DEST="$REF_DIR/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz"
    if command -v wget &> /dev/null; then
        wget -q "$NCBI_URL" -O "$DEST"
    elif command -v curl &> /dev/null; then
        curl -s -L "$NCBI_URL" -o "$DEST"
    fi
    if [ -s "$DEST" ]; then
        GENE_ANNOT="$DEST"
        echo "Downloaded: $GENE_ANNOT"
    else
        echo "HATA: indirme de basarisiz (bu compute node'un internet erisimi olmayabilir)."
        echo "Cozum: login node'dan (arf-ui1) asagidaki komutla dosyayi manuel indirip"
        echo "  ayni yola koy, sonra job'u tekrar gonder:"
        echo "  wget '$NCBI_URL' -O '$DEST'"
        exit 1
    fi
fi

echo "Using annotation file: $GENE_ANNOT"

# Gene features: column 3 == "gene" in both GTF and GFF3
if [[ "$GENE_ANNOT" == *.gz ]]; then
    zcat "$GENE_ANNOT" | awk -F'\t' '$3=="gene"' > $OUT_DIR/genes_raw.gtf
else
    awk -F'\t' '$3=="gene"' "$GENE_ANNOT" > $OUT_DIR/genes_raw.gtf
fi
wc -l $OUT_DIR/genes_raw.gtf
if [ ! -s "$OUT_DIR/genes_raw.gtf" ]; then
    echo "HATA: genes_raw.gtf bos -- $GENE_ANNOT icinde 3. kolon 'gene' olan satir yok."
    echo "Dosyanin ilk birkac satirini kontrol et:"
    if [[ "$GENE_ANNOT" == *.gz ]]; then zcat "$GENE_ANNOT" | head -20; else head -20 "$GENE_ANNOT"; fi
    exit 1
fi

# -----------------------------------------------------------------------
# 2. Background SNP counts per gene (ALL tested SNPs, from full annotated VCF)
#    and year-significant SNP counts per gene (from year_sig_annotated.vcf.gz)
#    Both deduplicated by unique chr+pos+gene (multi-transcript inflation fix)
# -----------------------------------------------------------------------
FULL_VCF="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/mega_all_years_annotated.vcf"
if [ ! -f "$FULL_VCF" ]; then
    FULL_VCF="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/mega_all_years_annotated.vcf.gz"
fi
if [ ! -f "$FULL_VCF" ]; then
    echo "Neither .vcf nor .vcf.gz found at the expected path. Searching within all_years_26 ..."
    FULL_VCF=$(find /arf/scratch/ssenkal/selin/sekanslar/all_years_26 -maxdepth 1 \
               -iname "mega_all_years_annotated.vcf*" 2>/dev/null | head -n 1)
fi
if [ -z "$FULL_VCF" ] || [ ! -f "$FULL_VCF" ]; then
    echo "HATA: mega_all_years_annotated.vcf(.gz) hicbir yerde bulunamadi."
    echo "all_years_26 icindeki .vcf* dosyalarini listeliyorum:"
    find /arf/scratch/ssenkal/selin/sekanslar/all_years_26 -maxdepth 1 -iname "*.vcf*" 2>/dev/null
    exit 1
fi
echo "Using FULL_VCF: $FULL_VCF"

# --- Sanity check: is this VCF actually usably annotated? ---
# Known failure mode: if SnpEff was run BEFORE chromosomes were renamed
# from NCBI accessions (NT_033779.5 etc.) to 2L/2R/3L/3R/4, every record
# gets ANN="...|ERROR_CHROMOSOME_NOT_FOUND" with an EMPTY Gene_Name field.
# Detect that here instead of silently producing an empty background table.
echo ">>> Sanity-checking FULL_VCF annotation quality (first 20 records)..."
PROBE=$(bcftools query -f '%INFO/ANN\n' "$FULL_VCF" 2>/dev/null | head -20)
if echo "$PROBE" | grep -q "ERROR_CHROMOSOME_NOT_FOUND"; then
    echo "TESPIT: FULL_VCF, kromozomlar rename edilmeden once SnpEff'ten gecmis"
    echo "(ERROR_CHROMOSOME_NOT_FOUND -- Gene_Name alani bos). Bu dosya kullanilamaz."
    echo ">>> Ham (annotasyonsuz) VCF'i dogru veritabaniyla yeniden annotate ediyorum..."

    RAW_VCF="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/mega_all_years.vcf"
    if [ ! -f "$RAW_VCF" ]; then
        echo "HATA: ham VCF de bulunamadi: $RAW_VCF"
        echo "all_years_26 icindeki .vcf dosyalari:"
        find /arf/scratch/ssenkal/selin/sekanslar/all_years_26 -maxdepth 1 -iname "*.vcf*" 2>/dev/null
        exit 1
    fi

    # Confirm the raw VCF actually uses renamed chromosomes (2L/2R/...).
    # If not, rename it here using the EXACT SAME substitutions the
    # original pipeline's mpileup-rename step used, rather than failing --
    # this is a one-time text substitution, cheap even at this file size.
    RAW_CHR_CHECK=$(bcftools query -f '%CHROM\n' "$RAW_VCF" 2>/dev/null | head -1)
    echo "Ham VCF'teki ilk kromozom adi: $RAW_CHR_CHECK"
    case "$RAW_CHR_CHECK" in
        2L|2R|3L|3R|4|X|Y)
            echo "OK, isimlendirme SnpEff BDGP6.32.105 ile zaten uyumlu."
            ;;
        *)
            echo "Ham VCF rename edilmemis ($RAW_CHR_CHECK). Pipeline'in orijinal"
            echo "rename adimini (mpileup icin kullanilanin ayni) burada uyguluyorum..."
            RENAMED_VCF="$OUT_DIR/mega_all_years_renamed.vcf"
            sed -e 's/NT_033779\.5/2L/g' \
                -e 's/NT_033778\.4/2R/g' \
                -e 's/NT_037436\.4/3L/g' \
                -e 's/NT_033777\.3/3R/g' \
                -e 's/NC_004353\.4/4/g' \
                "$RAW_VCF" > "$RENAMED_VCF"

            RENAMED_CHECK=$(bcftools query -f '%CHROM\n' "$RENAMED_VCF" 2>/dev/null | head -1)
            echo "Rename sonrasi ilk kromozom adi: $RENAMED_CHECK"
            case "$RENAMED_CHECK" in
                2L|2R|3L|3R|4|X|Y) echo "OK, rename basarili." ;;
                *) echo "HATA: rename sonrasi da beklenen isim gelmedi ($RENAMED_CHECK)."
                   echo "Ham VCF'in kromozom isimlendirmesi beklenenden farkli olabilir."
                   exit 1 ;;
            esac
            RAW_VCF="$RENAMED_VCF"
            ;;
    esac

    # Locate the SAME SnpEff install that successfully produced
    # year_sig_annotated.vcf.gz (confirmed to have BDGP6.32.105 data).
    SNPEFF_JAR=""
    for CAND in \
        "/arf/scratch/ssenkal/selin/programs/miniconda3/envs/r_af/share/snpeff-5.2-3/snpEff.jar" \
        "/arf/scratch/ssenkal/selin/programs/envs/vcf2sync_env/share/snpeff-5.1-0/snpEff/snpEff.jar"
    do
        if [ -f "$CAND" ]; then SNPEFF_JAR="$CAND"; break; fi
    done
    if [ -z "$SNPEFF_JAR" ]; then
        echo "Known paths'te snpEff.jar bulunamadi, scoped arama yapiyorum..."
        SNPEFF_JAR=$(find /arf/scratch/ssenkal/selin/programs -maxdepth 6 -iname "snpEff.jar" 2>/dev/null | head -n 1)
    fi
    if [ -z "$SNPEFF_JAR" ]; then
        echo "HATA: snpEff.jar hicbir yerde bulunamadi."
        exit 1
    fi
    echo "Using SnpEff jar: $SNPEFF_JAR"

    REANNOT_VCF="$OUT_DIR/mega_all_years_reannotated.vcf"
    java -Xmx28g -jar "$SNPEFF_JAR" -v BDGP6.32.105 "$RAW_VCF" > "$REANNOT_VCF" 2> "$OUT_DIR/snpeff_reannotate.log"

    if [ ! -s "$REANNOT_VCF" ]; then
        echo "HATA: yeniden annotasyon basarisiz oldu, log:"
        tail -30 "$OUT_DIR/snpeff_reannotate.log"
        exit 1
    fi

    REPROBE=$(bcftools query -f '%INFO/ANN\n' "$REANNOT_VCF" 2>/dev/null | head -20)
    if echo "$REPROBE" | grep -q "ERROR_CHROMOSOME_NOT_FOUND"; then
        echo "HATA: yeniden annotasyondan sonra hala ERROR_CHROMOSOME_NOT_FOUND var. Manuel kontrol gerekli."
        exit 1
    fi
    echo "Yeniden annotasyon basarili. Gene_Name alanlari artik dolu."
    FULL_VCF="$REANNOT_VCF"
else
    echo "OK, FULL_VCF annotasyonu gecerli gorunuyor (ERROR_CHROMOSOME_NOT_FOUND yok)."
fi
echo "Final FULL_VCF: $FULL_VCF"

SIG_VCF="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/year_sig_annotated.vcf.gz"

echo ">>> Extracting background (all tested) SNP-gene pairs..."
# NOTE: stderr is NOT suppressed here on purpose -- a silent bcftools failure
# (wrong path, missing index, etc.) previously produced an empty file with
# no visible error. If this errors, the message will now show in the .err log.
# Restricted to autosomes only (2L,2R,3L,3R,4) -- matches the original
# mpileup filtering step in the main pipeline (X/Y were never included
# there), and keeps this analysis on exactly the same SNP universe as the
# rest of the manuscript ("SNPs on all autosomes").
AUTOSOMES='"2L","2R","3L","3R","4"'
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$FULL_VCF" | \
  awk -F'\t' -v autos="2L 2R 3L 3R 4" '
    BEGIN { split(autos, a, " "); for (i in a) keep[a[i]] = 1 }
    {
      if (!($1 in keep)) next;
      n=split($3, anns, ",");
      for(i=1;i<=n;i++){
        split(anns[i], f, "|");
        gene=f[4];
        if(gene!="") print $1"\t"$2"\t"gene;
      }
    }' | sort -u > $OUT_DIR/background_chr_pos_gene.txt

if [ ! -s "$OUT_DIR/background_chr_pos_gene.txt" ]; then
    echo "HATA: background_chr_pos_gene.txt bos kaldi -- bcftools basarisiz olmus olabilir."
    echo "Elle test et:"
    echo "  bcftools query -f '%CHROM\\t%POS\\t%INFO/ANN\\n' $FULL_VCF | head -5"
    exit 1
fi

echo ">>> Extracting year-significant SNP-gene pairs..."
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$SIG_VCF" | \
  awk -F'\t' -v autos="2L 2R 3L 3R 4" '
    BEGIN { split(autos, a, " "); for (i in a) keep[a[i]] = 1 }
    {
      if (!($1 in keep)) next;
      n=split($3, anns, ",");
      for(i=1;i<=n;i++){
        split(anns[i], f, "|");
        gene=f[4];
        if(gene!="") print $1"\t"$2"\t"gene;
      }
    }' | sort -u > $OUT_DIR/sig_chr_pos_gene.txt

echo "Background unique chr+pos+gene rows: $(wc -l < $OUT_DIR/background_chr_pos_gene.txt)"
echo "Year-sig unique chr+pos+gene rows:   $(wc -l < $OUT_DIR/sig_chr_pos_gene.txt)"

# -----------------------------------------------------------------------
# 3. R: merge gene length + background count + sig count, test the confound
# -----------------------------------------------------------------------
cat > $OUT_DIR/gene_length_analysis.R << 'RSCRIPT'
library(data.table)
library(ggplot2)

OUT_DIR <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL"

# --- Gene lengths from GTF/GFF3 ---
# NOTE: the reference genome here is NCBI RefSeq (GCF_000001215.4), so the
# attribute column may be GFF3-style ("gene=Ank2;...", "Name=Ank2;...")
# rather than GTF-style ('gene_name "Ank2";'). Try both, in order of
# preference, instead of assuming GTF quoting.
gtf <- fread(file.path(OUT_DIR, "genes_raw.gtf"), header = FALSE, sep = "\t", quote = "")
setnames(gtf, c("chr","source","feature","start","end","score","strand","frame","attr"))
gtf[, length_bp := end - start + 1]

extract_attr <- function(attr, key) {
  gtf_pat <- paste0(key, ' "([^"]+)"')          # GTF:  key "value";
  gff_pat <- paste0('(^|;)\\s*', key, '=([^;]+)') # GFF3: key=value;
  out <- rep(NA_character_, length(attr))
  hit_gtf <- grepl(gtf_pat, attr)
  out[hit_gtf] <- sub(paste0('.*', gtf_pat, '.*'), '\\1', attr[hit_gtf])
  hit_gff <- is.na(out) & grepl(gff_pat, attr)
  out[hit_gff] <- sub(paste0('.*', gff_pat, '.*'), '\\2', attr[hit_gff])
  out
}

# Try, in order: gene_name (GTF) -> gene (GFF3, most common on RefSeq) ->
# Name (GFF3 fallback) -> gene_id / ID (last resort, may be FBgn/LOC IDs)
gtf[, gene_name := extract_attr(attr, "gene_name")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "gene")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "Name")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "gene_id")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "ID")]

n_unparsed <- sum(is.na(gtf$gene_name))
cat("Gene records with NO parsable name (dropped):", n_unparsed, "of", nrow(gtf), "\n")
if (n_unparsed > 0) {
  cat("Example unparsed attribute strings:\n")
  print(head(gtf[is.na(gene_name), attr], 3))
}
gtf <- gtf[!is.na(gene_name)]

# --- Restrict to autosomes only ---
# Uses the SAME 5 RefSeq accessions the main pipeline's mpileup filtering
# step used (2L=NT_033779.5, 2R=NT_033778.4, 3L=NT_037436.4, 3R=NT_033777.3,
# 4=NC_004353.4). X/Y/mitochondrial genes are excluded so this analysis
# covers exactly the same SNP universe as the rest of the manuscript
# ("SNPs on all autosomes").
AUTOSOME_ACCESSIONS <- c("NT_033779.5","NT_033778.4","NT_037436.4","NT_033777.3","NC_004353.4")
n_before <- nrow(gtf)
gtf <- gtf[chr %in% AUTOSOME_ACCESSIONS]
cat("Gene records restricted to autosomes:", nrow(gtf), "of", n_before,
    "(dropped X/Y/mito/unplaced scaffolds)\n")
if (nrow(gtf) == 0) {
  cat("WARNING: 0 rows after autosome filter -- the GFF's chromosome/seqid naming\n")
  cat("does not match the expected accessions. Example chr values seen in GTF:\n")
  print(unique(fread(file.path(OUT_DIR, "genes_raw.gtf"), header = FALSE, sep = "\t", quote = "")$V1))
}

gene_lengths <- gtf[, .(length_bp = max(length_bp)), by = gene_name]  # collapse duplicate records per gene
cat("Unique autosomal genes with length info:", nrow(gene_lengths), "\n")

# --- Background & significant SNP-gene tables ---
bg  <- fread(file.path(OUT_DIR, "background_chr_pos_gene.txt"), header = FALSE,
             col.names = c("chr","pos","gene_name"))
sig <- fread(file.path(OUT_DIR, "sig_chr_pos_gene.txt"), header = FALSE,
             col.names = c("chr","pos","gene_name"))

if (nrow(bg) == 0) {
  stop("background_chr_pos_gene.txt loaded with 0 rows -- the bcftools extraction ",
       "step failed upstream (see the .sh script's own check, which should have ",
       "already caught this before R even started). Re-run the bash script and ",
       "inspect its output/error log.")
}

# --- Sanity check BEFORE merging: do gene symbols actually overlap? ---
# This is the check that would have caught a gene_name/gene_id mismatch
# (e.g. VCF giving "Ank2" while the GTF only yielded "FBgn0000083") before
# silently producing a near-empty merged table.
vcf_genes <- unique(bg$gene_name)
overlap_n <- length(intersect(vcf_genes, gene_lengths$gene_name))
overlap_pct <- round(100 * overlap_n / length(vcf_genes), 1)
cat(sprintf("\nGene symbol overlap check: %d / %d (%.1f%%) of VCF-annotated genes found in GTF length table\n",
            overlap_n, length(vcf_genes), overlap_pct))
if (overlap_pct < 50) {
  cat("WARNING: overlap below 50% -- gene_name/gene_id mismatch is likely.\n")
  cat("Example VCF gene symbols NOT found in GTF table:\n")
  print(head(setdiff(vcf_genes, gene_lengths$gene_name), 10))
  cat("Example GTF gene_name values (for comparison):\n")
  print(head(gene_lengths$gene_name, 10))
  cat("Fix extract_attr()/attribute key choice above before trusting downstream results.\n")
}

bg_counts  <- bg[,  .(n_tested = .N), by = gene_name]
sig_counts <- sig[, .(n_sig    = .N), by = gene_name]

merged <- merge(bg_counts, sig_counts, by = "gene_name", all.x = TRUE)
merged[is.na(n_sig), n_sig := 0]
merged <- merge(merged, gene_lengths, by = "gene_name")
merged <- merged[length_bp > 0 & n_tested > 0]
merged[, prop_sig := n_sig / n_tested]

cat("\nGenes in final merged table:", nrow(merged), "\n")

# --- Test 1: raw sig-SNP count vs gene length (expected to correlate trivially) ---
ct1 <- cor.test(merged$length_bp, merged$n_sig, method = "spearman")
cat("\n== Raw significant-SNP count vs gene length (Spearman) ==\n")
print(ct1)

# --- Test 2: proportion significant (normalized by tested SNPs) vs gene length ---
#     This is the key test: if length no longer predicts the PROPORTION,
#     the raw-count correlation is just a density artifact, not a real bias.
ct2 <- cor.test(merged$length_bp, merged$prop_sig, method = "spearman")
cat("\n== Proportion significant (n_sig/n_tested) vs gene length (Spearman) ==\n")
print(ct2)

# --- Test 3: Poisson regression, n_sig ~ log(length), offset = log(n_tested) ---
#     Directly asks: controlling for how many SNPs were even testable in the
#     gene (which itself scales with length), does length still matter?
merged[, log_length := log(length_bp)]
pois_fit <- glm(n_sig ~ log_length + offset(log(n_tested)), data = merged, family = poisson())
cat("\n== Poisson regression: n_sig ~ log(gene length), offset = log(n_tested) ==\n")
print(summary(pois_fit))

# --- Save table and plots ---
fwrite(merged, file.path(OUT_DIR, "gene_length_vs_significance.tsv"), sep = "\t")

p1 <- ggplot(merged, aes(x = length_bp, y = n_sig)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_x_log10() +
  labs(title = "Raw significant SNP count vs gene length",
       x = "Gene length (bp, log10 scale)", y = "Year-significant SNPs per gene") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "raw_count_vs_length.pdf"), p1, width = 7, height = 5)

p2 <- ggplot(merged, aes(x = length_bp, y = prop_sig)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_x_log10() +
  labs(title = "Proportion significant (normalized) vs gene length",
       x = "Gene length (bp, log10 scale)", y = "Proportion of tested SNPs that are year-significant") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "proportion_vs_length.pdf"), p2, width = 7, height = 5)

sink(file.path(OUT_DIR, "gene_length_control_summary.txt"))
cat("Gene-length confound check\n===========================\n\n")
cat("Genes analyzed:", nrow(merged), "\n\n")
cat("Test 1 - raw count vs length (Spearman rho, p):\n")
cat(sprintf("  rho = %.4f, p = %.3g\n\n", ct1$estimate, ct1$p.value))
cat("Test 2 - proportion significant vs length (Spearman rho, p):\n")
cat(sprintf("  rho = %.4f, p = %.3g\n\n", ct2$estimate, ct2$p.value))
cat("Test 3 - Poisson regression coefficient on log(length), offset=log(n_tested):\n")
print(coef(summary(pois_fit)))
sink()

cat("\n>>> Done. Outputs in", OUT_DIR, "\n")
RSCRIPT

Rscript $OUT_DIR/gene_length_analysis.R

echo ">>> Job finished: $(date)"
