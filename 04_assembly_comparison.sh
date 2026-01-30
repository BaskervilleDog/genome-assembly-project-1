#!/bin/bash
#==============================================================================
# ASSEMBLY COMPARISON PIPELINE
#==============================================================================
# Autor: Gianlucca de Urzêda Alves
# Data: 30/01/2026
# Descrição: Compara todas as assemblies geradas e recomenda a melhor
#==============================================================================

set -e
set -u
set -o pipefail

#==============================================================================
# CONFIGURAÇÕES
#==============================================================================

# Assemblies para comparar
FLYE_ASSEMBLY="./flye_assembly/pilon_round2/pilon_round2.fasta"
FLYE_SCAFFOLDED="./scaffolded_assembly/scaffold/ragtag.scaffold.fasta"
RAVEN_ASSEMBLY="./raven_assembly/raven_assembly.fasta"
WTDBG2_ASSEMBLY="./wtdbg2_assembly/wtdbg2_assembly.fasta"

# Parâmetros
THREADS=12
SPECIES="Trichoderma_harzianum"

# Output
OUTPUT_DIR="./assembly_comparison"
LOG_DIR="${OUTPUT_DIR}/logs"

#==============================================================================
# FUNÇÕES
#==============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*"
}

log_section() {
    echo ""
    echo "========================================================================"
    echo "$*"
    echo "========================================================================"
}

#==============================================================================
# VERIFICAÇÕES INICIAIS
#==============================================================================

log_section "ASSEMBLY COMPARISON PIPELINE - START"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"

log_info "Checking available assemblies..."

ASSEMBLIES=()
LABELS=()

if [[ -f "${FLYE_ASSEMBLY}" ]]; then
    ASSEMBLIES+=("${FLYE_ASSEMBLY}")
    LABELS+=("Flye_Pilon")
    log_info "✓ Flye+Pilon assembly found"
fi

if [[ -f "${FLYE_SCAFFOLDED}" ]]; then
    ASSEMBLIES+=("${FLYE_SCAFFOLDED}")
    LABELS+=("Flye_Scaffolded")
    log_info "✓ Flye scaffolded assembly found"
fi

if [[ -f "${RAVEN_ASSEMBLY}" ]]; then
    ASSEMBLIES+=("${RAVEN_ASSEMBLY}")
    LABELS+=("Raven")
    log_info "✓ Raven assembly found"
fi

if [[ -f "${WTDBG2_ASSEMBLY}" ]]; then
    ASSEMBLIES+=("${WTDBG2_ASSEMBLY}")
    LABELS+=("wtdbg2")
    log_info "✓ wtdbg2 assembly found"
fi

if [[ ${#ASSEMBLIES[@]} -lt 2 ]]; then
    log_info "Found only ${#ASSEMBLIES[@]} assembly(ies)"
    log_info "Need at least 2 assemblies to compare!"
    log_info "Run assembly scripts first"
    exit 1
fi

log_info "Found ${#ASSEMBLIES[@]} assemblies to compare"

#==============================================================================
# ETAPA 1: QUAST COMPARATIVE
#==============================================================================

log_section "STEP 1/3: QUAST COMPARATIVE ANALYSIS"

log_info "Running QUAST on all assemblies..."
log_info "  Comparing: ${LABELS[@]}"

mamba run -n quast quast \
    "${ASSEMBLIES[@]}" \
    -o "${OUTPUT_DIR}/quast_comparison" \
    --labels "$(IFS=,; echo "${LABELS[*]}")" \
    --threads ${THREADS} \
    --fungus \
    2>&1 | tee "${LOG_DIR}/quast_comparison.log"

log_info "✓ QUAST comparison complete"
log_info "  Report: ${OUTPUT_DIR}/quast_comparison/report.html"

#==============================================================================
# ETAPA 2: COLLECT BUSCO RESULTS
#==============================================================================

log_section "STEP 2/3: COLLECTING BUSCO RESULTS"

BUSCO_DIR="${OUTPUT_DIR}/busco_comparison"
mkdir -p "${BUSCO_DIR}"

log_info "Collecting BUSCO results from individual assemblies..."

# Coletar BUSCO de cada assembly
for i in "${!ASSEMBLIES[@]}"; do
    label="${LABELS[$i]}"
    
    # Determinar caminho do BUSCO
    case ${label} in
        "Flye_Pilon")
            busco_path="./flye_assembly/quality_assessment/busco"
            ;;
        "Flye_Scaffolded")
            busco_path="./scaffolded_assembly/quality_assessment/busco"
            ;;
        "Raven")
            busco_path="./raven_assembly/quality_assessment/busco"
            ;;
        "wtdbg2")
            busco_path="./wtdbg2_assembly/quality_assessment/busco"
            ;;
    esac
    
    if [[ -d "${busco_path}" ]]; then
        summary=$(find "${busco_path}" -name "short_summary*.txt" | head -1)
        if [[ -f "${summary}" ]]; then
            cp "${summary}" "${BUSCO_DIR}/${label}_busco_summary.txt"
            log_info "✓ Collected BUSCO for ${label}"
        else
            log_info "⚠ BUSCO summary not found for ${label}"
            log_info "  Run BUSCO on this assembly first"
        fi
    else
        log_info "⚠ BUSCO results not found for ${label}"
    fi
done

#==============================================================================
# ETAPA 3: GENERATE FINAL COMPARISON REPORT
#==============================================================================

log_section "STEP 3/3: GENERATING FINAL COMPARISON REPORT"

REPORT="${OUTPUT_DIR}/FINAL_COMPARISON_REPORT.txt"

cat > "${REPORT}" << EOF
================================================================================
GENOME ASSEMBLY COMPARISON REPORT
================================================================================
Species: ${SPECIES}
Date: $(date)
Assemblies compared: ${#ASSEMBLIES[@]}

ASSEMBLIES:
-----------
EOF

for i in "${!ASSEMBLIES[@]}"; do
    echo "$((i+1)). ${LABELS[$i]}: ${ASSEMBLIES[$i]}" >> "${REPORT}"
done

cat >> "${REPORT}" << EOF

================================================================================
QUAST SUMMARY TABLE
================================================================================

EOF

# Extrair tabela do QUAST
if [[ -f "${OUTPUT_DIR}/quast_comparison/report.tsv" ]]; then
    # Cabeçalho
    head -1 "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
    
    # Métricas importantes
    grep "# contigs" "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
    grep "Total length" "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
    grep "^N50" "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
    grep "^L50" "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
    grep "GC (%)" "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
    grep "# N's per 100 kbp" "${OUTPUT_DIR}/quast_comparison/report.tsv" >> "${REPORT}"
fi

cat >> "${REPORT}" << EOF

================================================================================
BUSCO COMPARISON
================================================================================

EOF

# Adicionar BUSCO de cada assembly
for i in "${!LABELS[@]}"; do
    label="${LABELS[$i]}"
    busco_file="${BUSCO_DIR}/${label}_busco_summary.txt"
    
    if [[ -f "${busco_file}" ]]; then
        echo "${label}:" >> "${REPORT}"
        grep "C:" "${busco_file}" >> "${REPORT}"
        echo "" >> "${REPORT}"
    fi
done

cat >> "${REPORT}" << EOF

================================================================================
DETAILED BUSCO RESULTS
================================================================================

EOF

for i in "${!LABELS[@]}"; do
    label="${LABELS[$i]}"
    busco_file="${BUSCO_DIR}/${label}_busco_summary.txt"
    
    if [[ -f "${busco_file}" ]]; then
        echo "${label}:" >> "${REPORT}"
        echo "$(printf '%.0s-' {1..80})" >> "${REPORT}"
        grep -A 8 "Results:" "${busco_file}" >> "${REPORT}"
        echo "" >> "${REPORT}"
    fi
done

#==============================================================================
# RECOMENDAÇÃO AUTOMÁTICA
#==============================================================================

cat >> "${REPORT}" << EOF

================================================================================
RECOMMENDATION - HOW TO CHOOSE THE BEST ASSEMBLY
================================================================================

DECISION HIERARCHY:
-------------------
1. BUSCO Complete (%)           ← MOST IMPORTANT
2. N50 (Mb)                     ← If BUSCO similar
3. BUSCO Duplicated (%)         ← Lower is better
4. Number of contigs/scaffolds  ← Lower is better
5. Total length                 ← Should match expected size

QUALITY THRESHOLDS FOR FUNGI (~40 Mb):
---------------------------------------
BUSCO Complete:
  ✓ >99%      = Exceptional (reference-quality)
  ✓ 95-99%    = Excellent (publication-ready)
  ⚠ 90-95%    = Good (acceptable)
  ✗ <90%      = Poor (reprocess recommended)

N50:
  ✓ >2 Mb     = Excellent
  ✓ 500kb-2Mb = Good
  ⚠ 100-500kb = Moderate
  ✗ <100kb    = Fragmented

Contigs/Scaffolds:
  ✓ <30       = Excellent (chromosome-level)
  ✓ 30-100    = Good
  ⚠ 100-500   = Moderate
  ✗ >500      = Fragmented

DECISION RULES:
---------------
Rule 1: If one assembly has BUSCO >1% higher → Choose it
Rule 2: If BUSCO within 1% → Choose highest N50
Rule 3: If BUSCO and N50 similar → Choose lowest duplication
Rule 4: If all similar → Choose Flye (most reliable)

TYPICAL RESULTS FOR YOUR ORGANISM:
-----------------------------------
Expected ranking (from experience):
1. Flye+Pilon:        BUSCO 99-100%, N50 3-5 Mb    ← Usually best
2. Flye+Scaffolded:   BUSCO 99-100%, N50 4-7 Mb    ← If scaffolding worked
3. Raven:             BUSCO 98-99%,  N50 3-5 Mb    ← Good alternative
4. wtdbg2 (draft):    BUSCO 90-95%,  N50 2-4 Mb    ← Needs polishing

SPECIFIC RECOMMENDATIONS:
-------------------------
EOF

# Análise automática (simplificada)
best_assembly=""
best_busco=0

for i in "${!LABELS[@]}"; do
    label="${LABELS[$i]}"
    busco_file="${BUSCO_DIR}/${label}_busco_summary.txt"
    
    if [[ -f "${busco_file}" ]]; then
        busco_complete=$(grep "C:" "${busco_file}" | awk -F'[:%]' '{print $2}' | tr -d ' ')
        
        # Remover vírgulas se houver
        busco_complete=$(echo "${busco_complete}" | sed 's/,/./g')
        
        # Comparação simplificada
        if (( $(echo "${busco_complete} > ${best_busco}" | bc -l 2>/dev/null || echo 0) )); then
            best_busco=${busco_complete}
            best_assembly=${label}
        fi
    fi
done

if [[ -n "${best_assembly}" ]]; then
    cat >> "${REPORT}" << EOF
Based on BUSCO completeness, the recommended assembly is:

  → ${best_assembly} (BUSCO: ${best_busco}%)

However, please verify:
1. Open QUAST report: ${OUTPUT_DIR}/quast_comparison/report.html
2. Compare N50 values
3. Check for excessive duplication
4. Verify total length is reasonable (~40 Mb)

EOF
else
    cat >> "${REPORT}" << EOF
Unable to determine best assembly automatically.
Please review:
1. QUAST report: ${OUTPUT_DIR}/quast_comparison/report.html
2. BUSCO results above
3. Choose based on decision hierarchy

EOF
fi

cat >> "${REPORT}" << EOF

================================================================================
OUTPUT FILES
================================================================================
QUAST interactive report: ${OUTPUT_DIR}/quast_comparison/report.html
QUAST text report:        ${OUTPUT_DIR}/quast_comparison/report.txt
BUSCO summaries:          ${BUSCO_DIR}/
This report:              ${REPORT}

TO VIEW INTERACTIVE REPORT:
   firefox ${OUTPUT_DIR}/quast_comparison/report.html

================================================================================
NEXT STEPS
================================================================================
1. ✓ Review this report
2. ✓ Open QUAST HTML report (interactive charts and tables)
3. ✓ Choose best assembly based on metrics
4. ✓ Optional: Run Kraken2 contamination check on chosen assembly
5. ✓ Proceed with genome annotation (Funannotate, Augustus, or MAKER)

================================================================================
EOF

cat "${REPORT}"
log_info "Report saved: ${REPORT}"

#==============================================================================
# CRIAR SCRIPT DE VISUALIZAÇÃO
#==============================================================================

VIEW_SCRIPT="${OUTPUT_DIR}/view_results.sh"

cat > "${VIEW_SCRIPT}" << 'EOFSCRIPT'
#!/bin/bash
# Script para visualizar resultados da comparação

echo "Opening QUAST comparison report..."

if command -v firefox &> /dev/null; then
    firefox quast_comparison/report.html &
elif command -v google-chrome &> /dev/null; then
    google-chrome quast_comparison/report.html &
elif command -v chromium &> /dev/null; then
    chromium quast_comparison/report.html &
else
    echo "Browser not found. Please open manually:"
    echo "  $(pwd)/quast_comparison/report.html"
fi

echo ""
echo "Displaying comparison report:"
cat FINAL_COMPARISON_REPORT.txt
EOFSCRIPT

chmod +x "${VIEW_SCRIPT}"

log_info "View script created: ${VIEW_SCRIPT}"

#==============================================================================
# FINALIZAÇÃO
#==============================================================================

log_section "ASSEMBLY COMPARISON COMPLETE!"

log_info "Summary:"
log_info "  Assemblies compared: ${#ASSEMBLIES[@]}"
log_info "  QUAST report: ${OUTPUT_DIR}/quast_comparison/report.html"
log_info "  Final report: ${REPORT}"
echo ""
log_info "To view results:"
log_info "  cd ${OUTPUT_DIR} && ./view_results.sh"
echo ""
log_info "Or manually:"
log_info "  firefox ${OUTPUT_DIR}/quast_comparison/report.html"
log_info "  cat ${REPORT}"
echo ""
log_info "Pipeline completed at: $(date)"

exit 0