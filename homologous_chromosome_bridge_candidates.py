# Script combinado Bash + Python para priorizar proteínas que faciliten interacciones entre cromosomas homólogos en trigo
# ==========================================
# PARTE 1: QC y eliminación de redundancia
# ==========================================
# Filtrar secuencias menores de 50 aa
# seqkit seq -m 50 wheat_proteome.fasta -o wheat_filtered.fasta

# Agrupar isoformas/redundancias
# cd-hit -i wheat_filtered.fasta -o wheat_nr.fasta -c 0.98 -n 5 -T 8 -M 16000

# ==========================================
# PARTE 2: Anotación de dominios (InterProScan y HMMER)
# ==========================================
# InterProScan (local)
# interproscan.sh -i wheat_nr.fasta -f tsv -dp -appl Pfam,SMART,CDD,Prosite -o interpro.tsv

# HMMER vs Pfam-A
# hmmscan --cpu 8 --domtblout hmmscan.tbl Pfam-A.hmm wheat_nr.fasta
# Dominios a priorizar: Myb_DNA-binding, Myb_telobox, HMG_box, OB_fold, coiled-coil, Leucine zipper, SAM, TPR, ARM/HEAT

# ==========================================
# PARTE 3: Predicción de coiled-coil y dimerización
# ==========================================
# deepcoil.py wheat_nr.fasta > coils.tsv
# Interpretación: regiones con Coil_score > 0.5 indicarán capacidad de dimerización/oligomerización

# ==========================================
# PARTE 4: Búsqueda de homólogos remotos
# ==========================================
# JackHMMER / PSI-BLAST
# jackhmmer -N 5 -E 1e-3 --tblout jack_trb.tbl TRB_seed.fasta wheat_nr.fasta
# psiblast -query TRB_seed.fasta -db wheat_nr.fasta -num_iterations 5 -evalue 1e-3 -out psiblast_trb.out
# HHpred para similitudes estructurales remotas

# ==========================================
# PARTE 5: Predicción de unión a ADN
# ==========================================
# DP-Bind / DNABind / TransBind
# mapear residuos que podrían interactuar con telómeros

# ==========================================
# PARTE 6: Localización nuclear
# ==========================================
# cNLS Mapper / NucPred / DeepLoc
# marcar candidatos nucleares para asegurar potencial contacto con cromosomas homólogos

# ==========================================
# PARTE 7: Coexpresión y datos meiosis
# ==========================================
# expr = pd.read_csv('meiosis_expression.csv')  # columnas: Protein, TPM_meiosis
# Priorizar proteínas expresadas durante leptotene/zygotene, fase de formación del bouquet

# ==========================================
# PARTE 8: Integración y scoring (Python)
# ==========================================
import pandas as pd

# Cargar datos
interpro = pd.read_csv('interpro.tsv', sep='\t', header=None, names=['Protein', 'MD5', 'SeqType','Start','End','Score','Status','Date','Domain'])
coils = pd.read_csv('coils.tsv')  # ajustar columnas según output
expr = pd.read_csv('meiosis_expression.csv')  # columnas: Protein, TPM_meiosis

candidates = pd.DataFrame()
candidates['Protein'] = interpro['Protein'].unique()

# Función de scoring para puentes homólogos
def score_homolog_bridge(prot):
    # Evidencia de unión a ADN/telómero
    S_DNA = 1.0 if any(d in interpro[interpro['Protein']==prot]['Domain'].values for d in ['Myb_DNA-binding','Myb_telobox','HMG_box','OB_fold']) else 0.0
    # Evidencia de dimerización/oligomerización
    coil_score = coils[coils['Protein']==prot]['Coil_score'].max() if prot in coils['Protein'].values else 0.0
    S_DIMER = 1.0 if coil_score > 0.5 else 0.0
    # Expresión en meiosis (normalizada 0-1)
    S_MEIO = expr[expr['Protein']==prot]['TPM_meiosis'].values[0] / expr['TPM_meiosis'].max() if prot in expr['Protein'].values else 0.0
    # Localización nuclear proxy (si está en expresión meiosis, asumimos nuclear)
    S_LOC = 1.0 if prot in expr['Protein'].values else 0.0
    # Conservación / duplicación placeholder
    S_CONS = 1.0  # se puede actualizar si se tiene información filogenética
    # Score final ponderado
    score = 0.35*S_DNA + 0.35*S_DIMER + 0.15*S_MEIO + 0.10*S_LOC + 0.05*S_CONS
    return score

candidates['Score_final'] = candidates['Protein'].apply(score_homolog_bridge)
candidates['Prioridad'] = pd.cut(candidates['Score_final'], bins=[0,0.45,0.7,1.0], labels=['Bajo','Medio','Alto'])

# Guardar tabla final
candidates.to_csv('homolog_bridge_candidates_table.csv', index=False)

# ==========================================
# PARTE 9: Extraer FASTA de candidatos TOP para AF2-Multimer
# ==========================================
top_candidates = candidates[candidates['Prioridad']=='Alto']['Protein']
with open('top_homolog_bridge_candidates.fasta','w') as out_f:
    with open('wheat_nr.fasta','r') as fasta:
        write_flag = False
        for line in fasta:
            if line.startswith('>'):
                prot_id = line[1:].split()[0]
                write_flag = prot_id in top_candidates.values
            if write_flag:
                out_f.write(line)

# ==========================================
# NOTAS PARA EL USUARIO
# - Ajustar thresholds de Coil_score según la salida de DeepCoil.
# - Revisar dominios en interpro.tsv para incluir otros posibles dominios de unión a telómero.
# - Para AF2-Multimer: preparar combinaciones de candidatos top (homodímeros y heterodímeros) para modelar puentes.
# - Se recomienda filtrar duplicados y revisar expresión subgenómica A/B/D para trigo.
# - La columna Prioridad permite seleccionar candidatos Top para modelado y experimentos posteriores.