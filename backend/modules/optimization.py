from Bio.Seq import Seq
import math

# E. coli Class II (Highly Expressed Genes) Codon Usage Table (Approximate weights)
# Weight (w) = Frequency / Max Frequency for that AA
CODON_WEIGHTS = {
    'GCA': 0.586, 'GCC': 0.122, 'GCG': 1.000, 'GCT': 0.407, # Ala
    'TGC': 1.000, 'TGT': 0.444, # Cys
    'GAC': 1.000, 'GAT': 0.434, # Asp
    'GAA': 1.000, 'GAG': 0.259, # Glu
    'TTC': 1.000, 'TTT': 0.296, # Phe
    'GGA': 0.060, 'GGC': 0.603, 'GGG': 0.207, 'GGT': 1.000, # Gly
    'CAC': 1.000, 'CAT': 0.291, # His
    'ATA': 0.003, 'ATC': 1.000, 'ATT': 0.185, # Ile
    'AAA': 1.000, 'AAG': 0.253, # Lys
    'CTA': 0.007, 'CTC': 0.037, 'CTG': 1.000, 'CTT': 0.042, 'TTA': 0.020, 'TTG': 0.020, # Leu
    'ATG': 1.000, # Met (Start)
    'AAC': 1.000, 'AAT': 0.051, # Asn
    'CCA': 0.135, 'CCC': 0.012, 'CCG': 1.000, 'CCT': 0.070, # Pro
    'CAA': 0.163, 'CAG': 1.000, # Gln
    'AGA': 0.004, 'AGG': 0.002, 'CGA': 0.006, 'CGC': 0.356, 'CGG': 0.008, 'CGT': 1.000, # Arg
    'AGC': 0.410, 'AGT': 0.085, 'TCA': 0.077, 'TCC': 0.179, 'TCG': 0.068, 'TCT': 1.000, # Ser
    'ACA': 0.076, 'ACC': 1.000, 'ACG': 0.297, 'ACT': 0.347, # Thr
    'GTA': 0.495, 'GTC': 0.066, 'GTG': 0.221, 'GTT': 1.000, # Val
    'TGG': 1.000, # Trp
    'TAC': 1.000, 'TAT': 0.239, # Tyr
    'TAA': 1.000, 'TAG': 0.000, 'TGA': 0.290 # Stop
}

# Optimal Codon Map (The "Best" codon for each AA)
OPTIMAL_CODONS = {
    'A': 'GCG', 'C': 'TGC', 'D': 'GAC', 'E': 'GAA', 'F': 'TTC',
    'G': 'GGT', 'H': 'CAC', 'I': 'ATC', 'K': 'AAA', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT',
    'S': 'TCT', 'T': 'ACC', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAC', '*': 'TAA'
}

def get_codon_weight(codon):
    return CODON_WEIGHTS.get(codon, 0.1)

def calculate_cai(sequence):
    """Calculates Codon Adaptation Index (Geometric Mean of Weights)"""
    if len(sequence) < 3: return 0.0
    
    weights = []
    # Process in chunks of 3
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3]
        w = get_codon_weight(codon)
        weights.append(w)
        
    if not weights: return 0.0
    
    # Geometric Mean = exp( (1/N) * sum(ln(w)) )
    log_sum = sum(math.log(w) if w > 0 else -9.0 for w in weights)
    cai = math.exp(log_sum / len(weights))
    return round(cai, 2)

def generate_velocity_profile(sequence):
    """Generates a list of weights representing ribosome speed/efficiency"""
    return [get_codon_weight(sequence[i:i+3]) for i in range(0, len(sequence)-2, 3)]

def optimize_sequence(sequence):
    """
    Module 3: Optimization
    Returns Optimized DNA, Comparison Metrics, and Velocity Data.
    """
    # 1. Translate
    try:
        protein = Seq(sequence).translate()
    except:
        return {"error": "Invalid Sequence"}

    # 2. Optimize (Simple: Replace with Optimal Codon)
    optimized_dna = ""
    for aa in protein:
        optimized_dna += OPTIMAL_CODONS.get(aa, 'NNN')

    # 3. Calculate Metrics (Before vs After)
    cai_original = calculate_cai(sequence)
    cai_optimized = calculate_cai(optimized_dna)
    
    # 4. Generate Velocity Data (For Graph)
    # We limit data points to 100 max for chart performance
    vel_orig = generate_velocity_profile(sequence)
    vel_opt = generate_velocity_profile(optimized_dna)
    
    # Downsample if too long (simple skip)
    step = max(1, len(vel_orig) // 100)
    
    graph_data = {
        "labels": [i for i in range(1, len(vel_orig), step)],
        "original": vel_orig[::step],
        "optimized": vel_opt[::step]
    }
    
    # 5. GC Balance Check
    from Bio.SeqUtils import gc_fraction
    gc_orig = round(gc_fraction(sequence) * 100, 1)
    gc_opt = round(gc_fraction(optimized_dna) * 100, 1)

    return {
        "original_dna": sequence,
        "optimized_dna": optimized_dna,
        "cai_before": cai_original,
        "cai_after": cai_optimized,
        "gc_before": gc_orig,
        "gc_after": gc_opt,
        "velocity_graph": graph_data
    }