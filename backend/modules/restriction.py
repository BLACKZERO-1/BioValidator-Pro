from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes, Analysis
import math

def find_restriction_sites(sequence):
    """
    Module 2: Restriction Mapping (Pro Edition)
    Identifies cuts and classifies enzymes by frequency (Single vs Multi).
    Calculates Virtual Gel migration.
    """
    if not sequence:
        return {"status": "FAIL", "count": 0, "sites": [], "message": "No sequence."}

    seq_obj = Seq(sequence)
    
    # We use a standard analysis
    # Note: For speed in a demo, we might usually filter this list, 
    # but BioPython handles it reasonably well.
    rb = Analysis(AllEnzymes, seq_obj)
    full_result = rb.full()
    
    sites_data = []
    single_cutters = [] # The "Golden" enzymes
    double_cutters = []
    
    for enzyme, cuts in full_result.items():
        count = len(cuts)
        if count == 0:
            continue
            
        enz_name = str(enzyme)
        
        # Add to main list
        for cut in cuts:
            sites_data.append({
                "enzyme": enz_name,
                "position": cut,
                "frequency": count
            })
            
        # Classify
        if count == 1:
            single_cutters.append(enz_name)
        elif count == 2:
            double_cutters.append(enz_name)

    # Sort linear map by position
    sites_data.sort(key=lambda x: x['position'])

    # --- VIRTUAL GEL SIMULATION MATH ---
    # Calculates where the band appears on a 1% Agarose Gel.
    # Logic: Smaller DNA moves faster (further down).
    # We map logarithmic scale to 0-100% UI height.
    total_len = len(sequence)
    
    try:
        # Standard Ladder reference points: 10kb (top) to 0.5kb (bottom)
        # Formula: Inverse Logarithmic interpolation
        # High value (10000) -> Low Y (Top)
        # Low value (500) -> High Y (Bottom)
        
        min_log = math.log(500)   # Bottom of gel
        max_log = math.log(15000) # Top of gel
        val_log = math.log(max(total_len, 1))
        
        # Normalize to 0-1
        normalized = (val_log - min_log) / (max_log - min_log)
        
        # Invert (Big molecules stay at top/0%)
        gel_y = (1 - normalized) * 100
        
        # Clamp between 5% and 95% to stay inside the div
        gel_y = max(5, min(95, gel_y))
        
    except:
        gel_y = 50.0

    return {
        "status": "CLONING READY" if len(single_cutters) > 0 else "LIMITED",
        "message": f"{len(single_cutters)} unique cut sites found.",
        "count": len(sites_data),
        "sites": sites_data,
        "single_cutters": single_cutters,
        "double_cutters": double_cutters,
        "gel_pos": round(gel_y, 1),
        "total_len": total_len
    }