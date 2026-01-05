from Bio.SeqUtils import gc_fraction
import collections

def check_gc_content(sequence):
    """Check 1.2: GC Content Analysis"""
    gc_percent = gc_fraction(sequence) * 100
    status = "PASS"
    risk = "LOW"
    message = "Optimal"
    
    if gc_percent < 25 or gc_percent > 75:
        status = "FAIL"
        risk = "HIGH"
        message = "Critical Extreme"
    elif gc_percent < 35 or gc_percent > 65:
        status = "WARNING"
        risk = "MODERATE"
        message = "Suboptimal"
        
    return {"value": round(gc_percent, 2), "status": status, "risk": risk, "message": message}

def check_homopolymers(sequence, threshold=6):
    """Check 1.3: Find runs of 6+ identical nucleotides."""
    homopolymers = []
    current_char = ""
    current_count = 0
    start_pos = 0

    for i, char in enumerate(sequence):
        if char == current_char:
            current_count += 1
        else:
            if current_count >= threshold:
                homopolymers.append({
                    "start": start_pos + 1, "end": i, "length": current_count, "base": current_char
                })
            current_char = char
            current_count = 1
            start_pos = i
            
    if current_count >= threshold:
        homopolymers.append({
            "start": start_pos + 1, "end": len(sequence), "length": current_count, "base": current_char
        })

    status = "FAIL" if homopolymers else "PASS"
    risk = "HIGH" if homopolymers else "LOW"
    
    return {"status": status, "risk": risk, "count": len(homopolymers), "details": homopolymers, "message": "Detected" if homopolymers else "None"}

def check_length(sequence):
    """
    Check 1.1: Length Validation
    UPDATED: Raised limit to 10,000bp.
    """
    length = len(sequence)
    status = "PASS"
    risk = "LOW"
    message = "Optimal"

    if length < 100:
        status = "FAIL"
        risk = "HIGH"
        message = "Too Short (<100bp)"
    elif length > 10000:
        status = "WARNING"
        risk = "MODERATE"
        message = "Long (>10kb)"
    
    return {"value": length, "status": status, "risk": risk, "message": message}

def check_repeats(sequence):
    """Check 1.4: Repetitive Region Scanner"""
    length = len(sequence)
    repeats = []
    window_size = 20
    seen_chunks = {}

    for i in range(length - window_size + 1):
        chunk = sequence[i : i + window_size]
        if chunk in seen_chunks:
            repeats.append({
                "start": seen_chunks[chunk] + 1,
                "end": seen_chunks[chunk] + window_size,
                "sequence": chunk
            })
            repeats.append({
                "start": i + 1,
                "end": i + window_size,
                "sequence": chunk
            })
        else:
            seen_chunks[chunk] = i

    status = "WARNING" if repeats else "PASS"
    risk = "MODERATE" if repeats else "LOW"

    return {"status": status, "risk": risk, "count": len(repeats), "details": repeats, "message": "Detected" if repeats else "None"}

def get_gc_plot_data(sequence, window_size=50):
    """Generates Line Graph Data"""
    if len(sequence) < window_size:
        return {"labels": [1], "values": [round(gc_fraction(sequence) * 100, 2)]}

    data_points = []
    labels = []
    step = max(1, len(sequence) // 200) # Increased resolution for smoother graph
    
    for i in range(0, len(sequence) - window_size, step):
        subseq = sequence[i : i + window_size]
        gc = gc_fraction(subseq) * 100
        data_points.append(round(gc, 2))
        labels.append(i + 1)
        
    return {"labels": labels, "values": data_points}

def get_gc_histogram_data(sequence, window_size=50):
    """Generates Histogram Data"""
    bins = [0] * 10 
    step = max(1, len(sequence) // 200)
    
    for i in range(0, len(sequence) - window_size, step):
        subseq = sequence[i : i + window_size]
        gc = int(gc_fraction(subseq) * 100)
        bin_idx = min(gc // 10, 9)
        bins[bin_idx] += 1
        
    return bins