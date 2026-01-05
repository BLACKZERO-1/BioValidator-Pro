from Bio.Seq import Seq

def predict_problems(sequence):
    """
    Module 4: Safety & Stability Scanner (Pro Console Edition)
    """
    issues = []
    checks = {
        "stops": {"status": "PASS", "count": 0, "label": "Premature Stops", "risk_val": 5},
        "hairpins": {"status": "PASS", "count": 0, "label": "Strong Hairpins", "risk_val": 2},
        "chi": {"status": "PASS", "count": 0, "label": "Chi Sites", "risk_val": 3}
    }
    
    # Helper
    def find_all(sub, seq):
        start = 0
        while True:
            start = seq.find(sub, start)
            if start == -1: return
            yield start
            start += 1

    # 1. Chi Sites
    for pos in find_all("GCTGGTGG", sequence):
        issues.append({ "type": "Chi Site", "risk": "High", "start": pos+1, "end": pos+8, "color": "orange", "val": 3 })
        checks["chi"]["count"] += 1

    # 2. Hairpins
    window = 10
    for i in range(0, len(sequence) - window - 10, 50): 
        chunk = sequence[i:i+window]
        rev_comp = str(Seq(chunk).reverse_complement())
        search_area = sequence[i+window+5 : i+window+40]
        if rev_comp in search_area:
            issues.append({ "type": "Hairpin", "risk": "Medium", "start": i+1, "end": i+40, "color": "blue", "val": 2 })
            checks["hairpins"]["count"] += 1
            if checks["hairpins"]["count"] >= 5: break 

    # 3. Stops
    try:
        protein = Seq(sequence).translate(to_stop=False)
        for i, aa in enumerate(protein[:-1]): 
            if aa == "*":
                pos = i * 3
                issues.append({ "type": "Premature Stop", "risk": "Critical", "start": pos+1, "end": pos+3, "color": "red", "val": 5 })
                checks["stops"]["count"] += 1
    except: pass

    # --- RISK DENSITY CALCULATION (THIS IS THE MISSING PART CAUSING THE ERROR) ---
    bin_count = 20
    bin_size = max(1, len(sequence) // bin_count)
    risk_density = [0] * bin_count
    
    for issue in issues:
        bin_idx = min(bin_count - 1, (issue['start'] // bin_size))
        # Ensure issue has 'val', default to 1 if not found
        risk_val = issue.get('val', 1) 
        risk_density[bin_idx] += risk_val

    # Scoring
    score = 100
    score -= (checks["stops"]["count"] * 50)
    score -= (checks["chi"]["count"] * 20)
    score -= (checks["hairpins"]["count"] * 5)
    score = max(0, score)
    
    for key in checks:
        if checks[key]["count"] > 0: checks[key]["status"] = "FAIL"

    issues.sort(key=lambda x: x['start'])

    return {
        "score": score,
        "status": "SECURE" if score == 100 else ("UNSAFE" if score < 50 else "STABLE"),
        "issues": issues,
        "checks": checks,
        "risk_density": risk_density, 
        "bin_labels": [f"{i*bin_size}-{(i+1)*bin_size}" for i in range(bin_count)],
        "message": f"{len(issues)} anomalies detected."
    }