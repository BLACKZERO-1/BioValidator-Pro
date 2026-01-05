import datetime

def generate_report_text(sequence, results):
    """
    Generates a formatted text file content for the analysis.
    """
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    report = []
    report.append("==================================================")
    report.append("              BIOVALIDATOR LAB REPORT             ")
    report.append("==================================================")
    report.append(f"Date: {timestamp}")
    report.append(f"Sequence Length: {len(sequence)} bp")
    report.append("-" * 50)
    report.append("")
    
    # --- MODULE 1: SYNTHESIS ---
    report.append("[MODULE 1: SYNTHESIS FEASIBILITY]")
    s = results['synthesis']
    report.append(f"Status: {s['length']['status']} | GC Content: {s['gc']['value']}% ({s['gc']['status']})")
    if s['homopolymers']['count'] > 0:
        report.append(f"WARNING: {s['homopolymers']['count']} Homopolymer runs detected.")
    else:
        report.append("Homopolymers: Clean")
    report.append("-" * 50)
    report.append("")

    # --- MODULE 2: RESTRICTION ---
    report.append("[MODULE 2: RESTRICTION MAP]")
    r = results['restriction']
    report.append(f"Status: {r['message']}")
    if r['count'] > 0:
        report.append(f"{'Enzyme':<15} | {'Position':<10} | {'Pattern'}")
        report.append("-" * 45)
        for site in r['sites']:
            report.append(f"{site['enzyme']:<15} | {site['position']:<10} | {site['pattern']}")
    report.append("-" * 50)
    report.append("")

    # --- MODULE 3: OPTIMIZATION ---
    report.append("[MODULE 3: CODON OPTIMIZATION]")
    o = results['optimization']
    report.append(f"Original CAI Score: Poor (Estimated)")
    report.append(f"Optimized Score: {o['percent_improved']}% Improved")
    report.append(f"Codons Changed: {o['changes']}")
    report.append("\nOPTIMIZED SEQUENCE (E. COLI K12):")
    report.append(o['optimized_dna'])
    report.append("-" * 50)
    report.append("")

    # --- MODULE 4: SAFETY ---
    report.append("[MODULE 4: BIOLOGICAL SAFETY]")
    p = results['prediction']
    if p['count'] == 0:
        report.append("Status: SAFE (No issues detected)")
    else:
        report.append(f"Status: CRITICAL ({p['count']} Issues)")
        for issue in p['issues']:
            report.append(f"CRITICAL: {issue['type']} - {issue['message']}")
            
    report.append("==================================================")
    report.append("             END OF REPORT - BIOVALIDATOR         ")
    report.append("==================================================")

    return "\n".join(report)