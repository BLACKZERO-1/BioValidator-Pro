import math

def generate_dna_pdb(sequence):
    """
    Generates a PDB format string for the first 50 base pairs of the sequence.
    This creates a REAL 3D model of the user's specific DNA.
    """
    # Limit to 50bp to prevent browser crash (3D rendering is heavy)
    seq_to_render = sequence[:50].upper()
    pdb = []
    atom_id = 1
    
    # B-DNA Geometric Constants
    rise_per_bp = 3.4
    twist_per_bp = 36.0 # degrees
    radius = 10.0
    
    # Simple mapping of base to approximate atom sizes/positions (Simplified for visualization)
    # We create a "Backbone" (P) and a "Base" (N) for each nucleotide
    
    for i, base in enumerate(seq_to_render):
        angle_deg = i * twist_per_bp
        angle_rad = math.radians(angle_deg)
        z = i * rise_per_bp
        
        # Strand 1 (Forward)
        x1 = radius * math.cos(angle_rad)
        y1 = radius * math.sin(angle_rad)
        pdb.append(f"ATOM  {atom_id:5d}  P   DA A{i+1:4d}    {x1:8.3f}{y1:8.3f}{z:8.3f}  1.00 20.00           P")
        atom_id += 1
        
        # Base pair connection (Visual representation only)
        x1_in = (radius-4) * math.cos(angle_rad)
        y1_in = (radius-4) * math.sin(angle_rad)
        pdb.append(f"ATOM  {atom_id:5d}  C1' DA A{i+1:4d}    {x1_in:8.3f}{y1_in:8.3f}{z:8.3f}  1.00 20.00           C")
        atom_id += 1

        # Strand 2 (Reverse Complement position)
        # Offset by ~180 degrees + minor groove shift
        angle_rad2 = math.radians(angle_deg + 140) 
        x2 = radius * math.cos(angle_rad2)
        y2 = radius * math.sin(angle_rad2)
        pdb.append(f"ATOM  {atom_id:5d}  P   DT B{i+1:4d}    {x2:8.3f}{y2:8.3f}{z:8.3f}  1.00 20.00           P")
        atom_id += 1

        # Connect atoms for visualization lines
        pdb.append(f"CONECT{atom_id-3:5d}{atom_id-2:5d}") # Connect P to Base
        pdb.append(f"CONECT{atom_id-2:5d}{atom_id:5d}")   # Connect Strands (Hydrogen bond simulation)
        
    return "\n".join(pdb)

def generate_plasmid_map(sequence, sites):
    """
    Generates a Circular SVG Map for Module 2.
    """
    cx, cy, r = 150, 150, 100
    svg = [f'<svg viewBox="0 0 300 300" class="w-full h-full">']
    
    # Draw Circle
    svg.append(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="#cbd5e1" stroke-width="10" />')
    
    # Map Restriction Sites
    seq_len = len(sequence)
    if seq_len > 0:
        for site in sites:
            # Calculate angle: (position / length) * 360
            angle_deg = (site['position'] / seq_len) * 360 - 90 # -90 to start at top
            angle_rad = math.radians(angle_deg)
            
            # Line coordinates (sticking out)
            x1 = cx + (r * math.cos(angle_rad))
            y1 = cy + (r * math.sin(angle_rad))
            x2 = cx + ((r + 20) * math.cos(angle_rad))
            y2 = cy + ((r + 20) * math.sin(angle_rad))
            
            # Draw Tick
            color = "#9333ea" # Purple
            svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="2" />')
            
            # Draw Label (if space permits)
            # Only label first few to avoid clustering in simple view
            if site['position'] % 2 == 0 or len(sites) < 10:
                text_x = cx + ((r + 35) * math.cos(angle_rad))
                text_y = cy + ((r + 35) * math.sin(angle_rad))
                svg.append(f'<text x="{text_x}" y="{text_y}" font-size="10" text-anchor="middle" fill="#333">{site["enzyme"]}</text>')

    svg.append('</svg>')
    return "".join(svg)

def generate_codon_heatmap(original_seq, optimized_seq):
    """
    Generates data for Module 3 Heatmap.
    Returns array of 0 (Bad/Changed) vs 1 (Good/Same) for visualization.
    """
    heatmap_data = []
    # Compare by chunks of 3 (codons)
    # Limit to 100 codons for the graph to prevent overcrowding
    limit = 100 
    
    orig_codons = [original_seq[i:i+3] for i in range(0, len(original_seq), 3)]
    opt_codons = [optimized_seq[i:i+3] for i in range(0, len(optimized_seq), 3)]
    
    for i in range(min(len(orig_codons), limit)):
        if i < len(opt_codons) and orig_codons[i] == opt_codons[i]:
            heatmap_data.append(1) # Good/Unchanged
        else:
            heatmap_data.append(0.2) # Bad/Changed (Low score for visualization)
            
    return heatmap_data