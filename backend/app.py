from flask import Flask, render_template, request, Response
from io import StringIO
from Bio import SeqIO

# IMPORTS
from modules.synthesis import check_gc_content, check_homopolymers, check_length, check_repeats, get_gc_plot_data, get_gc_histogram_data
from modules.restriction import find_restriction_sites
from modules.optimization import optimize_sequence
from modules.prediction import predict_problems
from modules.report import generate_report_text
# Visualizations
from modules.visualization import generate_dna_pdb, generate_plasmid_map, generate_codon_heatmap

app = Flask(__name__, 
            template_folder='../frontend/templates',
            static_folder='../frontend/static')

def run_analysis(sequence):
    """Helper function to run all checks + Generate Visuals"""
    
    # 1. Synthesis Check
    synthesis_res = {
        "length": check_length(sequence),
        "gc": check_gc_content(sequence),
        "gc_plot": get_gc_plot_data(sequence),
        "gc_hist": get_gc_histogram_data(sequence), # <--- NEW
        "homopolymers": check_homopolymers(sequence),
        "repeats": check_repeats(sequence)
    }
    
    # 2. Restriction Check
    restriction_res = find_restriction_sites(sequence)
    restriction_res['map_svg'] = generate_plasmid_map(sequence, restriction_res['sites'])
    
    # 3. Optimization Check
    opt_res = optimize_sequence(sequence)
    opt_res['heatmap'] = generate_codon_heatmap(sequence, opt_res['optimized_dna'])
    
    # 4. Prediction Check
    pred_res = predict_problems(sequence)
    # Note: We removed pdb_structure since we switched to 2D Safety Map
    
    return {
        "synthesis": synthesis_res,
        "restriction": restriction_res,
        "optimization": opt_res,
        "prediction": pred_res
    }

@app.route('/', methods=['GET', 'POST'])
def index():
    results = None
    sequence = ""
    
    if request.method == 'POST':
        uploaded_file = request.files.get('file')
        if uploaded_file and uploaded_file.filename != '':
            try:
                file_content = uploaded_file.read().decode("utf-8")
                fasta_io = StringIO(file_content)
                record = next(SeqIO.parse(fasta_io, "fasta"))
                sequence = str(record.seq).upper()
            except:
                pass
        else:
            raw_text = request.form.get('sequence', '')
            if ">" in raw_text:
                lines = raw_text.splitlines()
                sequence = "".join([line.strip() for line in lines if not line.startswith(">")]).upper()
            else:
                sequence = "".join(raw_text.split()).upper()

        if sequence:
            results = run_analysis(sequence)
            
    return render_template('index.html', results=results, sequence=sequence)

@app.route('/download', methods=['POST'])
def download():
    sequence = request.form.get('sequence', '')
    if not sequence: return "No sequence", 400
    results = run_analysis(sequence)
    report_content = generate_report_text(sequence, results)
    return Response(report_content, mimetype="text/plain", headers={"Content-disposition": "attachment; filename=biovalidator_report.txt"})

if __name__ == '__main__':
    app.run(debug=True, port=5000)