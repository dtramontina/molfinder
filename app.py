import os
from flask import Flask, render_template, request, redirect, url_for, flash
from werkzeug.utils import secure_filename
from config import Config
from database.db import db
from database.models import AnalysisRun, StoichiometryResult, MoleculeInstance
from molfinder.processor import LammpsProcessor
from datetime import datetime

app = Flask(__name__)
app.config.from_object(Config)

db.init_app(app)

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

@app.route('/')
def index():
    run_id = request.args.get('run_id', None)
    run = None
    if run_id:
        run = db.session.get(AnalysisRun, run_id)
    return render_template('index.html', title="Home", run=run)

@app.route('/analyze', methods=['GET', 'POST'])
def analyze():
    if request.method == 'GET':
        return redirect(url_for('index'))

    # --- Form Data Retrieval ---
    if 'lammps_file' not in request.files:
        flash('No file part in the request.')
        return redirect(request.url)

    file = request.files['lammps_file']
    analysis_name = request.form.get('analysis_name', 'Untitled Analysis')
    atom_mapping_str = request.form.get('atom_types', '')

    if file.filename == '':
        flash('No file selected for uploading.')
        return redirect(url_for('index'))

    if file:
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], f"{datetime.utcnow().timestamp()}_{filename}")
        file.save(filepath)

        # --- Database Run Creation ---
        run = AnalysisRun(
            name=analysis_name,
            lammps_file=filepath,
            atom_mappings=atom_mapping_str,
            status='Processing'
        )
        db.session.add(run)
        db.session.commit()

        # --- Processing ---
        try:
            processor = LammpsProcessor(filepath, atom_mapping_str)
            stoichiometry_groups = processor.process()

            # --- Store Results in DB ---
            image_dir = os.path.join(app.static_folder, 'images', 'results', str(run.id))
            os.makedirs(image_dir, exist_ok=True)

            for stoich_formula, mol_group in stoichiometry_groups.items():
                stoich_result = StoichiometryResult(
                    formula=stoich_formula,
                    count=mol_group['count'],
                    run_id=run.id
                )
                db.session.add(stoich_result)
                db.session.flush() # To get the stoich_result.id

                for i, molecule in enumerate(mol_group['molecules']):
                    img_filename = f"{stoich_formula.replace('/','-')}_{i}.png"
                    img_path = os.path.join(image_dir, img_filename)
                    molecule.draw_to_file(img_path)

                    mol_instance = MoleculeInstance(
                        smiles=molecule.smiles,
                        name=molecule.name,
                        image_path=os.path.join('images', 'results', str(run.id), img_filename),
                        stoichiometry_id=stoich_result.id
                    )
                    db.session.add(mol_instance)

            run.status = 'Completed'
            db.session.commit()
            return redirect(url_for('index', run_id=run.id))

        except Exception as e:
            run.status = 'Failed'
            db.session.commit()
            flash(f"An error occurred during analysis: {e}")
            return redirect(url_for('index', run_id=run.id))

@app.route('/results/<int:run_id>')
def results(run_id):
    run = db.session.get(AnalysisRun, run_id)
    if not run:
        flash('Analysis run not found.')
        return redirect(url_for('index'))
    return render_template('results.html', run=run, title=f"Results for {run.name}")

@app.route('/history')
def history():
    page = request.args.get('page', 1, type=int)
    runs_pagination = db.paginate(
        db.select(AnalysisRun).order_by(AnalysisRun.timestamp.desc()),
        page=page,
        per_page=app.config['RUNS_PER_PAGE']
    )
    return render_template('run_history.html', runs_pagination=runs_pagination, title="Analysis History")

if __name__ == '__main__':
    app.run(debug=True)
