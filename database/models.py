from .db import db
from datetime import datetime

class AnalysisRun(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(128), nullable=False)
    timestamp = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    lammps_file = db.Column(db.String(256), nullable=False)
    atom_mappings = db.Column(db.String(256), nullable=False)
    status = db.Column(db.String(64), default='Pending')
    
    stoichiometries = db.relationship('StoichiometryResult', backref='run', lazy=True, cascade="all, delete-orphan")

    def __repr__(self):
        return f'<AnalysisRun {self.name}>'

class StoichiometryResult(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    formula = db.Column(db.String(128), index=True, nullable=False)
    count = db.Column(db.Integer, nullable=False)
    run_id = db.Column(db.Integer, db.ForeignKey('analysis_run.id'), nullable=False)

    molecules = db.relationship('MoleculeInstance', backref='stoichiometry', lazy=True, cascade="all, delete-orphan")

    def __repr__(self):
        return f'<Stoichiometry {self.formula} (Count: {self.count})>'

class MoleculeInstance(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    smiles = db.Column(db.String(512), index=True, nullable=False)
    name = db.Column(db.String(256), nullable=True)
    image_path = db.Column(db.String(256), nullable=False)
    stoichiometry_id = db.Column(db.Integer, db.ForeignKey('stoichiometry_result.id'), nullable=False)

    def __repr__(self):
        return f'<Molecule {self.name} ({self.smiles})>'