# LAMMPS Molecule Finder Web App

This is a Flask-based web application that analyzes LAMMPS data or dump files to identify molecular structures. It uses OVITO for the heavy-lifting of particle and bond analysis, and RDKit for chemical informatics.

## Features

-   Upload LAMMPS data or dump files.
-   Specify atom type mappings to chemical symbols.
-   Processes the file to identify molecular structures within clusters of atoms.
-   Displays a hierarchical report of found molecules:
    -   Grouped by stoichiometry (e.g., C₂H₄).
    -   Lists all isomers for each stoichiometry.
    -   Shows SMILES code, common name, and a 2D plot for each isomer.
-   Stores all analysis results in a database (SQLite by default).
-   Provides a history page to navigate and review past analyses.

## Prerequisites

This application relies on external scientific libraries which need to be installed in your environment.

1.  **Python 3.8+**

2.  **OVITO**: The OVITO Python package is required. Please follow the official installation instructions from ovito.org. It is essential that the Python environment running the app has access to the `ovito` package.

3.  **RDKit**: This is used for chemical informatics. The recommended way to install it is via Conda:
    ```bash
    conda install -c conda-forge rdkit
    ```

## Installation

1.  Clone this repository.
2.  Create a Python virtual environment and activate it.
3.  Install the required Python packages:
    ```bash
    pip install -r requirements.txt
    ```
4.  **Database Setup**: This application uses Flask-Migrate to manage database schemas.
    ```bash
    # Initialize the database (only needs to be run once)
    flask db init

    # Create the initial migration
    flask db migrate -m "Initial migration."

    # Apply the migration to the database
    flask db upgrade
    ```

## Usage

1.  Run the Flask application:
    ```bash
    flask run
    ```
2.  Open your web browser and navigate to `http://127.0.0.1:5000`.
3.  Fill out the form:
    -   **Analysis Name:** A descriptive name for your run.
    -   **Atom Type Mapping:** Provide a comma-separated list mapping LAMMPS atom type IDs to chemical symbols (e.g., `1:C, 2:H, 3:O`).
    -   **LAMMPS File:** Upload your LAMMPS data or dump file.
4.  Click "Start Analysis".
5.  You will be redirected to the results page for that run. You can view all past runs from the "History" link.

## Project Structure

molfinder-web-app/
├── uploads/
│   └── .gitkeep
├── static/
│   ├── css/
│   │   └── style.css
│   └── images/
│       └── results/
│           └── .gitkeep
├── templates/
│   ├── base.html
│   ├── index.html
│   ├── results.html
│   └── run_history.html
├── molfinder
│   ├── __init__.py
│   ├── molecule.py
│   └── processor.py
├── database/
│   ├── __init__.py
│   ├── models.py
│   └── db.py
├── migrations/
├── app.py
├── config.py
├── requirements.txt
└── README.md
