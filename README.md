# MolFinder: LAMMPS Molecule Analyzer

MolFinder is a web application that analyzes LAMMPS data or dump files to identify molecular structures. It uses OVITO for particle analysis and RDKit for chemical informatics, providing a detailed, interactive report of the molecules found in your simulation.

All analysis results are stored in a database, and can be reviewed at any time on the "History" page.

## Features

-   **File Upload:** Supports LAMMPS data and dump files.
-   **Automatic Atom Typing:** Intelligently determines element types by:
    1.  Parsing comments in the `Masses` section of a LAMMPS data file.
    2.  Using a user-provided, comma-separated list of elements (e.g., `C,H,O`).
    3.  Reading the `element` column directly from a LAMMPS dump file.
-   **Interactive Results:** Displays results in a collapsible tree view (`Run > Stoichiometry > Molecule`).
-   **Molecule Details:** Shows the 2D image, SMILES code, and formula for each unique molecule in a pop-up modal.
-   **Persistent History:** All analysis runs are saved to a database for later review.

---

## Local Setup and Development

Follow these steps to run the application on your local machine.

### 1. Prerequisites

-   Python 3.10+
-   Git

### 2. Installation

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/dtramontina/molfinder
    cd molfinder
    ```

2.  **Create a Virtual Environment:**
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    ```

3.  **Install Dependencies:**
    The project's dependencies are listed in `requirements.txt`.
    ```bash
    pip install -r requirements.txt
    ```

### 3. Running the Application

1.  **Start the Flask Server:**
    ```bash
    flask run
    ```

2.  **Access the Application:**
    Open your web browser and navigate to `http://127.0.0.1:5000`. The first time you access the app, it will automatically create the `app.db` SQLite database file in the project directory.

---

## Deployment to the Web

This application is configured for easy, automated deployment on **Render**, a cloud platform with a generous free tier. It uses the `Dockerfile` and `render.yaml` files in this repository.

Every time you push a change to the `main` branch, Render will automatically rebuild and deploy the new version of your app.

### One-Time Setup

1.  **Create a Render Account:**
    *   Sign up for a free account at [render.com](https://render.com/) using your GitHub account.

2.  **Create a New Blueprint Service:**
    *   On the Render Dashboard, click **New +** and select **Blueprint**.
    *   Select your `molfinder-web-app` repository from the list.
    *   Render will automatically detect the `render.yaml` file. Give your service a name (e.g., `molfinder`) and click **Apply**.

3.  **Initial Deployment:**
    *   Render will begin the first deployment. It will provision a PostgreSQL database and build the Docker container for your web service. This may take a few minutes.
    *   Once the deployment is complete, you will have a public URL (e.g., `https://molfinder.onrender.com`) where you can access your live application.

From now on, your application will be automatically updated every time you `git push` to your `main` branch.
