"""Capture the /annotate page workflow end-state as a PDF for supplementary figure.
Uploads a four-compound test CSV, processes it, captures the resulting page render."""
from playwright.sync_api import sync_playwright
import tempfile, os
TEST_CSV = (
    "id,smiles\n"
    "curcumin,COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O\n"
    "fisetin,Oc1ccc2c(c1)oc(c(c2=O)O)c1ccc(c(c1)O)O\n"
    "caffeine,Cn1cnc2c1c(=O)n(C)c(=O)n2C\n"
    "novel_synthetic,C1CN(CCN1CCO)CC2=CC=CC=N2\n"
)
URL = "https://theobroma.l3s.uni-hannover.de/annotate"
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={"width": 1400, "height": 1600})
    page.goto(URL)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(2000)  # CDN scripts (papaparse, xlsx) settle
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(TEST_CSV); csv_path = f.name
    page.set_input_files("#file-input", csv_path)
    page.wait_for_selector("#preview-section", state="visible")
    page.wait_for_timeout(500)
    page.click("#process-btn")
    page.wait_for_selector("#results-section", state="visible", timeout=30000)
    page.wait_for_timeout(1500)  # progress bar and results render settle
    page.pdf(
        path="fig_annotate_workflow.pdf",
        format="A4",
        print_background=True,
        margin={"top": "5mm", "bottom": "5mm", "left": "5mm", "right": "5mm"},
        scale=0.7,
    )
    os.unlink(csv_path)
    browser.close()
print("Saved: fig_annotate_workflow.pdf")
