from playwright.sync_api import sync_playwright
URL = "https://theobroma.l3s.uni-hannover.de/compound/THEO_0858442"
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={"width": 1400, "height": 900})
    page.goto(URL)
    page.wait_for_load_state("networkidle")
    page.wait_for_timeout(3000)  # buffer for async stereoisomer widget and structure images
    page.pdf(
        path="fig3_fisetin_detail.pdf",
        format="A4",
        print_background=True,
        margin={"top": "5mm", "bottom": "5mm", "left": "5mm", "right": "5mm"},
        scale=0.65,
    )
    browser.close()
# split the 2-page output into p1 and p2 for direct LaTeX include
import fitz
doc = fitz.open("fig3_fisetin_detail.pdf")
for i, page in enumerate(doc, 1):
    sub = fitz.open(); sub.insert_pdf(doc, from_page=i-1, to_page=i-1)
    sub.save(f"fig3_fisetin_p{i}.pdf"); sub.close()
doc.close()
print("Saved: fig3_fisetin_detail.pdf, fig3_fisetin_p1.pdf, fig3_fisetin_p2.pdf")
