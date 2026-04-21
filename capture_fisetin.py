from playwright.sync_api import sync_playwright

URL = "https://theobroma.l3s.uni-hannover.de/compound/THEO_0858442"

with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={"width": 1400, "height": 900})
    page.goto(URL)
    page.wait_for_load_state("networkidle")
    page.screenshot(path="fig3a_detail_page.png", full_page=True)
    browser.close()
print("Saved: fig3a_detail_page.png")
