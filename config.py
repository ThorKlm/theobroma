import os

class Config:
    DB_URI = os.environ.get("DATABASE_URL",
        "postgresql://theobroma:theobroma@localhost:5432/theobroma")
    DATA_DIR = os.environ.get("DATA_DIR", "/home/thorben.klamt/theobroma/data")
    PER_PAGE = 50

# ============================================================================
# Version
# ============================================================================
VERSION_INTERNAL = "1.34"
VERSION_EXTERNAL = "1.0"
VERSION_DISPLAY = f"v{VERSION_INTERNAL}"
