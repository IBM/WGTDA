import os
import sys

# Absolute path to the repo root (WGTDA)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Path to src/
SRC = os.path.join(ROOT, "src")

# Inject into sys.path
sys.path.insert(0, SRC)
sys.path.insert(0, ROOT)

print(">>> PYTHONPATH set for tests:")
print(sys.path[:5])
