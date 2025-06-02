import os
import zipfile
from datetime import datetime

curr_y = str(datetime.now().year)
curr_m = str(datetime.now().month)
curr_d = str(datetime.now().day)
curr_h = str(datetime.now().hour)
curr_s = str(datetime.now().minute)

# Configuration
ZIP_NAME = f"IFRS17_DATA_QUALITY_TOOL_{curr_y}_{curr_m}_{curr_d}_{curr_h}_{curr_s}.zip"
EXCLUDE = {"zip_app.bat", "zip_app.py", "1_input", "2_results", "3_download", ".idea", "__pycache__", "venv", "IFRS17_DATA_QUALITY_TOOL.zip"}

def should_exclude(path):
    top_level = os.path.normpath(path).split(os.sep)[0]
    return top_level in EXCLUDE

def zip_folder(base_path, zip_file):
    for foldername, subfolders, filenames in os.walk(base_path):
        for name in filenames:
            full_path = os.path.join(foldername, name)
            rel_path = os.path.relpath(full_path, base_path)
            if not should_exclude(rel_path):
                zip_file.write(full_path, rel_path)

if os.path.exists(ZIP_NAME):
    os.remove(ZIP_NAME)

with zipfile.ZipFile(ZIP_NAME, 'w', zipfile.ZIP_DEFLATED) as zipf:
    zip_folder('.', zipf)

print(f"Created: {ZIP_NAME}")