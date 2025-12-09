import os

file_content = '''"""utils.py
Shared utility functions
"""
import os

def safe_write(adata, path):
    try:
        adata.write_h5ad(path)
    except Exception as e:
        raise IOError(f'Failed to write AnnData to {path}: {e}')
'''
output_dir = "scripts"
os.makedirs(output_dir, exist_ok=True)

with open(os.path.join(output_dir, "utils.py"), "w") as f:
    f.write(file_content)
