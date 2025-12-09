import os

file_content = """def test_load_file_exists(tmp_path):
    from scripts.preprocessing import load_data
    import scanpy as sc
    ad = sc.AnnData(X=[[1,2],[3,4]])
    p = tmp_path / "test.h5ad"
    ad.write_h5ad(str(p))
    loaded = load_data(str(p))
    assert loaded.n_obs == 2
"""
output_dir = "tests"
os.makedirs(output_dir, exist_ok=True)

with open(os.path.join(output_dir, "test_loading.py"), "w") as f:
    f.write(file_content)
