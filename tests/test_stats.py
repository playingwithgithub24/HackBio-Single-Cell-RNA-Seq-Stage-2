import os

file_content = """def test_credible_interval():
    from scripts.stats import credible_interval
    low, high = credible_interval(5, 100)
    assert 0 <= low < high <= 1
"""
output_dir = "tests"
os.makedirs(output_dir, exist_ok=True)

with open(os.path.join(output_dir, "test_stats.py"), "w") as f:
    f.write(file_content)
