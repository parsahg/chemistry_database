from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"src"))
from kluyveromyces.chemistry.canonical_loader import build_canonical_dataset
RAW=ROOT/"data"/"raw"; OUT=ROOT/"data"/"processed"
df=build_canonical_dataset(RAW/"data_search1.xlsx",RAW/"data_search2.xlsx",RAW/"our_lab_data.xlsx")
OUT.mkdir(parents=True,exist_ok=True)
path=OUT/"canonical_experiments.csv"
df.to_csv(path,index=False,encoding="utf-8-sig")
print(f"Created: {path}")
print(f"Rows: {len(df)} | Columns: {len(df.columns)}")
