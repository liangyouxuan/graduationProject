import os

import pandas as pd
from rdkit import Chem

from mybackend import settings

EXCEL_PATH = os.path.join(settings.BASE_DIR, "data_Chi.csv")

def load_excel_data():
    """读取 Excel 文件并提取去重的 SMILES 前后缀部分"""
    try:
        df = pd.read_csv(EXCEL_PATH)  # 读取整个 CSV 文件
        if "ps_pair" not in df.columns:  # 确保 "ps_pair" 这一列存在
            return set(), set()

        unique_prefixes = set()
        unique_suffixes = set()

        for smiles_string in df["ps_pair"].dropna():
            parts = smiles_string.strip().split("_")  # 去掉首尾空格并分割
            if len(parts) == 2:
                prefix, suffix = parts
                prefix = normalize_smiles(prefix)
                suffix = normalize_smiles(suffix)

                if prefix:
                    unique_prefixes.add(prefix)
                if suffix:
                    unique_suffixes.add(suffix)

        return unique_prefixes, unique_suffixes
    except Exception as e:
        print(f"读取 Excel 失败: {e}")
        return set(), set()

def normalize_smiles(smiles):
    """规范化 SMILES，移除显式氢并转换为标准格式（不保留手性信息）"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)  # 禁用手性信息
    return None

def check_smiles(smiles):
    """检查 SMILES 字符串的有效性，并返回 RDKit 分子对象"""
    if not isinstance(smiles, str) or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        return mol
    except Exception:
        return None
# 加载 SMILES 数据
prefixes, suffixes = load_excel_data()

# 打印前缀部分
print("📌 **前缀部分 (prefixes):**")
for prefix in prefixes:
    print(prefix)

print("\n📌 **后缀部分 (suffixes):**")
# 打印后缀部分
for suffix in suffixes:
    print(suffix)

