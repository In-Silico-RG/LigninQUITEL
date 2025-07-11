import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import re

# Función para contar átomos a partir de la fórmula molecular
def atom_count(formula, atom):
    match = re.search(rf'{atom}(\d*)', formula)
    if match:
        count = match.group(1)
        return int(count) if count else 1
    return 0

# Cargar la base de datos
df = pd.read_csv('lignin_molecules_database.csv')

# Listas para almacenar las razones H/C y O/C
hc_ratios = []
oc_ratios = []

# Calcular la fórmula molecular y las razones H/C y O/C
for idx, row in df.iterrows():
    smiles = row['smiles']
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            hc_ratios.append(None)
            oc_ratios.append(None)
            continue
        formula = CalcMolFormula(mol)
        c = atom_count(formula, 'C')
        h = atom_count(formula, 'H')
        o = atom_count(formula, 'O')
        if c > 0:
            hc_ratios.append(h / c)
            oc_ratios.append(o / c)
        else:
            hc_ratios.append(None)
            oc_ratios.append(None)
    except Exception as e:
        hc_ratios.append(None)
        oc_ratios.append(None)

# Agregar las columnas al DataFrame
df['H/C'] = hc_ratios
df['O/C'] = oc_ratios

# Filtrar filas válidas
df_valid = df.dropna(subset=['H/C', 'O/C'])

# Graficar y guardar el diagrama de Van Krevelen
plt.figure(figsize=(7, 6))
plt.scatter(df_valid['O/C'], df_valid['H/C'], alpha=0.6, s=15)
plt.xlabel('O/C')
plt.ylabel('H/C')
plt.title('Diagrama de Van Krevelen - Lignin Structures')
plt.grid(True)
plt.tight_layout()
plt.savefig('van_krevelen.png', dpi=300)
plt.show()

print("Diagrama de Van Krevelen guardado como van_krevelen.png")
