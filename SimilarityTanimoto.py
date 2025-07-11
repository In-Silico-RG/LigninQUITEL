import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import numpy as np

# Cargar base de datos
df = pd.read_csv('lignin_molecules_database.csv')

# Asegúrate que la columna del grado de polimerización se llama 'dp'
# Si no, cambia 'dp' por el nombre correcto en tu CSV

# Agrupar por cada grado de polimerización
for dp_value, group in df.groupby('dp'):
    print(f'Calculando matriz de similaridad de Tanimoto para DP = {dp_value} ...')
    smiles_list = group['smiles'].tolist()
    ids = group['id'].tolist()

    # Generar huellas moleculares
    fps = []
    valid_ids = []
    for sm, id_ in zip(smiles_list, ids):
        mol = Chem.MolFromSmiles(sm)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(fp)
            valid_ids.append(id_)

    # Calcular matriz de similaridad
    n = len(fps)
    tanimoto_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            tanimoto_matrix[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])

    # Guardar matriz en un archivo CSV
    matrix_df = pd.DataFrame(tanimoto_matrix, index=valid_ids, columns=valid_ids)
    matrix_df.to_csv(f'tanimoto_matrix_dp_{dp_value}.csv')
    print(f'Matriz guardada en tanimoto_matrix_dp_{dp_value}.csv\n')
