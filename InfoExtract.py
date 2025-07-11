import json

# Ruta al archivo JSON
json_file = '/home/sebastiandvargas/lignin-structure-generator-main/Executable_Jar_and_Config/output/json/LigninStructs_8.json'

# Leer el archivo y extraer los SMILES
with open(json_file, 'r') as f:
    data = json.load(f)

# Suponiendo que es una lista de diccionarios y cada uno tiene el campo 'smiles'
smiles_list = [entry['smiles'] for entry in data if 'smiles' in entry]

# Mostrar o guardar los SMILES
for smiles in smiles_list:
    print(smiles)

# Si quieres guardarlos en un archivo de texto:
with open('smiles_LigninStructs_8.txt', 'w') as f:
    for smiles in smiles_list:
        f.write(smiles + '\n')
