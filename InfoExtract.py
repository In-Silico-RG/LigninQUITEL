import os
import json
import csv

# Ruta del directorio donde están los archivos JSON
json_dir = '/home/sebastiandvargas/lignin-structure-generator-main/Executable_Jar_and_Config/output/json'

# Archivo de salida CSV
output_csv = 'lignin_molecules_database.csv'

# Encabezados para el CSV
header = ['id', 'dp', 'smiles', 'molweight']

# Lista para almacenar la información de todas las moléculas
molecules = []

# Recorrer todos los archivos JSON en el directorio
for filename in os.listdir(json_dir):
    if filename.endswith('.json'):
        file_path = os.path.join(json_dir, filename)
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        # El grado de polimerización puede estar en "dp" o "evaluatedDP" (lo que esté disponible)
        dp = data.get('dp', data.get('evaluatedDP', 'NA'))
        
        # Recorrer todas las moléculas en ligninchains
        for entry in data['ligninchains']:
            mol_id = entry.get('lg_id', '')
            smiles = entry.get('smilestring', '')
            molweight = entry.get('molWeight', entry.get('molweight', ''))
            
            # Generar un ID personalizado si quieres: por ejemplo, archivo + lg_id
            unique_id = f"{filename.replace('.json','')}_{mol_id}"
            
            molecules.append([unique_id, dp, smiles, molweight])

# Escribir la base de datos en un archivo CSV
with open(output_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(molecules)

print(f"¡Base de datos creada en {output_csv} con {len(molecules)} moléculas!")
