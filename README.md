Estamos creando un codigo para extraer informacion de ligninas generadas con LIgnin Structure Generator a partir de archivos .json

Aquí tienes un script (InfoExtract.py) en Python que lee todos los archivos .json en el directorio que especifiques (por ejemplo, la carpeta donde tienes los LigninStructs_*.json), extrae para cada molécula:

Un identificador único (por ejemplo: el lg_id o, si quieres un ID personalizado, lo podemos generar),
El grado de polimerización (dp o evaluatedDP),
El código SMILES (smilestring),
El peso molecular (molWeight).
Luego, guarda todo en una base de datos sencilla: un archivo .csv (que puedes abrir en Excel, pandas, etc.).

Notas:

El script busca los archivos .json en la carpeta indicada.
Puedes personalizar el ID de la molécula como prefieras.
El archivo resultante lignin_molecules_database.csv contendrá una fila por cada molécula con los datos solicitados.
Si algún campo no existe, se deja vacío o como 'NA'.
