import uproot
import json
import sys

def root_to_json(root_file_path, json_file_name):
    # Öffnet die ROOT-Datei
    with uproot.open(root_file_path) as file:
        # data_dict = {}

        # Iteriere über alle Objekte im ROOT-File
        for key in file.keys():
            obj = file[key]
            print(f"Lese Objekt: {key} ({type(obj)})")

            # Für TTrees: lade alle Arrays
            if isinstance(obj, uproot.behaviors.TTree.TTree):
                data = obj.arrays(library="np")  # oder "ak" für awkward arrays
                # numpy arrays sind nicht direkt JSON-serialisierbar, also in Listen umwandeln
                # data_serializable = {k: v.tolist() for k, v in data.items()}
                data_dict = {k: v.tolist() for k, v in data.items()}
                # data_dict[key] = data_serializable
            else:
                # Nicht unterstützter Typ (z.B. TH1, TGraph), hier ggf. erweitern
                print(f"Überspringe nicht unterstütztes Objekt: {key}")
            json_file_path = json_file_name + '_' + key.replace(';','_') + '.json'
            # Schreibe die Daten als JSON
            with open(json_file_path, "w") as json_file:
                json.dump(data_dict, json_file, indent=1)

            print(f"Datei erfolgreich geschrieben: {json_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Verwendung: python root_to_json.py input.root output.json")
    else:
        root_to_json(sys.argv[1], sys.argv[2])
