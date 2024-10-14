from flask import Flask, render_template, request, redirect, url_for
import os
import json
import networkTrain
import networkInfer
import pprint
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import pickle


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads/'  # Folder where uploaded files will be stored
finalD = {}


# def read_actual_data(file_path):
#     actual_data = {}
#     with open(file_path, 'r') as file:
#         for line in file:
#             parts = list(map(int, line.split()))
#             key = parts[0]
#             actual_data[key] = set(parts[1:])
#     return actual_data

def read_actual_data(file_path, value_dict):
    actual_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()  # Split the line into words
            if not parts:
                continue  # Skip empty lines
            key = value_dict[parts[0]]  # Map the first item (key) using the dictionary
            actual_data[key] = set(value_dict[part] for part in parts[1:] if part in value_dict)  # Map the rest of the line
    return actual_data

def clean_file(filepath):
    cleaned_lines = []
    with open(filepath, 'r') as file:
        for line in file:
            # Split by space, filter out empty strings and rejoin
            cleaned_line = ' '.join([x for x in line.strip().split(' ') if x])
            cleaned_lines.append(cleaned_line)

    # Save the cleaned version of the file
    cleaned_filepath = filepath.replace('.txt', '_cleaned.txt')
    with open(cleaned_filepath, 'w') as cleaned_file:
        cleaned_file.write('\n'.join(cleaned_lines))
    
    return cleaned_filepath

def calculate_accuracy_per_key(actual_data, predicted_data):
    accuracy_per_key = {}

    for key in actual_data:
        actual_genes = actual_data[key]
        predicted_genes = predicted_data.get(key, set())

        # Union of both actual and predicted genes to cover all possible cases
        all_genes = actual_genes.union(predicted_genes)

        if len(all_genes) == 0:  # To avoid division by zero
            accuracy_per_key[key] = 1.0 if len(actual_genes) == len(predicted_genes) else 0.0
        else:
            # Calculate the accuracy for this key
            correct_predictions = len(actual_genes.intersection(predicted_genes))
            accuracy = correct_predictions / len(all_genes)
            accuracy_per_key[key] = accuracy*100

    return accuracy_per_key



def calculate_metrics(actual_data, predicted_data):
    y_true = []
    y_pred = []

    for key in actual_data:
        actual_genes = actual_data[key]
        predicted_genes = predicted_data.get(key, set())

        for gene in actual_genes.union(predicted_genes):  
            y_true.append(1 if gene in actual_genes else 0)
            y_pred.append(1 if gene in predicted_genes else 0)

    # accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    
    accuracy_per_key = calculate_accuracy_per_key(actual_data, predicted_data)
    average_accuracy = sum(accuracy_per_key.values()) / len(accuracy_per_key)
    average_accuracy
    return average_accuracy, precision, recall, f1, accuracy_per_key

@app.route('/', methods=['GET', 'POST'])
def index():
    global finalD
    arr = [[], [], [], [], [], [], [], [], [], [], [], []]
    accuracy = 0.0
    recall = 0.0
    f1 = 0.0
    accuracy_per_key = {}
    precision = 0.0
    

    if request.method == 'POST':
        if 'file1' not in request.files or 'file2' not in request.files:
            return "No file part in the request"
        print(request.files)
        file = request.files['file1']
        file2 = request.files['file2']
        if file.filename == '' or file2.filename == '':
            return "No selected file"
        
        if file:
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)
            filepath2 = os.path.join(app.config['UPLOAD_FOLDER'], file2.filename)
            file2.save(filepath2)

            cleaned_filepath = clean_file(filepath)

            networkTrain.trainModel(cleaned_filepath)
            finalD = networkInfer.getInferRes(cleaned_filepath)
            cell_cycle_dict = {
                "Cln3": 0,
                "MBF": 1,
                "SBF": 2,
                "Cln1": 3,
                "Cdh1": 4,
                "Swi5": 5,
                "Cdc20": 6,
                "Clb5": 7,
                "Sic1": 8,
                "Clb1": 9,
                "Mcm1": 10
            }
            file_path2 = filepath2
            actual_data = read_actual_data(file_path2, cell_cycle_dict)
            predicted_data = {k: set(v['targetGenes']) for k, v in finalD.items()}
            accuracy, precision, recall, f1, accuracy_per_key = calculate_metrics(actual_data, predicted_data)
            
            pprint.pprint(finalD)
            pprint.pprint(actual_data)
            print("\n")
            pprint.pprint(accuracy_per_key)
            
            
            print(f"Accuracy: {accuracy:.2f}")
            print(f"Precision: {precision:.2f}")
            print(f"Recall: {recall:.2f}")
            print(f"F1-Score: {f1:.2f}")
            
            

    for k, v in finalD.items():
        for c, i in enumerate(list(v['targetGenes'])):
            arr[k + 1].append(i + 1)
            # arr[i + 1].append(k + 1)

    sets = [list(set(lst)) for lst in arr]
    sets = sets[1:]
    # return render_template('index.html', listt=json.dumps(sets))
    # with open("helllllloooooo.pkl", 'wb') as pickle_file:
    #     pickle.dump(sets, pickle_file)
    print(sets)
    return render_template('index.html', 
                           listt=json.dumps(sets), 
                           accuracy=accuracy, 
                           precision=precision*100, 
                           recall=recall*100, 
                           f1=f1*100, 
                           accuracy_per_key=json.dumps(accuracy_per_key))


if __name__ == '__main__':

    if not os.path.exists(app.config['UPLOAD_FOLDER']):
        os.makedirs(app.config['UPLOAD_FOLDER'])
    
    app.run(debug=True)
