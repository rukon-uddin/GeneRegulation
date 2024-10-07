import pickle
import tensorflow as tf
from tensorflow.keras.models import Sequential, load_model
from tqdm import tqdm
import numpy as np
import os
import boolF

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

txtPath = "NetworkTransition.txt0_90.txt"

with open(txtPath) as nt:
    mainData = nt.readlines()
mainData = [i.strip(" \n").split(" ") for i in mainData]

with open('models_dict.pkl', 'rb') as file:
    loaded_models_dict = pickle.load(file)
    


# Directory where models are saved
# save_dir = 'models3'
save_dir = 'cnnModel1'

# Initialize the dictionary to store models
loaded_models_dict = {}

# Loop through all files in the directory
for file_name in os.listdir(save_dir):
    if file_name.endswith('.keras'):
        # Extract key and sub_key from the file name
        # Assuming file names are in the format "model_key_sub_key.keras"
        parts = file_name.split('_')
        key = int(parts[1])  # Convert to integer, assuming keys are integers
        sub_key = int(parts[2].split('.')[0])  # Extract sub_key and remove '.keras'

        # Load the model
        model_path = os.path.join(save_dir, file_name)
        model = load_model(model_path)

        # Store the model in the dictionary
        if key not in loaded_models_dict:
            loaded_models_dict[key] = {}
        loaded_models_dict[key][sub_key] = model


def generate_combinations_with_indices(arr, r):
    n = len(arr)
    indices = list(range(r))
    result = []

    while True:
        # Create a dictionary with index-value pairs
        combination_dict = {indices[i]: arr[indices[i]] for i in range(r)}
        result.append(combination_dict)

        # Find the rightmost index that can be incremented
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return result  # No more combinations can be generated

        # Increment this index
        indices[i] += 1

        # Update the indices that follow
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1

def generate_permutations_with_indices(arr, r):
    from itertools import permutations
    n = len(arr)
    result = []

    # Generate all possible permutations using indices
    for perm in permutations(range(n), r):
        # Create a dictionary with index-value pairs
        permutation_dict = {i: arr[i] for i in perm}
        result.append(permutation_dict)

    return result


def getnetwork():
    print("********* INFERENCE ***********")
    timeStamp = 20
    numOfGenes = 10
    tgI = 0
    finalResult = boolF.getBoolF(txtPath)
    
    
    while tgI != numOfGenes:
        if tgI in finalResult:
            print(f"Target Gene Index: {tgI} generated from boolF" )
            tgI+=1
            continue
        print("*"*100)
        print("Target Gene Index: ", tgI)
        yList = [i[tgI] for i in mainData]
        minMismatch = 1000000
        datasetDict = {}
        for t0 in range(timeStamp):
            if t0+1 == timeStamp or tgI == numOfGenes:
                break
            x = mainData[t0][:tgI] + mainData[t0][tgI+1:]
            x = [float(v) for v in x]
            for r in range(1, len(x) + 1):
                comb = generate_combinations_with_indices(x, r)
                for cc in comb:
                    keyList = []
                    valList = []
                    # print(i)
                    for k, v in cc.items():
                        if k >= tgI:
                            k = k+1
                        keyList.append(k)
                        valList.append(v)

                    if tuple(keyList) not in datasetDict:
                        datasetDict[tuple(keyList)] = []

                    datasetDict[tuple(keyList)].append(valList)

        print(len(datasetDict))
        for k, v in tqdm(datasetDict.items(), total=len(datasetDict)):
            x = v
            predList = []
            probPred = []
            for p in x:
                predData = np.array([p])
                predictions = loaded_models_dict[tgI][len(p)].predict(predData, verbose=0)
                # threshold = 0.13
                # prediction = [1.0 if pred > threshold else 0.0 for pred in predictions]
                # predList.append(prediction[0])
                probPred.append(list(predictions)[0][0])

            max_value = max(probPred) - 0.1

            predList = [1 if value >= max_value else 0 for value in probPred]
            yy = yList[1:]
            yy = [float(b) for b in yy]
            mismatches = sum([1 for pred, actual in zip(predList, yy) if pred != actual])

            if mismatches < minMismatch:
                minMismatch = mismatches
                finalResult[tgI] = {"targetGenes": k, "mismatch": mismatches}
                if mismatches == 0 or mismatches <= 1:
                    break

            # finalModels[tgI] = allModels
        tgI+=1
        print(finalResult)
    
    return finalResult
        # break