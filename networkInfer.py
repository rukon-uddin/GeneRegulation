import pickle
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tqdm import tqdm
import numpy as np
from boolF import getBoolF
import os
from tensorflow.keras.models import load_model
import pickle


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


def remove_elements_from_row(data, row, indexListToRemove):
    selected_row = data[row]
    result = [val for idx, val in enumerate(selected_row) if idx not in indexListToRemove]

    return result



def getInferRes(networkPath):
    size = 1
    ff = []
    with open(networkPath) as nt:
        mainData = nt.readlines()
    mainData = [i.strip(" \n").split(" ") for i in mainData]

    with open("counterG.pickle", "rb") as f:
        counterG = pickle.load(f)



    save_dir = 'M/myModels'

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

    print("********* BOOLF RES ***********")
    temp1 = getBoolF(networkPath)
    print(temp1)
    print("********* BOOLF RES ***********")
    print("********* NN INFERENCE ***********")
    timeStamp = len(mainData)
    numOfGenes = len(mainData[0]) 
    tgI = 0
    finalResult = getBoolF(networkPath)

    while tgI != numOfGenes:
        minMismatch = 1000000
        datasetDict = {}
        print("*"*50)
        print(f"Target gene Index: {tgI}")
        temp1 = getBoolF(networkPath)
        indexListToRemove = []
        for k, v in temp1.items():
            indexListToRemove.append(k)
        
        if tgI in indexListToRemove:
            tgI+=1
            continue
        else:
            indexListToRemove.append(tgI)
        
        print(indexListToRemove)
        
        yList = [i[tgI] for i in mainData]
        for t0 in range(timeStamp):
            if t0+1 == timeStamp or tgI == numOfGenes:
                break
            # x = mainData[t0][:tgI] + mainData[t0][tgI+1:]
            x = mainData[t0]
            # x = remove_elements_from_row(mainData, t0, indexListToRemove)
            alreadyTaken = {}

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
                    if len(keyList) == 0:
                        continue
                    # print(f"KEYS: {keyList}   VALS: {valList}")
                    if tuple(keyList) not in alreadyTaken:
                        alreadyTaken[tuple(keyList)] = 1
                        if tuple(keyList) not in datasetDict:
                            datasetDict[tuple(keyList)] = []

                        datasetDict[tuple(keyList)].append(valList)
            
        
        print(len(datasetDict))
        for k, v in tqdm(datasetDict.items(), total=len(datasetDict)):
            x = v
            predList = []
            probPred = []
            if len(x[0]) > 5:
                continue
            
            
            for p in x:
                predData = np.array([p])
                predictions = loaded_models_dict[tgI][len(p)].predict(predData, verbose=0)
                # threshold = 0.13
                # prediction = [1.0 if pred > threshold else 0.0 for pred in predictions]
                # predList.append(prediction[0])
                probPred.append(list(predictions)[0][0])

            max_value = max(probPred) - 0.1

            predList = [1 if value >= max_value else 0 for value in probPred]
            yy = []
            for p in x:
                counter0 = counterG[tgI][tuple(p)]["zeroCount"]
                counter1 = counterG[tgI][tuple(p)]["oneCount"]
                if counter0>=counter1:
                    yy.append(0.0)
                else:
                    yy.append(1.0)
            # yy = yList[1:]
            # yy = [float(b) for b in yy]
            mismatches = sum([1 for pred, actual in zip(predList, yy) if pred != actual])

            if mismatches < minMismatch:
                print(f"Previous minMismatch: {minMismatch} | Improved: {mismatches} | GeneCombination: {k}")
                minMismatch = mismatches
                finalResult[tgI] = {"targetGenes": k, "mismatch": mismatches}
                if mismatches == 0:
                    break

            # finalModels[tgI] = allModels
        tgI+=1
        print(finalResult)
    # break
    return finalResult

# if __name__ == "__main__":
#     getInferRes("NetworkTransition.txt0_1435.txt")