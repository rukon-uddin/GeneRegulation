import boolF
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import numpy as np
import os
from itertools import combinations, permutations
from tqdm import tqdm
from tensorflow.keras.layers import Conv1D, Flatten
import shutil


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

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

def trainModel(networkPath):
    with open(networkPath) as nt:
        mainData = nt.readlines()
    mainData = [i.strip(" \n").split(" ") for i in mainData]

    timeStamp = len(mainData)
    numOfGenes = len(mainData[0])

    print(boolF.getBoolF(networkPath))
    finalResult = {}
    tgI = 0
    finalModels = {}
    while tgI != numOfGenes:
        temp1 = boolF.getBoolF(networkPath)
        indexListToRemove = []
        for k, v in temp1.items():
            indexListToRemove.append(k)

        if tgI in indexListToRemove:
            tgI+=1
            continue
        else:
            indexListToRemove.append(tgI)

        print("*"*100)
        print("Target Genet Index: ", tgI)
        allModels = {}
        rgI = tgI + 1 # regulatory gene
        yList = [i[tgI] for i in mainData]

        finalDataset = []
        counterG = {}
        for i in range(timeStamp):
            if i+1 == timeStamp:
                break
            Y = yList[i+1]
            Y = float(Y)
            # for j in range(rgI, numOfGenes):
                # x = mainData[i][rgI:j+1]
            x = mainData[i][:tgI] + mainData[i][tgI+1:]
            # x = remove_elements_from_row(mainData, i, indexListToRemove)
            x = [float(v) for v in x]
            for r in range(1, len(x) + 1):
                comb = list(set(combinations(x, r))) # taking only unique values
                for cm in comb:
                    if cm not in counterG:
                        counterG[cm] = {"zeroCount":0, "oneCount": 0}
                    if Y == 0.0:
                        counterG[cm]["zeroCount"] += 1
                    else:
                        counterG[cm]["oneCount"] += 1
                        
                    # finalDataset.append([list(cm), Y])
 
        # print(counterG)
        for k, v in counterG.items():
            count_0 = v["zeroCount"]
            count_1 = v["oneCount"]
            if count_0 >= count_1:
                finalDataset.append([list(k), 0.0])
            else:
                finalDataset.append([list(k), 1.0])
        
        print(len(finalDataset))
        for ng in range(1, numOfGenes):
            size = ng
            ff = []
            for i in finalDataset:
                if len(i[0]) == size:
                    ff.append(i)

            if len(ff) == 0:
                continue

            data = ff

            X = np.array([item[0] for item in data])
            y = np.array([item[1] for item in data])

            model = Sequential([
                Dense(size*2, input_shape=(size,), activation='relu'),
                # Dense(20, activation='relu'),
                Dense(1, activation='sigmoid')
            ])

            model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

            # Train the model
            # print("*"*100)
            print("Gene Combination Size: ", size)
            model.fit(X, y, epochs=50, batch_size=4, verbose=0)
            allModels[size] = model

        finalModels[tgI] = allModels
        tgI+=1


    model_dir = "M/myModels"

    if os.path.exists(model_dir):
        shutil.rmtree(model_dir)

    os.makedirs(model_dir)

    for key, model_dict in finalModels.items():
        for sub_key, model in model_dict.items():
            model.save(f"{model_dir}/model_{key}_{sub_key}.keras")
        

if __name__ == "__main__":
    trainModel("NetworkTransition.txt0_1435.txt")