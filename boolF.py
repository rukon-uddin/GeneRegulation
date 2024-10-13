import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def getBoolF(txtPath):
    gene_text_file = txtPath

    df = pd.read_csv(gene_text_file, sep=' ', header=None)
    df.columns = columns = [i for i in range(df.shape[1])]

    regulation = {"Negation": {i: [] for i in columns}, "Activation": {i: [] for i in columns}}
    network = {i: [] for i in columns}
    finalNetworkList = {}

    for c in range(len(columns)):
        tg = columns[c]  # target gene
        target_gene = df[tg]
        for v in range(0, len(columns)):
            if v == c:
                continue
            rg = columns[v]
            reg_gene = df[rg]
            positive_regulate = 0
            negetion_regulate = 0
            bk = 0
            for idx in range(1, len(target_gene)):
                ri = idx - 1
                tv = target_gene[idx]  # target gene value
                rv = reg_gene[ri]  # regulatory gene value
                if tv != rv and positive_regulate == 0:
                    positive_regulate = 1
                if tv == rv and negetion_regulate == 0:
                    negetion_regulate = 1

                if positive_regulate == 1 and negetion_regulate == 1:
                    bk = 1
                    break

            if bk == 0:
                if positive_regulate == 0:
                    regulation["Activation"][tg].append(rg)
                elif negetion_regulate == 0:
                    regulation["Negation"][tg].append(rg)

            rev_positive_regulate = 0
            rev_negetion_regulate = 0
            bk = 0
            for idx in range(1, len(target_gene)):
                ri = idx - 1
                tv = target_gene[ri]  # target gene value
                rv = reg_gene[idx]  # regulatory gene value
                if tv != rv and rev_positive_regulate == 0:
                    rev_positive_regulate = 1
                if tv == rv and rev_negetion_regulate == 0:
                    rev_negetion_regulate = 1

                if rev_positive_regulate == 1 and rev_negetion_regulate == 1:
                    bk = 1
                    break
                    
            if bk == 0:
                if positive_regulate == 0:
                    regulation["Activation"][rg].append(tg)
                elif negetion_regulate == 0:
                    regulation["Negation"][rg].append(tg)

    # print("Negation: ")
    for k, v in regulation["Negation"].items():
        if v != []:
            # print(k, v)
            for ii in v:
                finalNetworkList[k] = {"targetGenes": (ii,), "mismatch": 0}
            

    # print("Activation")
    for k, v in regulation["Activation"].items():
        if v != []:
            # print(k, v)
            for ii in v:
                finalNetworkList[k] = {"targetGenes": (ii,), "mismatch": 0}
    
    return finalNetworkList
