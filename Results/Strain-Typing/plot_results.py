import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC
import csv
from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix
import random
from copy import deepcopy


def train_and_plot(kernel, X_train, y_train, X_test, y_test, title):
    svm = SVC(kernel=kernel, C=1).fit(X_train, y_train)
    svm_predictions = svm.predict(X_test)

    # model accuracy for X_test
    accuracy = svm.score(X_test, y_test)
    print('Accuracy: ', accuracy)
    print('Confusion Matrix :')
    # creating a confusion matrix
    cm = confusion_matrix(y_test, svm_predictions)

    # Plot non-normalized confusion matrix
    class_names = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7', 'IC8', 'None']
    disp = plot_confusion_matrix(svm, X_test, y_test,
                                 display_labels=class_names,
                                 cmap=plt.cm.Blues,
                                 normalize='true',
                                values_format='.3f')
    # Source:
    # https://stackoverflow.com/questions/57043260/how-change-the-color-of-boxes-in-confusion-matrix-using-sklearn

    disp.ax_.set_title(title)
    plt.show()


def read_files(file):
    """Reads csv file and extracts necessary data for plotting"""
    r = csv.reader(open(file))
    lines = list(r)
    return lines

def get_ICs():
    """reads in csv with labels and genomes, assigns them to dictionary"""
    # Only taking genomes that have been assigned by MLST
    r = csv.reader(open("genomes_and_labels.csv"))
    lines = list(r)
    ICs = {
        "IC1": [],
        "IC2": [],
        "IC3": [],
        "IC4": [],
        "IC5": [],
        "IC6": [],
        "IC7": [],
        "IC8": [],
        "None": []
    }

    for i in range(len(lines)):
        labels = [lines[i][2], lines[i][4]]
        if labels[0] == labels[1]:
            # Oxford and Pasteur MLST assigned same label
            ICs[labels[0]].append(lines[i][0])
        else:
            # Using the protocol, that classified an IC
            # Could be either Oxford or Pasteur
            if labels[0] == 'None':
                del labels[0]
            ICs[labels[0]].append(lines[i][0])

    return ICs


def get_scores(file, lst):
    r = csv.reader(open(file))
    lines = list(r)

    for i in range(len(lines)):
        for x in range(len(lst)):
            for y in range(len(lst[x])):

                if lines[i][0] == lst[x][y]:
                    lst[x][y] = lines[i][1:]
                else:
                    pass


def create_plots(vecs = 4):
    ICs = get_ICs()
    training_assembly = [[], [], [], [], [], [], [], [], []]
    training_core = [[], [], [], [], [], [], [], [], []]
    label = [[], [], [], [], [], [], [], [], []]
    # taking 4 genomes randomly for each possible classification
    keys = ["IC1", "IC2", "IC3", "IC4", "IC5", "IC6", "IC7", "IC8", "None"]

    for i in range(len(keys)):
        sample = random.sample(ICs[keys[i]], vecs)
        training_assembly[i] = sample
        training_core[i] = deepcopy(sample)
        label[i] = [keys[i]] * vecs

    get_scores("Assemblys as reference, all k-mers.csv", training_assembly)
    get_scores("Core-Genome as reference, all k-mers.csv", training_core)

    # Plotting 4 Trainingvectors per label, Polynomial Kernel, Assemblys as reference, all k-mers
    title = '4 Trainingvectors per label, Polynomial Kernel, \n Assemblys as reference, all k-mers'
    train_and_plot('poly',)
    # Plotting 4 Trainingvectors per label, Polynomial Kernel, Assemblys as reference, every 10th k-mer

    # Plotting 4 Trainingvectors per label, Polynomial Kernel, Core-genome as reference, all k-mers

    # Plotting 4 Trainingvectors per label, Polynomial Kernel, Core-genome as reference, every 10th k-mer

    # Plotting 4 Trainingvectors per label, Linear Kernel, Assemblys as reference, all k-mers

    # Plotting 4 Trainingvectors per label, Linear Kernel, Assemblys as reference, every 10th k-mer

    # Plotting 4 Trainingvectors per label, Linear Kernel, Core-genome as reference, all k-mers

    # Plotting 4 Trainingvectors per label, Linear Kernel, Core-genome as reference, every 10th k-mer

    # Plotting 2 Trainingvectors per label, Polynomial Kernel, Assemblys as reference, all k-mers

    # Plotting 2 Trainingvectors per label, Polynomial Kernel, Assemblys as reference, every 10th k-mer

    # Plotting 2 Trainingvectors per label, Polynomial Kernel, Core-genome as reference, all k-mers

    # Plotting 2 Trainingvectors per label, Polynomial Kernel, Core-genome as reference, every 10th k-mer

    # Plotting 2 Trainingvectors per label, Linear Kernel, Assemblys as reference, all k-mers

    # Plotting 2 Trainingvectors per label, Linear Kernel, Assemblys as reference, every 10th k-mer

    # Plotting 2 Trainingvectors per label, Linear Kernel, Core-genome as reference, all k-mers

    # Plotting 2 Trainingvectors per label, Linear Kernel, Core-genome as reference, every 10th k-mer

    # Plotting Reads with 4 vectors per label, polynomial Kernel, Assemblys as reference

    # Plotting Reads with 4 vectors per label, polynomial Kernel, Core-Genome as reference





create_plots(4)

#print("Anleitung")
#Kernel, Anzahl tranining, Art Training, reads oder Assembly,Quick oder normal

#OXA extra



