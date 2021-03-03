import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC
import csv
from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix
import random

def train_and_plot(kernel, X_train, y_train, X_test, y_test, title):
    svm = SVC(kernel=kernel, C=1).fit(X_train, y_train)
    svm_predictions = svm.predict(X_test)

    # model accuracy for X_test
    accuracy = str(round(svm.score(X_test, y_test), 4))
    title = title + accuracy
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


def get_scores(file, dic, keys):
    r = csv.reader(open(file))
    lines = list(r)

    for key in keys:
        for i in lines:
            if i[0] in dic[key]:
                pos = dic[key].index(i[0])
                # Converting data to float
                dic[key][pos] = [float(x) for x in i[1:]]
    return dic


def get_random_training_vectors(dic, keys, n):
    # Chooses n random vectors per class, returns them
    # as np.arrays and separates training/test-data
    X_train = []
    y_train = []
    X_test = []
    y_test = []

    # Choosing random training vectors and removing them from test-data
    for key in keys:
        pos = random.sample(range(0, len(dic[key])), n)
        pos.sort()
        for i in range(len(pos)-1, -1, -1):
            X_train.append(dic[key][pos[i]])
            del dic[key][pos[i]]
            y_train.append(key)
        # Creating test-data
        X_test = X_test + dic[key]
        y_test = y_test + [key] * len(dic[key])

    return np.array(X_train), np.array(y_train), np.array(X_test), np.array(y_test)


def get_read_data(file, dic, keys):
    """ Assigns reads to labels"""
    r = csv.reader(open(file))
    lines = list(r)
    vecs = []
    labels = []
    for key in keys:
        for i in dic[key]:
            for j in lines:
                if i in j[0]:
                    vecs.append(j[2:])
                    labels.append(key)

    return np.array(vecs), np.array(labels)


def create_plots():
    ICs_assembly = get_ICs()
    ICs_core = get_ICs()

    # Changing NCBI accession to Score-Vector
    keys = ["IC1", "IC2", "IC3", "IC4", "IC5", "IC6", "IC7", "IC8", "None"]

    # Full Run: all k-mers are tested
    ICs_assembly_full = get_scores("Assemblys as reference, all k-mers.csv", ICs_assembly, keys)
    ICs_core_full = get_scores("Core-Genome as reference, all k-mers.csv", ICs_core, keys)
    # Quick-run: only every 10th k-mer is used
    ICs_assembly_quick = get_scores("Assemblys as reference, quick run.csv", ICs_assembly, keys)
    ICs_core_quick = get_scores("Core-Genome as reference, quick run.csv", ICs_core, keys)

    # Note on Quick-run: The Genomes used for training are from the full run, the same genomes are not removed from
    # the testing data set in the quick-run, because only a subset of k-mers is evaluated
    quick_assembly = []
    quick_assembly_label = []
    quick_core = []
    quick_core_label = []
    for key in keys:
        # Assembly Reference
        quick_assembly = quick_assembly + ICs_assembly_quick[key]
        quick_assembly_label = quick_assembly_label + [key] * len(ICs_assembly_quick[key])
        # Core Reference
        quick_core = quick_core + ICs_core_quick[key]
        quick_core_label = quick_core_label + [key] * len(ICs_core_quick[key])

    # Assemblys

    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_assembly_full, keys, 4)

    title = 'SVM-Kernel: Polynomial, \n 4 Trainingvectors per class, \n Assembled Genomes as reference , \n all k-mers are tested, \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_test, y_test, title)

    title = 'SVM-Kernel: Linear, \n 4 Trainingvectors per class, \n Assembled Genomes as reference , \n all k-mers are tested, \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, X_test, y_test, title)

    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_assembly_full, keys, 2)
    title = 'SVM-Kernel: Polynomial, \n 2 Trainingvectors per class, \n Assembled Genomes as reference , \n all k-mers are tested, \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_test, y_test, title)

    title = 'SVM-Kernel: Linear, \n 2 Trainingvectors per class, \n Assembled Genomes as reference , \n all k-mers are tested, \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, X_test, y_test, title)
    
    title = 'SVM-Kernel: Polynomial, \n 4 Trainingvectors per class, \n Assembled Genomes as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, quick_assembly, quick_assembly_label, title)

    title = 'SVM-Kernel: Linear, \n 4 Trainingvectors per class, \n Assembled Genomes as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, quick_assembly, quick_assembly_label, title)

    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_assembly_full, keys, 2)
    title = 'SVM-Kernel: Polynomial, \n 2 Trainingvectors per class, \n Assembled Genomes as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, quick_assembly, quick_assembly_label, title)

    title = 'SVM-Kernel: Linear, \n 2 Trainingvectors per class, \n Assembled Genomes as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, quick_assembly, quick_assembly_label, title)
    """
    # Core Genome

    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_core_full, keys, 4)
    """
    title = 'SVM-Kernel: Polynomial, \n 4 Trainingvectors per class, \n Core-Genome as reference , \n all k-mers are tested , \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_test, y_test, title)

    title = 'SVM-Kernel: Linear, \n 4 Trainingvectors per class, \n Core-Genome as reference , \n all k-mers are tested , \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, X_test, y_test, title)

    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_core_full, keys, 2)
    title = 'SVM-Kernel: Polynomial, \n 2 Trainingvectors per class, \n Core-Genome as reference , \n all k-mers are tested , \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_test, y_test, title)

    title = 'SVM-Kernel: Linear, \n 2 Trainingvectors per class, \n Core-Genome as reference , \n all k-mers are tested , \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, X_test, y_test, title)
    
    title = 'SVM-Kernel: Polynomial, \n 4 Trainingvectors per class, \n Core-Genome as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, quick_core, quick_core_label, title)

    title = 'SVM-Kernel: Linear, \n 4 Trainingvectors per class, \n Core-Genome as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, quick_core, quick_core_label, title)

    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_core_full, keys, 2)
    title = 'SVM-Kernel: Polynomial, \n 2 Trainingvectors per class, \n Core-Genome as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, quick_core, quick_core_label, title)

    title = 'SVM-Kernel: Linear, \n 2 Trainingvectors per class, \n Core-Genome as reference , \n every 10th k-mer tested, \n Total Accuracy: '
    train_and_plot('linear', X_train, y_train, quick_core, quick_core_label, title)

    # Reads
    dic = get_ICs()
    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_assembly_full, keys, 4)
    # Assembly as Reference, unmodified reads
    X_reads, y_reads = get_read_data('Reads as Input, Assembly as Reference.csv', dic, keys)
    title = 'Sequence-Reads as Input, \n SVM-Kernel: Poly, \n 4 Trainingvectors per class, \n Assembled Genomes as reference,\n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_reads, y_reads, title)

    # Assembly as Reference, modified reads
    X_reads, y_reads = get_read_data('Reads as Input, Assembly as Reference, 15% of the nucleotides changed.csv', dic, keys)
    title = 'Reads with 15% changed Nucleotides as Input, \n SVM-Kernel: Poly, \n 4 Trainingvectors per class, \n Core-Genome as reference,\n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_reads, y_reads, title)

    # Core as Reference, unmodified Reads
    X_train, y_train, X_test, y_test = get_random_training_vectors(ICs_core_full, keys, 4)
    X_reads, y_reads = get_read_data('Reads as Input, Core-Genome as Reference.csv', dic, keys)
    title = 'Sequence-Reads as Input, \n SVM-Kernel: Poly, \n 4 Trainingvectors per class, \n Core-Genome as reference,\n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_reads, y_reads, title)

    # Core as Reference, modified reads
    X_reads, y_reads = get_read_data('Reads as Input, Core-Genome as Reference, 15% of the nucleotides changed.csv', dic, keys)
    title = 'Reads with 15% changed Nucleotides as Input, \n SVM-Kernel: Poly, \n 4 Trainingvectors per class, \n Core-Genome as reference,\n Total Accuracy: '
    train_and_plot('poly', X_train, y_train, X_reads, y_reads, title)

create_plots()





