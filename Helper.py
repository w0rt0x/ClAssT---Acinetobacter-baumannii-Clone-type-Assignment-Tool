import BF_v2
import pickle
from bitarray import bitarray
import os

def write_file():
    # https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python/899176
    itemlist = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5','IC6','IC7','IC8']
    with open(r'filter/FilterClonetypes.txt', 'wb') as fp:
        pickle.dump(itemlist, fp)

def write_file2():
    # https://stackoverflow.com/questions/899103/writing-a-list-to-a-file-with-python/899176
    itemlist = ['bla OXA-23', 'bla OXA-24', 'bla OXA-51', 'bla OXA-58', 'bla OXA-134',
                'bla OXA-143','bla OXA-211','bla OXA-213','bla OXA-214','bla OXA-229','bla OXA-286']
    with open(r'filter/OXAs/FilterOXA.txt', 'wb') as fp:
        pickle.dump(itemlist, fp)


def train():
    files = os.listdir(r'I:\OXA-Gene')
    for i in range(len(files) -1, -1, -1):
        if 'fasta' not in files[i]:
            del files[i]

    for i in range(len(files)):

        BF = BF_v2.AbaumanniiBloomfilter(80000)
        BF.set_arraysize(80000)
        BF.set_clonetypes(1)
        BF.set_hashes(7)
        BF.set_k(20)
        path = r'I:/OXA-Gene/' + files[i]
        name = files[i][:-6] + '.txt'
        print(name)
        result = r'C:/Users/SG/Desktop/' + name
        BF.train_sequence(path, 0)
        BF.save_clonetypes(result)
        BF.cleanup()


def train_Core():
    files = os.listdir(r'I:\Cores')
    for i in range(len(files) -1, -1, -1):
        if 'fna' not in files[i]:
            del files[i]

    for i in range(len(files)):

        BF = BF_v2.AbaumanniiBloomfilter(100000000)
        BF.set_arraysize(100000000)
        BF.set_clonetypes(1)
        BF.set_hashes(7)
        BF.set_k(20)
        path = r'I:/Cores/' + files[i]
        name = files[i][:-4] + '.txt'
        print(name)
        result = r'C:/Users/SG/Desktop/' + name
        BF.train_sequence(path, 0)
        BF.save_clonetypes(result)
        BF.cleanup()


def opene():
    with open(r'C:\Users\SG\Desktop\a.baumannii Filter\FilterOXA.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)
    print(clonetypes)


def pw():
    from flask_bcrypt import Bcrypt
    bcrypt = Bcrypt()
    print(bcrypt.generate_password_hash('user'))
    print(bcrypt.generate_password_hash('pwd'))

def Test():

    temp = bitarray(0)
    with open(r'filter\OXA51_IC1.txt', 'rb') as fh:
        temp.fromfile(fh)
        print(len(temp))


def Test_Core_for_OXA():
    with open(r'filter/FilterClonetypes.txt', 'rb') as fp:
        clonetypes = pickle.load(fp)

    BF = BF_v2.AbaumanniiBloomfilter(22000000)
    BF.set_arraysize(22000000)
    BF.set_hashes(7)
    BF.set_k(20)
    # User Options
    BF.set_reads(1000)

    paths = [r'filter/CoreIC1.txt',
             r'filter/CoreIC2.txt',
             r'filter/CoreIC3.txt',
             r'filter/CoreIC4.txt',
             r'filter/CoreIC5.txt',
             r'filter/CoreIC6.txt',
             r'filter/CoreIC7.txt',
             r'filter/CoreIC8.txt']

    BF.read_clonetypes(paths, clonetypes)

    Oxa_paths = [r'H:\bla-51-like\IC1\OXA69.fasta',
                 r'H:\bla-51-like\IC1\OXA92.fasta',
                 r'H:\bla-51-like\IC1\OXA107.fasta',
                 r'H:\bla-51-like\IC1\OXA110.fasta',
                 r'H:\bla-51-like\IC2\OXA66.fasta',
                 r'H:\bla-51-like\IC2\OXA82.fasta',
                 r'H:\bla-51-like\IC2\OXA172.fasta',
                 r'H:\bla-51-like\IC2\OXA201.fasta',
                 r'H:\bla-51-like\IC2\OXA202.fasta',
                 r'H:\bla-51-like\IC3\OXA71.fasta',
                 r'H:\bla-51-like\IC3\OXA113.fasta',
                 r'H:\bla-51-like\IC4\OXA51.fasta',
                 r'H:\bla-51-like\IC4\OXA219.fasta',
                 r'H:\bla-51-like\IC5\OXA65.fasta',
                 r'H:\bla-51-like\IC6\OXA90.fasta',
                 r'H:\bla-51-like\IC6\OXA200.fasta',
                 r'H:\bla-51-like\IC7\OXA64.fasta',
                 r'H:\bla-51-like\IC8\OXA68.fasta',
                 r'H:\bla-51-like\IC8\OXA128.fasta']

    for path in Oxa_paths:
        BF.lookup_sequence(path)
        score = BF.get_score()
        print(score)


def main():
    #Test_Core_for_OXA()
    #write_file()
    #opene()
    #write_file2()
    pw()
    #train_Core()


if __name__ == '__main__':
    main()
