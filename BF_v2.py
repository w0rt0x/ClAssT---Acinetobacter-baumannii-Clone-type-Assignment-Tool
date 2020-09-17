import mmh3
from bitarray import bitarray
from Bio import SeqIO
from copy import deepcopy
from OXA_Table import OXATable
from Bio.Seq import Seq


class AbaumanniiBloomfilter:
    """ Bloomfilter that can read FASTA and FASTQ files to assign the given file to a reference-genome"""
    # Implementation of the Bloomfilter Project for Acinetobacter baumannii
    # Code partly from https://github.com/Phelimb/BIGSI

    clonetypes = 1  # Number of IC's
    hits_per_filter = [0] * clonetypes  # Hit counter per IC
    array_size = 22000000  # Standard arraysize per IC is 22mio for Core-genome
    hashes = 7  # Number of used Hash-functions
    k = 20  # length of the k-meres
    names = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7', 'IC8']  # names if the IC's
    number_of_kmeres = 0  # counter of k-meres, will be used to calculate score
    reads = 1000  # standard read number

    def __init__(self, arraysize):
        """ creates empty matrix"""
        pass
        self.matrix = bitarray(arraysize)
        self.matrix.setall(False)
        self.array_size = arraysize

    # Setter

    def set_arraysize(self, new):
        """ changes Arraysize to new input-value, does not recreate matrix"""
        self.array_size = new

    def set_clonetypes(self, new):
        """ changes number of Clonetypes"""
        self.clonetypes = new
        self.hits_per_filter = [0] * self.clonetypes

    def set_hashes(self, new):
        """Changes number of used hash-functions"""
        self.hashes = new

    def set_k(self, new):
        """ Changes length of k-meres"""
        self.k = new

    def set_names(self, new):
        """ Changes Names of Filters, Input must be a List of names"""
        self.names = new

    def reset_counter(self):
        """resets counter"""
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes

    def set_reads(self, new):
        """ Changes number of reads to new value"""
        self.reads = new

    # Getter

    def get_score(self):
        """calculates score for all clonetypes
            Score is #hits / #kmeres"""

        score = []

        # calculates float for each value in [hits per filter]
        for i in range(self.clonetypes):
            if self.hits_per_filter[i] == 0:
                score.append(0.0)
            else:
                score.append(round(float(self.hits_per_filter[i]) / float(self.number_of_kmeres), 2))

        return score

    def get_norm(self):
        """ Divides each vector entry by sum of vector entrys"""
        s = sum(self.hits_per_filter)
        score = []

        # calculates float for each value in [hits per filter]
        for i in range(self.clonetypes):
            if self.hits_per_filter[i] == 0 or s == 0:
                score.append(0.0)
            else:
                score.append(round(float(self.hits_per_filter[i]) / s, 2))

        return score

    def get_reads(self):
        """ gets number of reads """
        return self.reads

    def get_hits_per_filter(self):
        """gets Hits per Filter"""
        return self.hits_per_filter

    def get_kmeres_per_sequence(self):
        """gets K-mer counter"""
        # returns number of k-meres per file
        return self.number_of_kmeres

    def get_names(self):
        """ gets names of filters"""
        return self.names

    # File management

    def save_clonetypes(self, path):
        """saves matrix as a binary file to the input-path"""
        # saving filters of clonetypes

        # creating file and saving matrix with the bitarray modul
        with open(path, 'wb') as fh:
            # writing to file with bitarray command
            self.matrix.tofile(fh)

    def read_clonetypes(self, paths, names):
        """ reads slices from files and concatsthem to a matrix,
        paths is list of paths and names is a string list"""

        # Updating parameters
        self.clonetypes = len(paths)
        self.names = names
        self.matrix = bitarray(0)
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes

        # creating matrix from single filters
        for path in paths:

            temp = bitarray()

            with open(path, 'rb') as fh:
                temp.fromfile(fh)
            self.matrix.extend(temp)

    # Bloomfilter

    def hash(self, kmer):
        """Hashes given string and returns Positions for the Array"""

        # Empty list for Array positions
        positions = []

        # Creating hashes for needed number of hash functions
        for i in range(self.hashes):
            # mmh3 takes that string and a seed,
            # each hash function takes an individual seed
            # after that, the hash-value will me divided by the array size until
            # a position in the array is guaranteed
            positions.append(mmh3.hash(kmer, i) % self.array_size)

        return positions

    def lookup(self, kmer, limit=False):
        """checks if an element is in the filters, returns list with True/False,
           takes kmer input string and checks all clonetypes if the k-mer is inside that set of kmers"""

        # getting positions
        positions = self.hash(str(kmer))

        # control if element is in filter
        hits = [True] * self.clonetypes

        for i in range(self.clonetypes):
            # all 7 Positions are hardcoded, the number of hashes is always(!) 7
            # if all positions  are True, then hits[i] will also stay True
            # (i*self.array_size) skips to the same position in the next filter
            hits[i] = (self.matrix[positions[0] + (i*self.array_size)] &
                       self.matrix[positions[1] + (i*self.array_size)] &
                       self.matrix[positions[2] + (i*self.array_size)] &
                       self.matrix[positions[3] + (i*self.array_size)] &
                       self.matrix[positions[4] + (i*self.array_size)] &
                       self.matrix[positions[5] + (i*self.array_size)] &
                       self.matrix[positions[6] + (i*self.array_size)])

            if hits[i]:
                if limit:
                    if self.table.lookup(self.names[i], kmer):
                        self.hits_per_filter[i] += 1
                else:
                    # Update hit counter
                    self.hits_per_filter[i] += 1

    def train(self, kmer, clonetype):
        """ trains specific filter for a k-mer, input is that kmer and the desired Filter"""

        # getting hash Values
        positions = self.hash(kmer)
        # changing 0s to 1 in filter
        for i in range(len(positions)):
            # getting position of cell
            self.matrix[self.array_size * clonetype + positions[i]] = True

    def train_sequence(self, filepath, clonetype):
        """trains whole sequence into filter, takes filepath to file and the desired filter as input"""
        # for each sequence (in multi-FASTA file)
        for sequence in SeqIO.parse(filepath, "fasta"):
            # for each k-mere
            for i in range(len(sequence.seq) - self.k + 1):
                # trains k-mere into filter
                self.train(str(sequence.seq[i: i + self.k]), clonetype)

    def train_lines(self, lines, ct):
        """ Trains Extracted lines of fasta/fna file, given as list of strings"""
        for j in range(len(lines)):
            for i in range(len(lines[j]) - self.k + 1):
                # trains k-mere into filter
                self.train(str(lines[j][i: i + self.k]), ct)

    def lookup_sequence(self, path):
        """uses lookup function for whole sequence, takes path to file: file must be FASTA"""

        # Counter of k-meres
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes

        # for each sequence (in multi-FASTA file)

        # counting hits in Intervals
        self.readset = [0] * self.clonetypes

        for sequence in SeqIO.parse(path, "fasta"):

            # for each k-mere
            for i in range(0, len(sequence.seq) - self.k + 1):

                # lookup for all k-meres in filter
                self.lookup(str(sequence.seq[i: i + self.k]))
                self.number_of_kmeres += 1

    def lookup_txt(self, reads, quick=False):
        """ Reading extracted fq-reads"""
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes

        if quick:
            # Quick: Non-overlapping k-mers
            for single_read in reads:
                # r is rest, so all kmers have size k
                for j in range(0, len(single_read) - self.k, 10):
                    self.number_of_kmeres += 1
                    self.lookup(single_read[j: j + self.k])
        else:
            for single_read in reads:
                for j in range(len(single_read) - self.k + 1):
                    # updating counter
                    self.number_of_kmeres += 1
                    # lookup for kmer
                    self.lookup(str(single_read[j: j + self.k]))

        reads = None

    def cleanup(self):
        """deletes matrix"""
        self.matrix = None

    def lookup_oxa(self, reads, ext):
        """ Looks for OXA Genes: Extension (ext) selects the fq-seach or fasta-search mode"""
        table = OXATable()
        table.read_dic(r'filter/OXAs_dict/oxa_dict.txt')
        if ext == 'fq':
            # fq mode
            self.table = table
            for i in range(len(reads)):
                # going through all reads, discarding those who don't get any hits with 3 test k-meres

                # Building 3 test-kmeres: first, last, and middle
                k1 = reads[i][0:self.k]  # first k-mer
                k2 = reads[i][len(reads[i]) - self.k:]  # last k-mer
                mid = len(reads[i])//2
                k3 = reads[i][mid:mid+self.k]  # k-mer in middle

                # Taking sum of list as reference, if sum has not increased after testing those 3 kmeres,
                # then the read won't be tested further
                hit_sum = sum(self.hits_per_filter)
                copy = deepcopy(self.hits_per_filter)
                self.lookup(k1, True)
                self.lookup(k2, True)
                self.lookup(k3, True)

                # needs at least 2 of 3 hits to continue with read
                if (sum(self.hits_per_filter) - hit_sum) > 1:

                    for j in range(1, len(reads[i]) - 1 - self.k + 1):
                        # Skipping first, last and middle k-mer
                        if j != mid:
                            self.lookup(reads[i][j:j + self.k], True)
                            self.number_of_kmeres += 1

                else:
                    # resetting hit counter
                    self.hits_per_filter = copy

                # same, but with reverse complement
                reads[i] = Seq(reads[i])
                reads[i] = reads[i].reverse_complement()
                k1 = reads[i][0:self.k]  # first k-mer
                k2 = reads[i][len(reads[i]) - self.k:]  # last k-mer
                mid = len(reads[i]) // 2
                k3 = reads[i][mid:mid + self.k]  # k-mer in middle

                # Taking sum of list as reference, if sum has not increased after testing those 3 kmeres,
                # then the read won't be tested further
                hit_sum = sum(self.hits_per_filter)
                copy = deepcopy(self.hits_per_filter)
                self.lookup(k1, True)
                self.lookup(k2, True)
                self.lookup(k3, True)

                # needs at least 2 of 3 hits to continue with read
                if (sum(self.hits_per_filter) - hit_sum) > 1:

                    for j in range(1, len(reads[i]) - 1 - self.k + 1):
                        # Skipping first, last and middle k-mer
                        if j != mid:
                            self.lookup(reads[i][j:j + self.k], True)
                            self.number_of_kmeres += 1

                else:
                    # resetting hit counter
                    self.hits_per_filter = copy

        else:
            # fasta mode

            for r in range(len(reads)):
                for j in range(len(reads[r]) - self.k + 1):
                    # New Lookup, Only k-meres that are found will be counted
                    positions = self.hash(reads[r][j:j+self.k])

                    # control if element is in filter
                    hits = [True] * self.clonetypes

                    for i in range(self.clonetypes):
                        # all 7 Positions are hardcoded, the number of hashes is always(!) 7
                        # if all positions  are True, then hits[i] will also stay True
                        # (i*self.array_size) skips to the same position in the next filter
                        hits[i] = (self.matrix[positions[0] + (i * self.array_size)] &
                                   self.matrix[positions[1] + (i * self.array_size)] &
                                   self.matrix[positions[2] + (i * self.array_size)] &
                                   self.matrix[positions[3] + (i * self.array_size)] &
                                   self.matrix[positions[4] + (i * self.array_size)] &
                                   self.matrix[positions[5] + (i * self.array_size)] &
                                   self.matrix[positions[6] + (i * self.array_size)])
                        if hits[i] and table.lookup(self.names[i], reads[r][j:j+self.k]):
                            # Update hit counter
                            self.hits_per_filter[i] += 1

                # same, but with reverse complement
                reads[r] = Seq(reads[r])
                reads[r] = reads[r].reverse_complement()
                for j in range(len(reads[r]) - self.k + 1):
                    # New Lookup, Only k-meres that are found will be counted
                    positions = self.hash(str(reads[r][j:j+self.k]))

                    # control if element is in filter
                    hits = [True] * self.clonetypes

                    for i in range(self.clonetypes):
                        # all 7 Positions are hardcoded, the number of hashes is always(!) 7
                        # if all positions  are True, then hits[i] will also stay True
                        # (i*self.array_size) skips to the same position in the next filter
                        hits[i] = (self.matrix[positions[0] + (i * self.array_size)] &
                                   self.matrix[positions[1] + (i * self.array_size)] &
                                   self.matrix[positions[2] + (i * self.array_size)] &
                                   self.matrix[positions[3] + (i * self.array_size)] &
                                   self.matrix[positions[4] + (i * self.array_size)] &
                                   self.matrix[positions[5] + (i * self.array_size)] &
                                   self.matrix[positions[6] + (i * self.array_size)])
                        if hits[i] and table.lookup(self.names[i], reads[r][j:j+self.k]):
                            # Update hit counter
                            self.hits_per_filter[i] += 1

        reads = None
        table.cleanup()

    def get_oxa_score(self):
        """ Returning hits per OXA/kmere in OXA-filter"""
        table = OXATable()
        counter = table.get_counter()
        score = []
        # calculates float for each value in [hits per filter]
        for i in range(self.clonetypes):
            if self.hits_per_filter[i] == 0:
                score.append(0.0)
            else:
                score.append(round(float(self.hits_per_filter[i]) / float(counter[self.names[i]]), 2))

        return score





