import PsiMod
import sys
from psiexceptions import *

class Table:

    def __init__(self, rows=(),
                 row_label_width=10,
                 row_label_precision=4,
                 cols=(),
                 width=16, precision=10):
        self.row_label_width = row_label_width
        self.row_label_precision = row_label_precision
        self.width = width
        self.precision = precision
        self.rows = rows

        if isinstance(cols, str):
            self.cols = (cols,)
        else:
            self.cols = cols

        self.labels = []
        self.data = []

    def format_label(self):
        #str = lambda x: (('%%%d.%df' % (self.row_label_width, self.row_label_precision)) % x)
        str = lambda x: (('%%%ds' % (self.row_label_width)) % x)
        return " ".join(map(str, self.labels))

    def format_values(self, values):
        str = lambda x: (('%%%d.%df' % (self.width, self.precision)) % x)
        return " ".join(map(str, values))

    def __getitem__(self, value):
        self.labels.append(value)
        return self

    def __setitem__(self, name, value):
        self.labels.append(name)
        label = self.format_label()
        self.labels = []

        if isinstance(value, list):
            self.data.append( (label, value ) )
        else:
            self.data.append( (label, [value] ) )

    def save(self, file):
        import pickle
        pickle_str = pickle.dumps(self)
        fileobj = open(file, "w")
        fileobj.write(str(self))
        fileobj.close()

    def __str__(self):
        rowstr = lambda x: '%%%ds' % self.row_label_width % x
        colstr = lambda x: '%%%ds' % self.width % x

        lines = []

        table_header = ""
        if isinstance(self.rows, str):
            table_header += "%%%ds" % self.row_label_width % self.rows
        else:
            table_header += " ".join(map(rowstr, self.rows))
        table_header += " ".join(map(colstr, self.cols))

        lines.append(table_header)

        for datarow in self.data:
            #print datarow
            row_data = datarow[0]
            row_data += self.format_values(datarow[1])
            lines.append(row_data)

        return "\n".join(lines) + "\n"

    def copy(self):
        import copy
        return copy.deepcopy(self)

    def absolute_to_relative(self, Factor = psi_hartree2kcalmol):
        import copy

        if len(self.data) == 0:
            return

        current_min = list(copy.deepcopy(self.data[0][1]))
        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                if current_min[col] > datarow[1][col]:
                    current_min[col] = datarow[1][col]

        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                #print datarow[1][col]
                datarow[1][col] = (datarow[1][col] - current_min[col]) * Factor

    def scale(self, Factor = psi_hartree2kcalmol):

        if len(self.data) == 0:
            return

        for datarow in self.data:
            for col in range(0, len(datarow[1])):
                #print datarow[1][col]
                datarow[1][col] = datarow[1][col] * Factor

def banner(text, type = 1, width = 35):
    lines = text.split('\n')
    max_length = 0
    for line in lines:
        if (len(line) > max_length):
            max_length = len(line)
    
    max_length = max([width, max_length]);
   
    null = '' 
    if type == 1:
        banner  = '  //' + null.center(max_length,'>') + '//\n'
        for line in lines:
            banner += '  //' + line.center(max_length) + '//\n'
        banner += '  //' + null.center(max_length,'<') + '//\n'
    
    if type == 2:
        banner = ''
        for line in lines:
            banner += (' ' + line + ' ').center(max_length,'=') 

    PsiMod.print_out(banner)

def print_stdout(stuff):
    print >> sys.stdout, stuff

def print_stderr(stuff):
    print >> sys.stderr, stuff
