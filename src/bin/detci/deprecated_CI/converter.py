import glob

rep_codes = """vectlen
num_blocks
icore
Ms0
Ia_code
Ib_code
Ia_size
Ib_size
offset
num_alpcodes
num_betcodes
nirreps
codes_per_irrep
buf_per_vect
buf_total
new_first_buf
maxvect
nvect
nunits
units
cur_vect
cur_buf
buffer_size
file_number
buf_size
buf2blk
buf_offdiag
first_ablk
last_ablk
decode
zero_blocks
blocks
buf_locked
buffer
in_file
extras
units_used
cur_unit
cur_size
first_unit"""

rep_codes = rep_codes.split()

#file_list = ['civect.cc']
file_list = ['opdm.cc', 'tpdm.cc', 'sigma.cc', 'mitrush_iter.cc']
for infile in file_list:
    f = open(infile)
    data = f.readlines()
    f.close()
    
    
    cnt = 0
    out = []
    for ln, line in enumerate(data):
        line = line[:-1] 
        other   = False
        #other  |= ('Param' in line)
        #other  |= ('/*' in line)
        #other  |= ('//' in line)
        #other  |= ('CalcIn' in line)
        #other  |= ('CIblks' in line)
    
        if ('**' in line) or ('::' in line) or other:
            out.append(line)
            continue
    
        for rep in rep_codes:
            pos = line.find(rep)
            if pos == -1: continue
            if line[pos-2:pos] == '->': continue
    
            cnt +=1
            for x in ['Ivec.', 'Jvec.', 'Dvec.', 'Hd.', 'C.', 'S.']:
                trep = x + rep
                line = line.replace(trep, trep + '_').replace(trep+'__', trep+'_')

            #line = line.replace(rep, rep + '_').replace(rep+'__', rep+'_')
        out.append(line)
    
    
    out = '\n'.join(out)
    f = open(infile, 'w')
    f.write(out)
    f.close()


