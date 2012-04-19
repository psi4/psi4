#!/usr/bin/python

# Extracts last git changelog item into html for psicode

gitlog = '../../../.git/logs/refs/heads/master'  # path to .git/logs/refs/heads/master
fgit = open(gitlog)
gitlines = fgit.readlines()
fgit.close()

fields = gitlines[-1].split()  # grab most recent changeset

changeset = '[' + str(fields[1][0:6]) + ']'  # abbreviated git id
    
indx = 0
while True:  # identify e-mail address
    if fields[indx].startswith('<') and fields[indx].endswith('>'): break
    if indx == 10: 
        print 'ERROR: Didn\'t find e-mail in changeset.\n'
        break
    indx += 1

email = fields[indx]  # don't want to post e-mail addresses
author = ' '.join(fields[2:indx])
comment = ' '.join(fields[indx+4:])
comment = comment[0:125]  # truncate comment string to fit in Trac Feed box

fhtml = open('latest_trac_changeset2.txt', 'w')  # write file for index.html to include
fhtml.write('<title>%s by %s:<br/>%s<br/></title>\n\n' % (changeset, author, comment))
fhtml.close()

