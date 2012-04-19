#!/usr/bin/python

import datetime

# Extracts last git changelog item into html for psicode

gitlog = '../../../.git/logs/refs/heads/master'  # path to .git/logs/refs/heads/master
gitlast = 'latest_trac_changeset2.txt'  # path to output file with latest changeset
githist = 'history_trac_changeset.txt'  # path to output file with last Nhist changesets
Nhist = 100

fgit = open(gitlog)
gitlines = fgit.readlines()
fgit.close()

flast = open(gitlast, 'w')
fhist = open(githist, 'w')

ichgst = -1
while ichgst >= (-1 * Nhist):

    try:
        fields = gitlines[ichgst].split()  # grab Nth most recent changeset
    except IndexError:
        break
    
    indx = 0
    while True:  # identify e-mail address
        if fields[indx].startswith('<') and fields[indx].endswith('>'): break
        if indx == 10: 
            print 'ERROR: Didn\'t find e-mail in changeset.\n'
            break
        indx += 1
    
    changeset = '[' + str(fields[1][0:6]) + ']'  # abbreviated git id
    email = fields[indx]  # don't want to post e-mail addresses
    humantime = datetime.datetime.fromtimestamp(int(fields[indx + 1])).strftime('%H:%M %m/%d/%Y')
    author = ' '.join(fields[2:indx])
    comment = ' '.join(fields[indx+4:])
    commentshort = comment[0:125]  # truncate comment string to fit in Trac Feed box
    
    fhist.write('<title>%s %s by %s:<br/>%s<br/></title>\n\n' % (humantime, changeset, author, comment))
    if ichgst == -1:
        flast.write('<title>%s %s by %s:<br/>%s<br/></title>\n\n' % (humantime, changeset, author, commentshort))

    ichgst -= 1

flast.close()
fhist.close()

