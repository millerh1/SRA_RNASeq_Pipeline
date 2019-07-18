import synapseclient
from synapseutils import walk
import csv
import sys


synID=sys.argv[1]
username=sys.argv[2]
password=sys.argv[3]
outpath=sys.argv[4]
outfile=outpath+"/synapse.csv"

syn = synapseclient.Synapse()
syn.login(username, password)

walkedPath=walk(syn, synID)
for filename in walkedPath:
    filename2=filename[2]
    with open(outfile,'w') as out:
      writer = csv.writer(out, delimiter=',', lineterminator='\n')
      writer.writerow(('filename','synID'))
      for row in filename2:
          writer.writerow(row)











