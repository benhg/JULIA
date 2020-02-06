import os
import csv

#create csv: 

csv_filename = "index_hopping_output.csv"
csv_fh = open(csv_filename, "write")
csv_writer = csv.DictWriter(csv_fh, fieldnames=["filename", "uniquely", "multi", "totalReads", "uniquelyAGAINST", "multiAGAINST", "totalReadsAGAINST", "percratio"])
csv_row = {"filename": "filename ", "uniquely": "uniquely ", "multi": "multi ", "totalReads": "totalreads ", "uniquelyAGAINST": "uniquelyC ", "multiAGAINST": "multiC ", "totalReadsAGAINST": "totalreadsC ", "percratio": "percRatio "}
csv_writer.writerow(csv_row)


#if its taxon labeling was correct -- source == label
csv_row = {"filename": str(filename), "uniquely": str(uniquely), "multi": str(multi), "totalReads": str(totalreads), "uniquelyAGAINST": str(uniquelyC), "multiAGAINST": str(multiC), "totalReadsAGAINST": str(totalreadsC), "percratio": str(percRatio)}

csv_writer.writerow(csv_row)

#if taxon label did not match, maybe a hopper: (entering the hopper loop)
os.system("echo " + filename + ">> ../maybehopper.txt") 


#if it is a hopper: 
os.system("echo " + str(filename) + "  " +  str(uniquely) + "  " + str(multi) + "  " + str(totalreads) + "  " +  str(uniquelyC) + "  " +   str(multiC) + "  " + str(totalreadsC) + "  " +  str(percRatio)   + "  >> ../hopper.txt")


#all are added here i believe (a text version of the csv?)
os.system("echo 'filename  uniquely  multi  totalReads  uniquelyAGAINST  multiA\
GAINST  totalreadsAGAINST  percRatio'  >> index_hopping_output.txt")

os.system("echo " + str(filename) + "  " +  str(uniquely) + "  " + str(multi) + "  " + str(totalreads) + "  " +  str(uniquelyC) + "  " +   str(multiC) + "  " + str(totalreadsC) + "  " +  str(percRatio)   + "  >> ../index_hopping_output.txt")

#if it never found another source: 
os.system("echo " + filename + " >> ../never.txt")
