#!/usr/bin/python

import sqlite3

conn = sqlite3.connect('hmdb.db')

lL = 150.10
uL = 150.15
condition = "MONOMASS >= %f AND MONOMASS < %f" % (lL, uL)
command = "SELECT * FROM HMDB WHERE %s" % condition
cursor = conn.execute(command)

for row in cursor:
    print (row[0], row[1], row[2], row[3], row[4], row[5], row[6])

conn.close()