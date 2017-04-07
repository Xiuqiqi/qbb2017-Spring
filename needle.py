#!/usr/bin/env python

"""
Perform Needleman-Wunsch global alignment of two nucleotide sequences.
usage: needle.py sequence_1 sequence_2
"""

import sys
import numpy

# HoxD70 matrix of Chiaromonte, Yap, Miller 2002,
#              A     C     G     T
sigma = [ [   91, -114,  -31, -123 ],
          [ -114,  100, -125,  -31 ],
          [  -31, -125,  100, -114 ],
          [ -123,  -31, -114,   91 ] ]

gap = 300

def substitution(char1, char2):
    if char1=='A':
        i=0
    if char1=='C':
        i=1
    if char1=='G':
        i=2
    if char1=='T':
        i=3
        
    if char2=='A':
        j=0
    if char2=='C':
        j=1
    if char2=='G':
        j=2
    if char2=='T':
        j=3
    
    value=sigma[i][j]
       
    return value
    

def compute_matrix( s, t ):
    """
    Fill in similarity score matrix, keeping traceback.
    """
    n = len( s )
    m = len( t )

    # n, vertical -- s 
    # m, horizontal -- t
    
    dp_matrix = numpy.zeros( (n+1,m+1), float )
    tb_matrix = numpy.zeros( (n+1,m+1), int )
    
    #Initilize dp
    for i in range(1,m+1):
        dp_matrix[0,i]=dp_matrix[0,i-1]-gap
    
    for i in range(1,n+1):
        dp_matrix[i,0]=dp_matrix[i-1,0]-gap
        
    # initilize tb
    # 1 - left, 2 - up, 3 - diagonal
    for i in range(1,m+1):
        tb_matrix[0,i]=1 # left
    
    for i in range(1,n+1):
        tb_matrix[i,0]=2 # up
        
        
    # fill in the matrix
    for i in range(1,n+1):   # row - s
        for j in range(1,m+1): # column - t
            v = dp_matrix[i-1,j]-gap
            h = dp_matrix[i,j-1]-gap
            d = dp_matrix[i-1,j-1]+substitution(s[i-1],t[j-1])
            dp_matrix[i,j]=max(v,h,d)
            if max(v,h,d)==v:
                tb_matrix[i,j]=2
            elif max(v,h,d)==h:
                tb_matrix[i,j]=1
            elif max(v,h,d)==d:
                tb_matrix[i,j]=3
                
    #print dp_matrix
    #print tb_matrix
    return dp_matrix, tb_matrix

def print_alignment( dp_matrix, tb_matrix, s, t ):
    line_length=80
    n = len( s )
    m = len( t )
    # n, vertical -- s 
    # m, horizontal -- t
    
    row1 = ""
    row2 = ""
    score=dp_matrix[n,m]
    
    while dp_matrix[n,m]:
        
    # 1 - left, 2 - up, 3 - diagonal
        if tb_matrix[n,m] == 1:
            row1= "-"+row1
            row2= t[m-1]+row2
            n=n
            m=m-1  
        elif tb_matrix[n,m] == 2:
            row1= s[n-1]+row1
            row2= "-"+row2
            n=n-1
            m=m
        elif tb_matrix[n,m] == 3:
            row1= s[n-1]+row1
            row2= t[m-1]+row2
            n=n-1
            m=m-1
        #print n
        #print m
    
    linenum=int(len(row1)/line_length)
    #print linenum
    
    for i in range(0,linenum+1):
        print row1[i*line_length:(i+1)*line_length]
        print row2[i*line_length:(i+1)*line_length]
        print
    #print row1[linenum*line_length:]
    #print row2[linenum*line_length:]
    
    print "Score: ", score

def parse(file):
    file_=open(file)
    txt=[]
    while True:
        line = file_.readline().rstrip("\r\n")
        if line == "":
            break
        if line.startswith("#"):
            continue
        
        #fields=line.split("\t")

        #if fields[0] == "chrX":
        txt.append(line)
    return ''.join(txt)

def main():
    s = parse(sys.argv[1])
    t = parse(sys.argv[2])
    dp_matrix, tb_matrix = compute_matrix( s, t )
    print_alignment( dp_matrix, tb_matrix, s, t )

if __name__ == "__main__":
    main()