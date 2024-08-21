#!/usr/bin/env python3
import subprocess
import numpy as np

N_min = -5

#while True would stop for an integer. This way the loop goes on until N_min > 0
while (N_min <= 0): 
	try:
		N_min = int( input( "Insert N_min, starting matrix dimension"+'\n'))
		#if the input is valid (integer), break the loop
	except ValueError:
		print("Invalid input: integers only, retry"+'\n')
		continue
	if (N_min <= 0):
		print( "Invalid input: dimension should be positive, retry"+'\n')
	#print("Invalid input: dimension should be a positive integer,retry")

N_max = -5

while (N_max <= N_min):
	try:
		N_max = int( input( "Insert N_max, maximum matrix dimension \n"))
	except:
		print("Invalid input: integers only, retry"+'\n')
	if (N_max <= N_min):
		print( "Invalid input: N_max should be bigger than N_min, retry"+'\n')

#https://numpy.org/doc/stable/reference/generated/numpy.logspace.html#numpy.logspace 
#maybe geomspace was better; from base**start to base**stop; num =how many steps

Dimension = np.int32(np.round(np.logspace(np.log10(N_min), np.log10(N_max), num = 50, base = 10)))

#pass 5, newline 5, newline two times for matrix dimensions. Text = True allows for input string; don't put stdin if there is input

for dim in Dimension:
	#watch out for numpy ruining integers
	string = str(int(dim)) + ' \n'
	proc = subprocess.run(['./matrixmult'], text = True, input = 4*string  )
	print(proc)


