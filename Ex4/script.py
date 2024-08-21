#!/usr/bin/env python3
import subprocess
import numpy as np


a = 10
omega = 1

for N in range(1000, 2501, 500):

	string =  str(a) + ' \n'  + str(int(N)) +  ' \n' + str(omega) + ' \n' 
	proc = subprocess.run(['./harm'], text = True, input = string  )
	print(proc)
