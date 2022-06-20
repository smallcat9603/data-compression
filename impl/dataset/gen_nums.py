'''
Created on 2022/06/20

@author: huyao
'''

num = 8192
fn = "/Users/smallcat/Documents/GitHub/data-compression/impl/dataset/float_" + str(num) + ".txt"
for i in range(num):
    s = "0.123456789\n"
    f = open(fn, "a")
    f.write(s)
    f.close()