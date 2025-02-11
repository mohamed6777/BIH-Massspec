import os
import sys

sample = sys.argv[1]
input_r1 = sys.argv[2]
input_r2 = sys.argv[3]
output_r1 = sys.argv[4]
output_r2 = sys.argv[5]

cmd = f"umi_tools dedup --stdin={input_r1} --stdout={output_r1} --paired --log={output_r1}.log"
os.system(cmd)

cmd = f"umi_tools dedup --stdin={input_r2} --stdout={output_r2} --paired --log={output_r2}.log"
os.system(cmd)
